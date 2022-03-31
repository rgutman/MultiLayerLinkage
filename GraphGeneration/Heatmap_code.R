rm(list = ls())
library(lubridate)
library(MCMCpack)
library(ggplot2)
library(reshape2)
library(clue)

########################################################################################################################
########Uses posterior parameter estimates from BRL, CIBRL, and MLBRL for one simulated dataset under 40% error in region, area income, and month of birth
########to generate heatmaps of the relative log-probabilities from the block-level and record-level confusion matrices.
########Day of birth is omitted as a linking variable.
########################################################################################################################


##########################################################################################################################
#The simulated data used to generate the following figures will derive from one of the simulated datasets with 
#40 percent error in Region, Area Income, and DOB month

set.seed(1234)
seed=round(runif(100,1,100000000))
set.seed(seed[1])
##############################################################################################################
################Generate Simulated Data

############################################################
##########Linking Variables
################################################
###Simulate Data for True Links
#DOB
V1 = rnorm(450,30,2)
V1 = 2018-V1
DOB = format(date_decimal(V1), "%d-%m-%Y")
DOB_year = format(date_decimal(V1), "%Y")
DOB_month = format(date_decimal(V1), "%m")
DOB_day = format(date_decimal(V1), "%d")

#Gender
Gender=round(runif(450,1,2))

#V4: Unique ID number for the linked individuals
MatchID = round(runif(450,1000000,10000000))

#Merge Linked Data together
LinkData.A = data.frame(MatchID, Gender, DOB_year, DOB_month, DOB_day)
LinkData.B = data.frame(MatchID, Gender, DOB_year, DOB_month, DOB_day)

############################################################
###Simulate additional unlinked records in File A
#Assume N_A=600, need to generate data for 150 unlinked records

#Distribution of X variables will be the same as the linked data
V1_A = rnorm(150,30,2)
V1_A = 2018-V1_A
DOB_A = format(date_decimal(V1_A), "%d-%m-%Y")
DOB_year_A = format(date_decimal(V1_A), "%Y")
DOB_month_A = format(date_decimal(V1_A), "%m")
DOB_day_A = format(date_decimal(V1_A), "%d")
Gender_A = round(runif(150,1,2))
MatchID_A = round(runif(150,1000000,10000000))

#Merge Unlinked data in file A together
N_A_unlinked = data.frame(MatchID_A, Gender_A, DOB_year_A, DOB_month_A, DOB_day_A)
names(N_A_unlinked)=names(LinkData.A)

#Combine Linked and Unlinked file A data to form dataset A
N_A = rbind(LinkData.A, N_A_unlinked)

############################################################
###Simulate additional unlinked records in File B
#Assume N_B=1200, need to generate data for 750 unlinked records

#Distribution of variables will be the same as the linked data
V1_B = rnorm(750,30,2)
V1_B = 2018-V1_B
DOB_B = format(date_decimal(V1_B), "%d-%m-%Y")
DOB_year_B = format(date_decimal(V1_B), "%Y")
DOB_month_B = format(date_decimal(V1_B), "%m")
DOB_day_B = format(date_decimal(V1_B), "%d")
Gender_B = round(runif(750,1,2))
MatchID_B = round(runif(750,1000000,10000000))

#Merge Unlinked data in file B together
N_B_unlinked = data.frame(MatchID_B, Gender_B, DOB_year_B, DOB_month_B, DOB_day_B)
names(N_B_unlinked) = names(LinkData.B)

#Combine Linked and Unlinked file B data to form dataset B
N_B = data.frame(rbind(LinkData.B, N_B_unlinked))
levels(N_A$DOB_year)=unique(c(levels(N_A$DOB_year),levels(N_B$DOB_year)))
levels(N_B$DOB_year)=unique(c(levels(N_A$DOB_year),levels(N_B$DOB_year)))

############################################################
##########Blocking Variables
############################################################
###Blocking Variables for True Blocks
Region_A = sample(c("N","W","MW","S"),size=30,replace=T)
Status_A = rbinom(30,1,.8)
Type_A = rbinom(30,1,.5)
SES_A=sample(c("L","M","H"),size=30,replace=T,prob=c(.25,.5,.25))
Income_A=rnorm(30,50000,10000)
Block_A_ID=round(runif(30,100,10000))
Block_A_S=seq_along(Block_A_ID)

###Blocking Variables for additional non-blocks in Data B
Region = sample(c("N","W","MW","S"),size=10,replace=T)
Status = rbinom(10,1,.8)
Type = rbinom(10,1,.5)
SES=sample(c("L","M","H"),size=10,replace=T,prob=c(.25,.5,.25))
Income=rnorm(10,50000,10000)

Block_B_ID=c(Block_A_ID,round(runif(10,100,10000)))
Block_B_T=seq_along(Block_B_ID)

Region_B=c(Region_A,Region)
Status_B=c(Status_A,Status)
Type_B=c(Type_A,Type)
SES_B=c(SES_A,SES)
Income_B=c(Income_A,Income)

###Append Blocking ID to each record pair
N_A$Block_ID=c(rep(Block_A_ID,each=15),rep(Block_A_ID,each=5))
N_B$Block_ID=c(rep(Block_B_ID,each=15),rep(Block_B_ID,each=15))
N_A$S=c(rep(Block_A_S,each=15),rep(Block_A_S,each=5))
N_B$T=c(rep(Block_B_T,each=15),rep(Block_B_T,each=15))

############################################################
###Introduce Error into entries from dataset A
############################################################
#Error in Region
Regionerrorprob=.4
Regionerror=rbinom(length(Region_A),1,Regionerrorprob)

for(i in 1:length(Region_A)){
  if(Regionerror[i]==1){
    Region_A[i]=sample(levels(as.factor(Region_A))[levels(as.factor(Region_A))!=Region_A[i]],1)
  }
}

#Error in Income
Income_A=Income_A+rnorm(length(Income_A),0,sd=500/qnorm(1-.4/2))

#Error in DOB
DOBerrorprob=.4
DOBerror=rbinom(length(N_A$DOB_month),1,DOBerrorprob)
for(i in 1:length(N_A$DOB_month)){
  if(DOBerror[i]==1){
    N_A$DOB_month[i]=sample(levels(N_A$DOB_month)[levels(N_A$DOB_month)!=N_A$DOB_month[i]],1)
  }
}


############################################################
##########Linking Variable Comparisons
############################################################
#Create Comparison matrix Gamma for each of the fields in A and B
DOB_year_comparison=t(matrix(N_B$DOB_year,nrow=length(N_B$DOB_year),ncol=nrow(N_A)))
DOB_month_comparison=t(matrix(N_B$DOB_month,nrow=length(N_B$DOB_month),ncol=nrow(N_A)))
DOB_day_comparison=t(matrix(N_B$DOB_day,nrow=length(N_B$DOB_day),ncol=nrow(N_A)))
Gender_comparison=t(matrix(N_B$Gender,nrow=length(N_B$Gender),ncol=nrow(N_A)))

#Compare each variable in dataset A with each element in dataset B element-wise
DOB_year_gamma=apply(DOB_year_comparison,2,'==',N_A$DOB_year)
DOB_month_gamma=apply(DOB_month_comparison,2,'==',N_A$DOB_month)
DOB_day_gamma=apply(DOB_day_comparison,2,'==',N_A$DOB_day)
Gender_gamma=apply(Gender_comparison,2,'==',N_A$Gender)

############################################################
##########Blocking Variable Comparisons
############################################################
#Create Comparison matrix Gamma for each of the fields in A and B
Region_comparison=t(matrix(Region_B,nrow=length(Region_B),ncol=length(Region_A)))
Status_comparison=t(matrix(Status_B,nrow=length(Status_B),ncol=length(Status_A)))
Type_comparison=t(matrix(Type_B,nrow=length(Type_B),ncol=length(Type_A)))
SES_comparison=t(matrix(SES_B,nrow=length(SES_B),ncol=length(SES_A)))
Income_comparison=t(matrix(Income_B,nrow=length(Income_B),ncol=length(Income_A)))

#Compare each variable in dataset A with each element in dataset B element-wise
Region_gamma=apply(Region_comparison,2,'==',Region_A)
Status_gamma=apply(Status_comparison,2,'==',Status_A)
Type_gamma=apply(Type_comparison,2,'==',Type_A)
SES_gamma=apply(SES_comparison,2,'==',SES_A)
Income_gamma=apply(Income_comparison,2,function(x) abs(x-Income_A)<=500)

############################################################
###Create Gamma matrices from element-wise comparisons
#Date of Birth
Gamma2=DOB_year_gamma-DOB_month_gamma
Gamma2[Gamma2!=1]=0
Gamma3=DOB_year_gamma+DOB_month_gamma-1
Gamma3[Gamma3!=1]=0
Gamma1=-(Gamma2+Gamma3)+1
Gamma1[Gamma1<=0]=0
#Gender
Gamma6=Gender_gamma
Gamma5=-(Gender_gamma)+1

#Region
Zeta1=-(Region_gamma)+1
Zeta2=Region_gamma
#Status
Zeta3=-(Status_gamma)+1
Zeta4=Status_gamma
#Type
Zeta5=-(Type_gamma)+1
Zeta6=Type_gamma
#SES
Zeta7=-(SES_gamma)+1
Zeta8=SES_gamma
#Income
Zeta9=-(Income_gamma)+1
Zeta10=Income_gamma


###Initialize blocking, linking, and blocking designation matrices
B_S=matrix(N_A$S,nrow=nrow(N_A),ncol=nrow(N_B))
B_T=t(matrix(N_B$T,nrow=nrow(N_B),ncol=nrow(N_A)))
B_ST=matrix(0,nrow=nrow(N_A),ncol=nrow(N_B))

B=matrix(0,nrow=length(Block_A_ID),ncol=length(Block_B_ID))

##############################################################################################################
###############Block Pair Log Probability Density Using CIBRL
###############Read in posterior estimates of Theta
theta_BM_Income=read.csv("theta_BM_Income_CIBRL_R4_I4_D4_noD_1.csv")
theta_BM_Region=read.csv("theta_BM_Region_CIBRL_R4_I4_D4_noD_1.csv")
theta_BM_Status=read.csv("theta_BM_Status_CIBRL_R4_I4_D4_noD_1.csv")
theta_BM_Type=read.csv("theta_BM_Type_CIBRL_R4_I4_D4_noD_1.csv")
theta_BM_SES=read.csv("theta_BM_SES_CIBRL_R4_I4_D4_noD_1.csv")
theta_BU_Income=read.csv("theta_BU_Income_CIBRL_R4_I4_D4_noD_1.csv")
theta_BU_Region=read.csv("theta_BU_Region_CIBRL_R4_I4_D4_noD_1.csv")
theta_BU_Status=read.csv("theta_BU_Status_CIBRL_R4_I4_D4_noD_1.csv")
theta_BU_Type=read.csv("theta_BU_Type_CIBRL_R4_I4_D4_noD_1.csv")
theta_BU_SES=read.csv("theta_BU_SES_CIBRL_R4_I4_D4_noD_1.csv")
theta_CM_DOB=read.csv("theta_CM_DOB_CIBRL_R4_I4_D4_noD_1.csv")
theta_CM_Gender=read.csv("theta_CM_Gender_CIBRL_R4_I4_D4_noD_1.csv")
theta_CU_DOB=read.csv("theta_CU_DOB_CIBRL_R4_I4_D4_noD_1.csv")
theta_CU_Gender=read.csv("theta_CU_Gender_CIBRL_R4_I4_D4_noD_1.csv")

################The estimate of each parameter will be the posterior mean of the last 100 iterations
theta_BM_Income_R=apply(theta_BM_Income[,c(400:500)],1,mean)
theta_BM_Region_R=apply(theta_BM_Region[,c(400:500)],1,mean)
theta_BM_Status_R=apply(theta_BM_Status[,c(400:500)],1,mean)
theta_BM_Type_R=apply(theta_BM_Type[,c(400:500)],1,mean)
theta_BM_SES_R=apply(theta_BM_SES[,c(400:500)],1,mean)

theta_BU_Income_R=apply(theta_BU_Income[,c(400:500)],1,mean)
theta_BU_Region_R=apply(theta_BU_Region[,c(400:500)],1,mean)
theta_BU_Status_R=apply(theta_BU_Status[,c(400:500)],1,mean)
theta_BU_Type_R=apply(theta_BU_Type[,c(400:500)],1,mean)
theta_BU_SES_R=apply(theta_BU_SES[,c(400:500)],1,mean)

theta_CM_DOB_R=apply(theta_CM_DOB[,c(400:500)],1,mean)
theta_CM_Gender_R=apply(theta_CM_Gender[,c(400:500)],1,mean)

theta_CU_DOB_R=apply(theta_CU_DOB[,c(400:500)],1,mean)
theta_CU_Gender_R=apply(theta_CU_Gender[,c(400:500)],1,mean)

###############Calculate the Ratio of P(Gamma_Bst|theta_BM)/P(Gamma_Bst|theta_BU) for the B matrix
num=(theta_BM_Income_R[1]^Zeta1)*(theta_BM_Income_R[2]^Zeta2)*(theta_BM_Status_R[1]^Zeta3)*(theta_BM_Status_R[2]^Zeta4)*
  (theta_BM_Type_R[1]^Zeta5)*(theta_BM_Type_R[2]^Zeta6)*(theta_BM_SES_R[1]^Zeta7)*(theta_BM_SES_R[2]^Zeta8)*(theta_BM_Income_R[1]^Zeta9)*(theta_BM_Income_R[2]^Zeta10)

den=(theta_BU_Income_R[1]^Zeta1)*(theta_BU_Income_R[2]^Zeta2)*(theta_BU_Status_R[1]^Zeta3)*(theta_BU_Status_R[2]^Zeta4)*
  (theta_BU_Type_R[1]^Zeta5)*(theta_BU_Type_R[2]^Zeta6)*(theta_BU_SES_R[1]^Zeta7)*(theta_BU_SES_R[2]^Zeta8)*(theta_BU_Income_R[1]^Zeta9)*(theta_BU_Income_R[2]^Zeta10)

#Generate Heatmap for Block Level probabilities for first 10 blocks using MLBRL
Block_Likelihood=t(apply((num/den),1,rev))
row_sum=apply(Block_Likelihood,1,sum)
Sum_mat=matrix(row_sum,nrow=30,ncol=40)
Block_LogLikelihood_CIBRL=log(Block_Likelihood/Sum_mat)
Block_LogLikelihood_CIBRL=Block_LogLikelihood_CIBRL[c(1:10),c(31:40)]
rownames(Block_LogLikelihood_CIBRL)=as.character(c(1:10))
colnames(Block_LogLikelihood_CIBRL)=as.character(sort(seq(1:10),decreasing=TRUE))
melted_Block_Likelihood_CIBRL=melt(Block_LogLikelihood_CIBRL)

png("CIBRL_block_small.png", width = 650, height = 550)
ggplot(data=melted_Block_Likelihood_CIBRL, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+coord_cartesian(xlim=c(1,10),ylim=c(1,10))+
  scale_x_continuous(breaks=seq(1:10))+scale_y_continuous(breaks=seq(1:10))+xlab("File 1 Block")+ylab("File 2 Block")+theme(legend.title=element_blank())+
  scale_fill_gradientn(limits = c(-30,0), colors = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(16)[6:16]))
dev.off()

###############Calculate the Ratio of P(Gamma_Cst|theta_CM)/P(Gamma_Cst|theta_CU) for the C matrix
num_C=(theta_CM_DOB_R[1]^Gamma1)*(theta_CM_DOB_R[2]^Gamma2)*(theta_CM_DOB_R[3]^Gamma3)*(theta_CM_Gender_R[1]^Gamma5)*(theta_CM_Gender_R^Gamma6)
den_C=(theta_CU_DOB_R[1]^Gamma1)*(theta_CU_DOB_R[2]^Gamma2)*(theta_CU_DOB_R[3]^Gamma3)*(theta_CU_Gender_R[1]^Gamma5)*(theta_CU_Gender_R^Gamma6)

#Generate Heatmap for Block Level probabilities for first 10 blocks using CIBRL
Block_ratio=num/den
Block_prob=Block_ratio[c(rep(1:nrow(Block_ratio), each=15), rep(1:nrow(Block_ratio), each=5)), c(rep(1:ncol(Block_ratio), each=15), rep(1:ncol(Block_ratio), each=15))]
Link_Likelihhood=t(apply((num_C/den_C*Block_prob),1,rev))
row_sum=apply(Link_Likelihhood,1,sum)
Sum_mat=matrix(row_sum,nrow=600,ncol=1200)
Link_LogLikelihood_CIBRL=log(Link_Likelihhood/Sum_mat)
Link_LogLikelihood_CIBRL=Link_LogLikelihood_CIBRL[c(1:150),c(1051:1200)]
rownames(Link_LogLikelihood_CIBRL)=as.character(c(1:150))
colnames(Link_LogLikelihood_CIBRL)=as.character(sort(seq(1:150),decreasing=TRUE))
melted_Link_Likelihood_CIBRL=melt(Link_LogLikelihood_CIBRL)

png("CIBRL_link_small.png", width = 650, height = 550)
ggplot(data=melted_Link_Likelihood_CIBRL, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+coord_cartesian(xlim=c(7,144),ylim=c(7,144))+
  scale_x_continuous(breaks=c(1,seq(15,150,length=10)))+scale_y_continuous(breaks=c(1,seq(15,150,length=10)))+xlab("File 1 Record")+ylab("File 2 Record")+ theme(legend.title=element_blank())+
  scale_fill_gradientn(limits = c(-30,0),
                       colours=c("navyblue", "darkmagenta", "darkorange1","red"))
dev.off()

##############################################################################################################
###############Block Pair Log Probability Density Using MLBRL
###############Read in posterior estimates of Theta
theta_BM_Income=read.csv("theta_BM_Income_MLBRL_R4_I4_D4_noD_1.csv")
theta_BM_Region=read.csv("theta_BM_Region_MLBRL_R4_I4_D4_noD_1.csv")
theta_BM_Status=read.csv("theta_BM_Status_MLBRL_R4_I4_D4_noD_1.csv")
theta_BM_Type=read.csv("theta_BM_Type_MLBRL_R4_I4_D4_noD_1.csv")
theta_BM_SES=read.csv("theta_BM_SES_MLBRL_R4_I4_D4_noD_1.csv")
theta_BU_Income=read.csv("theta_BU_Income_MLBRL_R4_I4_D4_noD_1.csv")
theta_BU_Region=read.csv("theta_BU_Region_MLBRL_R4_I4_D4_noD_1.csv")
theta_BU_Status=read.csv("theta_BU_Status_MLBRL_R4_I4_D4_noD_1.csv")
theta_BU_Type=read.csv("theta_BU_Type_MLBRL_R4_I4_D4_noD_1.csv")
theta_BU_SES=read.csv("theta_BU_SES_MLBRL_R4_I4_D4_noD_1.csv")
theta_CM_DOB=read.csv("theta_CM_DOB_MLBRL_R4_I4_D4_noD_1.csv")
theta_CM_Gender=read.csv("theta_CM_Gender_MLBRL_R4_I4_D4_noD_1.csv")
theta_CU_DOB=read.csv("theta_CU_DOB_MLBRL_R4_I4_D4_noD_1.csv")
theta_CU_Gender=read.csv("theta_CU_Gender_MLBRL_R4_I4_D4_noD_1.csv")
theta_CNB_DOB=read.csv("theta_CBU_DOB_MLBRL_R4_I4_D4_noD_1.csv")
theta_CNB_Gender=read.csv("theta_CBU_Gender_MLBRL_R4_I4_D4_noD_1.csv")

################The estimate of each parameter will be the posterior mean of the last 100 iterations
theta_BM_Income_R=apply(theta_BM_Income[,c(400:500)],1,mean)
theta_BM_Region_R=apply(theta_BM_Region[,c(400:500)],1,mean)
theta_BM_Status_R=apply(theta_BM_Status[,c(400:500)],1,mean)
theta_BM_Type_R=apply(theta_BM_Type[,c(400:500)],1,mean)
theta_BM_SES_R=apply(theta_BM_SES[,c(400:500)],1,mean)

theta_BU_Income_R=apply(theta_BU_Income[,c(400:500)],1,mean)
theta_BU_Region_R=apply(theta_BU_Region[,c(400:500)],1,mean)
theta_BU_Status_R=apply(theta_BU_Status[,c(400:500)],1,mean)
theta_BU_Type_R=apply(theta_BU_Type[,c(400:500)],1,mean)
theta_BU_SES_R=apply(theta_BU_SES[,c(400:500)],1,mean)

theta_CM_DOB_R=apply(theta_CM_DOB[,c(400:500)],1,mean)
theta_CM_Gender_R=apply(theta_CM_Gender[,c(400:500)],1,mean)

theta_CU_DOB_R=apply(theta_CU_DOB[,c(400:500)],1,mean)
theta_CU_Gender_R=apply(theta_CU_Gender[,c(400:500)],1,mean)

theta_CNB_DOB_R=apply(theta_CNB_DOB[,c(400:500)],1,mean)
theta_CNB_Gender_R=apply(theta_CNB_Gender[,c(400:500)],1,mean)

###############Calculate the Ratio of P(Gamma_Bst|theta_BM)/P(Gamma_Bst|theta_BU) for the B matrix
num=(theta_BM_Income_R[1]^Zeta1)*(theta_BM_Income_R[2]^Zeta2)*(theta_BM_Status_R[1]^Zeta3)*(theta_BM_Status_R[2]^Zeta4)*
  (theta_BM_Type_R[1]^Zeta5)*(theta_BM_Type_R[2]^Zeta6)*(theta_BM_SES_R[1]^Zeta7)*(theta_BM_SES_R[2]^Zeta8)*(theta_BM_Income_R[1]^Zeta9)*(theta_BM_Income_R[2]^Zeta10)

den=(theta_BU_Income_R[1]^Zeta1)*(theta_BU_Income_R[2]^Zeta2)*(theta_BU_Status_R[1]^Zeta3)*(theta_BU_Status_R[2]^Zeta4)*
  (theta_BU_Type_R[1]^Zeta5)*(theta_BU_Type_R[2]^Zeta6)*(theta_BU_SES_R[1]^Zeta7)*(theta_BU_SES_R[2]^Zeta8)*(theta_BU_Income_R[1]^Zeta9)*(theta_BU_Income_R[2]^Zeta10)

###############Calculate the Ratio of P(Gamma_Cst|theta_CM)/P(Gamma_Cst|theta_CU) for the C matrix
Gamma1_flat=as.vector(Gamma1)
Gamma2_flat=as.vector(Gamma2)
Gamma3_flat=as.vector(Gamma3)
Gamma5_flat=as.vector(Gamma5)
Gamma6_flat=as.vector(Gamma6)

Zeta1_flat=as.vector(Zeta1[c(rep(1:nrow(Zeta1), each=15), rep(1:nrow(Zeta1), each=5)), c(rep(1:ncol(Zeta1), each=15), rep(1:ncol(Zeta1), each=15))])
Zeta2_flat=as.vector(Zeta2[c(rep(1:nrow(Zeta2), each=15), rep(1:nrow(Zeta2), each=5)), c(rep(1:ncol(Zeta2), each=15), rep(1:ncol(Zeta2), each=15))])
Zeta3_flat=as.vector(Zeta3[c(rep(1:nrow(Zeta3), each=15), rep(1:nrow(Zeta3), each=5)), c(rep(1:ncol(Zeta3), each=15), rep(1:ncol(Zeta3), each=15))])
Zeta4_flat=as.vector(Zeta4[c(rep(1:nrow(Zeta4), each=15), rep(1:nrow(Zeta4), each=5)), c(rep(1:ncol(Zeta4), each=15), rep(1:ncol(Zeta4), each=15))])
Zeta5_flat=as.vector(Zeta5[c(rep(1:nrow(Zeta5), each=15), rep(1:nrow(Zeta5), each=5)), c(rep(1:ncol(Zeta5), each=15), rep(1:ncol(Zeta5), each=15))])
Zeta6_flat=as.vector(Zeta6[c(rep(1:nrow(Zeta6), each=15), rep(1:nrow(Zeta6), each=5)), c(rep(1:ncol(Zeta6), each=15), rep(1:ncol(Zeta6), each=15))])
Zeta7_flat=as.vector(Zeta7[c(rep(1:nrow(Zeta7), each=15), rep(1:nrow(Zeta7), each=5)), c(rep(1:ncol(Zeta7), each=15), rep(1:ncol(Zeta7), each=15))])
Zeta8_flat=as.vector(Zeta8[c(rep(1:nrow(Zeta8), each=15), rep(1:nrow(Zeta8), each=5)), c(rep(1:ncol(Zeta8), each=15), rep(1:ncol(Zeta8), each=15))])
Zeta9_flat=as.vector(Zeta9[c(rep(1:nrow(Zeta9), each=15), rep(1:nrow(Zeta9), each=5)), c(rep(1:ncol(Zeta9), each=15), rep(1:ncol(Zeta9), each=15))])
Zeta10_flat=as.vector(Zeta10[c(rep(1:nrow(Zeta10), each=15), rep(1:nrow(Zeta10), each=5)), c(rep(1:ncol(Zeta10), each=15), rep(1:ncol(Zeta10), each=15))])

#Obtain Unique Record-Level IDs and Append them to each record pair
Full_MatchID_A=N_A$MatchID#600
Full_MatchID_B=N_B$MatchID#1200

A_LinkID=rep(Full_MatchID_A,1200)
B_LinkID=rep(Full_MatchID_B,each=600)

#Obtain Unique Block-Level IDs and Append them to each record pair
Full_BlockID_A=rep(c(rep(Block_A_ID,each=15),rep(Block_A_ID,each=5)),30*40)
Full_BlockID_B=c(rep(Block_B_ID,each=30*20*15),rep(Block_B_ID,each=30*20*15))

GammaZetaFull=data.frame(A_LinkID,B_LinkID,Full_BlockID_A,Full_BlockID_B,Gamma1_flat,
                         Gamma2_flat,Gamma3_flat,Gamma5_flat,Gamma6_flat,
                         Zeta1_flat,Zeta2_flat,Zeta3_flat,Zeta4_flat,Zeta5_flat,
                         Zeta6_flat,Zeta7_flat,Zeta8_flat,Zeta9_flat,Zeta10_flat)

GammaZetaIDs=GammaZetaFull[,1:4]
GammaZeta=as.matrix(GammaZetaFull[,5:length(GammaZetaFull)])
GammaZeta_Link=GammaZeta[,c(1:5)]

#Extract paramters from posterior estimates
theta_M=c(theta_CM_DOB_R,theta_CM_Gender_R)
theta_U=c(theta_CU_DOB_R,theta_CU_Gender_R)
theta_NB=c(theta_CNB_DOB_R,theta_CNB_Gender_R)

mprod=apply(theta_M^t(GammaZeta_Link),2,prod)
uprod=apply(theta_U^t(GammaZeta_Link),2,prod)

w=mprod/uprod

IDweights=data.frame(GammaZetaFull,w)
IDweights$C=0

theta_C_likelihood=matrix(0,nrow=nrow(B),ncol=ncol(B))

#Apply Linear Sum Assignment algorithm to estimate a linking configuration within each block pair
for(r in 1:length(Block_B_ID)){
  for(q in 1:length(Block_A_ID)){
    block_qr=IDweights[(IDweights$Full_BlockID_B==Block_B_ID[r]) & (IDweights$Full_BlockID_A==Block_A_ID[q]),]
    w_qr=matrix(block_qr$w, nrow=length(unique(block_qr$A_LinkID)), ncol=length(unique(block_qr$B_LinkID)))
    C_qr=matrix(0, nrow=length(unique(block_qr$A_LinkID)), ncol=length(unique(block_qr$B_LinkID)))
    
    Locations_col=solve_LSAP(as.matrix(w_qr), maximum=TRUE)[]
    Locations_row=c(1:nrow(w_qr))
    
    #Identify links only if the linking weight is above the 75th percentile of all weights
    C_qr[cbind(Locations_row,Locations_col)]=1
    C_qr[w_qr<=quantile(w,0.75)]=0
    
    if(sum(apply(C_qr,2,sum)>1)>0 | sum(apply(C_qr,1,sum)>1)>0){print('LSAP error'); break}
    
    IDweights$C[(IDweights$Full_BlockID_B==Block_B_ID[r]) & (IDweights$Full_BlockID_A==Block_A_ID[q])]=C_qr
    
    Gamma_qr=IDweights[(IDweights$Full_BlockID_B==Block_B_ID[r]) & (IDweights$Full_BlockID_A==Block_A_ID[q]),c(5:9)]
    likelihood_qr=apply((theta_M^t(Gamma_qr*C_qr))*(theta_U^t(Gamma_qr*(1-C_qr))),2,prod)
    likelihood_nb_qr=apply(theta_NB^t(Gamma_qr),2,prod)
    theta_C_likelihood[q,r]=prod(likelihood_qr/likelihood_nb_qr)
  }
}

#Generate Heatmap for Block Level probabilities for first 10 blocks using MLBRL
Block_Likelihood_MLBRL=t(apply((num/den)*theta_C_likelihood,1,rev))
row_sum=apply(Block_Likelihood_MLBRL,1,sum)
Sum_mat=matrix(row_sum,nrow=30,ncol=40)
Block_LogLikelihood_MLBRL=log(Block_Likelihood_MLBRL/Sum_mat)
Block_LogLikelihood_MLBRL[Block_LogLikelihood_MLBRL<(-30)]=-30

Block_LogLikelihood_MLBRL=Block_LogLikelihood_MLBRL[c(1:10),c(31:40)]
rownames(Block_LogLikelihood_MLBRL)=as.character(c(1:10))
colnames(Block_LogLikelihood_MLBRL)=as.character(sort(seq(1:10),decreasing=TRUE))
melted_Block_Likelihood_MLBRL=melt(Block_LogLikelihood_MLBRL)

png("MLBRL_block_small.png", width = 650, height = 550)
ggplot(data=melted_Block_Likelihood_MLBRL, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+coord_cartesian(xlim=c(1,10),ylim=c(1,10))+
  scale_x_continuous(breaks=seq(1:10))+scale_y_continuous(breaks=seq(1:10))+theme(legend.title=element_blank())+labs(x="File 1 Block", y="File 2 Block")+ 
  scale_fill_gradientn(limits = c(-30,0), colors = rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(16)[6:16]))
dev.off()

###############Calculate the Ratio of P(Gamma_Cst|theta_CM)/P(Gamma_Cst|theta_CU) for the C matrix
num_C=(theta_CM_DOB_R[1]^Gamma1)*(theta_CM_DOB_R[2]^Gamma2)*(theta_CM_DOB_R[3]^Gamma3)*(theta_CM_Gender_R[1]^Gamma5)*(theta_CM_Gender_R^Gamma6)
den_C=(theta_CU_DOB_R[1]^Gamma1)*(theta_CU_DOB_R[2]^Gamma2)*(theta_CU_DOB_R[3]^Gamma3)*(theta_CU_Gender_R[1]^Gamma5)*(theta_CU_Gender_R^Gamma6)

#Generate Heatmap for Record Level probabilities for first 10 blocks using MLBRL
Block_ratio=(num/den)*theta_C_likelihood
Block_prob=Block_ratio[c(rep(1:nrow(Block_ratio), each=15), rep(1:nrow(Block_ratio), each=5)), c(rep(1:ncol(Block_ratio), each=15), rep(1:ncol(Block_ratio), each=15))]
Link_Likelihhood=t(apply((num_C/den_C*Block_prob),1,rev))
row_sum=apply(Link_Likelihhood,1,sum)
Sum_mat=matrix(row_sum,nrow=600,ncol=1200)
Link_LogLikelihood_MLBRL=log(Link_Likelihhood/Sum_mat)
Link_LogLikelihood_MLBRL=Link_LogLikelihood_MLBRL[c(1:150),c(1051:1200)]
rownames(Link_LogLikelihood_MLBRL)=as.character(c(1:150))
colnames(Link_LogLikelihood_MLBRL)=as.character(sort(seq(1:150),decreasing=TRUE))

Link_LogLikelihood_MLBRL[Link_LogLikelihood_MLBRL< (-40)]=-40
melted_Link_Likelihood_MLBRL=melt(Link_LogLikelihood_MLBRL)

png("MLBRL_link_small.png", width = 650, height = 550)
ggplot(data=melted_Link_Likelihood_MLBRL, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+coord_cartesian(xlim=c(7,144),ylim=c(7,144))+
  scale_x_continuous(breaks=c(1,seq(15,150,length=10)))+scale_y_continuous(breaks=c(1,seq(15,150,length=10)))+xlab("File 1 Record")+ylab("File 2 Record")+ theme(legend.title=element_blank())+
  scale_fill_gradientn(limits = c(-40,0),
                       colours=c("navyblue", "darkmagenta", "darkorange1","red"))
dev.off()

##############################################################################################################
###############Record Pair Log Probability Density Using BRL
###############Read in posterior estimates of Theta
theta_BM_Income=read.csv("theta_BM_Income_BRL_R4_I4_D4_noD_1.csv")
theta_BM_Region=read.csv("theta_BM_Region_BRL_R4_I4_D4_noD_1.csv")
theta_BM_Status=read.csv("theta_BM_Status_BRL_R4_I4_D4_noD_1.csv")
theta_BM_Type=read.csv("theta_BM_Type_BRL_R4_I4_D4_noD_1.csv")
theta_BM_SES=read.csv("theta_BM_SES_BRL_R4_I4_D4_noD_1.csv")
theta_BU_Income=read.csv("theta_BU_Income_BRL_R4_I4_D4_noD_1.csv")
theta_BU_Region=read.csv("theta_BU_Region_BRL_R4_I4_D4_noD_1.csv")
theta_BU_Status=read.csv("theta_BU_Status_BRL_R4_I4_D4_noD_1.csv")
theta_BU_Type=read.csv("theta_BU_Type_BRL_R4_I4_D4_noD_1.csv")
theta_BU_SES=read.csv("theta_BU_SES_BRL_R4_I4_D4_noD_1.csv")
theta_CM_DOB=read.csv("theta_CM_DOB_BRL_R4_I4_D4_noD_1.csv")
theta_CM_Gender=read.csv("theta_CM_Gender_BRL_R4_I4_D4_noD_1.csv")
theta_CU_DOB=read.csv("theta_CU_DOB_BRL_R4_I4_D4_noD_1.csv")
theta_CU_Gender=read.csv("theta_CU_Gender_BRL_R4_I4_D4_noD_1.csv")

################The estimate of each parameter will be the posterior mean of the last 100 iterations
theta_BM_Income_R=apply(theta_BM_Income[,c(400:500)],1,mean)
theta_BM_Region_R=apply(theta_BM_Region[,c(400:500)],1,mean)
theta_BM_Status_R=apply(theta_BM_Status[,c(400:500)],1,mean)
theta_BM_Type_R=apply(theta_BM_Type[,c(400:500)],1,mean)
theta_BM_SES_R=apply(theta_BM_SES[,c(400:500)],1,mean)

theta_BU_Income_R=apply(theta_BU_Income[,c(400:500)],1,mean)
theta_BU_Region_R=apply(theta_BU_Region[,c(400:500)],1,mean)
theta_BU_Status_R=apply(theta_BU_Status[,c(400:500)],1,mean)
theta_BU_Type_R=apply(theta_BU_Type[,c(400:500)],1,mean)
theta_BU_SES_R=apply(theta_BU_SES[,c(400:500)],1,mean)

theta_CM_DOB_R=apply(theta_CM_DOB[,c(400:500)],1,mean)
theta_CM_Gender_R=apply(theta_CM_Gender[,c(400:500)],1,mean)

theta_CU_DOB_R=apply(theta_CU_DOB[,c(400:500)],1,mean)
theta_CU_Gender_R=apply(theta_CU_Gender[,c(400:500)],1,mean)

#Incorporate the blocking variables as linking variables within RL algorithm
N_A$Region=c(rep(Region_A,each=15),rep(Region_A,each=5))
N_B$Region=c(rep(Region_B,each=15),rep(Region_B,each=15))
N_A$Status=c(rep(Status_A,each=15),rep(Status_A,each=5))
N_B$Status=c(rep(Status_B,each=15),rep(Status_B,each=15))
N_A$Type=c(rep(Type_A,each=15),rep(Type_A,each=5))
N_B$Type=c(rep(Type_B,each=15),rep(Type_B,each=15))
N_A$SES=c(rep(SES_A,each=15),rep(SES_A,each=5))
N_B$SES=c(rep(SES_B,each=15),rep(SES_B,each=15))
N_A$Income=c(rep(Income_A,each=15),rep(Income_A,each=5))
N_B$Income=c(rep(Income_B,each=15),rep(Income_B,each=15))

#Create Comparison matrix Gamma for each of the fields in A and B
Region_comparison=t(matrix(N_B$Region,nrow=length(N_B$Region),ncol=nrow(N_A)))
Status_comparison=t(matrix(N_B$Status,nrow=length(N_B$Status),ncol=nrow(N_A)))
Type_comparison=t(matrix(N_B$Type,nrow=length(N_B$Type),ncol=nrow(N_A)))
SES_comparison=t(matrix(N_B$SES,nrow=length(N_B$SES),ncol=nrow(N_A)))
Income_comparison=t(matrix(N_B$Income,nrow=length(N_B$Income),ncol=nrow(N_A)))

#Compare each variable in dataset A with each element in dataset B element-wise
Region_gamma=apply(Region_comparison,2,'==',N_A$Region)
Status_gamma=apply(Status_comparison,2,'==',N_A$Status)
Type_gamma=apply(Type_comparison,2,'==',N_A$Type)
SES_gamma=apply(SES_comparison,2,'==',N_A$SES)
Income_gamma=apply(Income_comparison,2,function(x) abs(x-Income_A)<=500)

#Region
Zeta1=-(Region_gamma)+1
Zeta2=Region_gamma
#Status
Zeta3=-(Status_gamma)+1
Zeta4=Status_gamma
#Type
Zeta5=-(Type_gamma)+1
Zeta6=Type_gamma
#SES
Zeta7=-(SES_gamma)+1
Zeta8=SES_gamma
#Income
Zeta9=-(Income_gamma)+1
Zeta10=Income_gamma


###############Calculate the Ratio of P(Gamma_Cst|theta_CM)/P(Gamma_Cst|theta_CU) for the C matrix
num=(theta_BM_Income_R[1]^Zeta1)*(theta_BM_Income_R[2]^Zeta2)*(theta_BM_Status_R[1]^Zeta3)*(theta_BM_Status_R[2]^Zeta4)*
  (theta_BM_Type_R[1]^Zeta5)*(theta_BM_Type_R[2]^Zeta6)*(theta_BM_SES_R[1]^Zeta7)*(theta_BM_SES_R[2]^Zeta8)*(theta_BM_Income_R[1]^Zeta9)*(theta_BM_Income_R[2]^Zeta10)*
  (theta_CM_DOB_R[1]^Gamma1)*(theta_CM_DOB_R[2]^Gamma2)*(theta_CM_DOB_R[3]^Gamma3)*(theta_CM_Gender_R[1]^Gamma5)*(theta_CM_Gender_R^Gamma6)

den=(theta_BU_Income_R[1]^Zeta1)*(theta_BU_Income_R[2]^Zeta2)*(theta_BU_Status_R[1]^Zeta3)*(theta_BU_Status_R[2]^Zeta4)*
  (theta_BU_Type_R[1]^Zeta5)*(theta_BU_Type_R[2]^Zeta6)*(theta_BU_SES_R[1]^Zeta7)*(theta_BU_SES_R[2]^Zeta8)*(theta_BU_Income_R[1]^Zeta9)*(theta_BU_Income_R[2]^Zeta10)*
  (theta_CU_DOB_R[1]^Gamma1)*(theta_CU_DOB_R[2]^Gamma2)*(theta_CU_DOB_R[3]^Gamma3)*(theta_CU_Gender_R[1]^Gamma5)*(theta_CU_Gender_R^Gamma6)

#Generate Heatmap for Record Level probabilities for first 10 blocks using BRL
Link_Likelihhood=t(apply((num/den),1,rev))
row_sum=apply(Link_Likelihhood,1,sum)
Sum_mat=matrix(row_sum,nrow=600,ncol=1200)
Link_LogLikelihood_BRL=log(Link_Likelihhood/Sum_mat)
Link_LogLikelihood_BRL=Link_LogLikelihood_BRL[c(1:150),c(1051:1200)]
rownames(Link_LogLikelihood_BRL)=as.character(c(1:150))
colnames(Link_LogLikelihood_BRL)=as.character(sort(seq(1:150),decreasing=TRUE))
melted_Link_Likelihood_BRL=melt(Link_LogLikelihood_BRL)

png("BRL_link_small.png", width = 650, height = 550)
ggplot(data=melted_Link_Likelihood_BRL, aes(x=Var1, y=Var2, fill=value))+ geom_tile()+coord_cartesian(xlim=c(7,144),ylim=c(7,144))+
  scale_x_continuous(breaks=c(1,seq(15,150,length=10)))+scale_y_continuous(breaks=c(1,seq(15,150,length=10)))+xlab("File 1 Record")+ylab("File 2 Record")+ theme(legend.title=element_blank())+
  scale_fill_gradientn(limits = c(-40,0),
                       colours=c("navyblue", "darkmagenta", "darkorange1","red"))
dev.off()
