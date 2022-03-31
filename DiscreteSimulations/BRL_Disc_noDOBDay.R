rm(list = ls())
library(coda)
library(lubridate)
library(MCMCpack)
library(crayon)
library(tidyr)
library(reshape2)
library(clue)
library(fastLink)
library(Rcpp)
library(Matrix)
library(stringdist)

##########################################################################
##Reading Arguments for the simulations#######################################################
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  Regionerrorprob=0.4
  IncomeError = 0.4
  DOBerrorprob=0.4
} else if (length(args)==3) {
  Regionerrorprob=as.numeric(args[1])/10
  IncomeError = as.numeric(args[2])/10
  DOBerrorprob =as.numeric(args[3])/10
}
print(sprintf("%f,%f,%f",Regionerrorprob, IncomeError, DOBerrorprob))

########################################################################################################################
########Implements BRL for the Simulation Setting ######################################################################
########################################################################################################################

set.seed(1234)
seed=round(runif(100,1,100000000))

###########################################################################################################################
#loading past iterations
###########################################################################################################################
nStartVal = 1
if(file.exists(sprintf("Sensitivity_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10)))
{
  Sensitivity=c(read.csv( sprintf("Sensitivity_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Sensitivity_sd=c(read.csv(sprintf("Sensitivity_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  PPV=c(read.csv(sprintf("PPV_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  PPV_sd=c(read.csv(sprintf("PPV_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  F1=c(read.csv(sprintf("F1_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  F1_sd=c(read.csv(sprintf("F1_sd_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Accuracy=c(read.csv(sprintf("Accuracy_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Accuracy_sd=c(read.csv(sprintf("Accuracy_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  nStartVal=which.min(Sensitivity[,2])
} else {
  #Overall simulation statistics to Export
  Sensitivity=rep(0,100); Sensitivity_sd=rep(0,100)
  PPV=rep(0,100); PPV_sd=rep(0,100)
  F1=rep(0,100); F1_sd=rep(0,100)
  Accuracy=rep(0,100); Accuracy_sd=rep(0,100)
}

#Overall simulation statistics to Export
Sensitivity=rep(0,100); Sensitivity_sd=rep(0,100)
PPV=rep(0,100); PPV_sd=rep(0,100)
F1=rep(0,100); F1_sd=rep(0,100)
Accuracy=rep(0,100); Accuracy_sd=rep(0,100)

#TrueLinkDesignation=c(1:450)
TrueLinkDesignation=c(1:450,rep(9999,150))
TrueBlockDesignation=c(1:30)

#Repeat the simulation 100 times
for(nrep in 1:100){
  set.seed(seed[nrep])
############################################################
####Basic Simulation Scenario
############################################################
###There will be two datasets: Data A and Data B.
###Blocking Structure: Data A will have 30 blocks, Data B will have 40 blocks. Each block in Data A will exist in Data B
###Linking Structure: Each block in Data A will have 20 records, each block in Data B will have 30 records
###Each block will have 15 true links 
###N_A = 600; N_B = 1200; S = 30; T = 40; N_Mst = 15, N_M = 450

##Linking Variables:
##Age ~ N(40, 2^2) then convert age to DOB
##Gender ~ Bernoulli(0.5)

##Blocking variables:
##Region ~ Multinomial 3 categories
##Status ~ Bernoulli(.8)
##Type ~ Bernoulli(.5)

###The variables MatchID_A and MatchID_B give the true unique identifier for records in file A and file B.
###Similarly, the variables BlockID_A and BlockID_B give the true block identifer in file A and B
############################################################
##########Linking Variables
############################################################
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

Regionerror=rbinom(length(Region_A),1,Regionerrorprob)

for(i in 1:length(Region_A)){
  if(Regionerror[i]==1){
    Region_A[i]=sample(levels(as.factor(Region_A))[levels(as.factor(Region_A))!=Region_A[i]],1)
  }
}

#Error in Income
Income_A=Income_A+rnorm(length(Income_A),0,sd=500/qnorm(1-IncomeError/2))

#Error in DOB
DOBerror=rbinom(length(N_A$DOB_month),1,DOBerrorprob)
for(i in 1:length(N_A$DOB_month)){
  if(DOBerror[i]==1){
    N_A$DOB_month[i]=sample(unique(N_A$DOB_month)[unique(N_A$DOB_month)!=N_A$DOB_month[i]],1)
  }
}

############################################################
##########Linking Variable Comparisons
############################################################
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

#Create Comparison matrix Gamma for each of the linking variables in A and B
DOB_year_comparison=t(matrix(N_B$DOB_year,nrow=length(N_B$DOB_year),ncol=nrow(N_A)))
DOB_month_comparison=t(matrix(N_B$DOB_month,nrow=length(N_B$DOB_month),ncol=nrow(N_A)))
DOB_day_comparison=t(matrix(N_B$DOB_day,nrow=length(N_B$DOB_day),ncol=nrow(N_A)))
Gender_comparison=t(matrix(N_B$Gender,nrow=length(N_B$Gender),ncol=nrow(N_A)))
Region_comparison=t(matrix(N_B$Region,nrow=length(N_B$Region),ncol=nrow(N_A)))
Status_comparison=t(matrix(N_B$Status,nrow=length(N_B$Status),ncol=nrow(N_A)))
Type_comparison=t(matrix(N_B$Type,nrow=length(N_B$Type),ncol=nrow(N_A)))
SES_comparison=t(matrix(N_B$SES,nrow=length(N_B$SES),ncol=nrow(N_A)))
Income_comparison=t(matrix(N_B$Income,nrow=length(N_B$Income),ncol=nrow(N_A)))

#Compare each variable in dataset A with each element in dataset B element-wise
DOB_year_gamma=apply(DOB_year_comparison,2,'==',N_A$DOB_year)
DOB_month_gamma=apply(DOB_month_comparison,2,'==',N_A$DOB_month)
#DOB_day_gamma=apply(DOB_day_comparison,2,'==',N_A$DOB_day)
Gender_gamma=apply(Gender_comparison,2,'==',N_A$Gender)
Region_gamma=apply(Region_comparison,2,'==',N_A$Region)
Status_gamma=apply(Status_comparison,2,'==',N_A$Status)
Type_gamma=apply(Type_comparison,2,'==',N_A$Type)
SES_gamma=apply(SES_comparison,2,'==',N_A$SES)
Income_gamma=apply(Income_comparison,2,function(x) abs(x-Income_A)<=500)

############################################################
###Create Comparison matrices from element-wise comparisons
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

#Specify values for hyperparameters of prior distributions
prior_CM_DOB=c(1,1,1); prior_CM_Gender=c(1,1); prior_CM_Region=c(1,1); prior_CM_Status=c(1,1); prior_CM_Type=c(1,1)
prior_CM_SES=c(1,1); prior_CM_Income=c(1,1)
prior_CU_DOB=c(1,1,1); prior_CU_Gender=c(1,1); prior_CU_Region=c(1,1); prior_CU_Status=c(1,1); prior_CU_Type=c(1,1)
prior_CU_SES=c(1,1); prior_CU_Income=c(1,1)
prior_pi=c(1,1)

#Specify number of iterations (including burn-in) and initialize parameter matrices
#nsim=500 lower number of replications, easier for testing
nsim=2000
theta_CM_DOB=matrix(0,nrow=length(prior_CM_DOB),ncol=nsim)
theta_CM_Gender=matrix(0,nrow=length(prior_CM_Gender),ncol=nsim)
theta_CU_DOB=matrix(0,nrow=length(prior_CU_DOB),ncol=nsim)
theta_CU_Gender=matrix(0,nrow=length(prior_CU_Gender),ncol=nsim)
theta_CM_Region=matrix(0,nrow=length(prior_CM_Region),ncol=nsim)
theta_CU_Region=matrix(0,nrow=length(prior_CU_Region),ncol=nsim)
theta_CM_Status=matrix(0,nrow=length(prior_CM_Status),ncol=nsim)
theta_CU_Status=matrix(0,nrow=length(prior_CU_Status),ncol=nsim)
theta_CM_Type=matrix(0,nrow=length(prior_CM_Type),ncol=nsim)
theta_CU_Type=matrix(0,nrow=length(prior_CU_Type),ncol=nsim)
theta_CM_SES=matrix(0,nrow=length(prior_CM_SES),ncol=nsim)
theta_CU_SES=matrix(0,nrow=length(prior_CU_SES),ncol=nsim)
theta_CM_Income=matrix(0,nrow=length(prior_CM_Income),ncol=nsim)
theta_CU_Income=matrix(0,nrow=length(prior_CU_Income),ncol=nsim)

#Initialize Vectors to store links and other linking statistics
LinkDesignation=matrix(0,nrow=nrow(N_A),ncol=nsim)
LinkProbability=matrix(0,nrow=nrow(N_A),ncol=nsim)

#Start with an empty linking matrix
C=matrix(0,nrow=nrow(N_A),ncol=nrow(N_B))

#Specify starting values for linking and blocking parameters 
theta_CM_DOB[,1]=c(.05,25,.7)
theta_CU_DOB[,1]=c(.7,.20,.10)
theta_CM_Gender[,1]=c(.25,.75)
theta_CU_Gender[,1]=c(.75,.25)
theta_CM_Region[,1]=c(.3,.7)
theta_CU_Region[,1]=c(.66,.34)
theta_CM_Status[,1]=c(.25,.75)
theta_CU_Status[,1]=c(.75,.25)
theta_CM_Type[,1]=c(.4,.6)
theta_CU_Type[,1]=c(.6,.4)
theta_CM_SES[,1]=c(.35,.65)
theta_CU_SES[,1]=c(.65,.35)
theta_CM_Income[,1]=c(.2,.8)
theta_CU_Income[,1]=c(.8,.2)

for(t in nStartVal:nsim){
  ########################################################################################################################
  #Sample the Posterior Distribution of linking and blocking parameters
  theta_CM_DOB[,t]=rdirichlet(1,c(prior_CM_DOB[1]+sum(Gamma1*C),prior_CM_DOB[2]+sum(Gamma2*C),prior_CM_DOB[3]+sum(Gamma3*C)))
  theta_CU_DOB[,t]=rdirichlet(1,c(prior_CU_DOB[1]+sum(Gamma1*(1-C)),prior_CU_DOB[2]+sum(Gamma2*(1-C)),prior_CU_DOB[3]+sum(Gamma3*(1-C))))
  theta_CM_Gender[,t]=rdirichlet(1,c(prior_CM_Gender[1]+sum(Gamma5*C),prior_CM_Gender[2]+sum(Gamma6*C)))
  theta_CU_Gender[,t]=rdirichlet(1,c(prior_CU_Gender[1]+sum(Gamma5*(1-C)),prior_CU_Gender[2]+sum(Gamma6*(1-C))))
  theta_CM_Region[,t]=rdirichlet(1,c(prior_CM_Region[1]+sum(Zeta1*C),prior_CM_Region[2]+sum(Zeta2*C)))
  theta_CU_Region[,t]=rdirichlet(1,c(prior_CU_Region[1]+sum(Zeta1*(1-C)),prior_CU_Region[2]+sum(Zeta2*(1-C))))
  theta_CM_Status[,t]=rdirichlet(1,c(prior_CM_Status[1]+sum(Zeta3*C),prior_CM_Status[2]+sum(Zeta4*C)))
  theta_CU_Status[,t]=rdirichlet(1,c(prior_CU_Status[1]+sum(Zeta3*(1-C)),prior_CU_Status[2]+sum(Zeta4*(1-C))))
  theta_CM_Type[,t]=rdirichlet(1,c(prior_CM_Type[1]+sum(Zeta5*C),prior_CM_Type[2]+sum(Zeta6*C)))
  theta_CU_Type[,t]=rdirichlet(1,c(prior_CU_Type[1]+sum(Zeta5*(1-C)),prior_CU_Type[2]+sum(Zeta6*(1-C))))
  theta_CM_SES[,t]=rdirichlet(1,c(prior_CM_SES[1]+sum(Zeta7*C),prior_CM_SES[2]+sum(Zeta8*C)))
  theta_CU_SES[,t]=rdirichlet(1,c(prior_CU_SES[1]+sum(Zeta7*(1-C)),prior_CU_SES[2]+sum(Zeta8*(1-C))))
  theta_CM_Income[,t]=rdirichlet(1,c(prior_CM_Income[1]+sum(Zeta9*C),prior_CM_Income[2]+sum(Zeta10*C)))
  theta_CU_Income[,t]=rdirichlet(1,c(prior_CU_Income[1]+sum(Zeta9*(1-C)),prior_CU_Income[2]+sum(Zeta10*(1-C))))
  
  ########################################################################################################################
  #Iterate through the rows of C
  for(i in 1:nrow(C)){
    #Resetting the result for row i
    C[i,]=0
    #Extracting individuals in file B that have links, not including the individual currently linked to i
    B_linked=c(1:nrow(N_B))[apply(C,2,sum)==1]
    #If there are no links, create empty link vector. Then extract non-linked individuals in dataset B
    if(length(B_linked)==0){
      B_unlinked=c(1:nrow(N_B))
    }else{B_unlinked=setdiff(c(1:nrow(N_B)),B_linked)}
    
    #Extract matrix of DOB and Gender Gamma
    Gamma_DOB=cbind(Gamma1[i,B_unlinked],Gamma2[i,B_unlinked],Gamma3[i,B_unlinked])
    Gamma_Gender=cbind(Gamma5[i,B_unlinked],Gamma6[i,B_unlinked])
    Zeta_Region=cbind(Zeta1[i,B_unlinked],Zeta2[i,B_unlinked])
    Zeta_Status=cbind(Zeta3[i,B_unlinked],Zeta4[i,B_unlinked])
    Zeta_Type=cbind(Zeta5[i,B_unlinked],Zeta6[i,B_unlinked])
    Zeta_SES=cbind(Zeta7[i,B_unlinked],Zeta8[i,B_unlinked])
    Zeta_Income=cbind(Zeta9[i,B_unlinked],Zeta10[i,B_unlinked])
    
    #Calculate ratio of likelihoods for Gamma comparisons
    num=apply(theta_CM_DOB[,t]^t(Gamma_DOB),2,prod)*apply(theta_CM_Gender[,t]^t(Gamma_Gender),2,prod)*
      apply(theta_CM_Region[,t]^t(Zeta_Region),2,prod)*apply(theta_CM_Status[,t]^t(Zeta_Status),2,prod)*
      apply(theta_CM_Type[,t]^t(Zeta_Type),2,prod)*apply(theta_CM_SES[,t]^t(Zeta_SES),2,prod)*
      apply(theta_CM_Income[,t]^t(Zeta_Income),2,prod)
    den=apply(theta_CU_DOB[,t]^t(Gamma_DOB),2,prod)*apply(theta_CU_Gender[,t]^t(Gamma_Gender),2,prod)*
      apply(theta_CU_Region[,t]^t(Zeta_Region),2,prod)*apply(theta_CU_Status[,t]^t(Zeta_Status),2,prod)*
      apply(theta_CU_Type[,t]^t(Zeta_Type),2,prod)*apply(theta_CU_SES[,t]^t(Zeta_SES),2,prod)*
      apply(theta_CU_Income[,t]^t(Zeta_Income),2,prod)
    Likelihood=num/den
    
    #Calculate the probability of individual i not linking
    p_nolink=(nrow(N_B)-length(B_linked))*(nrow(N_A)-length(B_linked)+prior_pi[2]-1)/(length(B_linked)+prior_pi[1])
    
    #Parse together possible moves and move probability
    B_unlinked=c(B_unlinked,nrow(N_B)+i)
    B_prob=c(Likelihood,p_nolink)/sum(Likelihood,p_nolink)
    
    #Sample new bipartite link for individual i_st
    link_designation=sample(B_unlinked,size=1,prob=B_prob)
    if(link_designation<=nrow(N_B)){C[i,link_designation]=1}
    
    #Export Posterior Probability of the given Link Designation
    LinkProbability[i,t]=B_prob[B_unlinked==link_designation]
    }

  #Export Designations of Linking Matrix after every iteration of updating the linking matrices
  a=which(C==1,arr.ind=TRUE)
  LinkDesignation[,t][apply(C,1,sum)>=1]=matrix(a[order(a[,1]),],ncol=2)[,2]
  
  # write.csv(theta_CM_Region,paste0("theta_BM_Region_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_Status,paste0("theta_BM_Status_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_Type,paste0("theta_BM_Type_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_SES,paste0("theta_BM_SES_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_Income,paste0("theta_BM_Income_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_Region,paste0("theta_BU_Region_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_Status,paste0("theta_BU_Status_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_Type,paste0("theta_BU_Type_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_SES,paste0("theta_BU_SES_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_Income,paste0("theta_BU_Income_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_DOB,paste0("theta_CM_DOB_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CM_Gender,paste0("theta_CM_Gender_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_DOB,paste0("theta_CU_DOB_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(theta_CU_Gender,paste0("theta_CU_Gender_BRL_R0_I0_D0_noD_",nrep,".csv"))
  # write.csv(LinkDesignation,paste0("LinkDesignation_BRL_R0_I0_D0_noD_",nrep,".csv"))
}

#n=apply(LinkDesignation[,401:500],2,function(x) sum(x!="0"))
#TP=apply(LinkDesignation[,401:500],2,function(x) sum(x==TrueLinkDesignation))
n=apply(LinkDesignation[,(nsim/2+1):nsim],2,function(x) sum(x!="0"))
TP=apply(LinkDesignation[,(nsim/2+1):nsim],2,function(x) sum(x==TrueLinkDesignation))
FP=n-TP

Sensitivity[nrep]=mean(TP/450)
Sensitivity_sd[nrep]=sd(TP/450)
PPV[nrep]=mean(TP/n)
PPV_sd[nrep]=sd(TP/n)
F1[nrep]=mean(2*((TP/450)*(TP/n))/(TP/450+TP/n))
F1_sd[nrep]=sd(2*((TP/450)*(TP/n))/(TP/450+TP/n))

write.csv(Sensitivity, sprintf("Sensitivity_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(Sensitivity_sd, sprintf("Sensitivity_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(PPV, sprintf("PPV_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(PPV_sd, sprintf("PPV_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(F1, sprintf("F1_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(F1_sd, sprintf("F1_sd_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(Accuracy, sprintf("Accuracy_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
write.csv(Accuracy_sd, sprintf("Accuracy_sd_BRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))

print(nrep)
Sys.sleep(0.01)
flush.console()
}