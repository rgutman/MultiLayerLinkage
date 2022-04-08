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
##Reading Arguments#######################################################
##########################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  Regionerrorprob=0
  IncomeError = 0
  DOBerrorprob=0.2
} else if (length(args)==3) {
  Regionerrorprob=as.numeric(args[1])/10
  IncomeError = as.numeric(args[2])/10
  DOBerrorprob =as.numeric(args[3])/10
}
print(sprintf("%f,%f,%f",Regionerrorprob, IncomeError, DOBerrorprob))
########################################################################################################################
########Implements MLBRL for the Simulation Setting with 40% Error in Region, 40% Error in Area Income, 40% Error in Day of Birth
########################################################################################################################

set.seed(1234)
seed=round(runif(100,1,100000000))

####loading past iterations
nStartVal = 1
if(file.exists(sprintf("Sensitivity_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10)))
{
  Sensitivity=c(read.csv( sprintf("Sensitivity_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Sensitivity_sd=c(read.csv(sprintf("Sensitivity_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  PPV=c(read.csv(sprintf("PPV_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  PPV_sd=c(read.csv(sprintf("PPV_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  F1=c(read.csv(sprintf("F1_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  F1_sd=c(read.csv(sprintf("F1_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Accuracy=c(read.csv(sprintf("Accuracy_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  Accuracy_sd=c(read.csv(sprintf("Accuracy_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10),header=TRUE)[,2])
  nStartVal=which.min(Sensitivity)
} else {
  #Overall simulation statistics to Export
  Sensitivity=rep(0,100); Sensitivity_sd=rep(0,100)
  PPV=rep(0,100); PPV_sd=rep(0,100)
  F1=rep(0,100); F1_sd=rep(0,100)
  Accuracy=rep(0,100); Accuracy_sd=rep(0,100)
}

#Overall simulation statistics to Export
print(nStartVal)
TrueLinkDesignation=c(1:450,rep(9999,150))
TrueBlockDesignation=c(1:30)

roundUp <- function(x) 10^ceiling(log10(x))

#Repeat the simulation 100 times
for(nrep in nStartVal:100){
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
  ##Income ~ Normal(30,50)
  ##SES ~multinomial 3 categories
  
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
  IncomeVar = 10
  Income_A=rnorm(30,50,IncomeVar)
  randNum = sample(100:10000,40,replace=FALSE)

  Block_A_ID=randNum[1:30]
  Block_A_S=seq_along(Block_A_ID)
  ###Blocking Variables for additional non-blocks in Data B
  Region = sample(c("N","W","MW","S"),size=10,replace=T)
  Status = rbinom(10,1,.8)
  Type = rbinom(10,1,.5)
  SES=sample(c("L","M","H"),size=10,replace=T,prob=c(.25,.5,.25))
  Income=rnorm(10,50,IncomeVar)
  
  Block_B_ID=randNum
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
  numMaxBlocks = max(length(unique((N_A$S))),length(unique(N_B$T)))
  numRound = roundUp(numMaxBlocks)
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
  #IncomeError = 0.4
  Income_A=Income_A+rnorm(length(Income_A),0,sd=2/qnorm(1-IncomeError/2))

  #Error in DOB
  #DOBerrorprob=.4
  DOBerror=rbinom(length(N_A$DOB_month),1,DOBerrorprob)
  for(i in 1:length(N_A$DOB_month)){
    if(DOBerror[i]==1){
      N_A$DOB_month[i]=sample(unique(N_A$DOB_month)[unique(N_A$DOB_month)!=N_A$DOB_month[i]],1)
    }
  }
  
  ############################################################
  ##########Linking Variable Comparisons
  ############################################################
  #Create Comparison matrix Gamma for each of the linking variables in A and B
  DOB_year_comparison=t(matrix(N_B$DOB_year,nrow=length(N_B$DOB_year),ncol=nrow(N_A)))
  DOB_month_comparison=t(matrix(N_B$DOB_month,nrow=length(N_B$DOB_month),ncol=nrow(N_A)))
  DOB_day_comparison=t(matrix(N_B$DOB_day,nrow=length(N_B$DOB_day),ncol=nrow(N_A)))
  Gender_comparison=t(matrix(N_B$Gender,nrow=length(N_B$Gender),ncol=nrow(N_A)))
  
  #Incorporate the blocking variables as linking variables in the orginal files
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
  
  
  #Compare each linking variable in dataset A with each element in dataset B element-wise
  DOB_year_gamma=apply(DOB_year_comparison,2,'==',N_A$DOB_year)
  DOB_month_gamma=apply(DOB_month_comparison,2,'==',N_A$DOB_month)
  #DOB_day_gamma=apply(DOB_day_comparison,2,'==',N_A$DOB_day)
  Gender_gamma=apply(Gender_comparison,2,'==',N_A$Gender)
  
  ############################################################
  ##########Blocking Variable Comparisons
  ############################################################
  #Create Comparison matrix Zeta for each of the blocking in A and B
  Region_comparison=t(matrix(Region_B,nrow=length(Region_B),ncol=length(Region_A)))
  Status_comparison=t(matrix(Status_B,nrow=length(Status_B),ncol=length(Status_A)))
  Type_comparison=t(matrix(Type_B,nrow=length(Type_B),ncol=length(Type_A)))
  SES_comparison=t(matrix(SES_B,nrow=length(SES_B),ncol=length(SES_A)))
  Income_comparison=t(matrix(Income_B,nrow=length(Income_B),ncol=length(Income_A)))
  
  RegionP_comparison=t(matrix(N_B$Region,nrow=length(N_B$Region),ncol=nrow(N_A)))
  StatusP_comparison=t(matrix(N_B$Status,nrow=length(N_B$Status),ncol=nrow(N_A)))
  TypeP_comparison=t(matrix(N_B$Type,nrow=length(N_B$Type),ncol=nrow(N_A)))
  SESP_comparison=t(matrix(N_B$SES,nrow=length(N_B$SES),ncol=nrow(N_A)))
  IncomeP_comparison=t(matrix(N_B$Income,nrow=length(N_B$Income),ncol=nrow(N_A)))
  
  #Compare each blocking variable in dataset A with each element in dataset B element-wise
  Region_gamma=apply(Region_comparison,2,'==',Region_A)
  Region_gamma[Region_gamma==TRUE]=1
  
  Status_gamma=apply(Status_comparison,2,'==',Status_A)
  Status_gamma[Status_gamma==TRUE]=1
  
  Type_gamma=apply(Type_comparison,2,'==',Type_A)
  Type_gamma[Type_gamma==TRUE]=1
  
  SES_gamma=apply(SES_comparison,2,'==',SES_A)
  SES_gamma[SES_gamma==TRUE]=1
  
  #Income_gamma=apply(Income_comparison,2,function(x) abs(x-Income_A)<=500)
  #Income_gamma[Income_gamma==TRUE]=1
  Income_Diff = apply(Income_comparison,2,function(x) {(x-Income_A)})
  
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
  Gamma6=as.numeric(Gender_gamma)
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
  #Zeta9=-(Income_gamma)+1
  #Zeta10=Income_gamma
  
  #Convert Comparison Matrices to Column Vectors
  Gamma1_flat=as.vector(Gamma1)
  Gamma2_flat=as.vector(Gamma2)
  Gamma3_flat=as.vector(Gamma3)
  Gamma5_flat=as.vector(Gamma5)
  Gamma6_flat=as.vector(Gamma6)
  GammaAll_flat = cbind(Gamma1_flat,Gamma2_flat,Gamma3_flat,Gamma5_flat,Gamma6_flat)
  compMatRows = rep(1:nrow(Gamma1), ncol(Gamma1))
  compMatCols = rep(1:ncol(Gamma1), each = nrow(Gamma1))
  
  #Obtain Unique Record-Level IDs and Append them to each record pair
  Full_MatchID_A=N_A$MatchID#600
  Full_MatchID_B=N_B$MatchID#1200
  
  A_LinkID=rep(Full_MatchID_A,1200)
  B_LinkID=rep(Full_MatchID_B,each=600)
  
  #Obtain Unique Block-Level IDs and Append them to each record pair
  Full_BlockID_A=rep(c(rep(Block_A_ID,each=15),rep(Block_A_ID,each=5)),30*40)
  Full_BlockID_B=c(rep(Block_B_ID,each=30*20*15),rep(Block_B_ID,each=30*20*15))

  DOB_yearP_gamma = as.vector(DOB_year_gamma*2)
  DOB_monthP_gamma=as.vector(DOB_month_gamma*2)
  GenderP_gamma = as.vector(Gender_gamma*2)
  RegionP_gamma=as.vector(apply(RegionP_comparison,2,'==',N_A$Region)*2)
  StatusP_gamma=as.vector(apply(StatusP_comparison,2,'==',N_A$Status)*2)
  TypeP_gamma=as.vector(apply(TypeP_comparison,2,'==',N_A$Type)*2)
  SESP_gamma=as.vector(apply(SESP_comparison,2,'==',N_A$SES)*2)
  IncomeP_gamma=as.vector(apply(IncomeP_comparison,2,function(x,Income_AP) abs(x-Income_AP)<=4,N_A$Income)*2)
  
  GammaZetaFullP=data.frame(A_LinkID,B_LinkID,Full_BlockID_A,Full_BlockID_B,DOB_yearP_gamma,
                            DOB_monthP_gamma,GenderP_gamma,RegionP_gamma,StatusP_gamma,
                            TypeP_gamma,SESP_gamma,IncomeP_gamma)
  
  numContVar = 1

  #running EM algorithm to get a good starting point
  g1 <- gammaKpar(N_A$Gender, N_B$Gender)
  g2 <- gammaKpar(N_A$DOB_year, N_B$DOB_year)
  g3 <- gammaKpar(N_A$DOB_month, N_B$DOB_month)
  g4 <- gammaKpar(N_A$Region, N_B$Region)
  g5 <- gammaKpar(N_A$Status, N_B$Status)
  g6 <- gammaKpar(N_A$Type, N_B$Type)
  g7 <- gammaKpar(N_A$SES, N_B$SES)
  g8 <- gammaNUMCK2par(N_A$Income, N_B$Income, cut.a=4)
  
  tc <- tableCounts(list(g1, g2, g3, g4, g5, g6, g7, g8), nobs.a = nrow(N_A), nobs.b = nrow(N_B))
  
  em.log <- emlinklog(tc, nobs.a = nrow(N_A), nobs.b = nrow(N_B))
  
  GammaZetaFullP$w = NA
  apply(em.log$patterns.w,1, function(x,wMat) { tempMat = sweep(wMat[,5:12],2,x[1:8]); whichRow = which(rowSums(abs(tempMat))==0); 
            GammaZetaFullP$w[whichRow] <<- x[10]},GammaZetaFullP)		
  
  GammaZetaFullP$w[which(is.na(GammaZetaFullP$w))] = min(GammaZetaFullP$w,na.rm=TRUE)
  GammaZetaFullP$w = exp(GammaZetaFullP$w)
  IDweights=data.frame(GammaZetaFullP[,1:4],GammaZetaFullP$w)
  names(IDweights) = c("A_LinkID", "B_LinkID","Full_BlockID_A","Full_BlockID_B","w")
  IDweights$C=0

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
      C_qr[w_qr<=quantile(w_qr,0.75)]=0
      
      if(sum(apply(C_qr,2,sum)>1)>0 | sum(apply(C_qr,1,sum)>1)>0){print('LSAP error'); break}
      
      IDweights$C[(IDweights$Full_BlockID_B==Block_B_ID[r]) & (IDweights$Full_BlockID_A==Block_A_ID[q])]=C_qr
    }
  }
  
  #Specify values for hyperparameters of prior distributions
  prior_BM_Region=c(1,1); prior_BM_Status=c(1,1); prior_BM_Type=c(1,1); prior_BM_SES=c(1,1); 
  prior_BU_Region=c(1,1); prior_BU_Status=c(1,1); prior_BU_Type=c(1,1); prior_BU_SES=c(1,1); 
  prior_CM_DOB=c(1,1,1); prior_CM_Gender=c(1,1)
  prior_CU_DOB=c(1,1,1); prior_CU_Gender=c(1,1)
  prior_CBU_DOB=c(1,1,1); prior_CBU_Gender=c(1,1)
  prior_pi=c(1,1)
  prior_BM_Income=c(0,0,1,5^2)
  prior_BU_Income=c(0,0,1,5^2)
  
  
  #Specify number of iterations (including burn-in) and initialize parameter matrices
  #nsim=500 lower number of replications, easier for testing
  nsim=1000
  
  theta_BM_Region=matrix(0,nrow=length(prior_BM_Region),ncol=nsim)
  theta_BM_Status=matrix(0,nrow=length(prior_BM_Status),ncol=nsim)
  theta_BM_Type=matrix(0,nrow=length(prior_BM_Type),ncol=nsim)
  theta_BM_SES=matrix(0,nrow=length(prior_BM_SES),ncol=nsim)
  theta_BM_Income=matrix(0,nrow=numContVar*2,ncol=nsim)
  theta_BU_Region=matrix(0,nrow=length(prior_BU_Region),ncol=nsim)
  theta_BU_Status=matrix(0,nrow=length(prior_BU_Status),ncol=nsim)
  theta_BU_Type=matrix(0,nrow=length(prior_BU_Type),ncol=nsim)
  theta_BU_SES=matrix(0,nrow=length(prior_BU_SES),ncol=nsim)
  theta_BU_Income=matrix(0,nrow=numContVar*2,ncol=nsim)
  theta_CM_DOB=matrix(0,nrow=length(prior_CM_DOB),ncol=nsim)
  theta_CM_Gender=matrix(0,nrow=length(prior_CM_Gender),ncol=nsim)
  theta_CU_DOB=matrix(0,nrow=length(prior_CU_DOB),ncol=nsim)
  theta_CU_Gender=matrix(0,nrow=length(prior_CU_Gender),ncol=nsim)
  theta_CBU_DOB=matrix(0,nrow=length(prior_CBU_DOB),ncol=nsim)
  theta_CBU_Gender=matrix(0,nrow=length(prior_CBU_Gender),ncol=nsim)
  
  #Initialize Vectors to store links and other linking statistics
  BlockDesignation=matrix(0,nrow=length(Block_A_ID),ncol=nsim)
  LinkDesignation=matrix(0,nrow=nrow(N_A),ncol=nsim)
  BlockMove=matrix(0,nrow=length(Block_A_ID),ncol=nsim)
  BlockMoveType=matrix(0,nrow=length(Block_A_ID),ncol=nsim)
  MoveProbability=matrix(0,nrow=length(Block_A_ID),ncol=nsim)
  
  ###Initialize blocking and blocking designation matrices
  allPairs = expand.grid(1:length(N_A$S),1:length((N_B$T)))
  vecAllBlock = apply(allPairs, 1, function(x) return(N_A$S[x[1]]*numRound + N_B$T[x[2]]))
  vecLinked = rep(0, length(vecAllBlock))
  numInA = nrow(N_A)
  numInB = nrow(N_B)
  
  B_S=matrix(N_A$S,nrow=nrow(N_A),ncol=nrow(N_B))
  B_T=t(matrix(N_B$T,nrow=nrow(N_B),ncol=nrow(N_A)))
  B_ST=matrix(0,nrow=nrow(N_A),ncol=nrow(N_B))
  
  B=matrix(0,nrow=length(Block_A_ID),ncol=length(Block_B_ID))
  
  ###Write function to project partitions of B onto the linking space B_ST
  ###Inputs for function will be B, numRow, numCol. The output will be B_ST
  block_projection=function(B, numRow, numCol){
    ST=which(B==1,arr.ind=T)
    vecLinked[1:length(vecLinked)] = 0
    linkBlock = c(ST[,1]*100 + ST[,2])
    vecLinked[which(vecAllBlock %in% linkBlock)] = 1
    return(matrix(vecLinked,nrow=numRow,ncol=numCol))
  }  
  ###Specify starting values for parameters, blocking, and linking matrix
  #Starting value for blocks will be a random permutation of the indexes
  BlockDesignation[,1]=sample(Block_B_T,size=30,replace=FALSE)
  
  ###Function to transform block row and column designations into a matrix of 0's and 1's
  create_block_matrix=function(Block_Designation,B){
    A=data.frame(seq_along(Block_Designation),BlockDesignation)
    B_mat=matrix(0,nrow=nrow(B),ncol=ncol(B))
    for (i in 1:dim(A)[1]){
      B_mat[A[i,1],A[i,2]]=1}
    return(B_mat)}
  
  #Need to transform BlockDesignation to a matrix B
  B=create_block_matrix(BlockDesignation[,1],B)
  #Project block designations onto B_ST
  B_ST=block_projection(B,numInA,numInB)
  
  #Start with an empty linking matrix
  C=matrix(0,nrow=nrow(N_A),ncol=nrow(N_B))
  
  
  #Specify starting values for blocking and linking parameters 
  theta_BM_Region[,1]=c(.3,.7)
  theta_BU_Region[,1]=c(.66,.34)
  theta_BM_Status[,1]=c(.25,.75)
  theta_BU_Status[,1]=c(.75,.25)
  theta_BM_Type[,1]=c(.4,.6)
  theta_BU_Type[,1]=c(.6,.4)
  theta_BM_SES[,1]=c(.35,.65)
  theta_BU_SES[,1]=c(.65,.35)
  #theta_BM_Income[,1]=c(.2,.8)
  #theta_BU_Income[,1]=c(.8,.2)
  theta_BM_Income[,1]=c(0,5^2)
  theta_BU_Income[,1]=c(0,10^2)
  theta_CM_DOB[,1]=c(.05,25,.7)
  theta_CU_DOB[,1]=c(.7,.20,.10)
  theta_CM_Gender[,1]=c(.25,.75)
  theta_CU_Gender[,1]=c(.75,.25)
  theta_CBU_DOB[,1]=c(.8,.125,.075)
  theta_CBU_Gender[,1]=c(.8,.2)

  nlink = 25
  
  for(t in 2:nsim){
    ########################################################################################################################
    #Sample the Posterior Distribution of linking and blocking parameters
    B_ST=block_projection(B, numInA, numInB)
    theta_BM_Region[,t]=rdirichlet(1,c(prior_BM_Region[1]+sum(Zeta1*B),prior_BM_Region[2]+sum(Zeta2*B)))
    theta_BU_Region[,t]=rdirichlet(1,c(prior_BU_Region[1]+sum(Zeta1*(1-B)),prior_BU_Region[2]+sum(Zeta2*(1-B))))
    theta_BM_Status[,t]=rdirichlet(1,c(prior_BM_Status[1]+sum(Zeta3*B),prior_BM_Status[2]+sum(Zeta4*B)))
    theta_BU_Status[,t]=rdirichlet(1,c(prior_BU_Status[1]+sum(Zeta3*(1-B)),prior_BU_Status[2]+sum(Zeta4*(1-B))))
    theta_BM_Type[,t]=rdirichlet(1,c(prior_BM_Type[1]+sum(Zeta5*B),prior_BM_Type[2]+sum(Zeta6*B)))
    theta_BU_Type[,t]=rdirichlet(1,c(prior_BU_Type[1]+sum(Zeta5*(1-B)),prior_BU_Type[2]+sum(Zeta6*(1-B))))
    theta_BM_SES[,t]=rdirichlet(1,c(prior_BM_SES[1]+sum(Zeta7*B),prior_BM_Type[2]+sum(Zeta8*B)))
    theta_BU_SES[,t]=rdirichlet(1,c(prior_BU_SES[1]+sum(Zeta7*(1-B)),prior_BU_Type[2]+sum(Zeta8*(1-B))))
    #theta_BM_Income[,t]=rdirichlet(1,c(prior_BM_Income[1]+sum(Zeta9*B),prior_BM_Income[2]+sum(Zeta10*B)))
    #theta_BU_Income[,t]=rdirichlet(1,c(prior_BU_Income[1]+sum(Zeta9*(1-B)),prior_BU_Income[2]+sum(Zeta10*(1-B))))
    n_M = sum(B)
    v_nM = prior_BM_Income[3] + n_M 
    sig_nM = (prior_BM_Income[3]*prior_BM_Income[4] + (n_M - 1)*var(Income_Diff[B==1]) + 
          prior_BM_Income[2]*n_M/(prior_BM_Income[2] + n_M)*(mean(Income_Diff[B==1])- prior_BM_Income[1])^2)/v_nM
    mu_nM = prior_BM_Income[2]/(prior_BM_Income[2] + n_M)*prior_BM_Income[1] + n_M/(prior_BM_Income[2] + n_M)*mean(Income_Diff[B==1])
    k_nM = prior_BM_Income[2]+n_M
    theta_BM_Income[2,t]=rinvgamma(1,v_nM/2,v_nM*sig_nM/2)
    theta_BM_Income[1,t]=rnorm(1,mean=mu_nM, sqrt(theta_BM_Income[2,t]/k_nM))
    
    
    n_U = sum(1-B)
    v_nU = prior_BU_Income[3] + n_U 
    sig_nU = (prior_BU_Income[3]*prior_BU_Income[4] + (n_U - 1)*var(Income_Diff[B==0]) + 
               prior_BU_Income[2]*n_U/(prior_BU_Income[2] + n_U)*(mean(Income_Diff[B==0])- prior_BU_Income[1])^2)/v_nU
    mu_nU = prior_BU_Income[2]/(prior_BU_Income[2] + n_U)*prior_BU_Income[1] + n_U/(prior_BU_Income[2] + n_U)*mean(Income_Diff[B==0])
    k_nU = prior_BU_Income[2]+n_U
    theta_BU_Income[2,t]=rinvgamma(1,v_nU/2,v_nU*sig_nU/2)
    theta_BU_Income[1,t]=rnorm(1,mean=mu_nU, sqrt(theta_BU_Income[2,t]/k_nU))
    
    theta_CM_DOB[,t]=rdirichlet(1,c(prior_CM_DOB[1]+sum(Gamma1*C*B_ST),prior_CM_DOB[2]+sum(Gamma2*C*B_ST),prior_CM_DOB[3]+sum(Gamma3*C*B_ST)))
    theta_CU_DOB[,t]=rdirichlet(1,c(prior_CU_DOB[1]+sum(Gamma1*(1-C)*B_ST),prior_CU_DOB[2]+sum(Gamma2*(1-C)*B_ST),prior_CU_DOB[3]+sum(Gamma3*(1-C)*B_ST)))
    theta_CM_Gender[,t]=rdirichlet(1,c(prior_CM_Gender[1]+sum(Gamma5*C*B_ST),prior_CM_Gender[2]+sum(Gamma6*C*B_ST)))
    theta_CU_Gender[,t]=rdirichlet(1,c(prior_CU_Gender[1]+sum(Gamma5*(1-C)*B_ST),prior_CU_Gender[2]+sum(Gamma6*(1-C)*B_ST)))
    theta_CBU_DOB[,t]=rdirichlet(1,c(prior_CBU_DOB[1]+sum(Gamma1*(1-B_ST)),prior_CBU_DOB[2]+sum(Gamma2*(1-B_ST)),prior_CBU_DOB[3]+sum(Gamma3*(1-B_ST))))
    theta_CBU_Gender[,t]=rdirichlet(1,c(prior_CBU_Gender[1]+sum(Gamma5*(1-B_ST)),prior_CBU_Gender[2]+sum(Gamma6*(1-B_ST))))
    
    ########################################################################################################################
    #Within every block, update the linking configurations C_ST
    b=which(B==1,arr.ind=TRUE); b=b[order(b[,1]),]
    C_flat = as.vector(C)
    for(s in 1:dim(b)[1]){
      blockST = b[s,1]*numRound + b[s,2]
      C_ST = C_flat[which(vecAllBlock==blockST)]
      rowVec = compMatRows[which(vecAllBlock==blockST)]
      colVec = compMatCols[which(vecAllBlock==blockST)]
      numRowUse = length(unique(rowVec))
      numColUse = length(unique(colVec))
      
      #Repeat the Gibbs Sampler nlink times to sufficiently sample the blocks
      for(j in 1:nlink){
        for(i_s in unique(rowVec)){
          #Reset the link result for row i_st
          C_ST[rowVec == i_s]=0
          #Extract units in C_ST that have links, not including the individual currently linked to i
          C_linked=unique(colVec[C_ST == 1])
          #If there are no links, create empty link vector. Then extract non-linked individuals in dataset B
          if(length(C_linked)==0){
            C_unlinked=unique(colVec)
          } else{C_unlinked=setdiff(unique(colVec),C_linked)}
          
          #Extract matrix of DOB and Gender Gamma
          GammaAll_flatU = GammaAll_flat[which(compMatRows == i_s & compMatCols %in% C_unlinked),]          
          
          #Calculate ratio of likelihoods for Gamma comparisons
          num=apply(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatU),2,prod)
          den=apply(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatU),2,prod)
          Likelihood=num/den
          
          #Calculate the probability of individual i not linking
          #This will be a function of the linkage overlap only in block ST
          p_nolink=abs((max(numRowUse, numColUse)-sum(C_ST))*(min(numRowUse, numColUse)-sum(C_ST)+prior_pi[2]-1)/(sum(C_ST)+prior_pi[1]))
          
          #Parse together possible moves and move probability
          C_unlinked=c(C_unlinked, max(compMatCols)+1)
          C_prob=c(Likelihood,p_nolink)/sum(Likelihood,p_nolink)
          
          #Sample new bipartite link for individual i_st
          link_designation=sample(C_unlinked,size=1,prob=C_prob)
          
          #Update linking matrix within block if a new link designation was selected
          if(link_designation<=max(compMatCols)){C_ST[which(rowVec == i_s & colVec == link_designation)]=1}
        }}
      
      #Update the linking matrix C with the new link configuration within block C_ST
      C_flat[which(vecAllBlock==blockST)] = C_ST
      
      #If the iteration is greater than 20, save the updated configuration in block ST as the proposal configuration for future MH steps
      if(t>20){
        IDweights$C[IDweights$Full_BlockID_A==Block_A_ID[b[s,1]] & IDweights$Full_BlockID_B==Block_B_ID[b[s,2]]]=C_ST
        }
    }
    
    ########################################################################################################################
    #Propose an update for each block s
    for(s in 1:dim(b)[1]){
      b=which(B==1,arr.ind=TRUE); b=b[order(b[,1]),]
      #Sample a new block designation 'r' for the proposal move. 
      #'r' will contain the set of indeces in B without 't' that 's' is currently partitioned with.
      r=sample(setdiff(c(1:dim(B)[2]),b[s,2]),size=1)
      #####################################################################
      #####There will be two possible moves depending on whether r is partitioned with another block in A or not
      ###Move 1: r is not currently partitioned with any block in A
      if(sum(b[,2]==r)==0){
        ###First, we need to calculate the potential links for block B_SR
        #Extract the elements that correspond to block SR
        blockSR = b[s,1]*numRound + r
        indBlockSR = which(vecAllBlock==blockSR)
        
        blockST = b[s,1]*numRound + b[s,2]
        indBlockST = which(vecAllBlock==blockST)
        
        blockSR = b[s,1]*numRound + r
        indBlockSR = which(vecAllBlock==blockSR)
        
        blockST = b[s,1]*numRound + b[s,2]
        indBlockST = which(vecAllBlock==blockST)
        
        C_SR=IDweights$C[IDweights$Full_BlockID_A==Block_A_ID[b[s,1]] & IDweights$Full_BlockID_B==Block_B_ID[r]]
        C_ST = C_flat[indBlockST]
        
        ##MH acceptance probability will have 4 parts
        ##Ratio of theta_BM and theta_BU
        theta_B_num=prod(theta_BM_Region[,t]^c(Zeta1[b[s,1],r],Zeta2[b[s,1],r]))*prod(theta_BM_Status[,t]^c(Zeta3[b[s,1],r],Zeta4[b[s,1],r]))*
          prod(theta_BM_Type[,t]^c(Zeta5[b[s,1],r],Zeta6[b[s,1],r]))*prod(theta_BM_SES[,t]^c(Zeta7[b[s,1],r],Zeta8[b[s,1],r]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[b[s,1],r],Zeta10[b[s,1],r]))*
          dnorm(Income_Diff[b[s,1],r],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[b[s,1],b[s,2]],Zeta2[b[s,1],b[s,2]]))*prod(theta_BU_Status[,t]^c(Zeta3[b[s,1],b[s,2]],Zeta4[b[s,1],b[s,2]]))*
          prod(theta_BU_Type[,t]^c(Zeta5[b[s,1],b[s,2]],Zeta6[b[s,1],b[s,2]]))*prod(theta_BU_SES[,t]^c(Zeta7[b[s,1],b[s,2]],Zeta8[b[s,1],b[s,2]]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[b[s,1],b[s,2]],Zeta10[b[s,1],b[s,2]]))
          dnorm(Income_Diff[b[s,1],b[s,2]],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))
        theta_B_den=prod(theta_BM_Region[,t]^c(Zeta1[b[s,1],b[s,2]],Zeta2[b[s,1],b[s,2]]))*prod(theta_BM_Status[,t]^c(Zeta3[b[s,1],b[s,2]],Zeta4[b[s,1],b[s,2]]))*
          prod(theta_BM_Type[,t]^c(Zeta5[b[s,1],b[s,2]],Zeta6[b[s,1],b[s,2]]))*prod(theta_BM_SES[,t]^c(Zeta7[b[s,1],b[s,2]],Zeta8[b[s,1],b[s,2]]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[b[s,1],b[s,2]],Zeta10[b[s,1],b[s,2]]))*
          dnorm(Income_Diff[b[s,1],b[s,2]],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[b[s,1],r],Zeta2[b[s,1],r]))*prod(theta_BU_Status[,t]^c(Zeta3[b[s,1],r],Zeta4[b[s,1],r]))*
          prod(theta_BU_Type[,t]^c(Zeta5[b[s,1],r],Zeta6[b[s,1],r]))*prod(theta_BU_SES[,t]^c(Zeta7[b[s,1],r],Zeta8[b[s,1],r]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[b[s,1],r],Zeta10[b[s,1],r]))
          dnorm(Income_Diff[b[s,1],r],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))
        
        
        GammaAll_flatSR = GammaAll_flat[indBlockSR,]
        theta_C_num=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatSR*C_SR)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatSR*c(1-C_SR))))

        GammaAll_flatST = GammaAll_flat[indBlockST,]
        theta_C_den=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatST*C_ST)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatST*c(1-C_ST))))
        
        theta_CB_num=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatST)))
        
        theta_CB_den=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatSR)))
        
        
        acceptance=sum(log(theta_B_num/theta_B_den), (theta_C_num-theta_C_den),(theta_CB_num-theta_CB_den))
        
        
        #If the proposal move is accepted, update both the linking blocking and linking matrices
        if(log(runif(1))<acceptance){
          B[b[s,1],r]=1; B[b[s,1],b[s,2]]=0
          C_flat[indBlockST]=0; C_flat[indBlockSR]=C_SR
        }
        
        BlockMove[b[s,1],t]=r
        BlockMoveType[b[s,1],t]=1
        MoveProbability[b[s,1],t]=acceptance
      } else{
        #Move2
        ###First, we need to calculate the potential links for blocks B_SR and B_QT
        #Extract the elements that correspond to block SR and QT
        blockSR = b[s,1]*numRound + r
        indBlockSR = which(vecAllBlock==blockSR)
        
        blockST = b[s,1]*numRound + b[s,2]
        indBlockST = which(vecAllBlock==blockST)
        
        q=b[b[,2]==r][1]
        blockQT = b[q,1]*numRound + b[s,2]
        indBlockQT = which(vecAllBlock==blockQT)
        
        blockQR = b[q,1]*numRound + r
        indBlockQR = which(vecAllBlock==blockQR)
        
        C_SR=IDweights$C[indBlockSR]
        C_QT=IDweights$C[indBlockQT]
        C_QR=C_flat[indBlockQR]
        C_ST=C_flat[indBlockST]
           
        ##MH acceptance probability will have 4 parts
        ##Ratio of theta_BM and theta_BU
        theta_B_num=prod(theta_BM_Region[,t]^c(Zeta1[b[s,1],r],Zeta2[b[s,1],r]))*prod(theta_BM_Status[,t]^c(Zeta3[b[s,1],r],Zeta4[b[s,1],r]))*
          prod(theta_BM_Type[,t]^c(Zeta5[b[s,1],r],Zeta6[b[s,1],r]))*prod(theta_BM_SES[,t]^c(Zeta7[b[s,1],r],Zeta8[b[s,1],r]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[b[s,1],r],Zeta10[b[s,1],r]))*
          dnorm(Income_Diff[b[s,1],r],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[b[s,1],b[s,2]],Zeta2[b[s,1],b[s,2]]))*prod(theta_BU_Status[,t]^c(Zeta3[b[s,1],b[s,2]],Zeta4[b[s,1],b[s,2]]))*
          prod(theta_BU_Type[,t]^c(Zeta5[b[s,1],b[s,2]],Zeta6[b[s,1],b[s,2]]))*prod(theta_BU_SES[,t]^c(Zeta7[b[s,1],b[s,2]],Zeta8[b[s,1],b[s,2]]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[b[s,1],b[s,2]],Zeta10[b[s,1],b[s,2]]))*
          dnorm(Income_Diff[b[s,1],b[s,2]],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))*
          prod(theta_BM_Region[,t]^c(Zeta1[q,b[s,2]],Zeta2[q,b[s,2]]))*prod(theta_BM_Status[,t]^c(Zeta3[q,b[s,2]],Zeta4[q,b[s,2]]))*
          prod(theta_BM_Type[,t]^c(Zeta5[q,b[s,2]],Zeta6[q,b[s,2]]))*prod(theta_BM_SES[,t]^c(Zeta7[q,b[s,2]],Zeta8[q,b[s,2]]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[q,b[s,2]],Zeta10[q,b[s,2]]))*
          dnorm(Income_Diff[q,b[s,2]],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[q,r],Zeta2[q,r]))*prod(theta_BU_Status[,t]^c(Zeta3[q,r],Zeta4[q,r]))*
          prod(theta_BU_Type[,t]^c(Zeta5[q,r],Zeta6[q,r]))*prod(theta_BU_Region[,t]^c(Zeta7[q,r],Zeta8[q,r]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[q,r],Zeta10[q,r]))
          dnorm(Income_Diff[q,r],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))
        
        theta_B_den=prod(theta_BM_Region[,t]^c(Zeta1[b[s,1],b[s,2]],Zeta2[b[s,1],b[s,2]]))*prod(theta_BM_Status[,t]^c(Zeta3[b[s,1],b[s,2]],Zeta4[b[s,1],b[s,2]]))*
          prod(theta_BM_Type[,t]^c(Zeta5[b[s,1],b[s,2]],Zeta6[b[s,1],b[s,2]]))*prod(theta_BM_SES[,t]^c(Zeta7[b[s,1],b[s,2]],Zeta8[b[s,1],b[s,2]]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[b[s,1],b[s,2]],Zeta10[b[s,1],b[s,2]]))*
          dnorm(Income_Diff[b[s,1],b[s,2]],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[b[s,1],r],Zeta2[b[s,1],r]))*prod(theta_BU_Status[,t]^c(Zeta3[b[s,1],r],Zeta4[b[s,1],r]))*
          prod(theta_BU_Type[,t]^c(Zeta5[b[s,1],r],Zeta6[b[s,1],r]))*prod(theta_BU_SES[,t]^c(Zeta7[b[s,1],r],Zeta8[b[s,1],r]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[b[s,1],r],Zeta10[b[s,1],r]))*
          dnorm(Income_Diff[b[s,1],r],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))*
          prod(theta_BM_Region[,t]^c(Zeta1[q,r],Zeta2[q,r]))*prod(theta_BM_Status[,t]^c(Zeta3[q,r],Zeta4[q,r]))*
          prod(theta_BM_Type[,t]^c(Zeta5[q,r],Zeta6[q,r]))*prod(theta_BM_Status[,t]^c(Zeta7[q,r],Zeta8[q,r]))*
          #prod(theta_BM_Income[,t]^c(Zeta9[q,r],Zeta10[q,r]))*
          dnorm(Income_Diff[q,r],mean=theta_BM_Income[1,t],sqrt(theta_BM_Income[2,t]))*
          prod(theta_BU_Region[,t]^c(Zeta1[q,b[s,2]],Zeta2[q,b[s,2]]))*prod(theta_BU_Status[,t]^c(Zeta3[q,b[s,2]],Zeta4[q,b[s,2]]))*
          prod(theta_BU_Type[,t]^c(Zeta5[q,b[s,2]],Zeta6[q,b[s,2]]))*prod(theta_BU_SES[,t]^c(Zeta7[q,b[s,2]],Zeta8[q,b[s,2]]))*
          #prod(theta_BU_Income[,t]^c(Zeta9[q,b[s,2]],Zeta10[q,b[s,2]]))
          dnorm(Income_Diff[q,b[s,2]],mean=theta_BU_Income[1,t],sqrt(theta_BU_Income[2,t]))
        
        ##Ratio of theta_CM and theta_CU
        GammaAll_flatSR = GammaAll_flat[indBlockSR,]
        theta_C_num_1=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatSR*C_SR)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatSR*(1-C_SR))))
        
        GammaAll_flatQT = GammaAll_flat[indBlockQT,]
        theta_C_num_2=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatQT*C_QT)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatQT*(1-C_QT))))
        
        GammaAll_flatST = GammaAll_flat[indBlockST,]
        theta_C_den_1=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatST*C_ST)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatST*(1-C_ST))))
        
        GammaAll_flatQR = GammaAll_flat[indBlockQR,]
        theta_C_den_2=sum(log(c(theta_CM_DOB[,t],theta_CM_Gender[,t])^t(GammaAll_flatQR*C_QR)))+
          sum(log(c(theta_CU_DOB[,t],theta_CU_Gender[,t])^t(GammaAll_flatQR*(1-C_QR))))
        
        
        theta_CB_num_1=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatST)))
        
        theta_CB_num_2=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatQR)))        
        
        theta_CB_den_1=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatSR)))
        
        theta_CB_den_2=sum(log(c(theta_CBU_DOB[,t],theta_CBU_Gender[,t])^t(GammaAll_flatQT)))
        
        acceptance=sum(log(theta_B_num/theta_B_den), (theta_C_num_1-theta_C_den_1), (theta_C_num_2-theta_C_den_2),
                              (theta_CB_num_1-theta_CB_den_1), (theta_CB_num_2-theta_CB_den_2))
        
        #If the proposal move is accepted, update both the linking blocking and linking matrices
        if(log(runif(1))<acceptance){
          B[b[s,1],r]=1; B[b[s,1],b[s,2]]=0; B[q,b[s,2]]=1; B[q,r]=0
          C_flat[indBlockST] = 0; C_flat[indBlockSR]=C_SR; C_flat[indBlockQR]=0; C_flat[indBlockQT]=C_QT;
        }
        BlockMove[b[s,1],t]=r
        BlockMoveType[b[s,1],t]=2
        MoveProbability[b[s,1],t]=acceptance
      }
    }
    C = matrix(C_flat, numInA, numInB)
    #Export Designations of Linking Matrix after every iteration of updating the linking matrices
    a=which(C==1,arr.ind=TRUE)
    LinkDesignation[,t][apply(C,1,sum)>=1]=matrix(a[order(a[,1]),],ncol=2)[,2]
    #Export Designations of Blocking Matrix after every iteration of the MH within Gibbs sampler
    BlockDesignation[,t]=b[order(b[,1]),][,2]
    print(paste("Iter, ",t))
    
    ##################################################################################################
    ##Export Simulation Parameters
    ##################################################################################################
    # write.csv(theta_BM_Region,paste0("theta_BM_Region_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BM_Status,paste0("theta_BM_Status_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BM_Type,paste0("theta_BM_Type_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BM_Income,paste0("theta_BM_Income_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BM_SES,paste0("theta_BM_SES_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BU_Region,paste0("theta_BU_Region_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BU_Status,paste0("theta_BU_Status_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BU_Type,paste0("theta_BU_Type_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BU_Income,paste0("theta_BU_Income_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_BU_SES,paste0("theta_BU_SES_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CM_DOB,paste0("theta_CM_DOB_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CM_Gender,paste0("theta_CM_Gender_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CU_DOB,paste0("theta_CU_DOB_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CU_Gender,paste0("theta_CU_Gender_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CBU_DOB,paste0("theta_CBU_DOB_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(theta_CBU_Gender,paste0("theta_CBU_Gender_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(BlockDesignation,paste0("BlockDesignation_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
    # write.csv(LinkDesignation,paste0("LinkDesignation_MLBRL_R4_I4_D4_noD_",nrep,".csv"))
  }
  
  n=apply(LinkDesignation[,(nsim/2+1):nsim],2,function(x) sum(x!="0"))
  TP=apply(LinkDesignation[,(nsim/2+1):nsim],2,function(x) sum(x==TrueLinkDesignation))
  FP=n-TP
  
  Sensitivity[nrep]=mean(TP/450)
  Sensitivity_sd[nrep]=sd(TP/450)
  PPV[nrep]=mean(TP/n)
  PPV_sd[nrep]=sd(TP/n)
  F1[nrep]=mean(2*((TP/450)*(TP/n))/(TP/450+TP/n))
  F1_sd[nrep]=sd(2*((TP/450)*(TP/n))/(TP/450+TP/n))
  
  Accuracy[nrep]=mean(apply(BlockDesignation[,(nsim/2+1):nsim],2,function(x) sum(x==TrueBlockDesignation))/30)
  Accuracy_sd[nrep]=sd(apply(BlockDesignation[,(nsim/2+1):nsim],2,function(x) sum(x==TrueBlockDesignation))/30)
  
  
  write.csv(Sensitivity, sprintf("Sensitivity_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(Sensitivity_sd, sprintf("Sensitivity_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(PPV, sprintf("PPV_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(PPV_sd, sprintf("PPV_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(F1, sprintf("F1_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(F1_sd, sprintf("F1_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(Accuracy, sprintf("Accuracy_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))
  write.csv(Accuracy_sd, sprintf("Accuracy_sd_Joint_MLBRL_R%d_I%d_D%d_noD.csv", Regionerrorprob*10, IncomeError*10, DOBerrorprob*10))

  print(paste("Simulation ",nrep))
  Sys.sleep(0.01)
  flush.console()
}


