## MMSDM and LQM Bias Correction Functions




## Setting the working directory
setwd("D:/Dissertation Stuff/Dissertation Work")




## Loading necessary libraries
library(tidyverse)
library(corrplot)
library(GGally)
library(modQR)
library(MASS)
library(lubridate)
library(qtlmt)
library(stringr)
library(ggcorrplot)
library(ggpubr)
library(data.table)











##################### FUNCTIONS NEEDED LATER #####################

## Spread function with multiple values
myspread <- function(df, key, value){
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

## Gather function with multiple values
mygather <- function(Dataframe,key,value){
  ## Gathering
  Dataframe2 = Dataframe
  Dataframe2 = gather(Dataframe2,key='VARIABLE',value='VALUE',-DATE,factor_key=T)
  Dataframe2$STATION = sub("\\_.*", "",Dataframe2$VARIABLE)
  Dataframe2$VarName = sub(".*_", "",Dataframe2$VARIABLE)
  ## Reordering columns
  Dataframe2 = Dataframe2[,c(1,4,5,3)]
  ## Spreading
  Dataframe3 = spread(Dataframe2,VarName,VALUE)
  return(Dataframe3)
}

## Indicator column creator
Indicator_Creator = function(Dataframe,key){
  ## Extracting unique levels in key
  UniqueLevels = unique(Dataframe[,key])
  ## Creating new columns
  for(i in 2:length(UniqueLevels)){
    ## Unique levels
    NewVar = paste(key,UniqueLevels[i],sep='.')
    ## Adding columns to dataframe
    Dataframe[,NewVar] = ifelse(Dataframe[,key]==UniqueLevels[i],1,0)
  }
  return(Dataframe)
}

## Bias correction function
Bias_Correction = function(Day,Vec,Bias_Dat,Step,Decimals){
  ## Variables to be used later
  Bias_Dat$Blank = 900000000
  Dec3 = (as.numeric(Vec*10^(Decimals)) %% 10)/10
  Column = match(paste('Bias_',floor(Vec*10^(Decimals-1))/10^(Decimals-1),sep=''),names(Bias_Dat))
  ## Bias Magnitude
  Bias_Magnitude = ifelse(Column==ncol(Bias_Dat),
                          Bias_Dat[Day,ncol(Bias_Dat)],NA)
  Bias_Magnitude = ifelse(Column<ncol(Bias_Dat),
                          Bias_Dat[cbind(Day,Column)]+Dec3*(Bias_Dat[cbind(Day,(Column+1))]-Bias_Dat[cbind(Day,Column)]),
                          Bias_Magnitude)
  return(Bias_Magnitude)
}

## Bias correction function for standard QM
Bias_Correction.QM = function(Vec,Bias_Dat,Step,Decimals){
  ## Variables to be used later
  LastDec = (as.numeric(Vec*10^Decimals) %% 10)/10
  LeadingDec = floor(10^(Decimals-1) * Vec) / 10^(Decimals-1)
  Row1 = LeadingDec*10^(Decimals-1)+1
  Row2 = ifelse(Row1<(10^(Decimals-1)+1),Row1+1,NA)
  
  ## Bias if exact step
  Bias_Magnitude = ifelse(LastDec==0,Bias_Dat[Row1,'Bias'],NA)
  
  ## Bias if not exact step
  Bias_Magnitude = ifelse(LastDec!=0,
                          Bias_Dat[Row1,'Bias']+LastDec*(Bias_Dat[Row2,'Bias']-Bias_Dat[Row1,'Bias']),
                          Bias_Magnitude)
  return(Bias_Magnitude)
}

##################### FUNCTIONS NEEDED LATER #####################














##################### MULTIVARIATE DOWNSCALING #####################

## Specifying model
#Model_Name = 'CCCMA-CanESM2'
Model_Name = 'GFDL-ESM2G'
#Model_Name = 'HadGEM2-ES'
#Model_Name = 'MRI-ESM1'

## Reading in the data
A.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''))
X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_X.mat.csv',sep=''))
RCP45.X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/RCP45_X.mat.csv',sep=''))
RCP85.X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/RCP85_X.mat.csv',sep=''))
ColNames = colnames(RCP85.X.mat)[2:5]
#ColNames = colnames(RCP85.X.mat)[6:13] ## USED FOR ESM1




## This function will take in historical predictands (GHCN)
## and historcal predictors (GCM). The function will also 
## take in the predictors for RCP45 and RCP85 scenarios (IA).
## It will ouptut a list containing several datasets:
## 1. B_hat (historical)
## 1. A matrix of fitted values (historical)
## 2. A matrix of predicted values for RCP45
## 3. A matrix of predicted values for RCP85

#TMAX_Values = cbind.data.frame(X.mat$tasmax,X.mat$tasmin,A.mat[,seq(2,116,5)])
#r = cor(as.matrix(TMAX_Values),use='complete.obs')
#ggcorrplot(r,hc.order = F,type = "lower",lab = TRUE)


## Function for multivariate downscaling
MMSDM_Train = function(A.mat,X.mat,RCP45.X.mat,RCP85.X.mat,ColNames){
  
  ## Storing list
  MMSDM_Output = list()
  
  ## Parameter estimation
  Y = as.matrix(A.mat[,2:ncol(A.mat)])
  X = as.matrix(cbind(Int=1,X.mat[,ColNames]))
  B.hat = solve(t(X) %*% X) %*% t(X) %*% Y
  
  ## Predicted
  A.hat = X %*% B.hat
  
  ## Error matrix
  E.mat = A.mat
  E.mat[,2:ncol(A.mat)] = E.mat[,2:ncol(A.mat)] - A.hat
  
  ## Calculating sigma matrix
  S.mat = diag(apply(E.mat[,2:ncol(A.mat)],2,sd))
  C.mat = cor(E.mat[,2:ncol(A.mat)])
  Sigma.mat = S.mat %*% C.mat %*% S.mat
  
  ## Adding in multivariate normal noise
  A.hat.2 = data.frame(A.hat + mvrnorm(n=nrow(A.hat),mu=rep(0,(ncol(A.mat)-1)),Sigma=Sigma.mat))
  A.hat.2$DATE = A.mat$DATE
  A.hat.2 = A.hat.2[,c(ncol(A.hat.2),1:(ncol(A.hat.2)-1))]
  
  ## Gathering and back-standardizing
  A.hat.2.gather = mygather(A.hat.2)
  A.hat.2.gather$Day = yday(as.Date(A.hat.2.gather$DATE))
  A.hat.2.gather$Day = ifelse((leap_year(as.Date(A.hat.2.gather$DATE)) & 
                                 A.hat.2.gather$Day>59),(A.hat.2.gather$Day-1),A.hat.2.gather$Day)
  A.hat.2.full = A.hat.2.gather
  setnames(A.hat.2.full,'Tmax','TMAX')
  setnames(A.hat.2.full,'Tmin','TMIN')
  
  ## Changing to matrices
  RCP45.X = as.matrix(cbind(Int=1,RCP45.X.mat[,ColNames]))
  RCP85.X = as.matrix(cbind(Int=1,RCP85.X.mat[,ColNames]))
  
  ## Predicting the output for RCP45 and RCP85
  RCP45.A.hat = RCP45.X %*% B.hat
  RCP85.A.hat = RCP85.X %*% B.hat
  
  ## Adding in the multivariate noise
  RCP45.A.hat2 = data.frame(RCP45.A.hat + mvrnorm(n=nrow(RCP45.A.hat),
                                                  mu=rep(0,ncol(RCP45.A.hat)),
                                                  Sigma=Sigma.mat))
  RCP45.A.hat2$DATE = RCP45.X.mat$Date
  RCP45.A.hat2 = RCP45.A.hat2[,c(ncol(RCP45.A.hat2),1:(ncol(RCP45.A.hat2)-1))]
  RCP85.A.hat2 = data.frame(RCP85.A.hat + mvrnorm(n=nrow(RCP85.A.hat),
                                                  mu=rep(0,ncol(RCP85.A.hat)),
                                                  Sigma=Sigma.mat))
  RCP85.A.hat2$DATE = RCP85.X.mat$Date
  RCP85.A.hat2 = RCP85.A.hat2[,c(ncol(RCP85.A.hat2),1:(ncol(RCP85.A.hat2)-1))]
  
  ## Gathering the future GCM data
  RCP45.gather = mygather(RCP45.A.hat2)
  RCP85.gather = mygather(RCP85.A.hat2)
  
  ## Back standardizing future GCM data RCP4.5
  RCP45.gather$Day = yday(as.Date(RCP45.gather$DATE))
  RCP45.gather$Day = ifelse((leap_year(as.Date(RCP45.gather$DATE)) & 
                               RCP45.gather$Day>59),(RCP45.gather$Day-1),RCP45.gather$Day)
  RCP45.full = RCP45.gather
  setnames(RCP45.full,'Tmax','TMAX')
  setnames(RCP45.full,'Tmin','TMIN')
  
  ## Back standardizing future GCM data RCP8.5
  RCP85.gather$Day = yday(as.Date(RCP85.gather$DATE))
  RCP85.gather$Day = ifelse((leap_year(as.Date(RCP85.gather$DATE)) & 
                               RCP85.gather$Day>59),(RCP85.gather$Day-1),RCP85.gather$Day)
  RCP85.full = RCP85.gather
  setnames(RCP85.full,'Tmax','TMAX')
  setnames(RCP85.full,'Tmin','TMIN')
  
  
  
  ### Temperature bias correction ###
  
  ## Observed historical data
  A.mat.gather = mygather(A.mat)
  A.mat.gather$Day = yday(as.Date(A.mat.gather$DATE))
  A.mat.gather$Day = ifelse((leap_year(as.Date(A.mat.gather$DATE)) & 
                               A.mat.gather$Day>59),(A.mat.gather$Day-1),A.mat.gather$Day)
  Day.MeansHist.Tmax = A.mat.gather %>%
    group_by(Day) %>%
    summarise(MeanTemp = mean(Tmax))
  Day.MeansHist.Tmin = A.mat.gather %>%
    group_by(Day) %>%
    summarise(MeanTemp = mean(Tmin))
  
  ## Fitted maximums
  Fitted = A.hat.2.full %>%
    group_by(Day) %>%
    summarise(MeanTemp = mean(TMAX))
  Fitted$Type = 'Fitted Values'
  
  ## Fitted minimums
  Fitted.min = A.hat.2.full %>%
    group_by(Day) %>%
    summarise(MeanTemp = mean(TMIN))
  Fitted.min$Type = 'Fitted Values'
  
  ## Difference between the two
  GCM_TMAX_Bias = Day.MeansHist.Tmax$MeanTemp - Fitted$MeanTemp
  GCM_TMIN_Bias = Day.MeansHist.Tmin$MeanTemp - Fitted.min$MeanTemp
  
  GCM_TMAX_Bias.lo = loess(GCM_TMAX_Bias ~ Day.MeansHist.Tmax$Day,span=.5)
  GCM_TMIN_Bias.lo = loess(GCM_TMIN_Bias ~ Day.MeansHist.Tmin$Day,span=.5)
  
  GCM_Bias = cbind.data.frame(Day=Day.MeansHist.Tmax$Day,
                              TMAX.Bias=GCM_TMAX_Bias.lo$fit,
                              TMIN.Bias=GCM_TMIN_Bias.lo$fit)
  
  ## Correcting
  RCP45.full = left_join(RCP45.full,GCM_Bias,by='Day')
  RCP45.full$TMAX.Corrected = RCP45.full$TMAX+RCP45.full$TMAX.Bias
  RCP45.full$TMIN.Corrected = RCP45.full$TMIN+RCP45.full$TMIN.Bias
  
  RCP85.full = left_join(RCP85.full,GCM_Bias,by='Day')
  RCP85.full$TMAX.Corrected = RCP85.full$TMAX+RCP85.full$TMAX.Bias
  RCP85.full$TMIN.Corrected = RCP85.full$TMIN+RCP85.full$TMIN.Bias
  
  ### Temperature bias correction ###
  
  
  
  ### Back transforming precipitation ###
  
  A.mat.gather = add_column(A.mat.gather,PRCP=(exp(A.mat.gather$R)-1),.after=4)
  A.hat.2.full = add_column(A.hat.2.full,PRCP=(exp(A.hat.2.full$R)-1),.after=4)
  RCP45.full = add_column(RCP45.full,PRCP=(exp(RCP45.full$R)-1),.after=4)
  RCP85.full = add_column(RCP85.full,PRCP=(exp(RCP85.full$R)-1),.after=4)
  
  ## Fixing Poc and PRCP
  A.hat.2.full$Poc = ifelse(A.hat.2.full$Poc>1,1,A.hat.2.full$Poc)
  A.hat.2.full$Poc = ifelse(A.hat.2.full$Poc<0,0,A.hat.2.full$Poc)
  
  ## Fixing Poc and PRCP
  RCP45.full$Poc = ifelse(RCP45.full$Poc>1,1,RCP45.full$Poc)
  RCP45.full$Poc = ifelse(RCP45.full$Poc<0,0,RCP45.full$Poc)
  
  ## Fixing Poc and PRCP
  RCP85.full$Poc = ifelse(RCP85.full$Poc>1,1,RCP85.full$Poc)
  RCP85.full$Poc = ifelse(RCP85.full$Poc<0,0,RCP85.full$Poc)
  
  
  ### Back transforming precipitation ###
  
  
  ## Adding to the output list
  MMSDM_Output[['A.mat']] = A.mat.gather
  MMSDM_Output[['X.mat']] = X.mat
  MMSDM_Output[['Y']] = Y
  MMSDM_Output[['X']] = X
  MMSDM_Output[['X.mat.RCP45']] = RCP45.X.mat
  MMSDM_Output[['X.mat.RCP85']] = RCP85.X.mat
  MMSDM_Output[['B.hat']] = B.hat
  MMSDM_Output[['E.mat']] = E.mat
  MMSDM_Output[['S.mat']] = S.mat
  MMSDM_Output[['C.mat']] = C.mat
  MMSDM_Output[['Sigma.mat']] = Sigma.mat
  MMSDM_Output[['A.hat.2']] = A.hat.2
  MMSDM_Output[['A.hat.2.full']] = A.hat.2.full
  MMSDM_Output[['RCP45.A.hat2']] = RCP45.A.hat2
  MMSDM_Output[['RCP85.A.hat2']] = RCP85.A.hat2
  MMSDM_Output[['RCP45.Downscaled']] = RCP45.full
  MMSDM_Output[['RCP85.Downscaled']] = RCP85.full
  
  ## Returning output list
  return(MMSDM_Output)
}

##################### MULTIVARIATE DOWNSCALING #####################








##################### LQM BIAS CORRECTION #####################

## Function for LOESS quantile mapping
LOESS_QM = function(MMSDM,Variable,Step,Decimals,Span){
  
  ## Initiating a storage list
  StorageList = list()
  StorageList[['Variable']] = Variable
  StorageList[['Step']] = Step
  StorageList[['Decimals']] = Decimals
  
  ## Creating step vector & new variable name
  Step_Vec = seq(0,1,Step)
  NewVar = paste(Variable,'.Corrected',sep='')
  
  ## For 360 day years
  DayVec = 1:365
  #DayVec = c(1:150,152:211,213:242,244:303,305:364)
  
  ## Looping through the step
  for(i in 1:length(Step_Vec)){
    
    ## Selecting step
    Step_Value = Step_Vec[i]
    
    ## Printing counter
    print(paste('Bias Data Step:',Step_Value))
    
    ## Historical
    Observed = MMSDM$A.mat %>%
      group_by(Day) %>%
      summarise(QuantileR = quantile(!!sym(Variable),Step_Value)) %>%
      as.data.frame()
    Observed$Type = 'Observed Historical'
    ## Fitted
    Fitted = MMSDM$A.hat.2.full %>%
      group_by(Day) %>%
      summarise(QuantileR = quantile(!!sym(Variable),Step_Value))
    Fitted$Type = 'Fitted Values'
    
    
    ## Creating dataframe containing biases
    if(i == 1){
      ## Estimating bias
      GCM_PRCP_Bias.lo = loess((Observed$QuantileR-Fitted$QuantileR) ~ 
                                 Observed$Day,span=Span)
      Bias_Dat = cbind.data.frame(Day=DayVec,value=predict(GCM_PRCP_Bias.lo))
      colnames(Bias_Dat)[i+1] = paste('Bias_',Step_Value,sep='')
    } else{
      ## Estimating bias
      GCM_PRCP_Bias.lo = loess((Observed$QuantileR-Fitted$QuantileR) ~ 
                                 Observed$Day,span=Span)
      Bias_Dat = cbind.data.frame(Bias_Dat,value=predict(GCM_PRCP_Bias.lo))
      colnames(Bias_Dat)[i+1] = paste('Bias_',Step_Value,sep='')
    }
  }
  
  ## For 360 day years
  FullLength = data.frame(Day=1:365)
  Bias_Dat = left_join(FullLength,Bias_Dat,by='Day')
  
  ## Storing the bias dataset
  StorageList[['Bias_Dat']] = Bias_Dat
  
  
  
  ## LQM bias correction conducted here
  
  ## Observed
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>%
    group_by(Day) %>%
    mutate(!!paste('Quantile.',Variable,sep='') := 
             round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
  
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>%
    mutate(!!paste('Bias_Factor',Variable,sep='') := 
             Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat))
  
  MMSDM$A.hat.2.full[,NewVar] = MMSDM$A.hat.2.full[,Variable] + 
    MMSDM$A.hat.2.full[,paste('Bias_Factor',Variable,sep='')]
  
  ## RCP4.5
  if('RCP45.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      group_by(Day) %>%
      mutate(!!paste('Quantile.',Variable,sep='') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      mutate(!!paste('Bias_Factor',Variable,sep='') :=
               Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat))
    
    MMSDM$RCP45.Downscaled[,NewVar] = MMSDM$RCP45.Downscaled[,Variable] + 
      MMSDM$RCP45.Downscaled[,paste('Bias_Factor',Variable,sep='')]
  }
  
  ## RCP8.5
  if('RCP85.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      group_by(Day) %>%
      mutate(!!paste('Quantile.',Variable,sep='') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      mutate(!!paste('Bias_Factor',Variable,sep='') := 
               Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat))
    
    MMSDM$RCP85.Downscaled[,NewVar] = MMSDM$RCP85.Downscaled[,Variable] + 
      MMSDM$RCP85.Downscaled[,paste('Bias_Factor',Variable,sep='')]
  }
  
  ## Adding to storage list
  StorageList[['MMSDM']] = MMSDM
  
  ## Returning data list
  return(StorageList)
}


## Specifying downscaled data
## CanESM2: set.seed(318)
## GFDL ESM2G: set.seed(31892)
## Hadley Center: set.seed(3189260)
## ESM1: set.seed(32608)
set.seed(31892)
MMSDM1 = MMSDM_Train(A.mat,X.mat,RCP45.X.mat,RCP85.X.mat,ColNames)
LQM1 = LOESS_QM(MMSDM1,'R',0.01,3)
LQM2 = LOESS_QM(LQM1$MMSDM,'Srad',0.01,3)


## Writing to csv files
PRCP_Bias = LQM1$Bias_Dat
SRAD_Bias = LQM2$Bias_Dat
Observed = LQM2$MMSDM$A.mat
Fitted = LQM2$MMSDM$A.hat.2.full
RCP45 = LQM2$MMSDM$RCP45.Downscaled
RCP85 = LQM2$MMSDM$RCP85.Downscaled

## Adding bias corrected precipitation
Fitted$PRCP.Corrected = Fitted$Poc*(exp(Fitted$R.Corrected) - 1)
Fitted$PRCP.Corrected = ifelse(Fitted$PRCP.Corrected<0,0,Fitted$PRCP.Corrected)
RCP45$PRCP.Corrected = RCP45$Poc*(exp(RCP45$R.Corrected) - 1)
RCP45$PRCP.Corrected = ifelse(RCP45$PRCP.Corrected<0,0,RCP45$PRCP.Corrected)
RCP85$PRCP.Corrected = RCP85$Poc*(exp(RCP85$R.Corrected) - 1)
RCP85$PRCP.Corrected = ifelse(RCP85$PRCP.Corrected<0,0,RCP85$PRCP.Corrected)


## Correcting negative solar radiation
Fitted$Srad.Corrected = ifelse(Fitted$Srad.Corrected<0,0,Fitted$Srad.Corrected)
RCP45$Srad.Corrected = ifelse(RCP45$Srad.Corrected<0,0,RCP45$Srad.Corrected)
RCP85$Srad.Corrected = ifelse(RCP85$Srad.Corrected<0,0,RCP85$Srad.Corrected)


## Writing to CSVs
#write.csv(PRCP_Bias,paste('Data - MMSDM/',Model_Name,'/Future Climate/PRCP_Bias.csv',sep=''),row.names=F)
#write.csv(SRAD_Bias,paste('Data - MMSDM/',Model_Name,'/Future Climate/SRAD_Bias.csv',sep=''),row.names=F)
#write.csv(Observed,paste('Data - MMSDM/',Model_Name,'/Future Climate/Observed.csv',sep=''),row.names=F)
#write.csv(Fitted,paste('Data - MMSDM/',Model_Name,'/Future Climate/Fitted.csv',sep=''),row.names=F)
#write.csv(RCP45,paste('Data - MMSDM/',Model_Name,'/Future Climate/RCP45.csv',sep=''),row.names=F)
#write.csv(RCP85,paste('Data - MMSDM/',Model_Name,'/Future Climate/RCP85.csv',sep=''),row.names=F)

##################### LQM BIAS CORRECTION #####################














##################### LQM VS QM BIAS CORRECTION #####################

Vec = MMSDM$A.hat.2.full$Quantile.R
Day = MMSDM$A.hat.2.full$Day
MMSDM = MMSDM1
Variable='R'
Step = 0.1
Span = 0.05
Decimals = 2


## Function for LOESS quantile mapping
LOESS_QM = function(MMSDM,Variable,Step,Decimals,Span){
  
  ## Initiating a storage list
  StorageList = list()
  StorageList[['Variable']] = Variable
  StorageList[['Step']] = Step
  StorageList[['Decimals']] = Decimals
  
  ## Creating step vector & new variable name
  Step_Vec = seq(0,1,Step)
  NewVar = paste(Variable,'.Corrected',sep='')
  
  ## For 360 day years
  DayVec = 1:365
  #DayVec = c(1:150,152:211,213:242,244:303,305:364)
  
  ## Looping through the step
  for(i in 1:length(Step_Vec)){
    
    ## Selecting step
    Step_Value = Step_Vec[i]
    
    ## Printing counter
    print(paste('Bias Data Step:',Step_Value))
    
    ## Historical
    Observed = MMSDM$A.mat %>%
      group_by(Day) %>%
      summarise(QuantileR = quantile(!!sym(Variable),Step_Value)) %>%
      as.data.frame()
    Observed$Type = 'Observed Historical'
    ## Fitted
    Fitted = MMSDM$A.hat.2.full %>%
      group_by(Day) %>%
      summarise(QuantileR = quantile(!!sym(Variable),Step_Value))
    Fitted$Type = 'Fitted Values'
    
    
    ## Creating dataframe containing biases
    if(i == 1){
      ## Estimating bias
      GCM_PRCP_Bias.lo = loess((Observed$QuantileR-Fitted$QuantileR) ~ 
                                 Observed$Day,span=Span)
      Bias_Dat = cbind.data.frame(Day=DayVec,value=predict(GCM_PRCP_Bias.lo))
      colnames(Bias_Dat)[i+1] = paste('Bias_',Step_Value,sep='')
    } else{
      ## Estimating bias
      GCM_PRCP_Bias.lo = loess((Observed$QuantileR-Fitted$QuantileR) ~ 
                                 Observed$Day,span=Span)
      Bias_Dat = cbind.data.frame(Bias_Dat,value=predict(GCM_PRCP_Bias.lo))
      colnames(Bias_Dat)[i+1] = paste('Bias_',Step_Value,sep='')
    }
  }
  
  ## For 360 day years
  FullLength = data.frame(Day=1:365)
  Bias_Dat = left_join(FullLength,Bias_Dat,by='Day')
  
  ## Storing the bias dataset
  StorageList[['Bias_Dat']] = Bias_Dat
  
  
  
  ## LQM bias correction conducted here
  
  ## Observed
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>% group_by(Day) %>%
    mutate(!!paste('Quantile.',Variable,sep='') := 
             round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
  
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>%
    mutate(!!paste('Bias_Factor',Variable,sep='') := 
             Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat,
                             Step,Decimals))
  
  MMSDM$A.hat.2.full[,NewVar] = MMSDM$A.hat.2.full[,Variable] + 
    MMSDM$A.hat.2.full[,paste('Bias_Factor',Variable,sep='')]
  
  ## RCP4.5
  if('RCP45.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      group_by(Day) %>%
      mutate(!!paste('Quantile.',Variable,sep='') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      mutate(!!paste('Bias_Factor',Variable,sep='') :=
               Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat,
                               Step,Decimals))
    
    MMSDM$RCP45.Downscaled[,NewVar] = MMSDM$RCP45.Downscaled[,Variable] + 
      MMSDM$RCP45.Downscaled[,paste('Bias_Factor',Variable,sep='')]
  }
  
  ## RCP8.5
  if('RCP85.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      group_by(Day) %>%
      mutate(!!paste('Quantile.',Variable,sep='') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      mutate(!!paste('Bias_Factor',Variable,sep='') := 
               Bias_Correction(Day,!!sym(paste('Quantile.',Variable,sep='')),Bias_Dat,
                               Step,Decimals))
    
    MMSDM$RCP85.Downscaled[,NewVar] = MMSDM$RCP85.Downscaled[,Variable] + 
      MMSDM$RCP85.Downscaled[,paste('Bias_Factor',Variable,sep='')]
  }
  
  ## Adding to storage list
  StorageList[['MMSDM']] = MMSDM
  
  ## Returning data list
  return(StorageList)
}



## Standard quantile mapping
QM = function(MMSDM,Variable,Step,Decimals){
  
  ## Initiating a storage list
  StorageList = list()
  StorageList[['Variable']] = Variable
  StorageList[['Step']] = Step
  StorageList[['Decimals']] = Decimals
  
  ## Creating step vector & new variable name
  Step_Vec = seq(0,1,Step)
  NewVar = paste(Variable,'.Corrected.QM',sep='')
  
  ## Creating bias dataframe
  Bias_Dat_QM = cbind.data.frame(Step_Vec)
  Bias_Dat_QM$Obs.Percentile = quantile(MMSDM$A.mat[,Variable],probs=Step_Vec)
  Bias_Dat_QM$Fit.Percentile = quantile(pull(MMSDM$A.hat.2.full, Variable),probs=Step_Vec)
  Bias_Dat_QM$Bias = Bias_Dat_QM$Obs.Percentile - Bias_Dat_QM$Fit.Percentile
  Bias_Dat = Bias_Dat_QM
  
  ## Storing the bias dataset
  StorageList[['Bias_Dat_QM']] = Bias_Dat
  
  
  
  ## LQM bias correction conducted here
  
  ## Fitted
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>%
    mutate(!!paste0('Quantile.',Variable,'.QM') := 
             round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
  MMSDM$A.hat.2.full = MMSDM$A.hat.2.full %>%
    mutate(!!paste0('Bias_Factor.',Variable,'.QM') := 
             Bias_Correction.QM(!!sym(paste0('Quantile.',Variable,'.QM')),
                                Bias_Dat,Step,Decimals))
  MMSDM$A.hat.2.full[,NewVar] = MMSDM$A.hat.2.full[,Variable] + 
    MMSDM$A.hat.2.full[,paste0('Bias_Factor.',Variable,'.QM')]
  
  ## RCP4.5
  if('RCP45.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      mutate(!!paste0('Quantile.',Variable,'.QM') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    MMSDM$RCP45.Downscaled = MMSDM$RCP45.Downscaled %>%
      mutate(!!paste0('Bias_Factor.',Variable,'.QM') :=
               Bias_Correction.QM(!!sym(paste0('Quantile.',Variable,'.QM')),
                                  Bias_Dat,Step,Decimals))
    MMSDM$RCP45.Downscaled[,NewVar] = MMSDM$RCP45.Downscaled[,Variable] + 
      MMSDM$RCP45.Downscaled[,paste0('Bias_Factor.',Variable,'.QM')]
  }
  
  ## RCP8.5
  if('RCP85.Downscaled' %in% names(MMSDM)){
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      mutate(!!paste0('Quantile.',Variable,'.QM') :=
               round(rank(!!sym(Variable))/length(!!sym(Variable)),digits=Decimals))
    MMSDM$RCP85.Downscaled = MMSDM$RCP85.Downscaled %>%
      mutate(!!paste0('Bias_Factor',Variable,'.QM') := 
               Bias_Correction.QM(!!sym(paste0('Quantile.',Variable,'.QM')),
                                  Bias_Dat,Step,Decimals))
    MMSDM$RCP85.Downscaled[,NewVar] = MMSDM$RCP85.Downscaled[,Variable] + 
      MMSDM$RCP85.Downscaled[,paste0('Bias_Factor',Variable,'.QM')]
  }
  
  ## Adding to storage list
  StorageList[['MMSDM']] = MMSDM
  #StorageList[[]] = A.hat.2.full
  #StorageList[[6]] = RCP45.Downscaled
  #StorageList[[7]] = RCP85.Downscaled
  
  ## Returning data list
  return(StorageList)
}



set.seed(789)
Span = .01
Step = 0.01
MMSDM1 = MMSDM_Train(A.mat,X.mat,RCP45.X.mat,RCP85.X.mat,ColNames)
LQM1 = LOESS_QM(MMSDM1,'R',Step,3,Span=Span)
LQM2 = QM(LQM1$MMSDM,'R',Step,3)


## Writing to csv files
Observed = LQM2$MMSDM$A.mat
Fitted = LQM2$MMSDM$A.hat.2.full
RCP45 = LQM2$MMSDM$RCP45.Downscaled
RCP85 = LQM2$MMSDM$RCP85.Downscaled

## Adding bias corrected precipitation
Fitted$PRCP.LQM = Fitted$Poc*(exp(Fitted$R.Corrected) - 1)
Fitted$PRCP.LQM = ifelse(Fitted$PRCP.LQM<0,0,Fitted$PRCP.LQM)
Fitted$PRCP.QM = Fitted$Poc*(exp(Fitted$R.Corrected.QM) - 1)
Fitted$PRCP.QM = ifelse(Fitted$PRCP.QM<0,0,Fitted$PRCP.QM)


## Plotting

### MEAN PRECIPITATION BY DAY ###

## Historical
Observed_Sum = Observed %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed'

## Fitted
Fitted$PRCP = ifelse(Fitted$PRCP<0,0,Fitted$PRCP)
Fitted_Sum = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP))
Fitted_Sum$Type = 'Fitted'

## LQM fitted
LQM_Fitted = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP.LQM))
LQM_Fitted$Type = 'LQM Corrected'

## QM corrected
QM_Fitted = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP.QM))
QM_Fitted$Type = 'QM Corrected'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,Fitted_Sum,LQM_Fitted,QM_Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
ggplot(Binded_Data,aes(x=Day,y=MeanPRCP,color=Type)) +
  ylab('Precipitation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line(lwd=0.5) + geom_point(size=1) + ylim(0,13) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14)) + 
  theme(plot.title = element_text(size=14)) +
  coord_cartesian(xlim=c(91,152)) + 
  scale_x_continuous(breaks=c(91,121,152),
                     labels=c('April 1',"May 1","June 1")) 
  #ggtitle(paste0('Smoothing parameter: ',Span,'; Step: ',Step))
  #ggtitle(paste0('Smoothing parameter: 0.10; Step: ',Step)) 

Plot1

ggarrange(Plot1,Plot2,Plot3,Plot4,Plot5,Plot6,
          common.legend=T,legend='bottom')

### MEAN PRECIPITATION BY DAY ###

##################### LQM VS QM BIAS CORRECTION #####################

