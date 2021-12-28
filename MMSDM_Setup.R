## Data formatting file




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










##################### READING IN DATA #####################

## Historical data
Historical = read.csv('Final.PR.Data.csv')

## Specifying model
#Model_Name = 'CCCMA-CanESM2'
Model_Name = 'GFDL-ESM2G'
#Model_Name = 'HadGEM2-ES'
#Model_Name = 'MRI-ESM1'

## GCM Data
GCM_Hist = read.csv(paste('Data - GCM/Processed Reduced Data/',Model_Name,'/Historical.csv',sep=''))
GCM_RCP45 = read.csv(paste('Data - GCM/Processed Reduced Data/',Model_Name,'/rcp45.csv',sep=''))
GCM_RCP85 = read.csv(paste('Data - GCM/Processed Reduced Data/',Model_Name,'/rcp85.csv',sep=''))

##################### READING IN DATA #####################















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
Bias_Correction = function(Day,Vec,Bias_Dat){
  ## Variables to be used later
  Bias_Dat$Blank = 900000000
  Dec3 = (as.integer(Vec*1000) %% 10)/10
  Column = match(paste('Bias_',floor(Vec*100)/100,sep=''),names(Bias_Dat))
  ## Bias Magnitude
  Bias_Magnitude = ifelse(Column==ncol(Bias_Dat),
                          Bias_Dat[Day,ncol(Bias_Dat)],NA)
  Bias_Magnitude = ifelse(Column<ncol(Bias_Dat),
                          Bias_Dat[cbind(Day,Column)]+Dec3*(Bias_Dat[cbind(Day,(Column+1))]-Bias_Dat[cbind(Day,Column)]),
                          Bias_Magnitude)
  return(Bias_Magnitude)
}

##################### FUNCTIONS NEEDED LATER #####################
















##################### FORMATTING DATA #####################

## Historical correlations
Historical_Reduced = Historical[,c(1,6,13,14,16,20)]

## Formatting historical
Historical_Reduced$Poc = ifelse(Historical_Reduced$PRCP>1,1,0)
Historical_Reduced$R = log(Historical_Reduced$PRCP+1)
Historical_Reduced$Day = yday(as.Date(Historical_Reduced$DATE,format='%m/%d/%Y'))
Historical_Reduced$Day = ifelse((leap_year(as.Date(Historical_Reduced$DATE,format='%m/%d/%Y')) & 
                                Historical_Reduced$Day>59),
                                (Historical_Reduced$Day-1),Historical_Reduced$Day)

## Removing leap day and trimming
Historical_Reduced = Historical_Reduced[as.Date(Historical_Reduced$DATE,format='%m/%d/%Y')>'1979-12-31' & 
                               as.Date(Historical_Reduced$DATE,format='%m/%d/%Y')<'2006-01-01' & 
                               format(as.Date(Historical_Reduced$DATE,format='%m/%d/%Y'),
                                      '%m%d')!='0229',]

## Just for the Hadley Center mode
#Historical_Reduced = Historical_Reduced[
#  !(as.numeric(format(as.Date(Historical_Reduced$DATE,format='%m/%d/%Y'),'%m'))>3 & 
#  as.numeric(format(as.Date(Historical_Reduced$DATE,format='%m/%d/%Y'),'%d'))>30),
#  ]


## Scaling temperatures
#Historical_Reduced = Historical_Reduced %>%
#  group_by(STATION,Day) %>%
#  mutate(TMAX.SD=sd(TMAX),TMIN.SD=sd(TMIN),
#         TMAX.Mean=mean(TMAX),TMIN.Mean=mean(TMIN),
#         TMAX.Scale=scale(TMAX),TMIN.Scale=scale(TMIN))

## Renaming columns
names(Historical_Reduced)[3:6] = c('Tmax','Tmin','PRCP','Srad')

## Standardized Hist_Spread
#Hist_Spread = myspread(Historical_Reduced[,c(1:2,14:15,6:8)],
#                       key=STATION,value=c(TMAX.Scale,TMIN.Scale,Poc,R,Srad))
#Hist_A.mat.Full = Hist_Spread[,c(1,rep(c(5:6,2:4),23)+rep(5*c(0:22),each=5))]
#Hist_A.mat = Hist_A.mat.Full[as.Date(Hist_A.mat.Full$DATE,format='%m/%d/%Y')>'1979-12-31' & 
#                             as.Date(Hist_A.mat.Full$DATE,format='%m/%d/%Y')<'2006-01-01' & 
#                             format(as.Date(Hist_A.mat.Full$DATE,format='%m/%d/%Y'),
#                                    '%m%d')!='0229',]
#Hist_A.mat$DATE = as.Date(Hist_A.mat$DATE,format='%m/%d/%Y')
#Hist_A.mat = Hist_A.mat[order(Hist_A.mat$DATE),]

## Unstandardized Hist_Spread
Hist_Spread_Unst = myspread(Historical_Reduced[,c(1:4,6:8)],
                            key=STATION,value=c(Tmax,Tmin,Poc,R,Srad))
Hist_A.mat.Full_Unst = Hist_Spread_Unst[,c(1,rep(c(5:6,2:4),23)+rep(5*c(0:22),each=5))]
Hist_A.mat_Unst = Hist_A.mat.Full_Unst[as.Date(Hist_A.mat.Full_Unst$DATE,format='%m/%d/%Y')>'1979-12-31' & 
                               as.Date(Hist_A.mat.Full_Unst$DATE,format='%m/%d/%Y')<'2006-01-01' & 
                               format(as.Date(Hist_A.mat.Full_Unst$DATE,format='%m/%d/%Y'),
                                      '%m%d')!='0229',]
Hist_A.mat_Unst$DATE = as.Date(Hist_A.mat_Unst$DATE,format='%m/%d/%Y')
Hist_A.mat_Unst = Hist_A.mat_Unst[order(Hist_A.mat_Unst$DATE),]


## Extracting the standardization key data
#Stand_Key = unique(Historical_Reduced[,c('STATION','Day','TMAX.SD','TMIN.SD','TMAX.Mean','TMIN.Mean')])






## Formatting historical GCM data
GCM_Hist$STATION = ifelse(GCM_Hist$Longitude<(-66),'STATION 1','STATION 3')
GCM_Hist$STATION = ifelse(GCM_Hist$Longitude<(-65.5) & GCM_Hist$Longitude>(-67),
                          'STATION 2',GCM_Hist$STATION)
Hist_X.mat = GCM_Hist[,c(1,9,6,7,4:5)]
Hist_X.mat$tasmax = Hist_X.mat$tasmax - 273.15
Hist_X.mat$tasmin = Hist_X.mat$tasmin - 273.15
Hist_X.mat = myspread(Hist_X.mat,key=STATION,value=c(tasmax,tasmin,pr,rsds))
Hist_X.mat = Hist_X.mat[as.Date(Hist_X.mat$Date)>'1979-12-31',]
Hist_X.mat$Month = as.numeric(format(as.Date(Hist_X.mat$Date),'%m'))
Hist_X.mat = Hist_X.mat[format(as.Date(Hist_X.mat$Date),'%m%d')!='0229',]
Hist_X.mat$Day = yday(as.Date(Hist_X.mat$Date))
Hist_X.mat$Day = ifelse((leap_year(as.Date(Hist_X.mat$Date)) & 
                          Hist_X.mat$Day>59),(Hist_X.mat$Day-1),Hist_X.mat$Day)
#Hist_X.mat = Hist_X.mat %>% group_by(Day) %>%
#  mutate(TMAX.SD=sd(tasmax),TMIN.SD=sd(tasmin),
#         TMAX.Mean=mean(tasmax),TMIN.Mean=mean(tasmin),
#         TMAX.Scale=scale(tasmax),TMIN.Scale=scale(tasmin))



## Future GCM data rcp45
GCM_RCP45$STATION = ifelse(GCM_RCP45$Latitude<18,'STATION 1','STATION 2')
GCM_RCP45_X.mat = GCM_RCP45[,c(1,9,6,7,4:5)]
GCM_RCP45_X.mat$tasmax = GCM_RCP45_X.mat$tasmax - 273.15
GCM_RCP45_X.mat$tasmin = GCM_RCP45_X.mat$tasmin - 273.15
GCM_RCP45_X.mat = myspread(GCM_RCP45_X.mat,key=STATION,value=c(tasmax,tasmin,pr,rsds))
#GCM_RCP45_X.mat$Date = as.Date(GCM_RCP45_X.mat$Date,format='%m/%d/%Y')
GCM_RCP45_X.mat$Month = as.numeric(format(as.Date(GCM_RCP45_X.mat$Date),'%m'))
GCM_RCP45_X.mat = GCM_RCP45_X.mat[format(as.Date(GCM_RCP45_X.mat$Date),'%m%d')!='0229',]
GCM_RCP45_X.mat$Day = yday(as.Date(GCM_RCP45_X.mat$Date))
GCM_RCP45_X.mat$Day = ifelse((leap_year(as.Date(GCM_RCP45_X.mat$Date)) & 
                           GCM_RCP45_X.mat$Day>59),(GCM_RCP45_X.mat$Day-1),GCM_RCP45_X.mat$Day)

#GCM_RCP45_X.mat = GCM_RCP45_X.mat %>% group_by(Day) %>%
#  mutate(TMAX.SD=sd(tasmax),TMIN.SD=sd(tasmin),
#         TMAX.Mean=mean(tasmax),TMIN.Mean=mean(tasmin),
#         TMAX.Scale=scale(tasmax),TMIN.Scale=scale(tasmin))



## Future GCM data rcp85
#GCM_RCP85$STATION = ifelse(GCM_RCP85$Longitude<(-66),'STATION 1','STATION 3')
#GCM_RCP85$STATION = ifelse(GCM_RCP85$Longitude<(-65.5) & GCM_RCP85$Longitude>(-67),
#                          'STATION 2',GCM_Hist$STATION)
GCM_RCP85$STATION = ifelse(GCM_RCP85$Latitude<18,'STATION 1','STATION 2')
GCM_RCP85_X.mat = GCM_RCP85[,c(1,9,6,7,4:5)]
GCM_RCP85_X.mat$tasmax = GCM_RCP85_X.mat$tasmax - 273.15
GCM_RCP85_X.mat$tasmin = GCM_RCP85_X.mat$tasmin - 273.15
GCM_RCP85_X.mat = myspread(GCM_RCP85_X.mat,key=STATION,value=c(tasmax,tasmin,pr,rsds))
#GCM_RCP85_X.mat$Date = as.Date(GCM_RCP85_X.mat$Date,format='%m/%d/%Y')
GCM_RCP85_X.mat$Month = as.numeric(format(as.Date(GCM_RCP85_X.mat$Date),'%m'))
GCM_RCP85_X.mat = GCM_RCP85_X.mat[format(as.Date(GCM_RCP85_X.mat$Date),'%m%d')!='0229',]
GCM_RCP85_X.mat$Day = yday(as.Date(GCM_RCP85_X.mat$Date))
GCM_RCP85_X.mat$Day = ifelse((leap_year(as.Date(GCM_RCP85_X.mat$Date)) & 
                           GCM_RCP85_X.mat$Day>59),(GCM_RCP85_X.mat$Day-1),GCM_RCP85_X.mat$Day)

#GCM_RCP85_X.mat = GCM_RCP85_X.mat %>% group_by(Day) %>%
#  mutate(TMAX.SD=sd(tasmax),TMIN.SD=sd(tasmin),
#         TMAX.Mean=mean(tasmax),TMIN.Mean=mean(tasmin),
#         TMAX.Scale=scale(tasmax),TMIN.Scale=scale(tasmin))



## Writing to CSVs
#write.csv(Hist_A.mat,paste('Data - MMSDM/',Model_Name,'/Hist_A.mat.csv',sep=''),row.names=F)
#write.csv(Hist_A.mat_Unst,paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''),row.names=F)
#write.csv(Hist_X.mat,paste('Data - MMSDM/',Model_Name,'/Hist_X.mat.csv',sep=''),row.names=F)
#write.csv(GCM_RCP45_X.mat,paste('Data - MMSDM/',Model_Name,'/RCP45_X.mat.csv',sep=''),row.names=F)
#write.csv(GCM_RCP85_X.mat,paste('Data - MMSDM/',Model_Name,'/RCP85_X.mat.csv',sep=''),row.names=F)
#write.csv(Stand_Key,paste('Data - MMSDM/',Model_Name,'/Standardization_Key.csv',sep=''),row.names=F)


##################### FORMATTING DATA #####################

