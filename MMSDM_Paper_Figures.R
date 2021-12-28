## MMSDM Organized Paper Figures





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





#PRCP_Bias2 = t(PRCP_Bias)
#write.csv(PRCP_Bias2,paste0('Data - MMSDM/',Model_Name,'/Future Climate/PRCP_Bias2.csv'),row.names=T)
#max(PRCP_Bias[,-1])
#min(PRCP_Bias[,-1])



##################### READING IN THE DATA #####################

## Specifying model
Model_Name = 'CCCMA-CanESM2'
Model_Name = 'MRI-ESM1'
Model_Name = 'HadGEM2-ES'
Model_Name = 'GFDL-ESM2G'


## Reading in the data
Observed_Full = read.csv('Final.PR.Data.csv')
#PRCP_Bias = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/PRCP_Bias.csv',sep=''))
#SRAD_Bias = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/SRAD_Bias.csv',sep=''))
Observed = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/Observed.csv',sep=''))
Fitted = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/Fitted.csv',sep=''))
RCP45 = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/RCP45.csv',sep=''))
RCP85 = read.csv(paste('Data - MMSDM/',Model_Name,'/Future Climate/RCP85.csv',sep=''))

## Adding column variable for year
Observed_Full$Year = year(as.Date(Observed_Full$DATE,format='%m/%d/%Y'))
Observed$Year = year(Observed$DATE)
RCP45$Year = year(RCP45$DATE)
RCP85$Year = year(RCP85$DATE)

## Adding column variable for month
Observed_Full$Month = month(as.Date(Observed_Full$DATE,format='%m/%d/%Y'))
Observed$Month = month(Observed$DATE)
RCP45$Month = month(RCP45$DATE)
RCP85$Month = month(RCP85$DATE)

## Adding column variable for season
Observed_Full$Season = ifelse(Observed_Full$Month<4 | Observed_Full$Month>11,'Dry Season','Late Rainfall Season')
Observed_Full$Season = ifelse(Observed_Full$Month<8 & Observed_Full$Month>3,'Early Rainfall Season',Observed_Full$Season)

Observed$Season = ifelse(Observed$Month<4 | Observed$Month>11,'Dry Season','Late Rainfall Season')
Observed$Season = ifelse(Observed$Month<8 & Observed$Month>3,'Early Rainfall Season',Observed$Season)

RCP45$Season = ifelse(RCP45$Month<4 | RCP45$Month>11,'Dry Season','Late Rainfall Season')
RCP45$Season = ifelse(RCP45$Month<8 & RCP45$Month>3,'Early Rainfall Season',RCP45$Season)

RCP85$Season = ifelse(RCP85$Month<4 | RCP85$Month>11,'Dry Season','Late Rainfall Season')
RCP85$Season = ifelse(RCP85$Month<8 & RCP85$Month>3,'Early Rainfall Season',RCP85$Season)

##################### READING IN THE DATA #####################



## Selecting a specific station
Fitted.pde = Fitted %>%
  filter(STATION == 'RQC00666992')
Observed_Full.pde = Observed_Full %>%
  filter(STATION == 'RQC00666992')
Observed.pde = Observed %>%
  filter(STATION == 'RQC00666992')
RCP45.pde = RCP45 %>%
  filter(STATION == 'RQC00666992')
RCP85.pde = RCP85 %>%
  filter(STATION == 'RQC00666992')












##################### SECTION 4.3 TEMPERATURE PLOT #####################

## Setting working directory
setwd("D:/Dissertation Stuff/Dissertation Work")

## Reading in data
Observed = read.csv('Final.PR.Data.csv')
Observed$Day2 = yday(as.Date(Observed$DATE,format='%m/%d/%Y'))
Observed$Day2 = ifelse((leap_year(as.Date(Observed$DATE,format='%m/%d/%Y')) & 
                          Observed$Day2>59),(Observed$Day2-1),Observed$Day2)
View(Observed)

## Downscaled data maximums
Day.Means = Observed %>%
  group_by(Year,Day2) %>%
  summarise(MeanTemp = mean(TMIN))
Day.Means$Year2 = ifelse(Day.Means$Year!=2017,'Other','2017')
Day.Means$Width = ifelse(Day.Means$Year!=2017,1,1.5)
Day.Means$Year2 = factor(Day.Means$Year2,levels=unique(Day.Means$Year2))


ggplot(Day.Means,aes(x=Day2,y=MeanTemp,color=Year2)) +
  ylab('Daily Minimum Temperature (C)') + 
  xlab('Julian Day') + labs(color="Year") + 
  geom_line(lwd=1.25) + 
  theme(legend.text=element_text(size=16)) +
  theme(text=element_text(size=16)) +
  guides(color=guide_legend(override.aes=list(size=1.5)))


##################### SECTION 4.3 TEMPERATURE PLOT #####################
























##################### SECTION 5.3 TEMPERATURE BIAS CORRECTION #####################


### STANDARDIZE PLOT 4 PANELS ###

## IMPORTANT NOTE HERE.
## In order to do the standardized 4 panel plots, you need to go back
## to the setup file and add back in the scaled variables.

## Specifying model
Model_Name = 'CCCMA-CanESM2'

## Reading in the data
#A.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat.csv',sep=''))
A.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''))
X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_X.mat.csv',sep=''))
RCP45.X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/RCP45_X.mat.csv',sep=''))
RCP85.X.mat = read.csv(paste('Data - MMSDM/',Model_Name,'/RCP85_X.mat.csv',sep=''))
Stand_Key = read.csv(paste('Data - MMSDM/',Model_Name,'/Standardization_Key.csv',sep=''))
#ColNames = c('TMAX.Scale','TMIN.Scale','pr','rsds')
ColNames = c('tasmax','tasmin','pr','rsds')

## Parameter estimation
Y = as.matrix(A.mat[,2:ncol(A.mat)])
X = as.matrix(X.mat[,ColNames])
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
A.hat.2 = data.frame(A.hat + mvrnorm(n=nrow(A.hat),mu=rep(0,ncol(A.hat)),Sigma=Sigma.mat))
A.hat.2$DATE = A.mat$DATE
A.hat.2 = A.hat.2[,c(ncol(A.hat.2),1:(ncol(A.hat.2)-1))]

## Gathering and back-standardizing
A.hat.2.gather = mygather(A.hat.2)
A.hat.2.gather$Day = yday(as.Date(A.hat.2.gather$DATE))
A.hat.2.gather$Day = ifelse((leap_year(as.Date(A.hat.2.gather$DATE)) & 
                               A.hat.2.gather$Day>59),(A.hat.2.gather$Day-1),A.hat.2.gather$Day)
A.hat.2.full = left_join(A.hat.2.gather,Stand_Key,by=c('STATION','Day'))
#A.hat.2.full$TMAX = A.hat.2.full$TMAX.Scale*A.hat.2.full$TMAX.SD + A.hat.2.full$TMAX.Mean
#A.hat.2.full$TMIN = A.hat.2.full$TMIN.Scale*A.hat.2.full$TMIN.SD + A.hat.2.full$TMIN.Mean
A.hat.2.full$TMAX = A.hat.2.full$Tmax
A.hat.2.full$TMIN = A.hat.2.full$Tmin

## Changing to matrices
RCP45.X = as.matrix(RCP45.X.mat[,ColNames])
RCP85.X = as.matrix(RCP85.X.mat[,ColNames])

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
RCP45.full = left_join(RCP45.gather,Stand_Key,by=c('STATION','Day'))
#RCP45.full$TMAX = RCP45.full$TMAX.Scale*RCP45.full$TMAX.SD + RCP45.full$TMAX.Mean
#RCP45.full$TMIN = RCP45.full$TMIN.Scale*RCP45.full$TMIN.SD + RCP45.full$TMIN.Mean
RCP45.full$TMAX = RCP45.full$Tmax
RCP45.full$TMIN = RCP45.full$Tmin

## Back standardizing future GCM data RCP8.5
RCP85.gather$Day = yday(as.Date(RCP85.gather$DATE))
RCP85.gather$Day = ifelse((leap_year(as.Date(RCP85.gather$DATE)) & 
                             RCP85.gather$Day>59),(RCP85.gather$Day-1),RCP85.gather$Day)
RCP85.full = left_join(RCP85.gather,Stand_Key,by=c('STATION','Day'))
#RCP85.full$TMAX = RCP85.full$TMAX.Scale*RCP85.full$TMAX.SD + RCP85.full$TMAX.Mean
#RCP85.full$TMIN = RCP85.full$TMIN.Scale*RCP85.full$TMIN.SD + RCP85.full$TMIN.Mean
RCP85.full$TMAX = RCP85.full$Tmax
RCP85.full$TMIN = RCP85.full$Tmin

## Downscaled data maximum
Day.Means45 = RCP45.full %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(TMAX))
Day.Means45$Type = 'Downscaled RCP4.5'

Day.Means85 = RCP85.full %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(TMAX))
Day.Means85$Type = 'Downscaled RCP8.5'

## Obseved maximums
#A.mat.unst = A.mat
A.mat.unst = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''))
A.mat.unst.gather = mygather(A.mat.unst)
A.mat.unst.gather$Day = yday(as.Date(A.mat.unst.gather$DATE))
A.mat.unst.gather$Day = ifelse((leap_year(as.Date(A.mat.unst.gather$DATE)) & 
                                  A.mat.unst.gather$Day>59),(A.mat.unst.gather$Day-1),A.mat.unst.gather$Day)
Day.MeansHist = A.mat.unst.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(Tmax))
Day.MeansHist$Type = 'Observed Historical'

## GCM historical maximum
X.mat.gather = X.mat[,1:3]
X.mat.gather$Day = yday(as.Date(X.mat.gather$Date))
X.mat.gather$Day = ifelse((leap_year(as.Date(X.mat.gather$Date)) & 
                             X.mat.gather$Day>59),(X.mat.gather$Day-1),X.mat.gather$Day)
Day.MeansHistGCM = X.mat.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmax))
Day.MeansHistGCM$Type = 'GCM Historical'

## Historical GCM efforts
X.mat.45gather = RCP45.X.mat[,1:3]
X.mat.45gather$Day = yday(as.Date(X.mat.45gather$Date))
X.mat.45gather$Day = ifelse((leap_year(as.Date(X.mat.45gather$Date)) & 
                               X.mat.45gather$Day>59),(X.mat.45gather$Day-1),X.mat.45gather$Day)
Day.Means45GCM = X.mat.45gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmax))
Day.Means45GCM$Type = 'GCM RCP4.5'

X.mat.85gather = RCP85.X.mat[,1:3]
X.mat.85gather$Day = yday(as.Date(X.mat.85gather$Date))
X.mat.85gather$Day = ifelse((leap_year(as.Date(X.mat.85gather$Date)) & 
                               X.mat.85gather$Day>59),(X.mat.85gather$Day-1),X.mat.85gather$Day)
Day.Means85GCM = X.mat.85gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmax))
Day.Means85GCM$Type = 'GCM RCP8.5'



## Binding the data
Binded_Data = rbind.data.frame(Day.MeansHist,Day.MeansHistGCM,
                               Day.Means45GCM,Day.Means85GCM,
                               Day.Means45,Day.Means85)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

Plot4 = ggplot(Binded_Data,aes(x=Day,y=MeanTemp,color=Type)) +
  geom_point(cex=.55) + ylab('Mean Temperature') + 
  xlab('Juian Day') + labs(color="Data Source") + 
  ylim(20.3,33.5) + 
  guides(colour=guide_legend(override.aes=list(size=2))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))


ggarrange(Plot1,Plot2,Plot3,Plot4,ncol=2,nrow=2,
          common.legend=T,legend="bottom")

### STANDARDIZE PLOT 4 PANELS ###





### MEAN TEMP BY YEAR ###

## Downscaled data maximums
RCP45.full$Year = year(as.Date(RCP45.full$DATE))
Day.Means45 = RCP45.full %>%
  group_by(Year) %>%
  summarise(MeanTemp = mean(TMAX))
Day.Means45$Type = 'Downscaled RCP4.5'

RCP85.full$Year = year(as.Date(RCP85.full$DATE))
Day.Means85 = RCP85.full %>%
  group_by(Year) %>%
  summarise(MeanTemp = mean(TMAX))
Day.Means85$Type = 'Downscaled RCP8.5'

## Binding the data
Binded_Data = rbind.data.frame(Day.Means45,Day.Means85)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

ggplot(Binded_Data,aes(x=Year,y=MeanTemp,color=Type)) +
  geom_line() + geom_point() + ylab('Mean Maximum Temperature (C)') + 
  xlab('Year') + labs(color="Data Source") + 
  guides(colour=guide_legend(override.aes=list(size=2,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

### MEAN TEMP BY YEAR ###





### TEMPERATURE MAX/MIN BIAS PLOTS ###

## Obseved maximums and minimums
A.mat.unst = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''))
A.mat.unst.gather = mygather(A.mat.unst)
A.mat.unst.gather$Day = yday(as.Date(A.mat.unst.gather$DATE))
A.mat.unst.gather$Day = ifelse((leap_year(as.Date(A.mat.unst.gather$DATE)) & 
                                  A.mat.unst.gather$Day>59),(A.mat.unst.gather$Day-1),A.mat.unst.gather$Day)
Day.MeansHist = A.mat.unst.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(Tmin))
Day.MeansHist$Type = 'Observed Historical'

## GCM historical maximum/minimum
X.mat.gather = X.mat[,1:3]
X.mat.gather$Day = yday(as.Date(X.mat.gather$Date))
X.mat.gather$Day = ifelse((leap_year(as.Date(X.mat.gather$Date)) & 
                             X.mat.gather$Day>59),(X.mat.gather$Day-1),
                          X.mat.gather$Day)
Day.MeansHistGCM = X.mat.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmin))
Day.MeansHistGCM$Type = 'GCM Historical'

## Fitted maximums
Fitted = A.hat.2.full %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(TMIN))
Fitted$Type = 'Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Day.MeansHist,Day.MeansHistGCM,Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plots
Plot3 = ggplot(Binded_Data,aes(x=Day,y=MeanTemp,color=Type)) +
  ylab('Minimum Temperature (C)') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + geom_point(cex=1) + ylim(17,32) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

## LOESS curve
GCM_PRCP_Bias.lo = loess((Day.MeansHist$MeanTemp-Fitted$MeanTemp) ~ 
                           Day.MeansHist$Day,span=.5)
DiffDat = cbind.data.frame(Day=1:365,Difference=(Day.MeansHist$MeanTemp-Fitted$MeanTemp),
                           LOESS=predict(GCM_PRCP_Bias.lo))

Plot4 = ggplot(DiffDat,aes(x=Day,y=Difference)) +
  geom_point(cex=1) + ylab('Bias Magnitude (C)') + 
  labs(color="Data Source") + ylim(-2,2) + 
  geom_line(aes(x=Day,y=LOESS),lwd=1.25,color='red') + 
  theme(text=element_text(size=14)) +
  xlab('Julian Day')

ggarrange(Plot1,Plot3,Plot2,Plot4,ncol=2,nrow=2,
          common.legend=T,legend="bottom")

### TEMPERATURE MAX/MIN BIAS PLOTS ###





### QUARTILE TEMP PLOTS TO SHOW LQM NOT NEEDED ###

## Obseved maximums and minimums
A.mat.unst = read.csv(paste('Data - MMSDM/',Model_Name,'/Hist_A.mat_Unst.csv',sep=''))
A.mat.unst.gather = mygather(A.mat.unst)
A.mat.unst.gather$Day = yday(as.Date(A.mat.unst.gather$DATE))
A.mat.unst.gather$Day = ifelse((leap_year(as.Date(A.mat.unst.gather$DATE)) & 
                                  A.mat.unst.gather$Day>59),(A.mat.unst.gather$Day-1),A.mat.unst.gather$Day)
Day.MeansHist = A.mat.unst.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = quantile(Tmax,.25))
Day.MeansHist$Type = 'Observed Historical'

Day.MeansHist.min = A.mat.unst.gather %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(Tmin))
Day.MeansHist.min$Type = 'Observed Historical'

## Fitted maximums
Fitted = A.hat.2.full %>%
  group_by(Day) %>%
  summarise(MeanTemp = quantile(TMAX,.75))
Fitted$Type = 'Fitted Values'

## Fitted minimums
Fitted.min = A.hat.2.full %>%
  group_by(Day) %>%
  summarise(MeanTemp = quantile(TMIN,.75))
Fitted.min$Type = 'Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Day.MeansHist,Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
ggplot(Binded_Data,aes(x=Day,y=MeanTemp,color=Type)) +
  ylab('Log Precipitation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + 
  ggtitle('Median') + geom_point() + 
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

## LOESS curve
GCM_PRCP_Bias.lo = loess((Day.MeansHist$MeanTemp-Fitted$MeanTemp) ~ 
                           Day.MeansHist$Day,span=.5)
DiffDat = cbind.data.frame(Day=1:365,Difference=(Day.MeansHist$MeanTemp-Fitted$MeanTemp),
                           LOESS=predict(GCM_PRCP_Bias.lo))

ggplot(DiffDat,aes(x=Day,y=Difference)) +
  geom_point(cex=1) + ylab('Bias Magnitude') + 
  xlab('Julian Day') + labs(color="Data Source") +
  geom_line(aes(x=Day,y=LOESS),lwd=1.5,color='red') + 
  theme(text=element_text(size=14)) +
  ggtitle('First Quartile')

### QUARTILE TEMP PLOTS TO SHOW LQM NOT NEEDED ###





### BIAS CORRECTION PLOTS ###

## Specifying downscaled data
set.seed(318)
MMSDM1 = MMSDM_Train(A.mat,X.mat,RCP45.X.mat,RCP85.X.mat,ColNames)
LQM1 = LOESS_QM(MMSDM1,'R',0.01,3)

## Observed historical
Day.MeansHist = LQM1$MMSDM$A.mat %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(Tmin))
Day.MeansHist$Type = 'Observed Historical'

## GCM historical RCP45
GCM.MeansHist = LQM1$MMSDM$X.mat %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmin))
GCM.MeansHist$Type = 'GCM Historical'

## GCM RCP45
GCM.45 = LQM1$MMSDM$X.mat.RCP45 %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmin))
GCM.45$Type = 'GCM RCP4.5'

## GCM RCP8.5
GCM.85 = LQM1$MMSDM$X.mat.RCP85 %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(tasmin))
GCM.85$Type = 'GCM RCP8.5'

## Downscaled RCP45
RCP45.Downscale = LQM1$MMSDM$RCP45.Downscaled %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(TMIN.Corrected))
RCP45.Downscale$Type = 'Downscaled RCP4.5'

## Downscaled RCP85
RCP85.Downscale = LQM1$MMSDM$RCP85.Downscaled %>%
  group_by(Day) %>%
  summarise(MeanTemp = mean(TMIN.Corrected))
RCP85.Downscale$Type = 'Downscale RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Day.MeansHist,GCM.MeansHist,
                               GCM.45,GCM.85,RCP45.Downscale,
                               RCP85.Downscale)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
Plot2 = ggplot(Binded_Data,aes(x=Day,y=MeanTemp,color=Type)) +
  ylab('Minimum Temperature (C)') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + geom_point() + ylim(17.5,33) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

ggarrange(Plot1,Plot2,ncol=2,nrow=1,
          common.legend=T,legend="bottom")

### BIAS CORRECTION PLOTS ###


##################### SECTION 5.3 TEMPERATURE BIAS CORRECTION #####################










##################### SECTION 5.4.2 LOESS QUANTILE MAPPING #####################


### Log PRCP CDF Plots ###

## Overlaid CDFs of raw, fitted, RCP4.5, and RCP8.5 precipitations
Raw = cbind.data.frame(PRCP=LQM1$MMSDM$A.mat$R,Type='Observed')
Fitted = cbind.data.frame(PRCP=LQM1$MMSDM$A.hat.2.full$R,Type='Fitted')
RCP45 = cbind.data.frame(PRCP=LQM1$MMSDM$RCP45.Downscaled$R,Type='RCP4.5')
RCP85 = cbind.data.frame(PRCP=LQM1$MMSDM$RCP85.Downscaled$R,Type='RCP8.5')
PRCP_Combined = rbind.data.frame(Raw,Fitted,RCP45,RCP85)
PRCP_Combined$Type = factor(PRCP_Combined$Type,levels=unique(PRCP_Combined$Type))

ggplot(PRCP_Combined,aes(x=PRCP,color=Type)) +
  stat_ecdf(lwd=1.5) + 
  coord_cartesian(xlim=c(-1,5)) +
  xlab('Log Precipitation') +
  ylab('Probability') + 
  labs(color="Data Source") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(override.aes=list(linetype=2.5))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

### Log PRCP CDF Plots ###





### Monthly means of specific variable ###

## Historical
Observed2 = Observed %>%
  group_by(Day) %>%
  summarise(MeanPRCP = quantile(R,0)) %>%
  as.data.frame()
Observed2$Type = 'Observed Historical'
## Fitted
Fitted2 = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = quantile(R,0))
Fitted2$Type = 'Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Observed2,Fitted2)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

### Monthly means of specific variable ###





### Quartile bias plots ###

## LOESS curve
GCM_PRCP_Bias.lo = loess((Observed2$MeanPRCP-Fitted2$MeanPRCP) ~ 
                           Observed2$Day,span=.5)
DiffDat = cbind.data.frame(Day=1:365,Difference=(Observed2$MeanPRCP-Fitted2$MeanPRCP),
                           LOESS=predict(GCM_PRCP_Bias.lo))

## Plotting
Plot3 = ggplot(Binded_Data,aes(x=Day,y=MeanPRCP,color=Type)) +
  ylab('Log Precipitation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + #ylim(-0.5,3) +
  ggtitle('Median') + geom_point() + 
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

## First bias magnitude plot
Plot2 = ggplot(DiffDat,aes(x=Day,y=Difference)) +
  geom_point(cex=1) + ylab('Bias Magnitude') + 
  xlab('Julian Day') + labs(color="Data Source") +
  geom_line(aes(x=Day,y=LOESS),lwd=1.5,color='red') + 
  ylim(-1.3,1.3) + theme(text=element_text(size=14)) +
  ggtitle('First Quartile')

ggarrange(Plot1,Plot3,Plot5,Plot2,Plot4,Plot6,
          nrow=2,ncol=3,legend='bottom',
          common.legend=T)

### Quartile bias plots ###





### Bias graphical example ###

## LOESS curve
GCM_PRCP_Bias.lo = loess((Observed2$MeanPRCP-Fitted2$MeanPRCP) ~ 
                           Observed2$Day,span=.5)
DiffDat0 = cbind.data.frame(Day=1:365,Difference=(Observed2$MeanPRCP-Fitted2$MeanPRCP),
                           LOESS=predict(GCM_PRCP_Bias.lo),
                           Type='81th Percentile Bias Magnitude')


## Binding 2 datasets
Diff = rbind.data.frame(DiffDat,DiffDat0)

## First bias magnitude plot
ggplot(Diff,aes(x=Day,y=Difference,color=Type)) +
  geom_point(cex=1,aes(color=Type)) + ylab('Bias Magnitude') + 
  xlab('Julian Day') + labs(color="Data Source") +
  geom_line(aes(x=Day,y=LOESS,color=Type),lwd=1) + 
  ylim(-1.25,1.25) + theme(legend.position="bottom") +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))


## Horizontal lines
Day136 = Diff[Diff$Day==139,]
Distance = Day136$LOESS[2] - Day136$LOESS[1]
Y1 = Day136$LOESS[1] + 1/10*Distance
Y2 = Day136$LOESS[1] + 2/10*Distance
Y3 = Day136$LOESS[1] + 3/10*Distance
Y4 = Day136$LOESS[1] + 4/10*Distance
Y5 = Day136$LOESS[1] + 5/10*Distance
Y6 = Day136$LOESS[1] + 6/10*Distance
Y7 = Day136$LOESS[1] + 7/10*Distance
Y8 = Day136$LOESS[1] + 8/10*Distance
Y9 = Day136$LOESS[1] + 9/10*Distance

## Zoomed bias magnitude plot
ggplot(Diff,aes(x=Day,y=Difference,color=Type)) +
  geom_point(cex=1,aes(color=Type)) + ylab('Bias Magnitude') + 
  xlab('Julian Day') + labs(color="Data Source") +
  geom_line(aes(x=Day,y=LOESS,color=Type),lwd=1) + 
  ylim(0,0.35) + theme(legend.position="bottom") +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14)) + xlim(101,200) +
  geom_vline(xintercept=139,linetype=1,color='black',size=1) +
  geom_hline(yintercept=Y1,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y2,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y3,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y4,linetype=1,color='purple',size=0.5) + 
  geom_hline(yintercept=Y5,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y6,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y7,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y8,linetype=2,color='gray10',size=0.5) + 
  geom_hline(yintercept=Y9,linetype=2,color='gray10',size=0.5)

### Bias graphical example ###





### MEAN PRECIPITATION BY DAY ###

## Historical
Observed_Sum = Observed %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## Fitted
Fitted$PRCP = ifelse(Fitted$PRCP<0,0,Fitted$PRCP)
Fitted_Sum = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP))
Fitted_Sum$Type = 'Fitted Values'

## Fitted
BC_Fitted = Fitted %>%
  group_by(Day) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected))
BC_Fitted$Type = 'Bias Corrected Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,Fitted_Sum,BC_Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
ggplot(Binded_Data,aes(x=Day,y=MeanPRCP,color=Type)) +
  ylab('Precipitation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + geom_point() + ylim(0,13.5) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

### MEAN PRECIPITATION BY DAY ###





### LQM CDF PLOT CHECK ###

## Creating dataframe
Rdat = cbind.data.frame(R = c(Observed$R,Fitted$R,Fitted$R.Corrected),
                        Type = c(rep('Observed R',218270),
                                 rep('Fitted R',218270),
                                 rep('Corrected R',218270)))
Rdat$Type = factor(Rdat$Type,levels=unique(Rdat$Type))


## Plotting
ggplot(Rdat,aes(x=R,color=Type)) +
  stat_ecdf(lwd=1.25) + 
  coord_cartesian(xlim=c(-2.5,5)) +
  xlab('Log Precipitation (mm)') +
  ylab('Probability') + 
  labs(color="Data Source") +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

### LQM CDF PLOT CHECK ###

##################### SECTION 5.4.2 LOESS QUANTILE MAPPING #####################










##################### SECTON 5.5: SRAD BIAS ASSESSMENT #####################


### Monthly means of specific variable ###

## Historical
Observed = LQM1$MMSDM$A.mat %>%
  group_by(Day) %>%
  summarise(MeanPRCP = quantile(Srad,.5)) %>%
  as.data.frame()
Observed$Type = 'Observed Historical'
## Fitted
Fitted = LQM1$MMSDM$A.hat.2.full %>%
  group_by(Day) %>%
  summarise(MeanPRCP = quantile(Srad,.5))
Fitted$Type = 'Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Observed,Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

### Monthly means of specific variable ###




### Quartile bias plots ###

## LOESS curve
GCM_PRCP_Bias.lo = loess((Observed$MeanPRCP-Fitted$MeanPRCP) ~ 
                           Observed$Day,span=.5)
DiffDat = cbind.data.frame(Day=1:365,Difference=(Observed$MeanPRCP-Fitted$MeanPRCP),
                           LOESS=predict(GCM_PRCP_Bias.lo))

## Plotting
Plot3 = ggplot(Binded_Data,aes(x=Day,y=MeanPRCP,color=Type)) +
  ylab('Solar Radiation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + ylim(225,500) +
  ggtitle('Median') + geom_point() + 
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

## First bias magnitude plot
Plot4 = ggplot(DiffDat,aes(x=Day,y=Difference)) +
  geom_point(cex=1) + ylab('Bias Magnitude') + 
  xlab('Julian Day') + labs(color="Data Source") +
  geom_line(aes(x=Day,y=LOESS),lwd=1.5,color='red') + 
  ylim(-85,85) + theme(text=element_text(size=14)) +
  ggtitle('Median')

ggarrange(Plot1,Plot3,Plot5,Plot2,Plot4,Plot6,
          nrow=2,ncol=3,legend='bottom',
          common.legend=T)

### Quartile bias plots ###




### EXAMINING FITTED VALUES PLOT ###

### Monthly means of specific variable ###

## Historical
Observed_Sum = Observed %>%
  group_by(Day) %>%
  summarise(MeanSRAD = mean(Srad)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## Fitted
Fitted_Sum = Fitted %>%
  group_by(Day) %>%
  summarise(MeanSRAD = mean(Srad))
Fitted_Sum$Type = 'Fitted Values'

## Fitted
BC_Fitted =Fitted %>%
  group_by(Day) %>%
  summarise(MeanSRAD = mean(Srad.Corrected))
BC_Fitted$Type = 'Bias Corrected Fitted Values'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,Fitted_Sum,BC_Fitted)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
ggplot(Binded_Data,aes(x=Day,y=MeanSRAD,color=Type)) +
  ylab('Solar Radiation') + 
  xlab('Julian Day') + labs(color="Data Source") + 
  geom_line() + geom_point() + ylim(275,450) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

### MEAN PRECIPITATION BY DAY ###

### EXAMINING FITTED VALUES PLOT ###


##################### SECTON 5.5: SRAD BIAS ASSESSMENT #####################


















##################### CHAPTER 6: FUTURE CLIMATE PLOTS TEMP/SRAD #####################


### TEMP PERCENTILE TABLE ###

quantile(Observed$Tmax,seq(0,1,.1))
quantile(RCP45$TMAX.Corrected,seq(0,1,.1))
quantile(RCP85$TMAX.Corrected,seq(0,1,.1))

quantile(Observed$Tmin,seq(0,1,.1))
quantile(RCP45$TMIN.Corrected,seq(0,1,.1))
quantile(RCP85$TMIN.Corrected,seq(0,1,.1))

### TEMP PERCENTILE TABLE ###





### SRAD PERCENTILE TABLE ###

quantile(Observed$Srad,seq(0,1,.1))
quantile(RCP45$Srad.Corrected,seq(0,1,.1))
quantile(RCP85$Srad.Corrected,seq(0,1,.1))

### SRAD PERCENTILE TABLE ###






### MEAN TEMPERATURE BY YEAR ###

## Historical
Observed_Sum2 = Observed_Full %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(TMIN)) %>%
  as.data.frame()
Observed_Sum2$Type = 'Observed Historical 2'

## Historical
Observed_Sum = Observed %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(Tmin)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## Historical
RCP45_Sum = RCP45 %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(TMIN.Corrected)) %>%
  as.data.frame()
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## Historical
RCP85_Sum = RCP85 %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(TMIN.Corrected)) %>%
  as.data.frame()
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,Observed_Sum2,RCP45_Sum,RCP85_Sum)
#Binded_Data = rbind.data.frame(Observed_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
Plot2 = ggplot(Binded_Data,aes(x=Year,y=MeanSRAD,color=Type)) +
  ylab('Minimum Temperature') + 
  xlab('Year') + labs(color="Data Source") + 
  geom_line() + geom_point() + #ylim(15.5,19) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14))

ggarrange(Plot1,Plot2,Plot3,Plot4,ncol=2,nrow=2,
          common.legend=T,legend='bottom')

### MEAN TEMPERATURE BY YEAR ###








### MEAN SRAD BY YEAR ###

## Historical
Observed_Sum = Observed %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(Srad)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## Historical
RCP45_Sum = RCP45 %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(Srad.Corrected)) %>%
  as.data.frame()
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## Historical
RCP85_Sum = RCP85 %>%
  group_by(Year) %>%
  summarise(MeanSRAD = mean(Srad.Corrected)) %>%
  as.data.frame()
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
Plot3 = ggplot(Binded_Data,aes(x=Year,y=MeanSRAD,color=Type)) +
  ylab('Solar Radiation') + 
  xlab('Year') + labs(color="Data Source") + 
  geom_line() + geom_point() + ylim(280,339) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=15)) +
  theme(text=element_text(size=15))

### MEAN SRAD BY YEAR ###








### SRAD VS TMAX PLOTS ###

## Historical
Observed_Sum = Observed %>%
  group_by(Year,Decade) %>%
  summarise_at(vars(Tmax,Srad),mean) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'
colnames(Observed_Sum)[3] = 'TMAX'

## RCP45
RCP45_Sum = RCP45 %>%
  group_by(Year,Decade) %>%
  summarise_at(vars(TMAX.Corrected,Srad),mean) %>%
  as.data.frame()
RCP45_Sum$Type = 'Downscaled RCP4.5'
colnames(RCP45_Sum)[3] = 'TMAX'

## RCP85
RCP85_Sum = RCP85 %>%
  group_by(Year,Decade) %>%
  summarise_at(vars(TMAX.Corrected,Srad),mean) %>%
  as.data.frame()
RCP85_Sum$Type = 'Downscaled RCP8.5'
colnames(RCP85_Sum)[3] = 'TMAX'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

ggplot(Binded_Data,aes(x=TMAX,y=Srad,color=Decade)) +
  geom_point()

### SRAD VS TMAX PLOTS ###


##################### CHAPTER 6: FUTURE CLIMATE PLOTS TEMP/SRAD #####################

















##################### CHAPTER 6: FUTURE CLIMATE PLOTS PRCP #####################


### PRCP PERCENTILE TABLE ###

quantile(Observed$PRCP,seq(0,1,.1))
quantile(RCP45$PRCP.Corrected,seq(0,1,.1))
quantile(RCP85$PRCP.Corrected,seq(0,1,.1))

### PRCP PERCENTILE TABLE ###








### MEAN PRCP BY YEAR ###

## Historical
Observed_Sum2 = Observed_Full.pde %>%
  group_by(Year) %>%
  summarise(MeanSRAD = var(PRCP)) %>%
  as.data.frame()
Observed_Sum2$Type = 'Observed Full Historical'

## Historical
Observed_Sum = Observed %>%
  group_by(Year) %>%
  summarise(MeanSRAD = var(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## Historical
RCP45_Sum = RCP45.pde %>%
  #filter(Year>2017) %>%
  group_by(Year) %>%
  summarise(MeanSRAD = var(PRCP.Corrected)) %>%
  as.data.frame()
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## Historical
RCP85_Sum = RCP85.pde %>%
  #filter(Year>2017) %>%
  group_by(Year) %>%
  summarise(MeanSRAD = var(PRCP.Corrected)) %>%
  as.data.frame()
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Plotting
Plot4 = ggplot(Binded_Data,aes(x=Year,y=MeanSRAD,color=Type)) +
  ylab('Precipitation') + 
  xlab('Year') + labs(color="Data Source") + 
  geom_line() + geom_point() + #ylim(9,27) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=15)) +
  theme(text=element_text(size=15))

### MEAN PRCP BY YEAR ###





### TOTAL PRCP BY YEAR AND MONTH ###

## Historical
Observed_Sum = Observed %>%
  group_by(Year,Season) %>%
  summarise(MeanSRAD = sum(PRCP)/23) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## RCP45
RCP45_Sum = RCP45 %>%
  group_by(Year,Season) %>%
  summarise(MeanSRAD = sum(PRCP.Corrected)/23) %>%
  as.data.frame()
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## RCP85
RCP85_Sum = RCP85 %>%
  group_by(Year,Season) %>%
  summarise(MeanSRAD = sum(PRCP.Corrected)/23) %>%
  as.data.frame()
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))
#Binded_Data$Season = factor(Binded_Data$Season,levels=c('Winter','Spring','Summer','Fall'))


## Plotting
Plot1=ggplot(Binded_Data,aes(x=Year,y=MeanSRAD,color=factor(Season))) +
  ylab('Total Precipitation') + 
  xlab('Year') + labs(color="Rainfall Season:") + 
  geom_line() + geom_point() +
  geom_smooth(aes(color=factor(Season))) + 
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=14)) +
  theme(text=element_text(size=14)) + 
  coord_cartesian(ylim=c(200,1200)) + 
  scale_y_continuous(breaks=seq(200,1200,200))

ggarrange(Plot1,Plot2,ncol=2,nrow=1,
          common.legend=T,
          legend='bottom')

### TOTAL PRCP BY YEAR AND MONTH ###







### PROPORTION OF DRY DAYS BY DECADE ###

## Creating decade variable
Observed$Decade = ifelse(Observed$Year<1990,'1980-1989','2000-2005')
Observed$Decade = ifelse(Observed$Year>1989 & Observed$Year<2000,
                         '1990-1999',Observed$Decade)

## RCP45
#RCP45$Decade = ifelse(as.Date(RCP45$DATE)>as.Date('2005-11-30') & RCP45$Year<2011,'Dec 2005-2010','2091-2099')
RCP45$Decade = ifelse(RCP45$Year>2005 & RCP45$Year<2011,'2006-2010','2091-2100')
RCP45$Decade = ifelse(RCP45$Year>2010 & RCP45$Year<2021,'2011-2020',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2020 & RCP45$Year<2031,'2021-2030',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2030 & RCP45$Year<2041,'2031-2040',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2040 & RCP45$Year<2051,'2041-2050',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2050 & RCP45$Year<2061,'2051-2060',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2060 & RCP45$Year<2071,'2061-2070',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2070 & RCP45$Year<2081,'2071-2080',RCP45$Decade)
RCP45$Decade = ifelse(RCP45$Year>2080 & RCP45$Year<2091,'2081-2090',RCP45$Decade)

## RCP85
RCP85$Decade = ifelse(RCP85$Year>2005 & RCP85$Year<2011,'2006-2010','2091-2100')
#RCP85$Decade = ifelse(as.Date(RCP85$DATE)>as.Date('2005-11-30') & RCP85$Year<2011,'Dec 2005-2010','2091-2099')
RCP85$Decade = ifelse(RCP85$Year>2010 & RCP85$Year<2021,'2011-2020',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2020 & RCP85$Year<2031,'2021-2030',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2030 & RCP85$Year<2041,'2031-2040',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2040 & RCP85$Year<2051,'2041-2050',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2050 & RCP85$Year<2061,'2051-2060',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2060 & RCP85$Year<2071,'2061-2070',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2070 & RCP85$Year<2081,'2071-2080',RCP85$Decade)
RCP85$Decade = ifelse(RCP85$Year>2080 & RCP85$Year<2091,'2081-2090',RCP85$Decade)

## Historical
Observed_Sum = Observed %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'
Observed_Sum$Ind = as.numeric(Observed_Sum$MeanPRCP < 3)
Observed_Sum2 = Observed_Sum %>%
  group_by(Decade,Month) %>%
  summarise(MeanInd = round(mean(Ind),digits=2)) %>%
  as.data.frame()
Observed_Sum3 = spread(Observed_Sum2,key=Month,value=MeanInd)
Observed_Sum3$YrMean = round(rowMeans(Observed_Sum3[,-1]),digits=2)

## RCP45
RCP45_Sum = RCP45 %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP45_Sum$Type = 'Observed Historical'
RCP45_Sum$Ind = as.numeric(RCP45_Sum$MeanPRCP < 3)
RCP45_Sum2 = RCP45_Sum %>%
  group_by(Decade,Month) %>%
  summarise(MeanInd = round(mean(Ind),digits=2)) %>%
  as.data.frame()
RCP45_Sum3 = spread(RCP45_Sum2,key=Month,value=MeanInd)
RCP45_Sum3$YrMean = round(rowMeans(RCP45_Sum3[,-1]),digits=2)

## RCP85
RCP85_Sum = RCP85 %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP85_Sum$Type = 'Observed Historical'
RCP85_Sum$Ind = as.numeric(RCP85_Sum$MeanPRCP < 3)
RCP85_Sum2 = RCP85_Sum %>%
  group_by(Decade,Month) %>%
  summarise(MeanInd = round(mean(Ind),digits=2)) %>%
  as.data.frame()
RCP85_Sum3 = spread(RCP85_Sum2,key=Month,value=MeanInd)
RCP85_Sum3$YrMean = round(rowMeans(RCP85_Sum3[,-1]),digits=2)

### PROPORTION OF DRY DAYS BY DECADE ###












### DISTRIBUTIONS OF CONSECUTIVE DRY DAYS ###

## Historical
Observed_Sum = Observed %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'
Observed_Sum$Ind = Observed_Sum$MeanPRCP<1
Consec_Dry = rle(Observed_Sum$Ind)$lengths[rle(Observed_Sum$Ind)$values]
Consec_Dry = cbind.data.frame(Days=Consec_Dry,Type='Observed Historical')

## RCP45
RCP45_Sum = RCP45 %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP45_Sum = RCP45_Sum[as.Date(RCP45_Sum$DATE)>as.Date('2005-12-31'),]
RCP45_Sum$Type = 'Downscaled RCP4.5'
RCP45_Sum$Ind = RCP45_Sum$MeanPRCP < 1
Consec_Dry45 = rle(RCP45_Sum$Ind)$lengths[rle(RCP45_Sum$Ind)$values]
Consec_Dry45 = cbind.data.frame(Days=Consec_Dry45,Type='Downscaled RCP4.5')

## RCP85
RCP85_Sum = RCP85 %>%
  group_by(DATE,Month,Decade) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP85_Sum = RCP85_Sum[as.Date(RCP85_Sum$DATE)>as.Date('2005-12-31'),]
RCP85_Sum$Type = 'Downscaled RCP8.5'
RCP85_Sum$Ind = RCP85_Sum$MeanPRCP < 1
Consec_Dry85 = rle(RCP85_Sum$Ind)$lengths[rle(RCP85_Sum$Ind)$values]
Consec_Dry85 = cbind.data.frame(Days=Consec_Dry85,Type='Downscaled RCP8.5')

## Overlaid CDFs of raw, fitted, RCP4.5, and RCP8.5 precipitations
PRCP_Combined = rbind.data.frame(Consec_Dry,Consec_Dry45,Consec_Dry85)
PRCP_Combined$Type = factor(PRCP_Combined$Type,levels=unique(PRCP_Combined$Type))

ggplot(PRCP_Combined,aes(x=Days,color=Type)) +
  stat_ecdf(lwd=1.5) + 
  coord_cartesian(xlim=c(0,20)) +
  xlab('Number of Consecutrive Dry Days') +
  ylab('Probability') + 
  labs(color="Data Source") +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(override.aes=list(linetype=2.5))) +
  theme(legend.text=element_text(size=15)) +
  theme(text=element_text(size=15))


### DISTRIBUTIONS OF CONSECUTIVE DRY DAYS ###








### NUMBER OF EXTREME PRECIPITATION EVENTS ###

## Defining Cutoff
Observed_Sum2 = Observed_Full %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
print(quantile(Observed_Sum2$MeanPRCP,0.5))
Observed_Sum2 = Observed_Sum2 %>%
  group_by(Year,Season) %>%
  summarise(ExtremePRCP = sum(MeanPRCP<5))
Observed_Sum2$Type = 'Observed Historical'

## Defining Cutoff
Observed_Sum = Observed.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
print(quantile(Observed_Sum$MeanPRCP,0.95))
Observed_Sum = Observed_Sum %>%
  group_by(Year,Season) %>%
  summarise(ExtremePRCP = sum(MeanPRCP>42.164))
Observed_Sum$Type = 'Observed Historical'

## RCP45
RCP45_Sum = RCP45.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame() %>%
  group_by(Year,Season)%>%
  summarise(ExtremePRCP = sum(MeanPRCP>42.164))
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## RCP45
RCP85_Sum = RCP85.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame() %>%
  group_by(Year,Season)%>%
  summarise(ExtremePRCP = sum(MeanPRCP>42.164))
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))
#Binded_Data$Season = factor(Binded_Data$Season,levels=c('Winter','Spring','Summer','Fall'))


## Plotting
Plot2 = ggplot(Binded_Data,aes(x=Year,y=ExtremePRCP,color=factor(Season))) +
  ylab('Pcio del Este Extreme Precipitation Events') + 
  xlab('Year') + labs(color="Data Source:") + 
  geom_line() + geom_point(cex=1) + #ylim(0,9) +
  geom_smooth(aes(color=factor(Season))) + 
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=13)) +
  theme(text=element_text(size=13)) + 
  coord_cartesian(ylim=c(0,35)) + 
  scale_y_continuous(breaks=seq(0,35,5))
  
## Gridded plots
ggarrange(Plot1,Plot2,ncol=2,nrow=1,
          common.legend=T,legend='bottom')


### NUMBER OF EXTREME PRECIPITATION EVENTS ###












### PROP OF RAIN IN WET VS DRY SEASON ###


## Rainy vs dry season indicator
Observed$Rainy = ifelse(Observed$Season=='Dry Season','Dry Season','Rainy Season')

## RCP45
RCP45$Rainy = ifelse(RCP45$Season=='Dry Season','Dry Season','Rainy Season')

## RCP85
RCP85$Rainy = ifelse(RCP85$Season=='Dry Season','Dry Season','Rainy Season')


## Historical
Observed_Sum = Observed %>%
  group_by(Year,Rainy) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
Observed_Sum$Type = 'Observed Historical'

## RCP45
RCP45_Sum = RCP45 %>%
  group_by(Year,Rainy) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'

## RCP85
RCP85_Sum = RCP85 %>%
  group_by(Year,Rainy) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'

## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Calculating ratio
RatioDat = cbind.data.frame(
  Binded_Data[1:(nrow(Binded_Data)/2)*2,],
  Ratio=Binded_Data[1:(nrow(Binded_Data)/2)*2,3]/
  Binded_Data[1:(nrow(Binded_Data)/2)*2-1,3])


## Plotting
Plot1 = ggplot(RatioDat,aes(x=Year,y=Ratio,color=factor(Type))) +
  ylab('Wet Season to Dry Season Total Rainfall Ratio') + 
  xlab('Year') + labs(color="Data Source") + 
  geom_line() + geom_point() + ylim(1,3) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=15)) +
  theme(text=element_text(size=15)) + 
  geom_smooth()


### PROP OF RAIN IN WET VS DRY SEASON ###















### PROP OF RAIN HIGH MAG VS LOW MAG ###

## Defining Cutoff
#Observed_Sum = Observed_Full.pde %>%
#  group_by(DATE,Year,Season) %>%
#  summarise(MeanPRCP = mean(PRCP)) %>%
#  as.data.frame()
#print(quantile(Observed_Sum$MeanPRCP,0.95))
#Observed_Sum$Extreme = ifelse(Observed_Sum$MeanPRCP>42.164,1,0)
#Observed_Sum = Observed_Sum %>%
#  group_by(Year,Extreme) %>%
#  summarise(ExtremePRCP_Total = sum(MeanPRCP))
#Observed_Sum$Type = 'Observed Historical'
#length(unique(Observed_Sum$Year))

## Defining Cutoff
Observed_Sum = Observed.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP)) %>%
  as.data.frame()
print(quantile(Observed_Sum$MeanPRCP,0.95))
Observed_Sum$Extreme = ifelse(Observed_Sum$MeanPRCP>42.164,1,0)
Observed_Sum = Observed_Sum %>%
  group_by(Year,Extreme) %>%
  summarise(ExtremePRCP_Total = sum(MeanPRCP))
Observed_Sum$Type = 'Observed Historical'
#length(unique(Observed_Sum$Year))


## RCP45
RCP45_Sum = RCP45.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP45_Sum$Extreme = ifelse(RCP45_Sum$MeanPRCP>42.164,1,0)
RCP45_Sum = RCP45_Sum %>%
  group_by(Year,Extreme) %>%
  summarise(ExtremePRCP_Total = sum(MeanPRCP))
#RCP45_Sum = RCP45_Sum[RCP45_Sum$Year>2005,]
RCP45_Sum$Type = 'Downscaled RCP4.5'
#MissYears = names(which(summary(factor(RCP45_Sum$Year))==1))
#MissYears = cbind.data.frame(Year = MissYears,
#                             Extreme = 1,
#                             ExtremePRCP_Total = 0,
#                             Type = 'Downscaled RCP4.5')
#RCP45_Sum = rbind.data.frame(RCP45_Sum,MissYears)
#RCP45_Sum = RCP45_Sum[order(RCP45_Sum$Year),]


## RCP85
RCP85_Sum = RCP85.pde %>%
  group_by(DATE,Year,Season) %>%
  summarise(MeanPRCP = mean(PRCP.Corrected)) %>%
  as.data.frame()
RCP85_Sum$Extreme = ifelse(RCP85_Sum$MeanPRCP>42.164,1,0)
RCP85_Sum = RCP85_Sum %>%
  group_by(Year,Extreme) %>%
  summarise(ExtremePRCP_Total = sum(MeanPRCP))
#RCP85_Sum = RCP85_Sum[RCP85_Sum$Year>2005,]
RCP85_Sum$Type = 'Downscaled RCP8.5'
#MissYears = names(which(summary(factor(RCP85_Sum$Year))==1))
#MissYears = cbind.data.frame(Year = MissYears,
#                             Extreme = 1,
#                             ExtremePRCP_Total = 0,
#                             Type = 'Downscaled RCP8.5')
#RCP85_Sum = rbind.data.frame(RCP85_Sum,MissYears)
#RCP85_Sum = RCP85_Sum[order(RCP85_Sum$Year),]


## Binding data
Binded_Data = rbind.data.frame(Observed_Sum,RCP45_Sum,RCP85_Sum)
Binded_Data$Type = factor(Binded_Data$Type,levels=unique(Binded_Data$Type))

## Calculating ratio
RatioDat = cbind.data.frame(
  Binded_Data[1:(nrow(Binded_Data)/2)*2-1,],
  Ratio=Binded_Data[1:(nrow(Binded_Data)/2)*2,3]/
    (Binded_Data[1:(nrow(Binded_Data)/2)*2-1,3]))
names(RatioDat)[5] = 'Ratio'
RatioDat$Year = as.numeric(RatioDat$Year)

## Plotting
Plot2 = ggplot(RatioDat,aes(x=Year,y=Ratio,color=factor(Type))) +
  ylab('Extreme Rainfall to Common Rainfall Ratio') + 
  xlab('Year') + labs(color="Data Source") + 
  geom_line() + geom_point() + #ylim(0,1.2) +
  guides(colour=guide_legend(override.aes=list(size=2.5,linetype=0))) +
  theme(legend.position="bottom") +
  theme(legend.text=element_text(size=15)) +
  theme(text=element_text(size=15)) + 
  geom_smooth()

## Gridded plots
ggarrange(Plot1,Plot2,ncol=2,nrow=1,
          common.legend=T,
          legend='bottom')

### PROP OF RAIN HIGH MAG VS LOW MAG ###






##################### CHAPTER 6: FUTURE CLIMATE PLOTS PRCP #####################



