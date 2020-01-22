# Respirometry Analysis
# Created 1/23/19
# Duncan Green

# Script calls on excel file "Respirometry Trials.xlsx"
setwd("~/Dropbox/Respirometry Data") #set working dir.
library(readxl)
library(ggplot2)
library(extrafont)
library(dplyr)
library(data.table)
library(tidyr)
library(nlstools)
# loadfonts()

# To skip all the processing steps and just retrieve oxygen consumption data, skip to ~ line 238

Resp.dat <- read_xlsx("Respirometry Trials.xlsx",1)  # Raw data from trials
Resp.meta <- read_xlsx("Respirometry Trials.xlsx",3) # Metadata from static trials
Resp.meta$FishWeight <- as.double(Resp.meta$FishWeight)
Resp.meta$FishLength <- as.double(Resp.meta$FishLength)
Resp.meta$Date <- as.Date(as.numeric(Resp.meta$Date),origin = "1899-12-30")
unique(Resp.dat$FishID) # examine all fish and background IDs

StatResp <- Resp.dat[is.na(Resp.dat$VFDfreq),] # Static trials only
unique(StatResp$FishID)

# Loop through raw data and manually select start and stop points for each trial

ConsumptionRate <- matrix(NA,100,12)
colnames(ConsumptionRate) <- c("FishID","Trial","S1","E1","S2","E2","S3","E3","S4","E4","S5","E5")

quartz() # Open Quartz graphic device
# Note: maximize quartz window size on mac screen before running code below for accuracy.
counter <- 0 # Start counter at 0 every time loop is run
for (f in unique(StatResp$FishID)){   # 'f' indicates current FishID (or background ID)
  for (t in unique(StatResp$Trial[StatResp$FishID==f])){    # likewise for trial number ('t')
counter <- counter+1
IDsub<- Resp.dat[Resp.dat$FishID==f & Resp.dat$Trial==t, ]
p <- IDsub$DOavg ~ IDsub$ElapsedTime
plot(p, pch=16,cex=0.8, main=c(f,t), xlab="Elapsed Time (min.)", ylab="D.O. (mg/L)")
segpts <- identify(p, plot=FALSE) # click on two points, hit escape to display row numbers
numNA <- ncol(ConsumptionRate) - (2 +length(segpts)) # subtract number of actual values, plus 2 for ID columns
ConsumptionRate[counter,] <- c(f,t,segpts,rep(NA,numNA)) #fill out consumption rate matrix
}}
dev.off() # Close quartz graphic device

# The resulting matrix indicates the point NUMBER (e.g., the xth pair of X,Y pts on the plot) of the 
# start (S) and end (E) of the linear section to be analyzed for rate of oxygen consumption
head(ConsumptionRate)
ConsumptionRate <- na.pass(as.data.frame(ConsumptionRate)) 
ConsumptionRate$Trial <- as.character(ConsumptionRate$Trial)
CRate2 <- ConsumptionRate 

# Change factor vectors to numeric
for (i in 2:12){
  CRate2[ ,i] <- as.numeric(as.character((CRate2[,i])))
}
CRate2$FishID <- as.character(CRate2$FishID)
len <- length(na.exclude(CRate2$FishID)) # number of non-NA rows
CRate2 <- CRate2[1:len,]
tail(CRate2)

# Need to use our new index (CRate2) of point numbers to query elapsed time and DO in StatResp
#for (i in 1:52) {# we have 52 rows of data in CRate2, i is row index
#  FTname <- c(as.character(CRate2[i,1]),as.character(CRate2[i,2])) # Fish & Trial name
    
#  for (j in )
#}

# Different approach
#CRate3 <- CRate2
#DOoutput <- matrix(NA,100,7)
#colnames(DOoutput) <- c("FishID", "Trial", "R1","R2","R3","R4","R5")

#counter <- 0 # Start counter at 0 every time loop is run
#for (f in unique(CRate3$FishID)){   # 'f' indicates current FishID (or background ID)
#  for (t in unique(CRate3$Trial[CRate3$FishID==f])){    # likewise for trial number ('t')
#    counter <- counter+1
    
#  }}

#Gave up at this point, decided to analyze manually
#write.csv(CRate2, "CRate2.csv") # export table of values 


# Lookup table method

CRate2get <- function(Data=CRate2, f, t) { # Function to return row index
  return(which(CRate2$FishID == f & CRate2$Trial == t))
}

# Name parameters
metpar2 <- matrix(NA,length(na.exclude(CRate2$FishID)),12) # empty matrix to hold metabolic equation parameters
colnames(metpar2) <- c("Slope1","Intercept1","AdjRsq1",
                      "Slope2","Intercept2","AdjRsq2",
                      "Slope3","Intercept3","AdjRsq3",
                      "Slope4","Intercept4","AdjRsq4")
rownames(metpar2) <- paste(na.exclude(CRate2$FishID),na.exclude(CRate2$Trial))
metpar2[,1:12] <- as.double(metpar2[,1:12]) #Coerce to numeric variables

# BEFORE you run this code!!! Change PDF save as name for each iteration/version (i.e. V2,V3)
counter <- 0 #reset our counter
for (f in unique(na.exclude(CRate2$FishID))){
  for (t in unique(na.exclude(CRate2$Trial[CRate2$FishID==f]))){
    b1<-NA;b2<-NA;b3<-NA;b4<-NA;m1<-NA;m2<-NA;m3<-NA;m4<-NA;
    r1<-NA;r2<-NA;r3<-NA;r4<-NA;fitsub1<-NA;fitsub2<-NA;
    fitsub3<-NA;fitsub4<-NA #reset all
    index <- CRate2get(f=f,t=t) # Row number
    IDsub <- Resp.dat[Resp.dat$FishID==f & Resp.dat$Trial==t, ]
    p <- IDsub$DOavg ~ IDsub$ElapsedTime
    pdf(paste("F",f,"T",t,"V5.pdf", sep =""))
    plot(p, pch=16,cex=2, main=paste(f,t), xlab="Elapsed Time (min.)", ylab="D.O. (mg/L)")
    tryCatch({
      fitsub1 <- lm(IDsub$DOavg[CRate2[index,3]:CRate2[index,4]] ~ IDsub$ElapsedTime[CRate2[index,3]:CRate2[index,4]])
      b1 <- fitsub1$coefficients[1] # Intercept
      m1 <- fitsub1$coefficients[2] # slope (MO2)
      r1 <- as.double(summary(fitsub1)[9]) # adjusted r-squared
      fitsub2 <- lm(IDsub$DOavg[CRate2[index,5]:CRate2[index,6]] ~ IDsub$ElapsedTime[CRate2[index,5]:CRate2[index,6]])
      b2 <- fitsub2$coefficients[1] # Intercept
      m2 <- fitsub2$coefficients[2] # slope (MO2)
      r2 <- as.double(summary(fitsub2)[9]) # adjusted r-squared
      fitsub3 <- lm(IDsub$DOavg[CRate2[index,7]:CRate2[index,8]] ~ IDsub$ElapsedTime[CRate2[index,7]:CRate2[index,8]])
      b3 <- fitsub3$coefficients[1] # Intercept
      m3 <- fitsub3$coefficients[2] # slope (MO2)
      r3 <- as.double(summary(fitsub3)[9]) # adjusted r-squared
      fitsub4 <- lm(IDsub$DOavg[CRate2[index,9]:CRate2[index,10]] ~ IDsub$ElapsedTime[CRate2[index,9]:CRate2[index,10]])
      b4 <- fitsub4$coefficients[1] # Intercept
      m4 <- fitsub4$coefficients[2] # slope (MO2)
      r4 <- as.double(summary(fitsub4)[9]) # adjusted r-squared
    }, error=function(e){}) # end of tryCatch section
    counter <- counter+1
    metpar2[counter, 1] <- m1 #slope1 parameter
    metpar2[counter, 2] <- b1 #intercept1 parameter
    metpar2[counter, 3] <- r1 #adjusted rsquared
    metpar2[counter, 4] <- m2 
    metpar2[counter, 5] <- b2 
    metpar2[counter, 6] <- r2
    metpar2[counter, 7] <- m3 
    metpar2[counter, 8] <- b3 
    metpar2[counter, 9] <- r3
    metpar2[counter, 10] <- m4 
    metpar2[counter, 11] <- b4 
    metpar2[counter, 12] <- r4
    tryCatch({
    points(IDsub$DOavg[CRate2[index,3]:CRate2[index,4]] ~ IDsub$ElapsedTime[CRate2[index,3]:CRate2[index,4]], pch=16,cex=2,col=4)
    fitsub1 <- lm(IDsub$DOavg[CRate2[index,3]:CRate2[index,4]] ~ IDsub$ElapsedTime[CRate2[index,3]:CRate2[index,4]])
    abline(fitsub1, lwd=2,col=2)
    fitsub1$coefficients[1] # Intercept
    fitsub1$coefficients[2] # slope (MO2)
    summary(fitsub1)[9] # adjusted r-squared
    points(IDsub$DOavg[CRate2[index,5]:CRate2[index,6]] ~ IDsub$ElapsedTime[CRate2[index,5]:CRate2[index,6]], pch=16,cex=2,col=4)
    fitsub2 <- lm(IDsub$DOavg[CRate2[index,5]:CRate2[index,6]] ~ IDsub$ElapsedTime[CRate2[index,5]:CRate2[index,6]])
    abline(fitsub2, lwd=2,col=2)
    fitsub2$coefficients[1] # Intercept
    fitsub2$coefficients[2] # slope (MO2)
    summary(fitsub2)[9] # adjusted r-squared
    points(IDsub$DOavg[CRate2[index,7]:CRate2[index,8]] ~ IDsub$ElapsedTime[CRate2[index,7]:CRate2[index,8]], pch=16,cex=2,col=4)
    fitsub3 <- lm(IDsub$DOavg[CRate2[index,7]:CRate2[index,8]] ~ IDsub$ElapsedTime[CRate2[index,7]:CRate2[index,8]])
    abline(fitsub3, lwd=2,col=2)
    fitsub3$coefficients[1] # Intercept
    fitsub3$coefficients[2] # slope (MO2)
    summary(fitsub3)[9] # adjusted r-squared
    points(IDsub$DOavg[CRate2[index,9]:CRate2[index,10]] ~ IDsub$ElapsedTime[CRate2[index,9]:CRate2[index,10]], pch=16,cex=2,col=4)
    fitsub4 <- lm(IDsub$DOavg[CRate2[index,9]:CRate2[index,10]] ~ IDsub$ElapsedTime[CRate2[index,9]:CRate2[index,10]])
    abline(fitsub4, lwd=2,col=2)
    fitsub4$coefficients[1] # Intercept
    fitsub4$coefficients[2] # slope (MO2)
    summary(fitsub4)[9] # adjusted r-squared
    }, error=function(e){}) # end of tryCatch section 
    dev.off() # Close PDF device to save plot
  }
}

metpar2 <- as.data.frame(metpar2)
setDT(metpar2, keep.rownames = "FishID")[]

# Mean slope for each fish/trial combo
for (i in 1:length(metpar2$FishID)){
  metpar2$MeanSlope[i] <- rowMeans(metpar2[i,c(2,5,8,11)],na.rm = T)
    #mean(c(metpar2[i,2],metpar2[i,5], metpar2[i,8],metpar2[i,11]),na.rm = TRUE)
  }

metpar2[,c(2,5,8,11,14)] # View slopes only

# Add columns to metadata

head(Resp.meta)
Resp.meta <- bind_cols(Resp.meta,metpar2[,2:14])
Resp.meta[51,1] <- "S31 T2" #Fixing redundant S31 trial labels (both were "S31")
head(Resp.meta)
tail(Resp.meta)

ggplot(Resp.meta, aes(as.factor(FishID),Slope1, col=as.factor(TempC))) + geom_point() + labs(title="mO2") +
  geom_point(aes(y=Slope2)) +
  geom_point(aes(y=Slope3)) +
  geom_point(aes(y=Slope4)) +
  geom_point(aes(y=MeanSlope),pch=3)

# Background
PerA <- mean(c(0,Resp.meta$MeanSlope[9])) # Background rate for period A; 0 used instead of an artificial pos. value for SRB2
PerB <- mean(c(Resp.meta$MeanSlope[9],Resp.meta$MeanSlope[13])) # Background rate for period B
PerC <- mean(c(Resp.meta$MeanSlope[13],Resp.meta$MeanSlope[20])) # Background rate for period C
PerD <- mean(c(Resp.meta$MeanSlope[20],Resp.meta$MeanSlope[35])) # Background rate for period D
PerE <- mean(c(Resp.meta$MeanSlope[36],Resp.meta$MeanSlope[52])) # Background rate for period E
Resp.meta$Background <- c(rep(0,5),rep(PerA,3),0,rep(PerB,3),0,rep(PerC,6),0,
                          rep(PerD,14),0,0,rep(PerE,15),0) # New variable for background resp.

Resp.meta$g.g.d <- rep(NA, length(Resp.meta$FishID)) # Set up new variable for code below:
Resp.meta$Cor.ggd <- rep(NA, length(Resp.meta$FishID)) # Same for background-corrected variable
#calculating SMR from above slopes. Clunky code for ease of comprehension.
for (i in 1:52){
s.avg <- (Resp.meta$MeanSlope[i]) # mean slope of oxygen decline for this fish/trial
mg.L.min <- -s.avg #just to keep track of what the number means (removed negative sign too)
rVolume <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume <- rVolume - fVolume # effective volume of water in chamber (liters)
mg.min <- mg.L.min*eVolume # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr <- mg.min*60 #mg/hr
mg.g.hr <- mg.hr/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d <- (mg.g.hr/1000)*24  # g O2 per g fish per day (as denoted in FB4) based on the MEAN of each fish/trial
Resp.meta$g.g.d[i] <- g.g.d
# Corrected variable:
s.avg2 <- (Resp.meta$MeanSlope[i] - Resp.meta$Background[i]) # mean slope of oxygen decline for this fish/trial
mg.L.min2 <- -s.avg2 #just to keep track of what the number means (removed negative sign too)
rVolume2 <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume2 <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume2 <- rVolume2 - fVolume2 # effective volume of water in chamber (liters)
mg.min2 <- mg.L.min2*eVolume2 # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr2 <- mg.min2*60 #mg/hr
mg.g.hr2 <- mg.hr2/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d2 <- (mg.g.hr2/1000)*24  # g O2 per g fish per day (as denoted in FB4) based on the MEAN of each fish/trial
Resp.meta$Cor.ggd[i] <- g.g.d2 # Correct for background respiration, calculated as average background respiration of previous and
# subsequent measurements
# Corrected variable, slope 1 only:
s.avg2 <- (Resp.meta$Slope1[i] - Resp.meta$Background[i]) #slope of oxygen decline for slope/period 1 of this fish/trial
mg.L.min2 <- -s.avg2 #just to keep track of what the number means (removed negative sign too)
rVolume2 <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume2 <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume2 <- rVolume2 - fVolume2 # effective volume of water in chamber (liters)
mg.min2 <- mg.L.min2*eVolume2 # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr2 <- mg.min2*60 #mg/hr
mg.g.hr2 <- mg.hr2/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d2 <- (mg.g.hr2/1000)*24  # g O2 per g fish per day (as denoted in FB4) 
Resp.meta$Cor.ggd1[i] <- g.g.d2 # Correct for background respiration, calculated as average background respiration of previous and
# Corrected variable, slope 2 only:
s.avg2 <- (Resp.meta$Slope2[i] - Resp.meta$Background[i]) #slope of oxygen decline for slope/period 2 of this fish/trial
mg.L.min2 <- -s.avg2 #just to keep track of what the number means (removed negative sign too)
rVolume2 <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume2 <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume2 <- rVolume2 - fVolume2 # effective volume of water in chamber (liters)
mg.min2 <- mg.L.min2*eVolume2 # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr2 <- mg.min2*60 #mg/hr
mg.g.hr2 <- mg.hr2/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d2 <- (mg.g.hr2/1000)*24  # g O2 per g fish per day (as denoted in FB4) 
Resp.meta$Cor.ggd2[i] <- g.g.d2 # Correct for background respiration, calculated as average background respiration of previous and
# Corrected variable, slope 3 only:
s.avg2 <- (Resp.meta$Slope3[i] - Resp.meta$Background[i]) #slope of oxygen decline for slope/period 3 of this fish/trial
mg.L.min2 <- -s.avg2 #just to keep track of what the number means (removed negative sign too)
rVolume2 <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume2 <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume2 <- rVolume2 - fVolume2 # effective volume of water in chamber (liters)
mg.min2 <- mg.L.min2*eVolume2 # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr2 <- mg.min2*60 #mg/hr
mg.g.hr2 <- mg.hr2/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d2 <- (mg.g.hr2/1000)*24  # g O2 per g fish per day (as denoted in FB4) 
Resp.meta$Cor.ggd3[i] <- g.g.d2 # Correct for background respiration, calculated as average background respiration of previous and
# Corrected variable, slope 4 only:
s.avg2 <- (Resp.meta$Slope4[i] - Resp.meta$Background[i]) #slope of oxygen decline for slope/period 4 of this fish/trial
mg.L.min2 <- -s.avg2 #just to keep track of what the number means (removed negative sign too)
rVolume2 <- Resp.meta$RespVolL[i] # volume of respirometer (in liters) for this trial
fVolume2 <- Resp.meta$FishWeight[i]/1000  # volume of fish (mL is ~mass in grams) in liters
eVolume2 <- rVolume2 - fVolume2 # effective volume of water in chamber (liters)
mg.min2 <- mg.L.min2*eVolume2 # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
mg.hr2 <- mg.min2*60 #mg/hr
mg.g.hr2 <- mg.hr2/Resp.meta$FishWeight[i] #rate in mg O2 per g fish per hour
g.g.d2 <- (mg.g.hr2/1000)*24  # g O2 per g fish per day (as denoted in FB4) 
Resp.meta$Cor.ggd4[i] <- g.g.d2 # Correct for background respiration, calculated as average background respiration of previous and

}

head(Resp.meta)

ggplot(Resp.meta, aes(as.factor(FishID),Cor.ggd1, col=as.factor(TempC))) + geom_point() + 
  labs(title="mO2",xlab="Fish ID",ylab="mO2 (g/g/day)") +
  geom_point(aes(y=Cor.ggd2)) +
  geom_point(aes(y=Cor.ggd3)) +
  geom_point(aes(y=Cor.ggd4)) +
  geom_point(aes(y=Cor.ggd),pch=3,col=1)

# Subset of only values for fish (no tests or BG trials)
fishsub <- Resp.meta[Resp.meta$Background !=0,]
badfish <- c("S2","S6","S12","S15","S15 T2","S15 T3","S16",
  "S25","S25 T2","S26","S26 T2", "S27","S31") # trials flagged in data (erratic swimming, air in chamber, etc.)

# Save below as a different number (fish3sub, fish4sub, etc.) if above analyses were redone
fish3sub <- fishsub[!fishsub$FishID %in% badfish,] #remove bad trials
fish3sub 

#Again, change below number!
write.csv(fish3sub, "Fish3sub.csv") # Save data externally 

# fish2sub <- read.csv("fish2sub.csv",1)     # Retrieve .csv if necessary
# fish2sub <- fish2sub[!fish2sub$FishID %in% badfish,] # double check that bad trials are removed

ggplot(fish3sub, aes(as.factor(FishID),Cor.ggd1, col=as.factor(TempC))) + geom_point() + 
  labs(title="mO2",xlab="Fish ID",ylab="mO2 (g/g/day)") +
  geom_point(aes(y=Cor.ggd2)) +
  geom_point(aes(y=Cor.ggd3)) +
  geom_point(aes(y=Cor.ggd4)) +
  geom_point(aes(y=Cor.ggd),pch=3,col=1)

sub5C <- fish3sub[fish3sub$TempC==5,] #5C trials only
sub10C <- fish3sub[fish3sub$TempC==10,] #10C trials only
sub15C <- fish3sub[fish3sub$TempC==15,] #15C trials only

#Graphical examination
plot(sub5C$Cor.ggd ~ sub5C$FishWeight)
lm5 <- lm(sub5C$Cor.ggd ~ sub5C$FishWeight)
summary(lm5)

plot(sub10C$Cor.ggd ~ sub10C$FishWeight)
lm10 <- lm(sub10C$Cor.ggd ~ sub10C$FishWeight)
summary(lm10)

plot(sub15C$Cor.ggd ~ sub15C$FishWeight)
lm15 <- lm(sub15C$Cor.ggd ~ sub15C$FishWeight)
summary(lm15)

# Reformat into long (untidy) data?

longfish <- rbind(fish3sub[,c(1:3,5:7,24)],
              setnames(fish3sub[,c(1:3,5:7,25)], names(fish3sub[,c(1:3,5:7,24)])),
              setnames(fish3sub[,c(1:3,5:7,26)], names(fish3sub[,c(1:3,5:7,24)])),
              setnames(fish3sub[,c(1:3,5:7,27)], names(fish3sub[,c(1:3,5:7,24)]))
              )
colnames(longfish)[7] <- "Cor.ggd"

ggplot(longfish, aes(as.factor(FishID),Cor.ggd, col=as.factor(TempC))) + geom_point() + 
  labs(title="mO2",xlab="Fish ID",ylab="mO2 (g/g/day)")

# Looking at minimum rate for each trial?
for (i in 1:length(fish3sub$FishID)){
  fish3sub$MinRate[i] <- min(fish3sub$Cor.ggd1[i],fish3sub$Cor.ggd2[i],
                          fish3sub$Cor.ggd3[i],fish3sub$Cor.ggd4[i], na.rm = TRUE)
}



#Graphical examination of all measurement periods instead of just means
sub5Cl <- longfish[longfish$TempC==5,] #5C trials only
sub10Cl <- longfish[longfish$TempC==10,] #10C trials only
sub15Cl <- longfish[longfish$TempC==15,] #15C trials only

plot(sub5Cl$Cor.ggd ~ sub5Cl$FishWeight)
lm5 <- lm(sub5Cl$Cor.ggd ~ sub5Cl$FishWeight)
summary(lm5)
abline(lm5)

plot(sub10Cl$Cor.ggd ~ sub10Cl$FishWeight)
lm10 <- lm(sub10Cl$Cor.ggd ~ sub10Cl$FishWeight)
summary(lm10)
abline(lm10)

plot(sub15Cl$Cor.ggd ~ sub15Cl$FishWeight)
lm15 <- lm(sub15Cl$Cor.ggd ~ sub15Cl$FishWeight)
summary(lm15)
abline(lm15)

plot(Cor.ggd ~ FishWeight, data = longfish)

longfish <- longfish[!is.na(longfish$Cor.ggd),]

longfish2 <- longfish[,-1] # Trying to get nlstools to work (need finite xlim)

f1 <- as.formula(Cor.ggd~RA*FishWeight^(RB)*exp(RQ*TempC))
preview(f1, data = longfish2, start = list(RA = 0.025, RB = -0.9, RQ= 0.05)) #requires package 'nlstools'
fit2 <- nls(f1, data = longfish2, start = list(RA = 0.025, RB = -0.9, RQ= 0.05))
overview(fit2)
preview(f1, data = longfish2, start = list(RA = coefficients(fit2)[1], 
                                           RB = coefficients(fit2)[2],
                                           RQ = coefficients(fit2)[3]))

# Fitting allometric function. Baseplot

plot(Cor.ggd ~ FishWeight, data = fish2sub) # all points (all temps)
lmlog <- lm(log(Cor.ggd) ~ log(FishWeight),data = fish2sub) #linear model of log-transformed data
summary(lmlog)
curve(exp(lmlog$coefficients[1])*x^(lmlog$coefficients[2]),add = TRUE,col=1) # coefficients [1] and [2] are 
# intercept and slope, respectively. This backtransforms linear coefs into our allometric (y~ax^b) form

lmlog5 <- lm(log(Cor.ggd) ~ log(FishWeight),data = sub5C) # 5C trials
summary(lmlog5)
curve(exp(lmlog5$coefficients[1])*x^(lmlog5$coefficients[2]),add = TRUE,col=4)
points(Cor.ggd ~ FishWeight, data = sub5C, col=4)

lmlog10 <- lm(log(Cor.ggd) ~ log(FishWeight),data = sub10C) # 10C trials
summary(lmlog10)
curve(exp(lmlog10$coefficients[1])*x^(lmlog10$coefficients[2]),add = TRUE,col=3)
points(Cor.ggd ~ FishWeight, data = sub10C, col=3)

lmlog15 <- lm(log(Cor.ggd) ~ log(FishWeight),data = sub15C) # 15C trials
summary(lmlog15)
curve(exp(lmlog15$coefficients[1])*x^(lmlog15$coefficients[2]),add = TRUE,col=2)
points(Cor.ggd ~ FishWeight, data = sub15C, col=2)

# Metabolic rate as a function of temperature
par(mfrow=c(2,2))
plot(Cor.ggd ~ TempC, data = fish2sub, xlim=c(0,15),ylim=c(-0.003,0.010),main="linear") # all points
lmtemp <- lm(Cor.ggd ~ TempC, data = fish2sub)
summary(lmtemp)
abline(h=0,lty=2)

lmlogtemp <- lm(log(Cor.ggd) ~ TempC, data = fish2sub) # log-linear relationship
summary(lmlogtemp)
plot(log(Cor.ggd) ~ TempC, data = fish2sub,xlim=c(0,15),main = "log-linear")
abline(lmlogtemp)


lmloglogtemp <- lm(log(Cor.ggd) ~ log(TempC),data = fish2sub) #log-log
summary(lmloglogtemp)
plot(log(Cor.ggd) ~ log(TempC),data = fish2sub, main="log-log",xlim=c(0,log(15)))
abline(lmloglogtemp)

fish2sub$TempK <- fish2sub$TempC + 273.15 #Temp Kelvin
fish2sub$invTempK <- 1000/fish2sub$TempK  # Inverse temp in Kelvin, x1000
lmArrtemp <- lm(log(Cor.ggd) ~ invTempK,data = fish2sub) #Arrhenius
summary(lmArrtemp)
plot(log(Cor.ggd) ~ invTempK,data = fish2sub ,main= "Arrhenius",xlim=c(3.46,3.66))
abline(lmArrtemp)
# Export pdf of plot here

# Fitting model using nls()

# Fit RQ parameter using RA and RB parameters from Generalized Coregonid (GC) model
f0 <- Cor.ggd ~ 0.0018*FishWeight^(-0.12)*exp(RQ*TempC)
fit.gc <- nls(f0, data = fish2sub, start = list(RQ= 0.047))
summary(fit.gc)
RQ.gc <- coefficients(fit.gc)[1]


plot(Cor.ggd ~ FishWeight, data = fish2sub)
points(Cor.ggd ~ FishWeight, data = sub5C, col=1) 
points(Cor.ggd ~ FishWeight, data = sub10C, col=3)
points(Cor.ggd ~ FishWeight, data = sub15C, col=2)
curve(0.0018*x^(-0.12)*exp(RQ.gc*5),add = TRUE,lty=3, col=1) # GC model with my RQ value
curve(0.0018*x^(-0.12)*exp(RQ.gc*10),add = TRUE,lty=3, col = 3)
curve(0.0018*x^(-0.12)*exp(RQ.gc*15),add = TRUE,lty=3, col=2)

# Madenjian (mj) et al. 2006 RA value
f0 <- Cor.ggd ~ 0.00085*FishWeight^(-0.12)*exp(RQ*TempC)
fit.mj <- nls(f0, data = fish2sub, start = list(RQ= 0.047))
summary(fit.mj)
RQ.mj <- coefficients(fit.mj)[1]

curve(0.00085*x^(-0.12)*exp(RQ.mj*5),add = TRUE,lty=1, col=1) # MJ model with my RQ value
curve(0.00085*x^(-0.12)*exp(RQ.mj*10),add = TRUE,lty=1, col = 3)
curve(0.00085*x^(-0.12)*exp(RQ.mj*15),add = TRUE,lty=1, col=2)

# Fitting all three parameters to my data. Note insignificance of Weight coefficients (not enough range in fish sizes).
f0 <- Cor.ggd ~ RA*FishWeight^(RB)*exp(RQ*TempC);
fit2 <- nls(f0, data = fish2sub, start = list(RA = 0.001, RB = -0.01, RQ= 0.01))
summary(fit2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*15),add = TRUE,lty=2, col=2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*10),add = TRUE,lty=2, col=3)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*5),add = TRUE,lty=2, col=1)

# function useful for plotting (needs work)
#f <- function(FishWeight,TempC,RA,RB,RQ) {RA*FishWeight*exp(RB)*exp(RQ*TempC)}
#co <- coefficients(fit2)
#curve(f(x, RA=co[1], RB=co[2], RQ=co[3], TempC=5), add = TRUE, col="green", lwd=2) 

# plotting oxygen consumption as a function of temperature (exponential) for a 1 g fish
plot(Cor.ggd ~ TempC, data = fish2sub, xlim=c(0,15),ylim=c(-0.003,0.010),main="exponential") # all points
f0 <- Cor.ggd ~ 0.00085*1^(-0.12)*exp(RQ*TempC)
fit.temp <- nls(f0, data = fish2sub, start = list(RQ= 0.047))
summary(fit.temp)
RQ.temp <- coefficients(fit.temp)[1]
curve(0.00085*1^(-0.12)*exp(RQ.temp*x),add = TRUE,lty=1, col=4) #using Madenjian RA/RB

f0 <- Cor.ggd ~ 0.0018*1^(-0.12)*exp(RQ*TempC)
fit.temp <- nls(f0, data = fish2sub, start = list(RQ= 0.047))
summary(fit.temp)
RQ.temp <- coefficients(fit.temp)[1]
curve(0.0018*1^(-0.12)*exp(RQ.temp*x),add = TRUE,lty=1, col=2) #using Gen. Cor. RA/RB

f0 <- Cor.ggd ~ 0.0018*FishWeight^(-0.12)*exp(RQ*TempC)
fit.temp <- nls(f0, data = fish2sub, start = list(RQ= 0.047))
summary(fit.temp)
RQ.temp <- coefficients(fit.temp)[1]
curve(0.0018*1^(-0.12)*exp(RQ.temp*x),add = TRUE,lty=1, col=3) #using Gen. Cor. RA/RB, for all fish weights (not just 1 g)

f0 <- Cor.ggd ~ RA*FishWeight^(RB)*exp(RQ*TempC);
fit.temp <- nls(f0, data = fish2sub, start = list(RA= 0.002, RB = -0.1, RQ= 0.05))
summary(fit.temp)
co <- coefficients(fit.temp)
curve(co[1]*1^(co[2])*exp(co[3]*x),add = TRUE,lty=2, col=1) #fitting all three pars with my data


diag(vcov(fit2)) #variances of model parameters
CI <- function(fitted.object, level=0.95) {
  # Function written by Franz Mueter. fitted.object has to be the result of lm, glm, gam, nls, or any other model
  # that has a method for the 'vcov' function
  cf <- coef(fitted.object) # extract coefficients
  alpha <- 1-level
  z <- qnorm(1-alpha/2) # Upper quantile corresponding to desired CI at 'level'
  SE <- sqrt(diag(vcov(fitted.object)))   # Standard errors of parameters
  cbind(lower = cf-z*SE, upper=cf+z*SE)   # Confidence intervals
}
CI(fit.gc)



#########
#### Again, using RA and RB values from generalised Coregonid (Rudstam et al. 1994) values
fo <- as.formula(Cor.ggd ~ RA*FishWeight^(RB)*exp(RQ*TempC))
preview(fo, data = fish2sub, start = list(RA = 0.00085, RB = -0.12, RQ= 0.047)) #requires package 'nlstools'
fit2 <- nls(fo, data = fish2sub, start = list(RA = 0.00085, RB = -0.12, RQ= 0.047))
overview(fit2)

# S3 method for nls - Code from Derek Ogle
library(FSA)
residPlot(
  fit2,
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "",
  pch = 16,
  col = "black",
  lty.ref = 3,
  lwd.ref = 1,
  col.ref = "black",
  resid.type = c("raw", "standardized", "studentized"),
  loess = FALSE,
  lty.loess = 2,
  lwd.loess = 1,
  col.loess = "black",
  trans.loess = 8,
  inclHist = TRUE
)

library(nleqslv)
fn <- function(x){
  c(-0.348475*x[1]*x[2] + x[1] + 0.25*x[2] - 0.6061,
    -0.339275*x[1]*x[2] + x[1] + 0.25*x[2] - 0.6429)
}

nleqslv(c(1, 1), fn)

# function taken from https://stackoverflow.com/questions/18305852/power-regression-in-r-similar-to-excel
power_eqn = function(df, start = list(a =0.001,b=0.1)){
  m = nls(Cor.ggd ~ a*FishWeight^b, start = start, data = df);
  eq <- substitute(italic(y) == a  ~italic(x)^b, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2)))
  as.character(as.expression(eq));                 
}

whichsub <- sub10C #change to subset of interest (sub15C, etc.)
ggplot(whichsub,aes(x = FishWeight,y = Cor.ggd)) +
  geom_point() + 
  stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start= c(a = 0.001,b=0.1)),se=FALSE) +  
  geom_text(x = 15, y = 0.003, label = power_eqn(whichsub), parse = TRUE)

# All temperatures combined
ggplot(fish2sub,aes(x = FishWeight,y = Cor.ggd,col=TempC)) +
  geom_point() + 
  stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start= c(a = 0.001,b=0.1)),se=FALSE) +  
  geom_text(x = 15, y = 0.003, label = power_eqn(fish2sub), parse = TRUE)

cv <- function(x) sd(x)/mean(x) # Coefficient of variation function




# Chabot p88: SMR linear with body mass when both are log transformed



# Example of TryCatch for continuing to run loop despite error message:
#for (i in 1:10) {
#  tryCatch({
#    print(i)
#    if (i==7) stop("Urgh, the iphone is in the blender !")
#  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#}


# fitsub <- lm(IDsub$DOavg[minsub:maxsub] ~ IDsub$ElapsedTime[minsub:maxsub])
# abline(fitsub, lwd=2,col=2)
# fitsub$coefficients[2] # slope (MO2)
# summary(fitsub)[9] # adjusted r-squared









# Old code:

#Pair each set of fish with the appropriate background trials:

fish <- "SRB2" # set to desired fishID (or background respiration ID)
trial <- "1" # likewise for trial number
IDsub <- Resp.dat[Resp.dat$FishID==fish & Resp.dat$Trial==trial, ]
p <- IDsub$DOavg ~ IDsub$ElapsedTime
plot(p, pch=16,cex=2, main="Dissolved Oxygen vs. Time", xlab="Elapsed Time (min.)", ylab="D.O. (mg/L)")
identify(p, plot=FALSE) # click on two points, hit escape to display row numbers
minsub <- 5 # point number from left indicating start of linear section
maxsub <- 24 # point number indicating end of linear section
points(IDsub$DOavg[minsub:maxsub] ~ IDsub$ElapsedTime[minsub:maxsub], pch=16,cex=2,col=4)
fitsub <- lm(IDsub$DOavg[minsub:maxsub] ~ IDsub$ElapsedTime[minsub:maxsub])
abline(fitsub, lwd=2,col=2)
fitsub$coefficients[2] # slope (MO2)
summary(fitsub)[9] # adjusted r-squared


# abline(v=c(12,52,80,127,153))

#calculating SMR from above slopes
s.avg <- (-0.02094847) #placeholder: average slope of linear portion of O2 decline across all trials of interest
mg.L.min <- -s.avg #just to keep track of what the number means (removed negative sign too)
rVolume <- 2.000 #placeholder: volume of respirometer in L
fVolume <- 0.015 #placeholder: volume of fish (~mass in grams) in L
eVolume <- rVolume - fVolume # effective volume of water in chamber
(mg.min <- mg.L.min*eVolume) # raw rate of oxygen decline (i.e., mg O2 per min: not volume specific)
(mg.hr <- mg.min*60) #mg/hr
(mg.g.hr <- mg.hr/(fVolume*1000)) #rate in mg O2 per g fish per hour
(g.g.d <- (mg.g.hr/1000)*24 ) # g O2 per g fish per day (as denoted in Deslauriers et al. 2017)

# Chabot p88: SMR linear with body mass when both are log transformed

# Same in ggplot2:
# needs work # Resp.dat.stat <- Resp.dat[Resp.dat$FishID==c("S1","S2","S3","S4","S5","S6","S7","S8"),] # only plotting static fish trials (no background trials)
p1 <- ggplot(Resp.dat, aes(ElapsedTime, DOavg, col=as.character(Trial))) + geom_point()
p1 + facet_wrap(~FishID)

# Obtain oxygen consumption rates for each trial
# ! Come back to this
MO2 Values<- matrix() # Output matrix

(allID <- unique(Resp.dat$FishID)) # save a list of all fish ID and background ID
for i in (1:length(allID)) {
  ID <- allID[i]
  IDsub <- Resp.dat[Resp.dat$FishID==ID, ]
  fitsub <- lm(IDsub$DOavg ~ IDsub$ElapsedTime)
  fitsub$coefficients[2]

MO2 <- 
}

#########
x <- c(1:11)
y <- c(slopepre,slope4.1.1, slope4.1.2, slope4.1.3,slope4.2.1, slope4.2.2, slope4.2.3,slope4.3.1, slope4.3.2,slopepost, slopeB14)
plot(x,y, pch=16)

##########

# Lookup table method

library(dplyr)
library(tidyr)
table %>%
  gather(key = "pet") %>%
  left_join(lookup, by = "pet") %>%
  spread(key = pet, value = class)