# Parameterization script: juvenile broad whitefish, Fish Bioenergetics 4.0
# Created by Duncan Green on 1/22/20

setwd("/Users/duncangreen/Documents/GreenThesis") # current working directory
library(readxl)
library(ggplot2)
library(extrafont)
library(dplyr)
library(data.table)
library(tidyr)
library(nlstools)

# RESPIRATION

load("RespDataProcessed.RData") # Processed respirometry data: same data in two forms, one tidy/long (RespDatR) and one wide (Fish3sub)
# NOTE: For raw data processing details, see 'RespirometryAnalysis.R', which calls on "Respirometry Trials.xlsx" and contains all processing steps. 

plot(MO2.corrected ~ FishWeight.g, data = RespDatP[RespDatP$Period!="Mean",]) # All measurements from respirometry trials at all three 
# temperatures (5,10,15C). Generally, each trial/fish had three measurement periods, from which I took the mean MO2 (oxygen consumption rate):

MeanResp <- RespDatP[RespDatP$Period == "Mean",] # Mean MO2 rates only
plot(MO2.corrected ~ FishWeight.g, data = MeanResp)

# Viewed by temperature treatment:
ggplot(MeanResp, aes(FishWeight.g, MO2.corrected, col=as.factor(TempC),pch=as.factor(TempC))) + geom_point(cex=3) +
  labs(xlab="Fish ID",ylab="mO2 (g/g/day)")


# Fitting respiration model (this is 'Model 1' in the Respiration section of FB4) using nls(). Using several of the "closest" species to borrow parameters from, 
# but there are more options within existing FB4 models. These are just some examples. 

# Model form: ' MO2 ~ RA*FishWeight^(RB) * exp(RQ*TempC) * ACT '
# Because my data were all recorded on sedentary fish, the ACT (activity multiplier) constant is simply 1, so I removed it from script below.

# 'Cor.ggd' is the same as MO2 in the code below, and is short for 'Corrected O2 consumption in grams oxygen per gram fish per day'. 'Corrected' in 
# the sense that I have subtracted background levels of bacterial oxygen consumption. 

# Create dataframe WITHOUT mean value observations
Resp2dat <- RespDatP[RespDatP$Period != "Mean",]
Resp2dat$Period <- as.factor(as.character(Resp2dat$Period)) # Clean up factor levels (remove "Mean" level)
Resp2dat <- Resp2dat[!is.na(Resp2dat$MO2.corrected),] # Most fish (all but 2) do not have a 4th trial period. Remove these missing Period 4s.

# Fit temperature (RQ) parameter only using weight (RA and RB) parameters from Generalized Coregonid (GC) model from Rudstam et al. 1994. 
f0 <- MO2.corrected ~ 0.0018*FishWeight.g^(-0.12)*exp(RQ*TempC) 
fit.gc <- nls(f0, data = Resp2dat, start = list(RQ= 0.047))
summary(fit.gc)
RQ.gc <- coefficients(fit.gc)[1]

plot(MO2.corrected ~ FishWeight.g, data = Resp2dat)
points(MO2.corrected ~ FishWeight.g, data = Resp2dat[Resp2dat$TempC==5,], col=4) 
points(MO2.corrected ~ FishWeight.g, data = Resp2dat[Resp2dat$TempC==10,], col=1)
points(MO2.corrected ~ FishWeight.g, data = Resp2dat[Resp2dat$TempC==15,], col=2)

curve(0.0018*x^(-0.12)*exp(RQ.gc*5),add = TRUE,lty=3, col=4) # GC model with my RQ value (5C)
curve(0.0018*x^(-0.12)*exp(RQ.gc*15),add = TRUE,lty=3, col=2) # GC model with my RQ value (15C)
curve(0.0018*x^(-0.12)*exp(RQ.gc*10),add = TRUE,lty=3, col = 1) # GC model with my RQ value (10C)


# Madenjian (mj) et al. 2006 RA value
f0 <- MO2.corrected ~ 0.00085*FishWeight.g^(-0.12)*exp(RQ*TempC)
fit.mj <- nls(f0, data = Resp2dat, start = list(RQ= 0.047))
summary(fit.mj)
RQ.mj <- coefficients(fit.mj)[1]

curve(0.00085*x^(-0.12)*exp(RQ.mj*5),add = TRUE,lty=2, col=4) # MJ model with my RQ value for the three temperatures
curve(0.00085*x^(-0.12)*exp(RQ.mj*15),add = TRUE,lty=2, col=2)
curve(0.00085*x^(-0.12)*exp(RQ.mj*10),add = TRUE,lty=2, col = 1)

# Karjalainen (kj) et al. 1997 RA value
f0 <- MO2.corrected ~0.006235*FishWeight.g^(-0.06258)*exp(RQ*TempC)
fit.kj <- nls(f0, data = Resp2dat, start = list(RQ= 0.0485))
summary(fit.kj)
RQ.kj <- coefficients(fit.kj)[1]

# Huuskonen (hk) et al. 1998 RA value
f0 <- MO2.corrected ~0.00584*FishWeight.g^(-0.05341)*exp(RQ*TempC)
fit.hk <- nls(f0, data = Resp2dat, start = list(RQ= 0.0506))
summary(fit.hk)
RQ.hk <- coefficients(fit.hk)[1]

# Fitting all three parameters to my data. Note lack of statistical significance of Weight coefficients (RA and RB).
f0 <- MO2.corrected ~ RA*FishWeight.g^(RB)*exp(RQ*TempC);
fit2 <- nls(f0, data = Resp2dat, start = list(RA = 0.001, RB = -0.01, RQ= 0.01))
summary(fit2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*5),add = TRUE,lty=1, col=4)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*15),add = TRUE,lty=1, col=2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*10),add = TRUE,lty=1, col=1)

# Mixed-effects approach ==============================

library(lme4)
library(car)

# Fitted to untransformed variables
f0 <- MO2.corrected ~ FishWeight.g + TempC + (1 | Period)
fit.lme <- lmer(f0, data=Resp2dat, REML = FALSE)
summary(fit.lme) # Variance of random effect is zero

# Again with log-transformed data
f0 <- log(MO2.corrected) ~ log(FishWeight.g) + TempC + (1 | Period)
fit.logme <- lmer(f0, data=Resp2dat, REML = FALSE)

fit.logme <- lmer(log(MO2.corrected) ~ log(FishWeight.g) + TempC + (log(FishWeight.g) | Period), data=Resp2dat, REML = TRUE)

Resp2dat %>% group_by(FishID) %>% summarize(number=n()) %>% View()
tapply(Resp2dat$Period,Resp2dat$FishID, FUN = length)


# Remove unnecessary 
f0 <- log(MO2.corrected) ~  (1 | FishID)
fit.logme <- lmer(f0, data=Resp2dat, REML = FALSE)


# REMOVE UNNECESSARY FACTOR LEVELS
Resp2dat.curry <- Resp2dat 
Resp2dat.curry$Period <- as.character(Resp2dat.curry$Period)
Resp2dat.curry$Period <- factor(Resp2dat.curry$Period)

f0 <- log(MO2.corrected) ~  (1 | Period)
fit.logme <- lmer(f0, data=Resp2dat.curry, REML = FALSE)



summary(fit.logme) # Variance of random effect is zero

ggplot(Resp2dat, aes(FishWeight.g, MO2.corrected, col=as.factor(TempC),pch=as.factor(TempC))) + geom_point(cex=3) +
  labs(xlab="Fish ID",ylab="mO2 (g/g/day)") + 
  geom_line(aes(y=predict(fit2), group=as.factor(TempC)),lty=1) +
  geom_line(aes(y=predict(fit.gc), group=as.factor(TempC)),lty=2) +
  geom_line(aes(y=predict(fit.lme), group=as.factor(TempC)),lty=3) +
  xlim(10,23)

#Another method, allowing extension of lines beyond points to see how curves behave for adult fish.
ggplot(Resp2dat, aes(FishWeight.g, MO2.corrected, xmax=500, col=as.factor(TempC),pch=as.factor(TempC))) + geom_point(cex=3) +
  labs(xlab="Fish ID",ylab="mO2 (g/g/day)") +
  stat_function(fun = function(x) coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*5),col=2) +
  stat_function(fun = function(x) coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*10),col=3) +
  stat_function(fun = function(x) coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*15),col=4) +
  stat_function(fun = function(x) 0.0018*x^(-0.12)*exp(coefficients(fit.gc)[1]*5),col=2,lty=2) +
  stat_function(fun = function(x) 0.0018*x^(-0.12)*exp(coefficients(fit.gc)[1]*10),col=3,lty=2) +
  stat_function(fun = function(x) 0.0018*x^(-0.12)*exp(coefficients(fit.gc)[1]*15),col=4,lty=2) +
  stat_function(fun = function(x) 0.00085*x^(-0.12)*exp(coefficients(fit.gc)[1]*5),col=2,lty=3) +
  stat_function(fun = function(x) 0.00085*x^(-0.12)*exp(coefficients(fit.gc)[1]*10),col=3,lty=3) +
  stat_function(fun = function(x) 0.00085*x^(-0.12)*exp(coefficients(fit.gc)[1]*15),col=4,lty=3)

# Creating a dataframe to streamline model coefficient selection. 'MyRQ' for each model fits an RQ parameter based on my data
#  to each existing model with RA and RB as constants:
parsResp <- data.frame("Model"=c("Madenjian","Rudstam","ThisStudy","Karjalainen","Huuskonen"),"RA"=NA,"RB"=NA,"RQ"=NA,"MyRQ"=NA)
parsResp[1,2:5] <- c(0.00085,-0.12,0.047,RQ.mj) # Add Madenjian et al.2006 coefficients
parsResp[2,2:5] <- c(0.0018,-0.12,0.047,RQ.gc) # Add Rudstam et al.1994 coefficients
parsResp[3,2:4] <- c(coefficients(fit2)[1],coefficients(fit2)[2],coefficients(fit2)[3]) # Add my coefficients (fit all three coefs to my data)
parsResp[4,2:5] <- c(0.006235,-0.06258,0.0485,RQ.kj) # Add Karjalainen et al.1997 coefficients
parsResp[5,2:5] <- c(0.00584,-0.05341,0.0506,RQ.hk) # Add Huuskonen et al.1998 coefficients

# all models (published RA,RB,RQ) at each temp.
temptemp <- 15 # pick one temperature treatment
ggplot(Resp2dat, aes(FishWeight.g, MO2.corrected, xmax=500,pch=as.factor(TempC))) + geom_point(cex=3) +
  labs(xlab="Fish ID",ylab="mO2 (g/g/day)") + 
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$RQ[1]*temptemp),col=2,lty=2) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$RQ[2]*temptemp),col=2,lty=3) +
  stat_function(fun = function(x) parsResp$RA[3]*x^parsResp$RB[3]*exp(parsResp$RQ[3]*temptemp),col=2,lty=1) +
  stat_function(fun = function(x) parsResp$RA[4]*x^parsResp$RB[4]*exp(parsResp$RQ[4]*temptemp),col=5,lty=4) +
  stat_function(fun = function(x) parsResp$RA[5]*x^parsResp$RB[5]*exp(parsResp$RQ[5]*temptemp),col=5,lty=5) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$MyRQ[1]*temptemp),col=1,lty=2) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$MyRQ[2]*temptemp),col=1,lty=3) +
  stat_function(fun = function(x) parsResp$RA[4]*x^parsResp$RB[4]*exp(parsResp$MyRQ[4]*temptemp),col=5,lty=4) +
  stat_function(fun = function(x) parsResp$RA[5]*x^parsResp$RB[5]*exp(parsResp$MyRQ[5]*temptemp),col=5,lty=5) +
  ggtitle(label=paste(temptemp, "degrees C model fits"))

#Huuskonen and Karjalainen models (turqouise) are obviously not good fits, so let's stop working with those.
# And since we only have 5 lines per temperature, let's plot them all at once again.

# CURRY: Black from from Madjengian biased low (at 5 and 10 deg C)

# Original Rudstam and 

ggplot(Resp2dat, aes(FishWeight.g, MO2.corrected,col=as.factor(TempC), xmax=100,pch=as.factor(TempC))) + geom_point(cex=3) +
  labs(xlab="Fish ID",ylab="mO2 (g/g/day)") + 
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$RQ[1]*5),col=1,lty=2) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$RQ[2]*5),col=2,lty=3) +
  stat_function(fun = function(x) parsResp$RA[3]*x^parsResp$RB[3]*exp(parsResp$RQ[3]*5),col=2,lty=1) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$MyRQ[1]*5),col=2,lty=4) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$MyRQ[2]*5),col=2,lty=5) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$RQ[1]*10),col=1,lty=2) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$RQ[2]*10),col=3,lty=3) +
  stat_function(fun = function(x) parsResp$RA[3]*x^parsResp$RB[3]*exp(parsResp$RQ[3]*10),col=3,lty=1) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$MyRQ[1]*10),col=3,lty=4) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$MyRQ[2]*10),col=3,lty=5) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$RQ[1]*15),col=4,lty=2) + #poor fit
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$RQ[2]*15),col=4,lty=3) + #poor fit
  stat_function(fun = function(x) parsResp$RA[3]*x^parsResp$RB[3]*exp(parsResp$RQ[3]*15),col=4,lty=1) +
  stat_function(fun = function(x) parsResp$RA[1]*x^parsResp$RB[1]*exp(parsResp$MyRQ[1]*15),col=4,lty=4) +
  stat_function(fun = function(x) parsResp$RA[2]*x^parsResp$RB[2]*exp(parsResp$MyRQ[2]*15),col=4,lty=5)
  

# At the 15 degree treatment, all fits look reasonable for my data except for the original gc(rudstam) and mj models, i.e., 
# all models with my fitted RQ value appear sensible.
# At the 10 degree treatment, all fits look reasonable except the original mj model.
# At the 5 degree treatment, all fits look somewhat reasonable, though Madenjian again appears to underestimate MO2.

AIC(fit2,fit.gc,fit.mj,fit.kj,fit.hk) #AIC confirms that Huuskonen and Karjalainen are inferior, and the other three
# are similar (Madenjian slightly inferior to the other two), though the model fitting all three coefficients to my 
# data (fit2) has 4 instead of 2 degrees of freedom.

# In all cases RA and RB parameters, and estimating third parameter



















# Experimental gibberish area. 2/19/20
f0 <- MO2.corrected ~ RA*FishWeight.g^(RB)*exp(RQ*5);
fit3 <- nls(f0, data = Resp2dat[Resp2dat$TempC==5,], start = list(RA = 0.001, RB = -0.01,RQ=0.01))
summary(fit3)

plot(MO2.corrected ~ FishWeight.g, data = Resp2dat[Resp2dat$TempC==5,])
new.weights <- data.frame(FishWeight.g = seq(1,28,by=1))
MO2.pred <- predict(fit3, newdata = new.weights$FishWeight.g)
par(new=T); plot(MO2.pred ~ new.weights$FishWeight.g, xlab="",ylab="",axes=F, xlim=c(10,30),ylim=c(0,0.007), type="l")
curve(coefficients(fit3)[1]*x^(coefficients(fit3)[2]),add = TRUE,lty=1, col=2)

