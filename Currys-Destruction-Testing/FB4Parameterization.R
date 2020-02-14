# Parameterization script: juvenile broad whitefish, Fish Bioenergetics 4.0
# Created by Duncan Green on 1/22/20

# setwd("/Users/duncangreen/Documents/GreenThesis") # current working directory
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
# but there are more options within existing FB4 models. These are just two examples. 

# Model form: ' MO2 ~ RA*FishWeight^(RB) * exp(RQ*TempC) * ACT '
# Because my data were all recorded on sedentary fish, the ACT (activity multiplier) constant is simply 1, so I removed it from script below.

# 'Cor.ggd' is the same as MO2 in the code below, and is short for 'Corrected O2 consumption in grams oxygen per gram fish per day'. 'Corrected' in 
# the sense that I have subtracted background levels of bacterial oxygen consumption. 

# Fit temperature (RQ) parameter only using weight (RA and RB) parameters from Generalized Coregonid (GC) model from Rudstam et al. 1994. 
f0 <- Cor.ggd ~ 0.0018*FishWeight^(-0.12)*exp(RQ*TempC) 
fit.gc <- nls(f0, data = fish3sub, start = list(RQ= 0.047))
summary(fit.gc)
RQ.gc <- coefficients(fit.gc)[1]


plot(Cor.ggd ~ FishWeight, data = fish3sub)
points(Cor.ggd ~ FishWeight, data = sub5C, col=1) 
points(Cor.ggd ~ FishWeight, data = sub10C, col=3)
points(Cor.ggd ~ FishWeight, data = sub15C, col=2)

curve(0.0018*x^(-0.12)*exp(RQ.gc*5),add = TRUE,lty=3, col=1) # GC model with my RQ value (5C)
curve(0.0018*x^(-0.12)*exp(RQ.gc*15),add = TRUE,lty=3, col=2) # GC model with my RQ value (15C)
curve(0.0018*x^(-0.12)*exp(RQ.gc*10),add = TRUE,lty=3, col = 3) # GC model with my RQ value (10C)


# Madenjian (mj) et al. 2006 RA value
f0 <- Cor.ggd ~ 0.00085*FishWeight^(-0.12)*exp(RQ*TempC)
fit.mj <- nls(f0, data = fish3sub, start = list(RQ= 0.047))
summary(fit.mj)
RQ.mj <- coefficients(fit.mj)[1]

curve(0.00085*x^(-0.12)*exp(RQ.mj*5),add = TRUE,lty=1, col=1) # MJ model with my RQ value for the three temperatures
curve(0.00085*x^(-0.12)*exp(RQ.mj*15),add = TRUE,lty=1, col=2)
curve(0.00085*x^(-0.12)*exp(RQ.mj*10),add = TRUE,lty=1, col = 3)


# Fitting all three parameters to my data. Note lack of statistical significance of Weight coefficients (RA and RB).
f0 <- Cor.ggd ~ RA*FishWeight^(RB)*exp(RQ*TempC);
fit2 <- nls(f0, data = fish3sub, start = list(RA = 0.001, RB = -0.01, RQ= 0.01))
summary(fit2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*5),add = TRUE,lty=2, col=1)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*15),add = TRUE,lty=2, col=2)
curve(coefficients(fit2)[1]*x^(coefficients(fit2)[2])*exp(coefficients(fit2)[3]*10),add = TRUE,lty=2, col=3)

# To sum up my issue again now that you've seen some data, I am trying to determine the 'best' way to parameterize these RA and RB coefficients
# despite having almost no variation in fish weights in my experimental data.


