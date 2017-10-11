## Hidden Markov Models for time series, Part 2

# Added to Github
# First load the requisite packages
rm(list=ls()) # Clear memory
graphics.off() # Clears graphics

# Install packages for zero-inflated hidden semi-Markov models
#install.packages("ziphsmm")
library(ziphsmm)

# Install packages for HMM models
#install.packages("depmixS4")
library(depmixS4)

# Install packages to calculate inverse logit
#install.packages("boot")
library(boot)

# Load data
fitbit <- read.csv("./Mike_Fitbit_data.csv", as.is = TRUE, header = TRUE)
fitbit$Date <- as.Date(fitbit$Date)

# Remove technical outliers where not wearing (either make = mean or zero)
# Create indicator based on time
fitbit$timetotal <- fitbit$Minutes.Sedentary + fitbit$Minutes.Lightly.Active + fitbit$Minutes.Fairly.Active + fitbit$Minutes.Very.Active + fitbit$Time.in.Bed
fitbit$perc.recorded <- fitbit$timetotal/(max(fitbit$timetotal, na.rm=TRUE))


#Remove fitbit recordings with less than 70% of time captured
fitbit <- subset(fitbit, fitbit$perc.recorded >= 0.7)

###############################################################################

#HMM Models for Fitbit Activity


# Graph and Histogram
graph_summary <- function(x, d, string) {
  fileName = paste0(string, "_summGraph.jpeg")
  jpeg(filename=fileName, width = 600, height = 400, quality=90)
  plot(d, x, type="c", col= "black",
       main= toupper(string), cex.main=2, #Size of Main font (cex)
       xlab="Date",
       ylab= string,
       cex.lab=1.4) #Size of label
  lines(lowess(d, x, f=1/6), col="green", lwd=2)
  segments(min(d), mean(x), max(d), mean(x), col="red", lwd=2)
  abline((linear=lm(x~d)), col="blue", lwd=2)
  legend(min(d), max(x), #Location
         c("Smoothed (Lowess)", paste("Average (=", round(mean(x),1), "min/day)"), paste("Linear Fit (Trend=", round(linear$coefficients[2], 1), "min/day)")), #Text
         col=c("green", "red", "blue"), #Line colors
         lty=c(1, 1, 1), #Line types
         lwd=c(2.5, 2.5, 2.5), #Line thickness
         bty="n", #No border ("o" if border)
         cex=1.3, #Text size
         y.intersp=0.85 #Spacing between text/lines
  )
  dev.off()
}
graph_summary(fitbit$Steps, fitbit$Date, "Steps")

fileName2 = "_histogram_Steps.jpeg"
jpeg(filename = fileName2, width = 600, height = 400, quality = 90)
hist(fitbit$Steps, n = 20, main = "Histogram of Daily Steps")
dev.off()


# 2-state model
model_dep2g <- depmix(response = Steps ~ 1, data = fitbit, nstates = 2, 
                      family = gaussian(),
                      trstart = runif(4))
fm_dep2g <- fit(model_dep2g, emc=em.control(rand=TRUE))
fm_dep2g
summary(fm_dep2g)
sumMat2g <- data.frame(summary(fm_dep2g))
colnames(sumMat2g) <- c("Mean", "SD")
sumMat2g

# Decode states
esttrans2g <- posterior(fm_dep2g)

ts.plot(esttrans2g$state)

# 3-state model
model_dep3g <- depmix(response = Steps ~ 1, data = fitbit, nstates = 3, 
                      family = gaussian(),
                      trstart = runif(9))
fm_dep3g <- fit(model_dep3g, emc=em.control(rand=TRUE))
fm_dep3g
summary(fm_dep3g)
sumMat3g <- data.frame(summary(fm_dep3g))
colnames(sumMat3g) <- c("Mean", "SD")
sumMat3g

# Decode states
esttrans3g <- posterior(fm_dep3g)

ts.plot(esttrans3g$state)

# 4-state model
model_dep4g <- depmix(response = Steps ~ 1, data = fitbit, nstates = 4, 
                      family = gaussian(),
                      trstart = runif(16))
fm_dep4g <- fit(model_dep4g, emc=em.control(rand=TRUE))
fm_dep4g
summary(fm_dep4g)
sumMat4g <- data.frame(summary(fm_dep4g))
colnames(sumMat4g) <- c("Mean", "SD")
sumMat4g

# Decode states
esttrans4g <- posterior(fm_dep4g)

ts.plot(esttrans4g$state)

# Compare with simple linear model
model_lin <- lm(fitbit$Steps~1)
AIC(model_lin)
BIC(model_lin)

fm_dep2g; fm_dep3g; fm_dep4g

# Plot best 2 models -- here 2 state
estrans2gMeans <- sumMat2g$Mean[esttrans2g$state]
graphData <- data.frame(cbind(as.Date(fitbit$Date), fitbit$Steps, estrans2gMeans, esttrans2g$state))
colnames(graphData) <- c("Date", "Steps", "StateMeans", "State")
graphData3$Date <- as.Date(fitbit$Date)

fileName = "./HMM_steps2.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
plot(graphData$Date, graphData$Steps, type = "l", main = "Hidden Markov Model for Activity",
     xlab = "Time", ylab = "Steps/day")
lines(graphData$Date, graphData$StateMeans, type = "o", col = "red")
abline(h = sumMat2g$Mean, col = "green", lty = 2)
legend("topleft", #Location
       c("Steps per day", "State transitions", "State means (= 12091 steps/day & 23252 steps/day)"), #Text
       col=c("black", "red", "green"), #Line colors
       lty=c(1, 1, 2), #Line types
       lwd=c(2.5, 2.5, 2.5), #Line thickness
       bty="n", #No border ("o" if border)
       cex=1.3, #Text size
       y.intersp=0.85 #Spacing between text/lines
)
dev.off()

estrans3gMeans <- sumMat3g$Mean[esttrans3g$state]
graphData3 <- data.frame(cbind(as.Date(fitbit$Date), fitbit$Steps, estrans3gMeans, esttrans3g$state))
colnames(graphData3) <- c("Date", "Steps", "StateMeans", "State")
graphData3$Date <- as.Date(fitbit$Date)

fileName = "./HMM_steps3.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
plot(graphData3$Date, graphData3$Steps, type = "l", main = "Hidden Markov Model for Activity",
     xlab = "Time", ylab = "Steps/day")
lines(graphData3$Date, graphData3$StateMeans, type = "o", col = "red")
abline(h = sumMat3g$Mean, col = "green", lty = 2)
legend("topleft", #Location
       c("Steps per day", "State transitions", "State means (= 9275 steps/day, 12817 steps/day, 14210 steps/day)"), #Text
       col=c("black", "red", "green"), #Line colors
       lty=c(1, 1, 2), #Line types
       lwd=c(2.5, 2.5, 2.5), #Line thickness
       bty="n", #No border ("o" if border)
       cex=1.3, #Text size
       y.intersp=0.85 #Spacing between text/lines
)
dev.off()

####################################################################
# Binary for 8-hour sleep (480 minutes a night)
fitbit$sleepGoal <- NA
fitbit$sleepGoal[fitbit$Minutes.Asleep >= 480] <- 1
fitbit$sleepGoal[fitbit$Minutes.Asleep < 480] <- 0

# Plot time series for sleep goal
fileName = "./TSPlot_sleepGoal.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
plot(fitbit$Date, fitbit$sleepGoal, type = "o", pch = 5,
     col = "blue",
     xlab = "Date", 
     ylab = "Goal reached",
     main = "Days reaching sleep goal",
     yaxt = "n")
axis(2, at = c(0,1))
dev.off()

# No states
model_bin1 <- glm(fitbit$sleepGoal~1, family = binomial())
summary(model_bin1)
inv.logit(model_bin1$coefficients)

# 2-state model
model_bin2g <- depmix(response = sleepGoal ~ 1, data = fitbit, nstates = 2, 
                      family = binomial(),
                      trstart = runif(4))
fm_bin2g <- fit(model_bin2g, emc=em.control(rand=TRUE))
fm_bin2g
summary(fm_bin2g)
sumMatbin2g <- data.frame(summary(fm_bin2g))
colnames(sumMatbin2g) <- c("logit")
sumMatbin2g$prob <- inv.logit(sumMatbin2g$logit)
sumMatbin2g


# 3-state model
model_bin3g <- depmix(response = sleepGoal ~ 1, data = fitbit, nstates = 3, 
                      family = binomial(),
                      trstart = runif(9))
fm_bin3g <- fit(model_bin3g, emc=em.control(rand=TRUE))
fm_bin3g
summary(fm_bin3g)
sumMatbin3g <- data.frame(summary(fm_bin3g))
colnames(sumMatbin3g) <- c("logit")
sumMatbin3g$prob <- inv.logit(sumMatbin3g$logit)
sumMatbin3g


# 4-state model
model_bin4g <- depmix(response = sleepGoal ~ 1, data = fitbit, nstates = 4, 
                      family = binomial(),
                      trstart = runif(16))
fm_bin4g <- fit(model_bin4g, emc=em.control(rand=TRUE))
fm_bin4g
summary(fm_bin4g)
sumMatbin4g <- data.frame(summary(fm_bin4g))
colnames(sumMatbin4g) <- c("logit")
sumMatbin4g$prob <- inv.logit(sumMatbin4g$logit)
sumMatbin4g

BIC(model_bin1); AIC(model_bin1); fm_bin2g; fm_bin3g; fm_bin4g

# Plot 2 state (best of HMMs)
# Decode
# Decode states
esttransbin2g <- posterior(fm_bin2g)

estransbin2gProbs <- sumMatbin2g$prob[esttransbin2g$state]
graphDataBin <- data.frame(cbind(as.Date(fitbit$Date), fitbit$sleepGoal, estransbin2gProbs, esttransbin2g$state))
colnames(graphDataBin) <- c("Date", "SleepGoal", "StateProb", "State")
graphDataBin$Date <- as.Date(fitbit$Date)

fileName = "./HMM_sleepGoal.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
plot(graphDataBin$Date, graphDataBin$SleepGoal, type = "o", main = "Hidden Markov Model for Sleep Goal",
     xlab = "Time", ylab = "Sleep Goal Reached")
lines(graphDataBin$Date, graphDataBin$StateProb, type = "o", col = "red")
abline(h = sumMatbin2g$prob, col = "green", lty = 2)
abline(h = inv.logit(model_bin1$coefficients), col = "purple", lty = 4)
legend("topleft", #Location
       c("Sleep Goal Y/N", "State transitions", "State mean probabilities", "Overall mean probability"), #Text
       col=c("black", "red", "green", "purple"), #Line colors
       lty=c(1, 1, 2, 4), #Line types
       lwd=c(2.5, 2.5, 2.5, 2.5), #Line thickness
       bty="n", #No border ("o" if border)
       cex=1.3, #Text size
       y.intersp=0.85 #Spacing between text/lines
)
dev.off()



#####################################################################
# HMM models
# Plot time series and histogram of exercise
fileName = "./TSPlot_exercise.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
plot(fitbit$Date, fitbit$Minutes.Very.Active, type = "l",
     col = "blue",
     xlab = "Date", 
     ylab = "Exericse minutes/day",
     main = "Daily Exercise")
dev.off()

fileName = "./Hist_exercise.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
hist(fitbit$Minutes.Very.Active, n=20, main= "Histogram: Exercise")
dev.off()


## 2 State zero inflated model
prior_init <- c(0.5, 0.5)
emit_init <- c(10, 100)
omega <- matrix(c(0.5,0.5,0.5,0.5),2,2,byrow=TRUE)
model_2s <- fasthmmfit(fitbit$Minutes.Very.Active, M = 2, 
                       prior_init = prior_init, 
                       tpm_init = omega, 
                       emit_init = emit_init,
                       zero_init = 0.5)
str(model_2s)
aic_2s = 2*length(model_2s$working_parm) + 2*model_2s$negloglik
bic_2s = log(length(fitbit$Minutes.Very.Active))*length(model_2s$working_parm) + 2*model_2s$negloglik


hmmviterbi(fitbit$Minutes.Very.Active, M = 2, 
           prior_init = prior_init, 
           tpm_init = omega, 
           emit_init = emit_init,
           zero_init = c(0.5, 0),
           plot = TRUE)

# 3-State zero inflated model
prior_init <- c(0.334, 0.333, 0.333)
emit_init <- c(10, 60, 200)
omega <- matrix(c(0.34, 0.33, 0.33, 0.34, 0.33, 0.33, 0.34, 0.33, 0.33),3,3,byrow=TRUE)
model_3s <- fasthmmfit(fitbit$Minutes.Very.Active, M = 3, 
                       prior_init = prior_init, 
                       tpm_init = omega, 
                       emit_init = emit_init,
                       zero_init = 0.5)
str(model_3s)
aic_3s = 2*length(model_3s$working_parm) + 2*model_3s$negloglik
bic_3s = log(length(fitbit$Minutes.Very.Active))*length(model_3s$working_parm) + 2*model_3s$negloglik


hmmviterbi(fitbit$Minutes.Very.Active, M = 3, 
           prior_init = prior_init, 
           tpm_init = omega, 
           emit_init = emit_init,
           zero_init = c(0.5, 0, 0),
           plot = TRUE)

# 4-State zero inflated model 
prior_init <- c(0.25, 0.25, 0.25, 0.25)
emit_init <- c(10, 20, 100, 200)
omega <- matrix(c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25),4,4,byrow=TRUE)
model_4s <- fasthmmfit(fitbit$Minutes.Very.Active, M = 4, 
                       prior_init = prior_init, 
                       tpm_init = omega, 
                       emit_init = emit_init,
                       zero_init = 0.5)
str(model_4s)
aic_4s = 2*length(model_4s$working_parm) + 2*model_4s$negloglik
bic_4s = log(length(fitbit$Minutes.Very.Active))*length(model_4s$working_parm) + 2*model_4s$negloglik

hmmviterbi(fitbit$Minutes.Very.Active, M = 4, 
           prior_init = prior_init, 
           tpm_init = omega, 
           emit_init = emit_init,
           zero_init = c(0.5, 0, 0, 0),
           plot = TRUE)

# 5-State zero inflated model
prior_init <- c(0.2, 0.2, 0.2, 0.2, 0.2)
emit_init <- c(10, 20, 50, 100, 200)
omega <- matrix(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),5,5,byrow=TRUE)
model_5s <- fasthmmfit(fitbit$Minutes.Very.Active, M = 5, 
                       prior_init = prior_init, 
                       tpm_init = omega, 
                       emit_init = emit_init,
                       zero_init = 0.5)
str(model_5s)
aic_5s = 2*length(model_5s$working_parm) + 2*model_5s$negloglik
bic_5s = log(length(fitbit$Minutes.Very.Active))*length(model_5s$working_parm) + 2*model_5s$negloglik

aic_2s; bic_2s; aic_3s; bic_3s; aic_4s; bic_4s; aic_5s; bic_5s

fileName = "./ZeroInfl_exercise_5state.jpeg"
jpeg(filename=fileName, width = 600, height = 400, quality=90)
hmmviterbi(fitbit$Minutes.Very.Active, M = 5, 
           prior_init = prior_init, 
           tpm_init = omega, 
           emit_init = emit_init,
           zero_init = c(0.5, 0, 0, 0, 0),
           plot = TRUE)
title("Zero-inflated 5-State HMM")
dev.off()

model_poisson <- glm(fitbit$Minutes.Very.Active ~ 1, family = poisson())
summary(model_poisson)
AIC(model_poisson)
BIC(model_poisson)
