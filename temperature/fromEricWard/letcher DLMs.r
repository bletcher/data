
load("e.RData")

par(mfrow=c(4,4),mai=c(0.3,0.3,0.1,0.1)) 


# For each year, look at ACF of differenced data 
acfWater = NA
acfAir = NA
for(i in min(e$year):max(e$year)) {
subE = subset(e, year==i)	
acfWater[i] = acf(diff(subE$waterTemp), na.action=na.pass)[1]
acfAir[i] = acf(diff(subE$airTemp), na.action=na.pass)[1]
}



#####################################################################################
# DLM # 1: this is just to illustrate time varying coefficients. Only applied to 
# water temperature series from 2002
#####################################################################################

# Create a bunch of variables for the DLM
subE = subset(e, year==2002)
airTemp = subE$airTemp
waterTemp = subE$waterTemp
N = length(airTemp)
D2 = diag(2)
Gt = diag(2)
Gt[1,2] = 1
Ft = matrix(c(1,0),1,2) 

dlmModel = cat("
model {
  # wishart prior on MVN for coefficient evolution
  coefTau ~ dwish(D2, 2);  	
  tauResid ~ dgamma(0.001,0.001);	
  # init priors on coef
  coef[1,1] ~ dnorm(0,1);
  coef[2,1] ~ dnorm(0,1);
  
  for(i in 2:N) {
  	predCoef[1:2,i] <- Gt%*%coef[1:2,i-1];
  	coef[1:2,i] ~ dmnorm(predCoef[1:2,i],coefTau[1:2,1:2]);
  	predY[i] <- inprod(Ft, coef[,i])
  	waterTemp[i] ~ dnorm(predY[i], tauResid);
  }
  
}

", file="model.txt")

library(runjags)
library(R2jags)
library(gtools)
library(gdata)

# These are all MCMC parameters...same as the ones we used
mcmc.chainLength <- as.integer(10000)  # burn-in plus post-burn
mcmc.burn <- as.integer(5000)
mcmc.thin = 1
mcmc.chains = 3      # needs to be at least 2 for DIC

jags.data = list("N","D2","Ft","Gt","airTemp","waterTemp")
jags.params=c("coef","tauResid","coefTau")
model.loc = paste("model.txt",sep="")

jags.model = jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.loc, n.chains = mcmc.chains, n.burnin = mcmc.burn, n.thin = mcmc.thin, n.iter = mcmc.chainLength, DIC = TRUE) 

attach.jags(jags.model)

par(mfrow = c(2,1), mai=c(0.7,0.7,0.3,0.2))
plot(apply(coef[,1,],2,mean), type="l",lwd=3)
lines(apply(coef[,1,],2,quantile,0.025), type="l",lwd=1)
lines(apply(coef[,1,],2,quantile,0.975), type="l",lwd=1)

plot(apply(coef[,2,],2,mean), type="l",lwd=3, ylim=c(-2,2))
lines(apply(coef[,2,],2,quantile,0.025), type="l",lwd=1)
lines(apply(coef[,2,],2,quantile,0.975), type="l",lwd=1)

# Make histogram of coef correlation - like regression, we expect to be negative
coefSigma = apply(coefTau,1,solve)
hist(coefSigma[2,]/(sqrt(coefSigma[1,]*coefSigma[4,])),100)


#####################################################################################
# DLM # 2: illustrate the covariate airTemp included. For simplicity, only the slope
# parameter is allowed to have a random walk, but the max() statement prevents predictions
# from dropping < 0
#####################################################################################

# Create a bunch of variables for the DLM
subE = subset(e, year==2002)
airTemp = subE$airTemp
waterTemp = subE$waterTemp
N = length(airTemp)
D2 = diag(2)
Gt = diag(2)
Gt[1,2] = 1
Ft = matrix(c(1,0),1,2) 

dlmModel2 = cat("
model {
  
  coefTau[1] ~ dgamma(0.001,0.001);
  coefTau[2] ~ dgamma(0.001,0.001);    	
  tauResid ~ dgamma(0.001,0.001);	
  # init priors on coef
  Intercept[1] ~ dnorm(0,1);
  slope[1] ~ dnorm(0,1);
  
  predY[1] <- Intercept[1];
  waterTemp[1] ~ dnorm(predY[1], tauResid);
  for(i in 2:N) {
  	predSlope[i] <- slope[i-1];
  	slope[i] ~ dnorm(predSlope[i], coefTau[1]);
  	#predIntercept[i] <- Intercept[i-1];
  	Intercept[i] <- Intercept[1]#~ dnorm(predIntercept[i], coefTau[2]);  	
  	predY[i] <- max(Intercept[i] + slope[i]*airTemp[i],0);
  	waterTemp[i] ~ dnorm(predY[i], tauResid);
  }
  residSS <- sum(pow(predY[2:N] - waterTemp[2:N],2))
  
}

", file="model2.txt")

library(runjags)
library(R2jags)
library(gtools)
library(gdata)

# These are all MCMC parameters...same as the ones we used
mcmc.chainLength <- as.integer(10000)  # burn-in plus post-burn
mcmc.burn <- as.integer(5000)
mcmc.thin = 1
mcmc.chains = 3      # needs to be at least 2 for DIC

jags.data = list("N","D2","Ft","Gt","airTemp","waterTemp")
jags.params=c("slope","tauResid","coefTau","predY","Intercept")
model.loc = paste("model2.txt",sep="")

jags.model = jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.loc, n.chains = mcmc.chains, n.burnin = mcmc.burn, n.thin = mcmc.thin, n.iter = mcmc.chainLength, DIC = TRUE) 

attach.jags(jags.model)

par(mfrow=c(2,1))
plot(airTemp, apply(slope,2,mean))

plot(apply(slope,2,mean), type="l",lwd=3, ylim=c(0,1),xlab="Day",ylab="Air Temp Effect")
lines(apply(slope,2,quantile,0.025), type="l",lwd=1)
lines(apply(slope,2,quantile,0.975), type="l",lwd=1)

plot(apply(Intercept,2,mean), type="l",lwd=3,xlab="Day",ylab="Air Temp Effect")
lines(apply(Intercept,2,quantile,0.025), type="l",lwd=1)
lines(apply(Intercept,2,quantile,0.975), type="l",lwd=1)

##
plot( apply(predY,2,mean))
lines(waterTemp/1)


#####################################################################################
# DLM # 3: unlike # 2, the random walks occur both in the intercept and in the 
# time varying coefficient of air temp. No max() is used.
#####################################################################################

# Create a bunch of variables for the DLM
subE = subset(e, year==2003)
airTemp = subE$airTemp
waterTemp = subE$waterTemp
N = length(airTemp)
D2 = diag(2)
Gt = diag(2)
Gt[1,2] = 1
Ft = matrix(c(1,0),1,2) 

dlmModel2 = cat("
model {
  
  coefTau[1] ~ dgamma(0.001,0.001);
  coefTau[2] ~ dgamma(0.001,0.001);    	
  tauResid ~ dgamma(0.001,0.001);	
  # init priors on coef
  Intercept[1] ~ dnorm(0,1);
  slope[1] ~ dnorm(0,1);
  
  predY[1] <- Intercept[1];
  waterTemp[1] ~ dnorm(predY[1], tauResid);
  for(i in 2:N) {
  	predSlope[i] <- slope[i-1];
  	slope[i] ~ dnorm(predSlope[i], coefTau[1]);
  	predIntercept[i] <- Intercept[i-1];
  	Intercept[i]~ dnorm(predIntercept[i], coefTau[2]);  	
  	predY[i] <- Intercept[i] + slope[i]*airTemp[i];
  	waterTemp[i] ~ dnorm(predY[i], tauResid);
  }
  residSS <- sum(pow(predY[2:N] - waterTemp[2:N],2))
  
}

", file="model3.txt")

library(runjags)
library(R2jags)
library(gtools)
library(gdata)

# These are all MCMC parameters...same as the ones we used
mcmc.chainLength <- as.integer(10000)  # burn-in plus post-burn
mcmc.burn <- as.integer(5000)
mcmc.thin = 1
mcmc.chains = 3      # needs to be at least 2 for DIC

jags.data = list("N","D2","Ft","Gt","airTemp","waterTemp")
jags.params=c("slope","tauResid","coefTau","predY","Intercept")
model.loc = paste("model3.txt",sep="")

jags.model = jags(jags.data, inits = NULL, parameters.to.save= jags.params, model.file=model.loc, n.chains = mcmc.chains, n.burnin = mcmc.burn, n.thin = mcmc.thin, n.iter = mcmc.chainLength, DIC = TRUE) 

attach.jags(jags.model)

par(mfrow=c(2,1))
plot(airTemp, apply(slope,2,mean))

plot(apply(slope,2,mean), type="l",lwd=3, ylim=c(0,1),xlab="Day",ylab="Air Temp Effect")
lines(apply(slope,2,quantile,0.025), type="l",lwd=1)
lines(apply(slope,2,quantile,0.975), type="l",lwd=1)
points(waterTemp/20)

plot(apply(Intercept,2,mean), type="l",lwd=3,xlab="Day",ylab="Air Temp Effect")
lines(apply(Intercept,2,quantile,0.025), type="l",lwd=1)
lines(apply(Intercept,2,quantile,0.975), type="l",lwd=1)
points(waterTemp/2)

##
plot( apply(predY,2,mean))
lines(waterTemp/2)


lines(waterTemp)
