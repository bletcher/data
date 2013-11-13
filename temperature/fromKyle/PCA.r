#This code predicts Mohseni Model parameters using a principal component analysis, mean parameters, and nearest neighbor parameters.
#The parameters are all then run through the model and plotted in a comparison with the calibration.

rm(list=ls())
library(pls)
library(MASS)

setwd("E:/USGS/Stream Temperature/BP_Analysis")
BasinData <- read.table("Basin_Chars_195.txt", header=TRUE)
	
#================================================================================================
#Post processing:
Mu <- BasinData$Mu
Alpha <- BasinData$Alpha
Theta <- BasinData$Theta
Beta <- BasinData$Beta

DA <- BasinData$Drainage_Area_Km2
Elev <- BasinData$Mean_Elevation_m
Slope <- BasinData$Mean_Slope_P
Aspect <- BasinData$Mean_Aspect
BFI <- BasinData$BFI_P
Clay <- BasinData$P_Clay
Sand <- BasinData$P_Sand
Silt <- BasinData$P_Silt
Biomass <- BasinData$Biomass_tons/BasinData$Drainage_Area_Km2
Stem <- BasinData$Stem_Count/BasinData$Drainage_Area_Km2
Urban <- BasinData$UrbanKm2/BasinData$Drainage_Area_Km2*100
Agro <- BasinData$AgroKm2/BasinData$Drainage_Area_Km2*100
Swamp <- BasinData$SwampsKm2/BasinData$Drainage_Area_Km2*100
Grass <- BasinData$GrassKm2/BasinData$Drainage_Area_Km2*100
Peat <- BasinData$PeatlandKm2/BasinData$Drainage_Area_Km2*100
Rock <- BasinData$RockKm2/BasinData$Drainage_Area_Km2*100
Decid <- BasinData$DeciduousKm2/BasinData$Drainage_Area_Km2*100
Oakpine <- BasinData$OakPineKm2/BasinData$Drainage_Area_Km2*100
Boreal <- BasinData$BorealKm2/BasinData$Drainage_Area_Km2*100
Water <- BasinData$OpenWaterKm2/BasinData$Drainage_Area_Km2*100
AugMeanAir <- BasinData$aug_mean_air_degC
AugMaxAir <- BasinData$aug_max_air_degC
AugPrcp <- BasinData$aug_mean_precip_mm	
FebMeanAir <- BasinData$feb_mean_air_degC	
FebMaxAir <- BasinData$feb_max_air_degC	
FebPrcp <- BasinData$feb_mean_precip_mm

#Bind em up:
chars <- data.frame(cbind(Mu, Alpha, Theta, Beta, DA, Elev, Slope, Aspect, BFI, Sand, Silt, Biomass, Stem, Urban, Agro, Swamp, Grass, Peat, Decid, Oakpine, Water, AugMeanAir, AugMaxAir, AugPrcp, FebPrcp, FebMeanAir, FebMaxAir))
#Need to remove values that sum to 1: (Excluded Clay, Boreal, Rock)

#================================================================================================
#                                 PRINCIPAL COMPONENT ANALYSIS
#================================================================================================
prin0.s <- prcomp(chars[,5:27],scale=TRUE,center=TRUE)
plot(prin0.s)
biplot(prin0.s, cex.lab = 2, cex.axis=2)
summary(prin0.s)
loadings(prin0.s)
scores(prin0.s)[,c(1:23)]

#-------------------------------
#Build linear models of the PCs:
#-------------------------------
MODEL_Mu <- lm(Mu ~ scores(prin0.s)) #PCs: 8
summary(MODEL_Mu)
mu.mod <- lm(Mu ~ scores(prin0.s)[,c(8)])
summary(mu.mod)

MODEL_Alpha <- lm(Alpha ~ scores(prin0.s)) #PCs: 1,3,4,5,10,12,18
summary(MODEL_Alpha)
alpha.mod <- lm(Alpha ~ scores(prin0.s)[,c(1,3,4)])#12 instead of 4?
summary(alpha.mod)

MODEL_Theta <- lm(Theta ~ scores(prin0.s)) #PCs: 1,2,3,5,8,9,11,12,14
summary(MODEL_Theta)
theta.mod <- lm(Theta ~ scores(prin0.s)[,c(1,2,3)])#,5,8,9,11,12,14)])
summary(theta.mod)

MODEL_Beta <- lm(Beta ~ scores(prin0.s)) #PCs: 1,3,4,7,14,19,20,22
summary(MODEL_Beta)
beta.mod <- lm(Beta ~ scores(prin0.s)[,c(1,3,4)])#,7,14,19,20,22)])
summary(beta.mod)



############################################################################################################################################################################################################
############################################################################################################################################################################################################
#====================================================================================================
#====================================================================================================
#====================================================================================================
#                         PREDICT PARAMETERS AT SITES IN THE PCA MODEL (CALIBRATION)
#====================================================================================================
#====================================================================================================
#====================================================================================================

#-------------------------------------------
#Set up an array for the prediction metrics:
#-------------------------------------------
pred_sites <-BasinData$Site
metrics.array <- array(NA,c(length(pred_sites),8))
colnames(metrics.array) <- c("Calibration_RMSE","Calibration_R2","PCA_RMSE", "PCA_R2","Mean_RMSE","MeanR2","Neighbor_RMSE","Neighbor_R2")
rownames(metrics.array) <-  pred_sites

#----------------------
#Load Temperature Data:
#----------------------
setwd("C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis")														#>>>>>>>> INPUT REQUIRED <<<<<<<<<
load("Breakpoint_Analysis_Input_Data_195.Rdata")	#temp_frame,site_list
temp_frame <- master.data


########################################## Add in if loop to reference MA or CT sites in BasinData

#----------------
#Mean Parameters:
#----------------
mu_mean <- mean(BasinData$Mu)#[1:46])
alpha_mean <- mean(BasinData$Alpha)#[1:46])
theta_mean <- mean(BasinData$Theta)#[1:46])
betas_mean <- mean(BasinData$Beta)#[1:46])

#------------------------------------------
#Borrow Some Parameters From The Neighbors:
#------------------------------------------
library(fossil)
	
Lat <- BasinData$Latitude#[1:46]
Lon <- BasinData$Longitude#[1:46]
sites <- BasinData$Site#[1:46]

coords <- array(NA,c(length(Lat),2))
coords[,1] <- Lon
coords[,2] <- Lat

distances <- earth.dist(coords, dist = TRUE)	#Distances are in km
dists <- as.matrix(distances)

mins <- apply(dists, 2, FUN = function(x) {min(x[x > 0])})

near.neighb <- array(NA,c(length(sites),2))
near.neighb[,1] <- as.character(sites)
colnames(near.neighb) <- c("Site of Interest","Nearest Site")

for (i in 1:length(pred_sites)){
near.neighb[i,2] <- as.character(sites[which(dists[i,] == mins[i])])
}


#Where in the MVR file are the sites of interest:
positions <- array(NA,c(length(pred_sites)))
for (i in 1:length(pred_sites)) {
positions[i] <- which(BasinData$Site == as.character(pred_sites[i]))
}

metrics.array <- array(NA,c(length(pred_sites),8))
colnames(metrics.array) <- c("Calibrated_RMSE", "Calibrated_R2", "PCA_RMSE", "PCA_R2", "Mean_RMSE", "Mean_R2", "Nearest_RMSE", "Nearest_R2")
rownames(metrics.array) <-  pred_sites

#====================================================================================================
#                                            PLOTTING
#====================================================================================================
#setwd("C:/Users/STUDENT/Documents/Stream Temp/MODEL_Mohseni et al")
setwd("C:/Users/STUDENT/Desktop")

#This loop automatically saves the 4x4 plot files. Enter the number of these files you wish to create for j (multiple of 15).
for(j in 1:11){		
png(filename=paste0("CT_Param_Compare_(",j,").png"),width=2175, height=1425, bg="white")															#>>>>>>>>> CHANGE NAME <<<<<<<<<<<

beg <- (j*15-14)
if((j*15) >= length(pred_sites)) fin <- length(pred_sites) else fin <- (j*15)

op<-par(mfrow=c(4,4),mar=c(4,5,4,2),oma=c(2,2,2,2))
for (i in beg:fin){		#length(pred_sites)

#Observed Air & Stream Temperatures:
Air.T <- temp_frame$airTemp[which(temp_frame$site == as.character(pred_sites[i]))]
Stream.T <- temp_frame$temp[which(temp_frame$site == as.character(pred_sites[i]))]
#Air.T <- met.data[,src_sites[i]+1,1]
#Stream.T <- temp.array[,src_sites[i]+1,1]

#Calibrated Mohseni Model Parameters:
mu <- Mu[positions[i]]
alpha <- Alpha[positions[i]]
theta <- Theta[positions[i]]
betas <- Beta[positions[i]]

stream.est <- mu+(alpha-mu)/(1+exp((4*tan(theta)/(alpha-mu))*(betas-Air.T)))


#Predicted PCA Parameters:
mu_pred <- predict(mu.mod)[positions[i]]
alpha_pred <- predict(alpha.mod)[positions[i]]
theta_pred <- predict(theta.mod)[positions[i]]
betas_pred <- predict(beta.mod)[positions[i]]

stream.pca <- mu_pred+(alpha_pred-mu_pred)/(1+exp((4*tan(theta_pred)/(alpha_pred-mu_pred))*(betas_pred-Air.T)))

#Mean Massachusetts Parameters:
stream.mean <- mu_mean+(alpha_mean-mu_mean)/(1+exp((4*tan(theta_mean)/(alpha_mean-mu_mean))*(betas_mean-Air.T)))


#Nearest Neighbor Parameters:

a <- near.neighb[which(near.neighb[,1] == pred_sites[i]),2]
near <- which(BasinData$Site == a)

mu_neighbor <- mean(BasinData$Mu[near])
alpha_neighbor <- mean(BasinData$Alpha[near])
theta_neighbor <- mean(BasinData$Theta[near])
betas_neighbor <- mean(BasinData$Beta[near])

stream.neighbor <- mu_neighbor+(alpha_neighbor-mu_neighbor)/(1+exp((4*tan(theta_neighbor)/(alpha_neighbor-mu_neighbor))*(betas_neighbor-Air.T)))

#Error Calcs:
notna <- which(is.na(Stream.T) == FALSE)


top1 = sum((Stream.T[notna]-stream.est[notna])^2)
RMSE1 <- sqrt(top1/(length(notna)))
R2_1 <- (cor(Stream.T[notna],stream.est[notna]))^2
metrics.array[i,1] <- RMSE1
metrics.array[i,2] <- R2_1

top2 = sum((Stream.T[notna]-stream.pca[notna])^2)
RMSE2 <- sqrt(top2/(length(notna)))
R2_2 <- (cor(Stream.T[notna],stream.pca[notna]))^2
metrics.array[i,3] <- RMSE2
metrics.array[i,4] <- R2_2

top3 = sum((Stream.T[notna]-stream.mean[notna])^2)
RMSE3 <- sqrt(top3/(length(notna)))
R2_3 <- (cor(Stream.T[notna],stream.mean[notna]))^2
metrics.array[i,5] <- RMSE3
metrics.array[i,6] <- R2_3

top4 = sum((Stream.T[notna]-stream.neighbor[notna])^2)
RMSE4 <- sqrt(top4/(length(notna)))
R2_4 <- (cor(Stream.T[notna],stream.neighbor[notna]))^2
metrics.array[i,7] <- RMSE4
metrics.array[i,8] <- R2_4

plot(Air.T,Stream.T, main = paste0(pred_sites[i]))
lines(sort(Air.T),sort(stream.pca),col="red", lwd=2)
lines(sort(Air.T),sort(stream.est),col="blue", lwd=2)
lines(sort(Air.T),sort(stream.mean),col="green", lwd=2)
lines(sort(Air.T),sort(stream.neighbor),col="yellow", lwd=2)
}

plot(0,0,pch=c(-1))
legend(-1, 1, legend=c("Observed","Mohseni Model","PCA Prediction","Mean MA Parameters","Neighbor's Parameters"), pch=c(1,-1,-1,-1,-1), pt.bg=c("white","white","white","white","white"), col=c("black","blue","red","green","yellow"),lty=c(-1,1,1,1,1), lwd=c(-1,2,2,2,2), ncol=1,cex=3)


dev.off()
}


#png(filename=paste0("Prediction Metrics Boxplot (Westbrook).png"),width=700, height=1050, bg="white")
a <- c("Calibrated", "PCA", "Mean ","Nearest Site")
op<-par(mfrow=c(2,1),mar=c(4,5,4,2),oma=c(2,2,2,2))
boxplot(metrics.array[,1],metrics.array[,3],metrics.array[,5],metrics.array[,7],ylab=expression(paste("RMSE ("^"o","C)",sep="")),names = a, main="160 Calibration Sites", cex=2,cex.lab=2,cex.axis=2)
boxplot(metrics.array[,2],metrics.array[,4],metrics.array[,6],metrics.array[,8],ylab=expression(paste("R"^"2",sep="")),names = a, cex=2,cex.lab=2,cex.axis=2)
#






cal.array <- metrics.array












#====================================================================================================
#                                            PLOTTING SELECT SITES
#====================================================================================================
setwd("C:/Users/STUDENT/Documents/Stream Temp/MODEL_Mohseni et al")


full_names <- c("","","", "Center Brook, MA", "Muddy Brook, MA")


op<-par(mfrow=c(1,2),mar=c(4,5,4,2),oma=c(2,2,2,2))
for (i in 4:5){		#length(pred_sites)

#Observed Air & Stream Temperatures:
Air.T <- met.data[,src_sites[i]+1,1]
Stream.T <- temp.array[,src_sites[i]+1,1]

#Calibrated Mohseni Model Parameters:
mu <- Mu[positions[i]]
alpha <- Alpha[positions[i]]
theta <- Theta[positions[i]]
betas <- Beta[positions[i]]

stream.est <- mu+(alpha-mu)/(1+exp((4*tan(theta)/(alpha-mu))*(betas-Air.T)))


#Predicted PCA Parameters:
mu_pred <- predict(mu.mod)[positions[i]]
alpha_pred <- predict(alpha.mod)[positions[i]]
theta_pred <- predict(theta.mod)[positions[i]]
betas_pred <- predict(beta.mod)[positions[i]]

stream.pca <- mu_pred+(alpha_pred-mu_pred)/(1+exp((4*tan(theta_pred)/(alpha_pred-mu_pred))*(betas_pred-Air.T)))

#Mean Massachusetts Parameters:
stream.mean <- mu_mean+(alpha_mean-mu_mean)/(1+exp((4*tan(theta_mean)/(alpha_mean-mu_mean))*(betas_mean-Air.T)))


#Nearest Neighbor Parameters:

a <- near.neighb[which(near.neighb[,1] == pred_sites[i]),2]
near <- which(BasinData$Site == a)

mu_neighbor <- mean(BasinData$Mu[near])
alpha_neighbor <- mean(BasinData$Alpha[near])
theta_neighbor <- mean(BasinData$Theta[near])
betas_neighbor <- mean(BasinData$Beta[near])

stream.neighbor <- mu_neighbor+(alpha_neighbor-mu_neighbor)/(1+exp((4*tan(theta_neighbor)/(alpha_neighbor-mu_neighbor))*(betas_neighbor-Air.T)))

#Error Calcs:
notna <- which(is.na(Stream.T) == FALSE)


top1 = sum((Stream.T[notna]-stream.est[notna])^2)
RMSE1 <- sqrt(top1/(length(notna)))
R2_1 <- (cor(Stream.T[notna],stream.est[notna]))^2
metrics.array[i,1] <- RMSE1
metrics.array[i,2] <- R2_1

top2 = sum((Stream.T[notna]-stream.pca[notna])^2)
RMSE2 <- sqrt(top2/(length(notna)))
R2_2 <- (cor(Stream.T[notna],stream.pca[notna]))^2
metrics.array[i,3] <- RMSE2
metrics.array[i,4] <- R2_2

top3 = sum((Stream.T[notna]-stream.mean[notna])^2)
RMSE3 <- sqrt(top3/(length(notna)))
R2_3 <- (cor(Stream.T[notna],stream.mean[notna]))^2
metrics.array[i,5] <- RMSE3
metrics.array[i,6] <- R2_3

top4 = sum((Stream.T[notna]-stream.neighbor[notna])^2)
RMSE4 <- sqrt(top4/(length(notna)))
R2_4 <- (cor(Stream.T[notna],stream.neighbor[notna]))^2
metrics.array[i,7] <- RMSE4
metrics.array[i,8] <- R2_4

plot(Air.T,Stream.T, main = paste0(full_names[i]), ylab="Stream Temperature (deg C)", xlab = "Stream Temperature (deg C)")
lines(sort(Air.T),sort(stream.pca),col="red", lwd=2)
lines(sort(Air.T),sort(stream.est),col="blue", lwd=2)
lines(sort(Air.T),sort(stream.mean),col="green", lwd=2)
lines(sort(Air.T),sort(stream.neighbor),col="yellow", lwd=2)


legend(-12, 20, legend=c("Observed","Calibrated Model","PCA Prediction","Mean MA Parameters","Neighbor's Parameters"), pch=c(1,-1,-1,-1,-1), pt.bg=c("white","white","white","white","white"), col=c("black","blue","red","green","yellow"),lty=c(-1,1,1,1,1), lwd=c(-1,2,2,2,2), ncol=1,cex=1)
}

plot(0,0,pch=c(-1))
legend(-1, 1, legend=c("Observed","Mohseni Model","PCA Prediction","Mean MA Parameters","Neighbor's Parameters"), pch=c(1,-1,-1,-1,-1), pt.bg=c("white","white","white","white","white"), col=c("black","blue","red","green","yellow"),lty=c(-1,1,1,1,1), lwd=c(-1,2,2,2,2), ncol=1,cex=3)


dev.off()
}













############################################################################################################################################################################################################
############################################################################################################################################################################################################
#====================================================================================================
#====================================================================================================
#====================================================================================================
#                               PREDICT PARAMETERS AT COMPLETELY NEW SITES (VERIFICATION)
#====================================================================================================
#====================================================================================================
#====================================================================================================
#Predict using this model at new sites:

setwd("C:/Users/STUDENT/Documents/Stream Temp/MODEL_Mohseni et al")
preds <- read.table("Pred_Sites_Phys_Chars.txt",header=TRUE)	#Run "Aggregate Basin Chars" for NH sites + any others to predict

#--------------------------------------
#Pre-processing of physical properties:
#--------------------------------------
pred_Mu <- preds $Mu
pred_Alpha <- preds $Alpha
pred_Theta <- preds $Theta
pred_Beta <- preds $Beta

pred_DA <- preds $Drainage_Area_Km2
pred_Elev <- preds $Mean_Elevation_m
pred_Slope <- preds $Mean_Slope_P
pred_Aspect <- preds $Mean_Aspect
pred_BFI <- preds $BFI_P
pred_Clay <- preds $P_Clay
pred_Sand <- preds $P_Sand
pred_Silt <- preds $P_Silt
pred_Biomass <- preds $Biomass_tons/preds $Drainage_Area_Km2
pred_Stem <- preds $Stem_Count/preds $Drainage_Area_Km2
pred_Urban <- preds $UrbanKm2/preds $Drainage_Area_Km2*100
pred_Agro <- preds $AgroKm2/preds $Drainage_Area_Km2*100
pred_Swamp <- preds $SwampsKm2/preds $Drainage_Area_Km2*100
pred_Grass <- preds $GrassKm2/preds $Drainage_Area_Km2*100
pred_Peat <- preds $PeatlandKm2/preds $Drainage_Area_Km2*100
pred_Rock <- preds $RockKm2/preds $Drainage_Area_Km2*100
pred_Decid <- preds $DeciduousKm2/preds $Drainage_Area_Km2*100
pred_Oakpine <- preds $OakPineKm2/preds $Drainage_Area_Km2*100
pred_Boreal <- preds $BorealKm2/preds $Drainage_Area_Km2*100
pred_Water <- preds $OpenWaterKm2/preds $Drainage_Area_Km2*100
pred_AugMeanAir <- preds $aug_mean_air_degC
pred_AugMaxAir <- preds $aug_max_air_degC
pred_AugPrcp <- preds $aug_mean_precip_mm	
pred_FebMeanAir <- preds $feb_mean_air_degC	
pred_FebMaxAir <- preds $feb_max_air_degC	
pred_FebPrcp <- preds $feb_mean_precip_mm

pred_chars <- data.frame(cbind(pred_DA, pred_Elev, pred_Slope, pred_Aspect, pred_BFI, pred_Sand, pred_Silt, pred_Biomass, pred_Stem, pred_Urban, pred_Agro, pred_Swamp, pred_Grass, pred_Peat, pred_Decid, pred_Oakpine, pred_Water, pred_AugMeanAir, pred_AugMaxAir, pred_AugPrcp, pred_FebMeanAir, pred_FebMaxAir, pred_FebPrcp))
#Excluded: Clay, Boreal, Rock

#--------------------------
#Set up for new parameters:
#--------------------------
num_sites <- length(pred_DA)

new_params <- array(NA, c(num_sites,4))
colnames(new_params) <- c("Alpha","Beta","Mu","Theta")

#-----------------------------------
#Standardizing the prediction chars:
#----------------------------------- 
#Get mean and sd:
pred.norm <- array(NA,c(2,23))
for (i in 1:23){
pred.norm[1,i] <- mean(pred_chars[,i])
pred.norm[2,i] <- sapply(pred_chars[i],sd)
}

#Loop through all sites:
pca_array <- array(NA, c(length(pred_DA), 23))
for (k in 1:23){
for (i in 1:length(pred_DA)){
a = 0
for (j in 1:23){
#Standardizes the basin characteristics:
x <- loadings(prin0.s)[j,k]*((pred_chars[i,j]-pred.norm[1,j])/pred.norm[2,j])
a = a + ifelse(is.nan(x),0,x)
}
pca_array[i,k] <- a
}
}

# j is to loop through the variables (basin chars)
# i is to loop through the number of sites
# k is to loop through the PCs

loadings(prin0.s)[,k]
pred_chars[i,]
pred.norm[1,]

#--------------------------------------------
#Predict each parameter using the PCA method: 
#--------------------------------------------

#Alpha
#-----
alpha.pcs <- c(1,3,4)
for (j in 1:num_sites){
for (i in 2:(length(alpha.mod$coef))){
a = 0 
a = a + pca_array[j,alpha.pcs[i-1]]*alpha.mod$coeff[i]
}
b <- alpha.mod$coeff[1] + a
new_params[j,1]<- b
}

#Beta
#----
beta.pcs <- c(1,3,4)
for (j in 1:num_sites){
for (i in 2:(length(beta.mod$coef))){
a = 0 
a = a + pca_array[j,beta.pcs[i-1]]*beta.mod$coeff[i]
}
b <- beta.mod$coeff[1] + a
new_params[j,2]<- b
}

#Mu
#--
mu.pcs <- c(8)
for (j in 1:num_sites){
for (i in 2:(length(mu.mod$coef))){
a = 0 
a = a + pca_array[j,mu.pcs[i-1]]*mu.mod$coeff[i]
}
b <- mu.mod$coeff[1] + a
new_params[j,3]<- b
}

#Theta
#-----
theta.pcs <- c(1,2,3)
for (j in 1:num_sites){
for (i in 2:(length(theta.mod$coef))){
a = 0 
a = a + pca_array[j,theta.pcs[i-1]]*theta.mod$coeff[i]
}
b <- theta.mod$coeff[1] + a
new_params[j,4]<- b
}

#-----------------------
#Mean Parameters Method:
#-----------------------
#Will need a bit of finagling to get mean of sites within a state:													# <-- ACTION NEEDED
mu_mean_pred <- mean(preds$Mu)
alpha_mean_pred <- mean(preds$Alpha)
theta_mean_pred <- mean(preds$Theta)
betas_mean_pred <- mean(preds$Beta)

#------------------------------------------
#Borrow Some Parameters From The Neighbors:
#------------------------------------------
library(fossil)
	
pred_Lat <- preds$Latitude
pred_Lon <- preds$Longitude
pred_near <- preds$Site

pred_coords <- array(NA,c(length(pred_Lat),2))
pred_coords[,1] <- pred_Lon
pred_coords[,2] <- pred_Lat

pred_distances <- earth.dist(pred_coords, dist = TRUE)	#Distances are in km
pred_dists <- as.matrix(pred_distances)

pred_mins <- apply(pred_dists, 2, FUN = function(x) {min(x[x > 0])})

pred_neighb <- array(NA,c(length(pred_near),2))
#pred_neighb[,1] <- pred_near
pred_neighb[,1] <- as.character(pred_near)

colnames(pred_neighb) <- c("Site of Interest","Nearest Site")

for (i in 1:length(pred_near)){
#pred_neighb[i,2] <- pred_near[which(pred_dists[i,] == pred_mins[i])]
pred_neighb[i,2] <- as.character(pred_near[which(pred_dists[i,] == pred_mins[i])])
}

#-------------------------------------------
#Set up an array for the prediction metrics:
#-------------------------------------------
pred_sites <- preds$Site
metrics.array <- array(NA,c(length(pred_sites),8))
colnames(metrics.array) <- c("Calibration_RMSE","Calibration_R2","PCA_RMSE", "PCA_R2","Mean_RMSE","MeanR2","Neighbor_RMSE","Neighbor_R2")
rownames(metrics.array) <-  pred_sites

#----------------------
#Load Temperature Data:
#----------------------
setwd("C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis")														#>>>>>>>> INPUT REQUIRED <<<<<<<<<
load("Breakpoint_Analysis_Input_Data_195.Rdata")	#temp_frame,site_list
temp_frame <- master.data

#setwd("C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data")
setwd("C:/Users/STUDENT/Desktop")


#This loop automatically saves the 4x4 plot files. Enter the number of these files you wish to create for j (multiple of 15).
for(j in 1:3){		
#png(filename=paste0("PCA_Predictions_(",j,").png"),width=2175, height=1425, bg="white")															#>>>>>>>>> CHANGE NAME <<<<<<<<<<<

beg <- (j*15-14)
if((j*15) >= length(pred_sites)) fin <- length(pred_sites) else fin <- (j*15)

#op<-par(mfrow=c(4,4),mar=c(4,5,4,2),oma=c(2,2,2,2))
for (i in beg:fin){		#length(pred_sites)

#i <- 3

Air.T <- temp_frame$airTemp[which(temp_frame$site == as.character(pred_sites[i]))]
Stream.T <- temp_frame$temp[which(temp_frame$site == as.character(pred_sites[i]))]


#Calibrated Mohseni Model Parameters:
#------------------------------------
mu_pred <- preds$Mu[i]
alpha_pred <- preds$Alpha[i]
theta_pred <- preds$Theta[i]
betas_pred <- preds$Beta[i]

stream.est_pred <- mu_pred+(alpha_pred-mu_pred)/(1+exp((4*tan(theta_pred)/(alpha_pred-mu_pred))*(betas_pred-Air.T)))


#Predicted PCA Parameters:
#-------------------------
mu_pca_pred <- new_params[i,3]
alpha_pca_pred <- new_params[i,1]
theta_pca_pred <- new_params[i,4]
betas_pca_pred <- new_params[i,2]

stream.pca_pred <- mu_pca_pred+(alpha_pca_pred-mu_pca_pred)/(1+exp((4*tan(theta_pca_pred)/(alpha_pca_pred-mu_pca_pred))*(betas_pca_pred-Air.T)))


#Mean Massachusetts Parameters:
#------------------------------
stream.mean_pred <- mu_mean_pred+(alpha_mean_pred-mu_mean_pred)/(1+exp((4*tan(theta_mean_pred)/(alpha_mean_pred-mu_mean_pred))*(betas_mean_pred-Air.T)))


#Nearest Neighbor Parameters:
#----------------------------
a <- pred_neighb[which(pred_neighb[,1] == pred_sites[i]),2]
near <- which(preds$Site == a)

mu_neighbor_pred <- mean(preds$Mu[near])
alpha_neighbor_pred <- mean(preds$Alpha[near])
theta_neighbor_pred <- mean(preds$Theta[near])
betas_neighbor_pred <- mean(preds$Beta[near])

stream.neighbor_pred <- mu_neighbor_pred+(alpha_neighbor_pred-mu_neighbor_pred)/(1+exp((4*tan(theta_neighbor_pred)/(alpha_neighbor_pred-mu_neighbor_pred))*(betas_neighbor_pred-Air.T)))

#Error Calcs:
#------------
notna <- which(is.na(Stream.T) == FALSE)

#Calibration Error:
#------------------
top1 = sum((Stream.T[notna]-stream.est_pred[notna])^2)
RMSE1 <- sqrt(top1/(length(notna)))
R2_1 <- (cor(Stream.T[notna],stream.est_pred[notna]))^2
metrics.array[i,1] <- RMSE1
metrics.array[i,2] <- R2_1

#PCA_Error:
#----------
top2 = sum((Stream.T[notna]-stream.pca_pred[notna])^2)
RMSE2 <- sqrt(top2/(length(notna)))
R2_2 <- (cor(Stream.T[notna],stream.pca_pred[notna]))^2
metrics.array[i,3] <- RMSE2
metrics.array[i,4] <- R2_2

#Mean Error:
#-----------
top3 = sum((Stream.T[notna]-stream.mean_pred[notna])^2)
RMSE3 <- sqrt(top3/(length(notna)))
R2_3 <- (cor(Stream.T[notna],stream.mean_pred[notna]))^2
metrics.array[i,5] <- RMSE3
metrics.array[i,6] <- R2_3

#Neighbor Error:
#---------------
top4 = sum((Stream.T[notna]-stream.neighbor_pred[notna])^2)
RMSE4 <- sqrt(top4/(length(notna)))
R2_4 <- (cor(Stream.T[notna],stream.neighbor_pred[notna]))^2
metrics.array[i,7] <- RMSE4
metrics.array[i,8] <- R2_4


#par(mfrow=c(1,1),mar=c(4,5,4,2),oma=c(2,2,2,2))
#plot(Air.T,Stream.T, main = paste0(pred_sites[i]), xlab=expression(paste("Stream Temperature ("^"o","C)",sep="")), ylab=expression(paste("Air Temperature ("^"o","C)",sep="")), cex.lab=2, cex.axis=2)
#lines(sort(Air.T),sort(stream.pca_pred),col="red", lwd=4)
#lines(sort(Air.T),sort(stream.est_pred),col="black", lwd=4)
#lines(sort(Air.T),sort(stream.mean_pred),col="green", lwd=4)
#lines(sort(Air.T),sort(stream.neighbor_pred),col="blue", lwd=4)
#legend(locator(1), legend=c("Observed","Mohseni Model","PCA Prediction","Mean Parameters","Nearest Site"), pch=c(1,-1,-1,-1,-1), pt.bg=c("white","white","white","white","white"), col=c("black","black","red","green","blue"),lty=c(-1,1,1,1,1), lwd=c(-1,4,4,4,4), ncol=1,cex=2)



}

#plot(0,0,pch=c(-1))
legend(-1, 1, legend=c("Observed","Mohseni Model","PCA Prediction","Mean Parameters","Neighbor's Parameters"), pch=c(1,-1,-1,-1,-1), pt.bg=c("white","white","white","white","white"), col=c("black","blue","red","green","yellow"),lty=c(-1,1,1,1,1), lwd=c(-1,2,2,2,2), ncol=1,cex=3)

#dev.off()
}




a <- c("Calibrated", "PCA", "Mean ","Nearest Site")

#setwd("C:/Users/STUDENT/Documents/Stream Temp/MODEL_Mohseni et al/PCA_Predictions")
#

#png(filename=paste0("Prediction Metrics Boxplot (Westbrook).png"),width=700, height=1050, bg="white")
op<-par(mfrow=c(2,1),mar=c(4,5,4,2),oma=c(2,2,2,2))
boxplot(metrics.array[,1],metrics.array[,3],metrics.array[,5],metrics.array[,7],ylab=expression(paste("RMSE ("^"o","C)",sep="")),names = a, main="35 Prediction Sites", cex=2,cex.lab=2,cex.axis=2)
boxplot(metrics.array[,2],metrics.array[,4],metrics.array[,6],metrics.array[,8],ylab=expression(paste("R"^"2",sep="")),names = a, cex=2,cex.lab=2,cex.axis=2)
#





#dev.off()





#png(filename=paste0("Prediction Metrics Boxplot (Westbrook).png"),width=700, height=1050, bg="white")
a <- c("Calibrated", "PCA", "Mean ","Nearest Site")
op<-par(mfrow=c(2,2),mar=c(1,5,4,2),oma=c(1,2,2,2))
boxplot(cal.array[,1],cal.array[,3],cal.array[,5],cal.array[,7],ylab=expression(paste("RMSE ("^"o","C)",sep="")), main="Calibration", cex=2, cex.lab=2, cex.axis=2, cex.main = 2)
boxplot(metrics.array[,1],metrics.array[,3],metrics.array[,5],metrics.array[,7], main="Validation", cex=2,cex.lab=2, cex.axis=2, cex.main = 2)
boxplot(cal.array[,2],cal.array[,4],cal.array[,6],cal.array[,8],ylab=expression(paste("R"^"2",sep="")), main="Calibration", names = a, cex=2, cex.lab=2, cex.axis=2, cex.main = 2)
boxplot(metrics.array[,2],metrics.array[,4],metrics.array[,6],metrics.array[,8], main="Validation", names = a, cex=2, cex.lab=2, cex.axis=2, cex.main = 2)

#



