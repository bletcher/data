# to do
# fill in all temp and airTemp data so we can get segmented to work
# compare reuslts with segements based on movingSD to previous segments (temp, season)

# Set up years as july-july.
# use of max vs min airTemp in addition to mean air temp
# need ot estimate summerBP for incomplete years
rm(list=ls())

library(ggplot2)
library(reshape2)
library(mgcv)
library(nlme)
library(plyr)
library(segmented)
library(zoo)

setwd('C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis')
load("Breakpoint_Analysis_Data_ALL_BEST_SITES.RData")	#master.data

site_list <- unique(master.data$site)

#For checking the order of e:
count <- seq(from=1, to=length(master.data$year), by=1)
e <- master.data
e$count <- count

names(e) <- c('date','dOY','year','site','temp','airTemp','rain','count')

e$season <- ifelse( e$dOY<80,1,
            ifelse( e$dOY<172,2,
            ifelse( e$dOY<264,3, 
            ifelse( e$dOY<355,4,1  ))))
e$season <- ifelse( e$dOY>355,1,e$season )

# this is the key metric for estimating the synchrony between air and water temp 
e$tempIndex <- (e$temp-e$airTemp)/e$temp

#############################################
#For temporarily looking at individual sites:
	#e <- e[which(e$site == 1805),]
#############################################

#===============================================================================
#Get moving mean and SD of temp index for each site and put into the data frame:
#===============================================================================
L.sites <- length(site_list)

movingMean <- -9999
movingSD <- -9999

for (sitex in 1:L.sites){

cur_site <- which(e$site == site_list[sitex])

window <- 10

cur_mean <-  rollapply(e$tempIndex[cur_site], width=window, fill=NA, mean)
cur_SD <-  rollapply(e$tempIndex[cur_site], width=window, fill=NA, sd)

movingMean <- append(movingMean, cur_mean)
movingSD <- append(movingSD, cur_SD)
}

movingMean <- movingMean[2:length(movingMean)]
movingSD <- movingSD[2:length(movingSD)]

e <- cbind(e,movingMean,movingSD)

e$meanSDDiff <- e$movingSD - e$movingMean

# just to make sure the merge doens't skrew up order
e <- e[order(e$count),]

#====================================================================================================================
#Segmented regression of water temperature to get breakpoint between ascending and descending temps by year and site:
#====================================================================================================================
start.yr <- min(e$year)
end.yr <- max(e$year)
por <- end.yr-start.yr+1

breaks <- data.frame(array(NA,c(por*L.sites,7)))
names(breaks) <- c('site','year','springBP','summerBP','fallBP','quantile_lo','quantile_hi')

for (j in 1:L.sites){

cur_site <- which(e$site == site_list[j])
e1 <- e[cur_site,]

#are data complete for a year?
complete <- as.data.frame(table( e1$year,is.na(e1$temp)))
incompleteYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq<365]))
completeYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq>=365]))

i=0
	for (year in c(start.yr:end.yr)){
		#specify location in stoarage array:
		i=i+1
		position <- i+(j-1)*por
		
		breaks$year[position] <- as.numeric(year)
		breaks$site[position] <- as.character(site_list[j])

		if(year %in% completeYears){

			out <- segmented( lm(temp~dOY, data=e1[e1$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
			breaks$summerBP[position] <- summary(out)$psi[2,2]
		}
	}
}

temp.breaks<-breaks

#===================================================================================
#Use segmented regression of the movingMeanSD to define fall and winter breakpoints:
#===================================================================================
for (j in 1:L.sites){

	cur_site <- which(e$site == site_list[j])
	e1 <- e[cur_site,]

	#are data complete for a year?
	complete <- as.data.frame(table( e1$year,is.na(e1$temp)))
	incompleteYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq<365]))
	completeYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq>=365]))

	#quantiles <- data.frame(year=c(min(e1$year):max(e1$year)))
	#quantiles$lo <- NA
	#quantiles$hi <- NA
	numForward <- 7 + 1  
	
	k <- 0
	for (year in c(start.yr:end.yr)){ 
		#print(year)

		k=k+1
		position <- k+(j-1)*por

		# Use predetermined breakpoints and add/subtract to find range to feed TIQ
		TIQ <- quantile(e1[e1$year  %in% year & e1$dOY %in% 125:275,'tempIndex'], probs=c(0.005,0.5,0.995),na.rm=T)
		movingSDQ <- quantile(e1[e1$year  %in% year & e1$dOY %in% 125:275,'movingSD'], probs=c(0.005,0.5,0.995),na.rm=T)

		#quantiles$lo[ quantiles$year == year ] <- TIQ[1]
		#quantiles$hi[ quantiles$year == year ] <- TIQ[3]
		breaks$quantile_lo[position] <- TIQ[1]
		breaks$quantile_hi[position] <- TIQ[3]
		
		runs <- data.frame(array(NA,c(1,numForward)))
		eYear <- e1[e1$year == year, ] 

		lmOut <- data.frame(array(NA,c(1,9)))
		winFactor <- 2


		if(year %in% completeYears){
			
			out <- segmented( lm(temp~dOY, data=e1[e1$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
			#tryCatch({out <- segmented( lm(temp~dOY, data=e1[e1$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

			# time forward until approx breakpoint ascending/descending water temp (in 'breaks')
			for (i in 1:(200)){ #-window*winFactor )){
				for (ii in 2:numForward ){
				runs[ i,ii ] <- 1*((eYear$movingMean[ i+ii-2 ] >= TIQ[1]) & (eYear$movingMean[ i+ii-2 ] <= TIQ[3]))
				runs[ i,1 ] <- prod( runs[ i, 2:numForward ] )
				}
			}
	
			# Make this not arbitrary, FOOL! use water temp break point
			# time backward until approx breakpoint ascending/descending water temp (in 'breaks')
			for (i in 366:(201)){ #-window*winFactor )){
				for (ii in 2:numForward ){
					runs[ i,ii ] <- 1*((eYear$movingMean[ i-ii+2 ] >= TIQ[1]) & (eYear$movingMean[ i-ii+2 ] <= TIQ[3]))
					runs[ i,1 ] <- prod( runs[ i, 2:numForward ] )
				}
			}
			
			#slopeQ <- quantile(lmOut[150:250,c(4)], probs=c(0.005,0.5,0.995))
			#diffQ <- quantile(lmOut[150:250,c(9)], probs=c(0.005,0.5,0.995))

			breaks$springBP[position] <- min(which(runs[,1] == 1))
			breaks$fallBP[position] <- max(which(runs[,1] == 1))
			
		} #completeYears loop
	} #year loop
} #site loop

#----------------------------------------------------------------
#Merge the breakpoints and quantiles with the existing dataframe:
#----------------------------------------------------------------
et <- merge( x=e, y=breaks, by.x=c('year','site'),all.x=T)

#Split up the year into 2 parts: (rising and falling temps)
et$segWaterTemp <- ifelse( et$dOY < et$summerBP,1,2 )

# check order:
et <- et[order(et$count),]

#-----------------
#Save the results:
#-----------------
temp_frame <- et
setwd('C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis')
save(temp_frame,site_list,file="Breakpoint_Data.Rdata")


###########################################################################################################################
# Look at some results:
###########################################################################################################################
library(ggplot2)

# Look at some results:
setwd('C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis')
load("Breakpoint_Data.Rdata")

#Distribution of breakpoints:
boxplot(temp_frame$springBP,temp_frame$summerBP,temp_frame$fallBP,names=c("Spring","Summer","Fall"))

#Select a Site to look at:
#-------------------------
#for (N in 1:195){			<--- Can't figure out why a simple for loop doesn't work here. Just produces blank files...

N <- 34 # <- which site to view.

cur <- which(temp_frame$site == site_list[N] & !is.na(temp_frame$temp))
temp <- temp_frame[cur,]

setwd('C:/Users/STUDENT/Documents/Stream Temp/Stream Temp Data/Breakpoint_Analysis')
# graph for breakpoints

png(filename=paste0("BP_Plot_Site_",temp$site[1],".png"),width=1000, height=400, bg="white")
ggplot( temp, aes(dOY,tempIndex)) +
  geom_point() +
  geom_point(aes(dOY,movingSD), colour='red') +
  geom_point(aes(dOY,movingMean), colour='blue') +
  geom_hline( aes(yintercept=as.numeric(quantile_lo)), colour='black') +
  geom_hline( aes(yintercept=as.numeric(quantile_hi)), colour='black') +
  geom_vline( aes(xintercept=as.numeric(springBP)), colour='green') +
  geom_vline( aes(xintercept=as.numeric(fallBP)), colour='orange') +
  geom_vline( aes(xintercept=as.numeric(summerBP)), colour='blue') +
  ylim(c(-10,10))+
 #  xlim(c(80,120))+
  ggtitle(paste(temp$site[1],sep=" ")) +
  facet_wrap(~year)
dev.off()

#}
 

 
 
 

#	"Error: No layers in plot"
ggplot(breaks, aes(year,springBP),colour='green' +
  geom_point(aes(year,fallBP),colour='red')+
  geom_point(aes(year,summerBP),colour='blue'))

 

 
 
 
 
 
 

 
 
 
 
 
 
 
 
 
#====================================================================================================================
#Get summer breakpoint across all data at a site (not split up by years):
#====================================================================================================================
start.yr <- min(e$year)
end.yr <- max(e$year)
por <- end.yr-start.yr+1

summer_BPs <- data.frame(array(NA,c(L.sites,2)))
names(summer_BPs) <- c('site','summerBP')

for (j in 1:L.sites){

cur_site <- which(e$site == site_list[j])
e1 <- e[cur_site,]

#are data complete for a year?
#complete <- as.data.frame(table( e1$year,is.na(e1$temp)))
#incompleteYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq<365]))
#completeYears <- as.numeric(as.character(complete$Var1[complete$Var2 == 'FALSE' & complete$Freq>=365]))

#num.yrs <- length(completeYears)

#i=0
	#for (year in c(start.yr:end.yr)){
	#for (year in c(completeYears[1]:completeYears[num.yrs])){
		#specify location in stoarage array:
		#i=i+1
		#position <- i+(j-1)*por
		
		#summer_BPs$site[position] <- as.numeric(year)
		summer_BPs$site[j] <- as.character(site_list[j])

		#if(year %in% completeYears){

			tryCatch({out <- segmented( lm(temp~dOY, data=e1), seg.Z=~dOY, psi=list(dOY=c(100,200)))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			summer_BPs$summerBP[j] <- summary(out)$psi[2,2]
		#}
	#}
}

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#                                                           END KO Stuff
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################  
#####################################################################################################################################
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
##############
# graphs for checking lm and breakpoints

ggplot( e[e$year  %in% year ,], aes(dOY,tempIndex)) +
  geom_point( aes(dOY,lmOut[dOY,3]), colour='red' ) +
  geom_line( aes(dOY,lmOut[dOY,3]), colour='red' ) +
  geom_point( aes(dOY,lmOut[dOY,4]), colour='green' ) +
 # geom_point( aes(dOY,lmOut[dOY,4]), colour='blue' ) +
  geom_point( aes(dOY,lmOut[dOY,8]), colour='orange' ) +
  geom_point( aes(dOY,lmOut[dOY,9]), colour='blue' ) +
  geom_point() +
  geom_line()+
  geom_hline( yintercept=slopeQ[1], colour='green') +
  geom_hline( yintercept=slopeQ[3], colour='green') +
  geom_hline( yintercept=diffQ[3], colour='blue') +
  geom_hline( yintercept=diffQ[1], colour='blue') +
  geom_hline( yintercept=TIQ[3], colour='black') +
  geom_hline( yintercept=TIQ[1], colour='black') +
  geom_vline( xintercept=firstBP, colour='blue') +
  geom_vline( xintercept=lastBP, colour='blue') +
  ylim(c(-10,10))+
 # xlim(c(80,220))+
  ggtitle(paste(year,window,winFactor,sep=" "))

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#############
# try strucchange. doesn't connect regression lines
#this doesn't seem to line up with visual inspection

year=2004
holdData <- ts(e[e$year == year,'movingSD'], start=1, end=366)
hold <- breakpoints( holdData ~ 1 )
bp <- data.frame(breakpoints(hold)[1] )

#ggplot( e[e$year  %in% year ,], aes(dOY,tempIndex)) +
ggplot( e[e$year  %in% year ,], aes(dOY,movingSD)) +
  geom_point() +
  geom_vline( aes( xintercept =  breakpoints  ),bp)
#


# loop over years
breaksSD <- data.frame(array(NA,c(20,3)))
i=0
for (year in c(1998:1999,2001:2007,2010)){
  i=i+1
  print(c(i,year))
  out <- segmented( lm(tempIndex~dOY, data=e[e$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
  # out <- segmented( lm(movingSD~dOY, data=e[e$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
  
  breaksSD[i,1] <- year
  breaksSD[i,2] <- summary(out)$psi[1,2]
  breaksSD[i,3] <- summary(out)$psi[2,2]
}


##################################################################

#################################################################
# segmented regression to find break points for waterTemp for winter, ascending, and descending across all years
# NOT USED
segE <- segmented( lm(temp~dOY, data=e[e$year %in% c(1998:2007,2010),]), seg.Z=~dOY, psi=list(dOY=c(100,200)))

plot(segE)
points(e$dOY,e$temp)
slope(segE)

break1E <- summary(segE)$psi[1,2]
break2E <- summary(segE)$psi[2,2]

e$segE <- ifelse( e$dOY < break1E,1,
                  ifelse( e$dOY < break2E,2,3 ))
###########################################################################################################
# by year
breaks <- data.frame(array(NA,c(20,3)))
i=0
for (year in c(1998:2007,2010)){
  i=i+1
  out <- segmented( lm(temp~dOY, data=e[e$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
  # out <- segmented( lm(movingSD~dOY, data=e[e$year == year,]), seg.Z=~dOY, psi=list(dOY=c(100,200)))
  
  breaks[i,1] <- year
  breaks[i,2] <- summary(out)$psi[1,2]
  breaks[i,3] <- summary(out)$psi[2,2]
}
# 'b2 is the break between the ascending and descending limbs
names(breaks) <- c('year','b1','b2')

e <- merge( x=e, y=breaks, by='year',all.x=T)
# just to make sure the merge doens't skrew up order
e <- e[order(e$date),]


e$segByYear <- ifelse( e$dOY < e$b1,1,
                       ifelse( e$dOY < e$b2,2,3 ))

ggplot(breaks, aes(year,b1)) +
  geom_point()+
  geom_point(aes(year,b2))

############################
# plot raw data
ggplot( e, aes(dOY,temp))+
  geom_point() +
  #geom_point( aes(dOY,airTemp), colour='red') +
  geom_point( aes(dOY,airTemp), colour='red') +
  facet_wrap(~year)

ggplot( e[e$year  %in% c(1998:2007,2010) ,], aes(dOY,temp/airTemp))+
  geom_point() +
  geom_line() +
  ylim(c(0,2)) +
  facet_wrap(~year)

ggplot( e, aes(studyDay,tempIndex,colour=factor(year)))+
  geom_point() +
  geom_line() +
  geom_line( aes( studyDay,movingMean ),colour='black')+
  geom_line( aes( studyDay,movingSD ),colour='blue')+
  ylim(c(-2,15)) +
  xlim(c(0,1000))
 
ggplot( e[e$year<2011,], aes(dOY,tempIndex,colour=factor(year)))+
 # geom_point() +
#  geom_line() +
 # geom_line( aes( dOY,movingMean ),colour='red')+
  geom_line( aes( dOY,movingSD ),colour='blue')+
  ylim(c(-2,50)) +
#  xlim(c(1500,2500)) +
#  xlim(c(0,4100)) +
  facet_wrap(~year)
 # scale_y_log10()

eTemp <- e[e$year %in% 2000:2001,]
eT <- eTemp$temp/eTemp$airTemp
eTS <- ts(data=eT)
bp <- breakpoints(eTS~1,breaks=10)

ggplot( e, aes(airTemp,temp) ) +
  geom_point() +
  facet_wrap(~year)

ggplot( e[e$year<2011,], aes(airTemp,temp, colour=factor(season)) ) +
  geom_point(  )+
  geom_smooth()+
  facet_wrap(~year)

ggplot( e[e$year<2011,], aes(airTemp,temp, colour=factor(segE)) ) +
  geom_point(  )+
  geom_smooth(method='lm')+
  facet_wrap(~year)

ggplot( e[e$year<2011,], aes(airTemp,temp, colour=factor(segByYear)) ) +
  geom_point(  )+
  geom_smooth(method='lm')+
  facet_wrap(~year)

ggplot( e[e$year<2011,], aes(tAirMax,temp, colour=factor(segE)) ) +
  geom_point(  )+
  geom_smooth(method='lm')+
  facet_wrap(~year)

ggplot( e[e$year<2011,], aes(dOY,temp))+#, colour=factor(season)) ) +
  geom_point(  )+
  geom_smooth(method='gam', formula= y ~ s( x, bs = "cc", k=15 ))+
  facet_wrap(~year)
########################################################################
# linear models by segByYear

# run lin mod with gls so we can compare other gls models
e1 <- gls( temp ~  airTemp, data=e )
e2 <- gls( temp ~  airTemp * factor(segByYear) * factor(year), na.action='na.exclude',data=e[e$year %in% c(1998:2007,2010),] )

acf(resid(e2))
AIC(ee1,ee2)
#plot(ee2)

ggplot( e[e$year  %in% c(1998:2007,2010) ,], aes(airTemp,temp, colour=factor(segByYear))) +
  geom_point() +
#  geom_line( aes(airTemp,predict(e2)) ) +
  facet_grid(year~seg)

########################################################################
#start with just one year
ee <- e[e$year==2003,]

# try segmented regression to identify breakpoints for winter, increasing, decreasing temps
lmEE <- lm(temp~dOY,data=ee)
summary(lmEE)

seg <- segmented( lmEE, seg.Z=~dOY, psi=list(dOY=c(100,200)))
plot(seg)
points(ee$dOY,ee$temp)
slope(seg)

break1 <- summary(seg)$psi[1,2]
break2 <- summary(seg)$psi[2,2]

ee$seg <- ifelse( ee$dOY < break1,1,
                  ifelse( ee$dOY < break2,2,3 ))

ee <- cbind(ee, ddply( ee, .(seg),summarize, tempScaled=scale(temp),airTempScaled=scale(airTemp))[,c(2,3)] )

##########################


ggplot( ee, aes(dOY,temp, colour=factor(season)))+
  geom_point()

ggplot( ee, aes(airTemp,temp, colour=(dOY)) ) +
  geom_point(  )+  
  geom_smooth()+
  facet_wrap(~season)

ggplot( ee, aes(airTemp,temp, colour=factor(season)) ) +
  geom_point(  )+
  geom_smooth()


ggplot( ee, aes(dOY,temp, colour=(dOY)) ) +
  geom_point(  )+
  geom_line(  )+
  geom_smooth()+
  geom_point( aes(dOY,airTemp), colour='red') +geom_line( aes(dOY,airTemp), colour='red') +
  facet_wrap(~season)

ggplot( ee, aes(dOY,temp-airTemp, colour=(dOY)) ) +
  geom_point(  )+
  geom_line(  )+
  geom_smooth() +
  ylim(c(-5,5)) +
  #geom_point( aes(dOY,airTemp), colour='red') +geom_line( aes(dOY,airTemp), colour='red') +
  #facet_wrap(~season)

ggplot( ee, aes(airTemp,temp, colour=factor(seg)) ) +
  geom_point()+
  geom_smooth(method='lm')

ggplot( ee, aes(airTempScaled,tempScaled, colour=factor(seg)) ) +
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~seg)

##############################################################################
# model for each seg, with and without temporal autocorrelation
ee0 <- lm( temp ~  airTemp * factor(seg), data=ee )
plot(ee1)

ggplot( ee, aes(airTemp,temp, colour=factor(seg))) +
  geom_point() +
  geom_line( aes(airTemp,predict(ee1)) ) +
  facet_wrap(~seg)

# run lin mod with gls so we can compare other gls models
ee1 <- gls( temp ~  airTemp, data=ee )
ee2 <- gls( temp ~  airTemp * factor(seg), data=ee )

acf(resid(ee2))
AIC(ee1,ee2)
#plot(ee2)

ggplot( ee, aes(airTemp,temp, colour=factor(seg))) +
  geom_point() +
  geom_line( aes(airTemp,predict(ee2)) ) +
  facet_wrap(~seg)

## add temporal autocorrelation
ee3 <- gls( temp ~  airTemp * factor(seg), correlation = corAR1(form=~dOY), data=ee)

ggplot( ee, aes(airTemp,temp, colour=(dOY))) +
  geom_point() +
  geom_line( aes(airTemp,predict(ee3)) ) +
  scale_colour_gradient(limits=c(0, 365), high='green',low="red")+
  facet_wrap(~seg)

acf(resid(ee3))

AIC(ee1,ee2,ee3)




pData <- cbind(min(ee$dOY):max(ee$dOY),predict(ee_3))
names(pData) <- c('dOY','t')

ggplot( ee, aes(dOY,tempScaled))+
  geom_point()+
  geom_point( aes(dOY,t),data=pData)


###################################
# just season 2
ee3 <- ee[ee$season == 2, ]

pairs(ee3[,c('dOY','temp','airTemp','tAirMax','tAirMin','flow','rain','precip','residT','residA')])

ggplot( ee3, aes(dOY,temp, colour=(dOY)) ) +
  geom_point(  )+
  geom_smooth()+
  geom_point( aes(dOY,airTemp), colour='red')

ggplot( ee3, aes(dOY,temp) ) +
  geom_point(  )+ geom_line( aes(dOY,temp)) +
  geom_point( aes(dOY,airTemp), colour='darkgreen') + geom_line( aes(dOY,airTemp), colour='darkgreen') +
  #geom_smooth()+
  geom_point( aes(dOY,tAirMax), colour='red') + geom_line( aes(dOY,tAirMax), colour='red') +
  geom_point( aes(dOY,tAirMin), colour='blue') +geom_line( aes(dOY,tAirMin), colour='blue')


ee3_1 <- lm( temp~airTemp, data=ee3)
plot(ee3_1)

ee3_2 <- gls( temp ~ airTemp, data=ee3)
plot(ee3_2)

ee3_3 <- gls( temp ~ airTemp, correlation = corAR1(form=~dOY), data=ee3)
points(ee3_3)
plot(ee3$temp, col='red')

pData <- cbind(min(ee3$dOY):max(ee3$dOY),predict(ee3_3))
names(pData) <- c('dOY','t')

ggplot( ee3, aes(dOY,temp))+
  geom_point()+
geom_point( aes(dOY,t),data=pData)


acf(ee3_3)
####################################

win.graph(); par(mfrow=c(1,1));

ggplot( ee, aes(dOY,temp))+
  geom_point() +
  geom_smooth(method='gam',formula=y ~ s(x, bs='cr', k=52))

ggplot( ee, aes(dOY,airTemp))+
  geom_point() +
  geom_smooth(method='gam',formula=y ~ s(x, bs='cr', k=52))

#  geom_point( aes(dOY,airTemp), colour='red') +
#geom_smooth(colour='red')
  


# simple gam with just day of Year
m1 <- gam( temp ~ s(dOY, bs='cr', k=15)
           , data=ee) 
gam.check(m1)
plot(m1, residuals=T)
acf(residuals(m1)) # strong pattern

# turn away from gam for now, try a sin/cos with autocorrelation
m2 <- lm(temp ~ cos(0.0172*dOY)+sin(0.0172*dOY), ,method="ML", data = ee )
plot(fitted(m2))
acf(resid(m2))

ggplot( ee, aes(dOY,temp)) +
  geom_point()+
  geom_point(aes(dOY,fitted(m2)))
# works well for some years, but not well for others (e.g. 2001)


# add in correlation
m3 <- gls(temp ~ cos(0.0172*dOY)+sin(0.0172*dOY),correlation=corAR1(form=~dOY) , data = ee )
plot(fitted(m3))
acf(resid(m3)) 
#still have strong autocorr

# autocorr of raw data - very high, as expected
acf(ee$temp)




# drop sin/cos, need to stick with spline
# add in correlation
m4 <- gls(temp ~ dOY,correlation=corAR1(form=~dOY) , data = ee )
plot(fitted(m4))
acf(resid(m4)) 
#still have strong autocorr
    
# 


m2 <- gamm( temp ~ s(dOY, bs='cr', k=15),
            correlation=corAR1(form = ~ 1 ),
            data=ee) 
acf(residuals(m2)) # strong pattern





#############################################
# spline across years
# didn't finish the function...
getSpline <- function( x,y,name, ... ){
  g <- gam( y ~ s( x, bs = "cc", k=15 ), data=e )
  j <- data.frame(dOY=1:366)
  p <- predict(g,j, type = "response")
  pred <- cbind(j,as.numeric(p))
  names(pred) <- c('dOY',name)
  return(pred)
}

predSpline <- getSpline( e$dOY, e$temp, 'temp', data=e )

gAir <- gam( airTemp  ~ s(dOY, bs = "cc", k=15), data=e )

p <- predict(g,j, type = "response")
pAir <- predict(gAir,j, type = "response")


e <- merge( x=e, y=predT, all.x=T )
e <- e[ order(e$date), ]

e$residT <- e$temp - e$splineT
e$residA <- e$airTemp - e$splineA


