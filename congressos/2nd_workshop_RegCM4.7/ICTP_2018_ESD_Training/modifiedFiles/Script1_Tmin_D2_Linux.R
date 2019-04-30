#2nd Workshop on Regional Climate Modeling and Extreme Events over South America
#Sao Paulo, Brazil

#Training Session Part 1
#Empirical Statistical Downscaling: Perfect Prognosis Activity
##############################################################

#In this activity we will use the open tool climate4R developed by the Santander Meteorology Group
#https://www.meteo.unican.es/en/climate4R
#https://github.com/SantanderMetGroup

rm(list=ls())

#The tool is formed by 4 R packages: loadeR, downscaleR, transformeR and visualizeR
#We first load the packages 
library(loadeR)
library(downscaleR)
library(transformeR)
library(visualizeR)

#You can find many examples of the use of the packages in the wiki
#https://github.com/SantanderMetGroup/loadeR/wiki
#https://github.com/SantanderMetGroup/downscaleR/wiki
#https://github.com/SantanderMetGroup/visualizeR/wiki
#https://github.com/SantanderMetGroup/transformeR/wiki

#You can also find many examples in the R help

####################################################################################

# Loading PREDICTANDS
# We first load the Station Data, which should have an specific format, see the format at
#https://github.com/SantanderMetGroup/loadeR/wiki/Standard-(ASCII)-format-for-station-data

#The CAM/SAM Station Data is in the folder GSN_CAM_SAM. We specify the path 
#and explore the location and information of stations.

CAM_SAM <- "/media/ICTP_2018_ESD_Training/GSN_CAM_SAM/"
setwd(CAM_SAM)
di <- dataInventory(CAM_SAM)
stationInfo(CAM_SAM)


#To load the data, we use the loadStationData function

#Loading station data using the station codes provided by the inventory
example1 <- loadStationData(dataset = CAM_SAM, 
                            var="tmin", 
                            stationID = c("082331", "087585"), #Manaus and Buenos Aires Stations
                            season = 6:8, 
                            years = 1981:2000)

#Loading station data from geographical coordinates
example2 <- loadStationData(dataset = CAM_SAM, 
                            var="tmax", 
                            lonLim = -70, 
                            latLim = 0, 
                            season = 6:8, 
                            years = 1981:2000)

#Loading station data within a given geographical bounding box
example3 <- loadStationData(dataset = CAM_SAM, 
                            var="precip", 
                            lonLim = c(-65,-50), 
                            latLim = c(0,20), 
                            season = 6:8, 
                            years = 1981:2000)

#Loading all stations
example4 <- loadStationData(dataset = CAM_SAM, 
                            var="tmax", 
                            season = c(12,1,2), 
                            years = 1981:2000)

#We can use the functions spatialPlot and climatology to visualize the data climatology
spatialPlot(climatology(example4), backdrop.theme = "countries", colorkey = T)


#We can plot the time series 
#time <- as.POSIXlt(example1$Dates$start)
#First station
#plot(time, example1$Data[,1], ty = 'l', col = "blue", xlab = "time", ylab = "T (?C)",ylim=c(-10, 30))
#Second station
#lines(time, example1$Data[,2], ty = 'l', col = "red")
#legend("bottomright", c("Manaus", "Buenos Aires"), col = c("blue", "red"), lty = 1)
#title("Tmin - JJA (1981-2000)")


######################################################################################
#For our ESD example, we choose tmin in Buenos Aires and Corrientes Stations in Argentina
CAM_SAM_Tmin<- loadStationData(dataset = CAM_SAM, 
                              var="tmin", 
                              stationID = c("087585","087166"), 
                              season = 6:8, 
                              years = 1981:2002)

y <- CAM_SAM_Tmin #our predictand data

#####################################################################################
#Loading PREDICTORS

#For this example we'll use NCEP Reanlaysis 1
# Definition of common time domain when loading the data
#seas=c(12,1,2)
#seas=c(3,4,5,6)
#seas=c(9,10,11)
seas=c(6,7,8)
yrs = 1990:2005

#We load the netcdf data using loadGridData
slp <- loadGridData(dataset = "/media/ICTP_2018_ESD_Training/R1/slp.R1.1978a2014.SA.nc",
                    var = "slp",season = seas, years = yrs)

spatialPlot(climatology(slp),backdrop.theme = "countries") #to visualize the slp mean field


#We choose the PREDICTORS in a subregion of Southeastern South America
slp <- loadGridData(dataset = "/media/42AF-A583/ICTP_2018_ESD_Training/R1/slp.R1.1978a2014.SA.nc",
                    var = "slp",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)

spatialPlot(climatology(slp),backdrop.theme = "countries")


T850 <- loadGridData(dataset = "/media/ICTP_2018_ESD_Training/R1/air.R1.1978a2014.SA.nc",
                     var = "air@850",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)

q850 <- loadGridData(dataset = "/media/rocio/42AF-A583/ICTP_2018_ESD_Training/R1/shum.R1.1978a2014.SA.nc",
                     var = "shum@850",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)


################################################################################
#Preparing the Predictors and Predictands to perform the downscaling
#Training period:1990:2000


xoriginal <- makeMultiGrid(slp,T850) #Creating the multigrid predictor data
getVarNames(xoriginal)# Displays the variables names in the multigrid

x <- subsetGrid(xoriginal, years = 1990:2000) #We choose a sub-period to train the model

#We get sure that we are working in the same predictor and predictand training period
y <- getTemporalIntersection(obs = y,prd = x, "obs" )
x <- getTemporalIntersection(obs = y,prd = x, "prd" )

#We should standarize predictors as they are different physical variables
x <- scaleGrid(x, type = "standardize")


# Prepare predictors and predictands to get the input for the downscale.train function
xyT<- prepareData(x = x, y = y)

# We can also reduce the predictors dimension using Principal Components
xyT.pc<- prepareData(x = x,y = y,spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))

################################################################################
#Downscaling Minimum Temperature
#Function: downscale.train

# Training the model

#model <- downscale.train(xyT, method = "analogs",
#                         sel.fun = "mean")

#Other specifications for the method
model <- downscale.train(xyT, method = "analogs",n.analogs=1,window=7,sel.fun = "mean")


BA.Obs.1995 <- subsetGrid(y,station.id = "087585",years = 1995)
BA.Pred.1995 <- subsetGrid(model$pred,station.id = "087585",years = 1995)
temporalPlot(BA.Obs.1995, BA.Pred.1995)




#Application of the model using Reanalysis Data
#Prediction: downscale.predict

#Option 1: in the same 1 period 1990-2000 
xyt     <- prepareNewData(newdata = x, data.structure = xyT)
pred1 <- downscale.predict(xyt, model = model)

BA.Obs.1995 <- subsetGrid(y,station.id = "087585",years = 1995)
BA.pred1.1995 <- subsetGrid(pred1,station.id = "087585",years = 1995)
temporalPlot(BA.Obs.1995,BA.pred1.1995)


#Option 2: in an independent period, for instance year 2002
xt <- subsetGrid(xoriginal,years = 2002)      
xt <- scaleGrid(xt, type = "standardize") 
xyt     <- prepareNewData(newdata = xt, data.structure = xyT)
pred2 <- downscale.predict(xyt, model = model)


BA.Obs.2002 <- subsetGrid(CAM_SAM_Tmin,station.id = "087585",years=2002)

#BA.Obs.2002 <- loadStationData(dataset = CAM_SAM, 
#                               var="tmin", 
#                               stationID = "087585",seas=c(6,7,8),years=2002)
BA.pred2.2002 <- subsetGrid(pred2,station.id = "087585",years = 2002)
temporalPlot(BA.Obs.2002,BA.pred2.2002)


#Option 3: K-folding cross-validation
#Function: downscale.cv
#analog.cv <- downscale.cv(x = x, y = y, method = "analogs",folds = 5, type = "chronological", n.analogs = 1)

analog.cv <- downscale.cv(x = x, y = y, method = "analogs",folds = 5, type = "chronological", n.analogs = 1, 
                          spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))

#New version 

analog.cv <- downscale.cv(x = x, y = y, method = "analogs",folds = 5,sampling.strategy = "kfold.chronological", n.analogs = 1, 
                          spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))


BA.Obs.1995 <- subsetGrid(y,station.id = "087585",years = 1995)
BA.pred3.1995 <- subsetGrid(analog.cv,station.id = "087585",years = 1995)

#temporalPlot(BA.Obs.1995,BA.pred3.1995)

temporalPlot(BA.Obs.1995,BA.pred1.1995,BA.pred3.1995)


######################################################################################################
#TESTING DIFFERENT PREDICTOR VARIABLES
#We try with the combination of T850 and q850

xoriginal2 <- makeMultiGrid(T850,q850) #Creating the multigrid predictor data
getVarNames(xoriginal2)# Displays the variables names

x2 <- subsetGrid(xoriginal2, years = 1990:2000) #We choose a sub-period to train the model

#We get sure that we are working in the same predictor and predictand training period
y <- getTemporalIntersection(obs = y,prd = x2, "obs" )
x2 <- getTemporalIntersection(obs = y,prd = x2, "prd" )
x2 <- scaleGrid(x2, type = "standardize")

#We perform the downscaling using x2 (T850-q850) as predictors and 
#y (tmin) as predictands
#We apply PCA to the predictor variables, and retain those that explain more than 90% of the variance
analog.cv2 <- downscale.cv(x = x2, y = y, method = "analogs",folds = 5, sampling.strategy = "kfold.chronological", n.analogs = 1, 
                          spatial.predictors = list(which.combine = getVarNames(x2),v.exp = 0.9))

#We select the year 1995 to visualize the time series
BA.Obs.1995 <- subsetGrid(y,station.id = "087585",years = 1995)
BA.pred3.1995 <- subsetGrid(analog.cv,station.id = "087585",years = 1995) #slp-T850
BA.pred4.1995 <- subsetGrid(analog.cv2,station.id = "087585",years = 1995) #T850-q850

temporalPlot(BA.Obs.1995,BA.pred3.1995)
temporalPlot(BA.Obs.1995,BA.pred3.1995,BA.pred4.1995)

#Validating the complete period 1990:2000
BA.Obs<-subsetGrid(y,station.id = "087585",years = 1990:2000)
BA.pred4 <- subsetGrid(analog.cv2,station.id = "087585",years = 1990:2000)

temporalPlot(BA.Obs,BA.pred4)


#The following lines make a basic comparison of the observed and downscaled time series
# The test period is coerced to a POSIXlt class object to define the temporal dimension
test.period <- as.POSIXlt(BA.Obs$Dates$start, tz = "GMT")

# Figure 1: Temporal series
plot(test.period, analog.cv2$Data[,1], ty = 'l', ylab = "tmin")
lines(test.period, BA.Obs$Data, col = "red", lty = 2)
legend("topleft", c("Observed", "Downscaled"), lty = c(1,2), col = c(1,2))

# Figure 2: Observed vs Downscaled daily values
plot(BA.Obs$Data, analog.cv2$Data[,1], asp = 1, ylab ="Observed" , xlab = "Downscaled", col = "blue")
lm1 <- lm(BA.Obs$Data ~ analog.cv2$Data[,1]) #Linear Regression, red line
abline(reg = lm1, col = "red")
r2 <- round(summary(lm1)$adj.r.squared, 2) #Variance explained
bias <- round(mean(analog.cv2$Data[,1] - BA.Obs$Data, na.rm = TRUE), 2) #BIAS 
rmse <- round(sqrt(mean((analog.cv2$Data[,1] - BA.Obs$Data)^2)),digits=2) #RMSE
sd.obs<-sd(BA.Obs$Data) #Observed Standard Deviation
ks <- ks.test(analog.cv2$Data[,1],BA.Obs$Data,alternative="two.sided") #Comparison of pdf
mtext(paste("Rsq =", r2, "  bias = ", bias, " rmse/sd.obs = ",round(rmse/sd.obs,digits=2), "  ks-pvalue = ",ks$p.value))

#Figure 3: qqPlot, comparison of pdf and percentiles
qqplot(BA.Obs$Data, analog.cv2$Data[,1], asp = 1, ylab = "Observed", xlab = "Downscaled", col = "blue")
abline(a=0,b=1, col = "red")
mtext(paste("P1.Obs =", quantile(BA.Obs$Data,probs=0.01), " / P1.Pred = ", quantile(analog.cv2$Data[,1],probs=0.01)))
title("qqPlot Buenos Aires")


