#Application on GCM outputs 
library(loadeR)
library(downscaleR)
library(transformeR)
library(visualizeR)

#HISTORICAL PERIOD
#We select the same season of the training period
seas=c(6,7,8) 
yrsH=1990:2000

#The CSIRO outputs for the historical run are in the folder HistoricalCSIRO
#########They are already regridded to the NCEP R1 grid#################

#Loading GCM outputs

slp.HIST <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/HistoricalCSIRO/slp.CSIRO.1970a2005.SA.nc",
                        var = "slp",lonLim=c(-65.0,-50.0),latLim=c(-37.5,-25.0),season = seas, years = yrsH)

T850.HIST <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/HistoricalCSIRO/air.CSIRO.1975a2005.SA.85000.nc",
                          var = "air",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrsH)

q850.HIST <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/HistoricalCSIRO/shum.CSIRO.1975a2005.SA.85000.nc",
                          var = "shum",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrsH)
#library(ncdf4)
#nc<-nc_open("C:/ML/ICTP_2018_ESD_Training/HistoricalCSIRO/slp.1970a2005.CSIRO-Mk3-6-0_historical_r1i1p1.NCEP.SA.nc")
#Type nc to get the information about the netcdf file

#PREPARING THE PREDICTOR DATA FOR DOWNSCALING

q850.HIST[["Variable"]][["level"]]<-850
T850.HIST[["Variable"]][["level"]]<-850

xGCMHIST <- makeMultiGrid(T850.HIST,q850.HIST) #Multigrid combination
#xGCMHIST <- subsetGrid(xGCMHIST, years = 1992:1993) 

xGCMHIST <- scaleGrid(xGCMHIST, type = "standardize") #Standarizing 

#We prepare the GCM outputs to have the same structure of the predictor variables
#during the training period, that is, xyT

xyt     <- prepareNewData(newdata = xGCMHIST, data.structure = xyT)

predH <- downscale.predict(xyt, model = model)



###############################################################################
#We can visualize the time series 
#However, we should keep in mind that there is no temporal synchronization
#between GCM downscaled outputs and observations


BA.Obs <- subsetGrid(y,station.id = "087585",years = 1990:2000)
#BA.NCEP.ESD <- subsetGrid(model$pred,station.id = "087585",years = 1990:2000)
BA.predH <- subsetGrid(predH,station.id = "087585",years = 1990:2000)

#temporalPlot(BA.Obs,BA.NCEP.ESD,BA.predH)

#In order to explore the Added Value of ESD, we can also plot the time series of 
#the raw GCM output for Tmin in the closest grid point to Buenos Aires.

Tmin.HIST <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/HistoricalCSIRO/tasmin.CSIRO.1990a2005.nc",
                          var = "tasmin",lonLim=c(-60.0,-56),latLim=c(-36.4,-34.5),
                          season = seas, years = 1990:2000)

dates1990.2000 <- as.POSIXlt(BA.Obs$Dates$start, tz = "GMT")

# Plot time series
plot(dates1990.2000, BA.Obs$Data, ty = 'l', ylab = "tmin") #OBSERVATIONS
lines(dates1990.2000, Tmin.HIST$Data[,2,2]-273, col = "red", lty = 2)#GCM RAW OUTPUTS
lines(dates1990.2000, BA.predH$Data, col="blue")#GCM.ESD
legend("topleft", c("Observed", "RawGCM","GCM.ESD"), lty = c(1,2,1), col = c("black","red","blue"))

#We add to the plot some basics statistics:Mean
mtext(paste("OBSmean =", round(mean(BA.Obs$Data),digits=2),
            "  GCMmean = ", round(mean(Tmin.HIST$Data[,2,2]-273),digits=2),
            " GCM.ESDmean = ",round(mean(BA.predH$Data),digits=2)))

#SD
#sd.obs<-round(sd(BA.Obs$Data),digits=2) #Observed Standard Deviation
#sd.GCM<-round(sd(Tmin.HIST$Data[,2,2]-273),digits=2)
#sd.GCM.ESD<-round(sd(BA.predH$Data),digits=2)


########################################################################################
#FUTURE CLIMATE PROJECTIONS
#We choose the RCP4.5 
#Period analyzed: 2086-2100

slpRCP45 <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/RCP4.5CSIRO/slp.CSIRO.2086a2100.SA.nc",
                    var = "slp",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = 2086:2100)

T850RCP45 <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/RCP4.5CSIRO/air.CSIRO.2086a2100.SA.850.nc",
                     var = "air",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = 2086:2100)

q850RCP45 <- loadGridData(dataset = "C:/ML/ICTP_2018_ESD_Training/RCP4.5CSIRO/shum.CSIRO.2086a2100.SA.850.nc",
                          var = "shum",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = 2086:2100)


T850RCP45[["Variable"]][["level"]]<-850
q850RCP45[["Variable"]][["level"]]<-850

xGCMRCP45 <- makeMultiGrid(T850RCP45,q850RCP45)
#xGCMRCP45 <- subsetGrid(xGCMRCP45, years = 2087:2090) 

#xtRCP45 <- xGCMRCP45  
xGCMRCP45 <- scaleGrid(xGCMRCP45, type = "standardize") 
xytRCP45     <- prepareNewData(newdata = xGCMRCP45, data.structure = xyT)
#We use the "model" trained with NCEP R1 to downscale RCP4.5 CSIRO outputs in the future
predRCP45 <- downscale.predict(xytRCP45, model = model)

OBSmean1990.2000<-round(mean(BA.Obs$Data),digits=2)
#GCMmean1990.2000<-round(mean(Tmin.HIST$Data[,2,2]-273),digits=2)
GCMH.CSIRO.mean1990.2000<-round(mean(BA.predH$Data),digits=2)
GCMRCP45.CSIRO.mean2086.2100<-round(mean(predRCP45$Data[,1]),digits=2)


BA.predRCP45.2086 <- subsetGrid(predRCP45,station.id = "087585",years = 2086)
temporalPlot(BA.predRCP45.2086)
