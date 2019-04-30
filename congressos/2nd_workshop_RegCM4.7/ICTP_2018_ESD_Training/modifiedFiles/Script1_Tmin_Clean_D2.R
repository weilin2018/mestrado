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

CAM_SAM <- "/home/aluno/ICTP_2018_ESD_Training/GSN_CAM_SAM/"

di <- dataInventory(CAM_SAM)
stationInfo(CAM_SAM)

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


#We choose the PREDICTORS in a subregion of Southeastern South America
slp <- loadGridData(dataset = "/home/aluno/ICTP_2018_ESD_Training/R1/slp.R1.1978a2014.SA.nc",
                   var = "slp",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)

T850 <- loadGridData(dataset = "/home/aluno/ICTP_2018_ESD_Training/R1/air.R1.1978a2014.SA.nc",
                     var = "air@850",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)

q850 <- loadGridData(dataset = "/home/aluno/ICTP_2018_ESD_Training/R1/shum.R1.1978a2014.SA.nc",
                     var = "shum@850",lonLim=c(-60.0,-47.5),latLim=c(-37.5,-30.0),season = seas, years = yrs)


################################################################################
#Preparing the Predictors and Predictands to perform the downscaling
#############using T850 and q850 as predictors###############################
#We will use the information of the period 1990:2000 to 
#downscale HISTORICAL and RCP outputs


xoriginal2 <- makeMultiGrid(T850,q850) #Creating the multigrid predictor data
getVarNames(xoriginal2)# Displays the variables names

x2 <- subsetGrid(xoriginal2, years = 1990:2000) #Selecting the 1990:2000 period

#We get sure that we are working in the same predictor and predictand calibration period
y <- getTemporalIntersection(obs = y,prd = x2, "obs" )
x2 <- getTemporalIntersection(obs = y,prd = x2, "prd" )
x2 <- scaleGrid(x2, type = "standardize")


# Prepare predictors and predictands to get the input for the downscale.train function
xyT<- prepareData(x = x2, y = y)

# We can also reduce the predictors dimension using Principal Components
#xyT.pc<- prepareData(x2 = x,y = y,spatial.predictors = list(which.combine = getVarNames(x),v.exp = 0.9))

model <- downscale.train(xyT, method = "analogs",n.analogs=1,window=7,sel.fun = "mean")


#WE WILL USE "model" TO DOWNSCALE THE HISTORICAL AND RCP4.5 GCM OUTPUTS 

################################################################################

