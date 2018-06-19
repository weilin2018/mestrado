#!/bin/bash

if [ "$1" == "-h" ]; then
	echo " 
Routine to extract temperature, salinity, elevation and current data from gcmplt.cdf,
created by a ECOM simulation.

Give, as 1st argument, the directory where the file is and, as 2nd argument, the filename.

4 new netcdf file will be creatad in the directory given.

Example:
	$ ./cut_netcdf_with_cdo.sh '/home/user/path/to/file/' 'exp03.cdf'
	"
	exit 0
fi
clear

# atribuindo nomes melhores as variaveis enviadas pelo terminal
BASE_DIR=$1
FILENAME=$2
SAVE_DIR=$1"exp03_variables/"

INFILE=$1$2

echo "Processing file $FILENAME"
echo "############"
echo "Extracting temperature data ..."
tmp=$SAVE_DIR'tmp.nc'
temp_fin=$SAVE_DIR'temp_exp03.cdf'
cdo -selname,temp $INFILE $tmp
cdo -seltimestep,112/352 $tmp $temp_fin
echo "Removing temporary files ..."
rm $tmp

echo "############"
echo "Extracting salinity data ..."
tmp=$SAVE_DIR'tmp.nc'
salt_fin=$SAVE_DIR'salt_exp03.cdf'
cdo -selname,salt $INFILE $tmp
cdo -seltimestep,112/352 $tmp $salt_fin
echo "Removing temporary files ..."
rm $tmp

echo "############"
echo "Extracting elevation data ..."
tmp=$SAVE_DIR'tmp.nc'
elev_fin=$SAVE_DIR'elev_exp03.cdf'
cdo -selname,elev $INFILE $tmp
cdo -seltimestep,112/352 $tmp $elev_fin
echo "Removing temporary files ..."
rm $tmp

echo "############"
echo "Extracting current data ..."
tmp=$SAVE_DIR'tmp.nc'
curr_fin=$SAVE_DIR'curr_exp03.cdf'
cdo -selname,u,v,w $INFILE $tmp
cdo -seltimestep,112/352 $tmp $curr_fin
echo "Removing temporary files ..."
rm $tmp

echo "Finish."