#!/usr/bin/env python

'''
script para utilizar o cdo (Climate Data Operators) para recortar os arquivos grb2 que tenho
para uma região específica


'''

import os
import glob


import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/FTP/')
SAVE_DIR = BASE_DIR.replace('github/', 'ventopcse/data/FTP/recortados/')

nfiles = glob.glob(DATA_DIR+'*.grib2')
nfiles.sort()

os.system('clear')

for f in nfiles:
	print("Cutting: %s"%(f.split('/')[-1]))
	os.system('cdo -sellonlatbox,-50,-40,-20,-30 %s %s'%(f, f.replace('FTP/','FTP/recortados/')))
	os.system('mv %s %s'%(f, f.replace('FTP/','FTP/processados/')))

# apos recortar, converter os arquivos para netcdf usando o proprio cdo

nfiles = glob.glob(DATA_DIR.replace('FTP/','FTP/recortados/*.grib2'))
nfiles.sort()

os.system('clear')
os.system('Converting grib2 files to netCDF files')

for f in nfiles:
	os.system('cdo -f nc copy %s %s' % (f, f.replace('.grib2', '.nc')))

os.system('mv %s %s' % (SAVE_DIR+'*.grib2', SAVE_DIR+'grb2/'))
