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
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR_FTP/')
SAVE_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR_FTP/recortados/')

nfiles = glob.glob(DATA_DIR+'*.grb2')
nfiles.sort()

for f in nfiles:
	print("Cutting: %s"%(f.split('/')[-1]))
	os.system('cdo -sellonlatbox,-50,-40,-20,-30 %s %s'%(f, f.replace('CFSR_FTP/','CFSR_FTP/recortados/')))

