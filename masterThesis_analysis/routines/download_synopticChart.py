import os
import wget
import requests
import numpy as np
import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

BASE_DIR = oceano.make_dir()
SAVE_DIR = BASE_DIR.replace('github/','ventopcse/data/synopticChart/')
FTP_BASE = 'http://img0.cptec.inpe.br/~rgptimg/Produtos-Pagina/Carta-Sinotica/Analise/Superficie/'

# creating vectors with years, months and day to download figures
years  = '2012 2013 2014 2015 2016 2017'.split(' ')
months = '01 02 03 04 05 06 07 08 09 10 11 12'.split(' ')


days = {
    '01': np.arange(1,32,1),
    '02': np.arange(1,29,1),
    '03': np.arange(1,32,1),
    '04': np.arange(1,32,1),
    '05': np.arange(1,32,1),
    '06': np.arange(1,31,1),
    '07': np.arange(1,32,1),
    '08': np.arange(1,32,1),
    '09': np.arange(1,31,1),
    '10': np.arange(1,32,1),
    '11': np.arange(1,31,1),
    '12': np.arange(1,32,1)
}

# creating names of files
fnames = []
for y in years:
    for m in months:
        for d in days[m]:
            d = str(d)
            base = 'superficie_'+y+m+d.zfill(2)+'00.gif'
            fnames.append(FTP_BASE+base)

fnames = np.asarray(fnames)

# entering the directory to save
os.chdir(SAVE_DIR)

# download
os.system('clear')
for f in fnames:
    wget.download(f)

##########################################################
### downloading a specific period
SAVE_DIR = BASE_DIR.replace('github/','ventopcse/data/synopticChart/201312/')

year = '2013'
mont = '12'
days = '19 20 21 22 23 24 25 26'.split(" ")
hour = '00 06 12 18'.split(" ")

fnames = []
for m in mont:
    for d in days:
        for h in hour:
            fnames.append(FTP_BASE+'superficie_'+year+m+d+h+'.gif')

fnames = np.asarray(fnames)

# entering the directory to save
os.chdir(SAVE_DIR)

# download
os.system('clear')
for f in fnames:
    wget.download(f)
