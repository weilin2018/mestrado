import os
import wget
import numpy as np

FTP_BASE = 'ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/'

SAVE_DIR = '/home/tparente/danilo/mestrado/ventopcse/data/CFSR_FTP/'

#os.chdir(SAVE_DIR)

YEARS = np.arange(1996,2014,1)

MONTHS= '11 12 01 02 03'.split(' ')

DATES = []

for year in YEARS:
    for month in MONTHS:
        DATES.append(str(year)+month)

os.system('clear')

for date in DATES:
	fname = FTP_BASE+date+'/wnd10m.gdas.'+date+'.grb2'
	print('Downloading: %s \n'%(fname))

	# os.system('wget %s' % (fname))
	wget.download(fname)

