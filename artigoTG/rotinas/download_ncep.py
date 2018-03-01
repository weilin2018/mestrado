import os

FTP_BASE = 'ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/'
SAVE_DIR = '/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/'

DATES = '199612 199701 199702 199703 199704 199705 199706'.split(' ')

os.chdir('/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR')

for date in DATES:
	fname = FTP_BASE+date+'/wnd10m.gdas.'+date+'.grb2'
	os.system('wget %s' % (fname))

	fname = FTP_BASE+date+'/prmsl.gdas.'+date+'.grb2'
	os.system('wget %s' % (fname))

