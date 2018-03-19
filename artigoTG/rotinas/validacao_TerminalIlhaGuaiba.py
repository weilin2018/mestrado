#-*-coding: utf-8-*-
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import pickle
#import scipy.spatial import cKDTree
import glob
import xray as xr 
import os
import string
from scipy.stats.stats import pearsonr 


import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('../rotinas/artigoTGpack/')

import artigoTGpack as oceano

# criar BASE_DIR

BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/','artigo_data/BNDO/')
SIMS_DIR = BASE_DIR.replace('github/','/artigo_data/simulacoes/')

os.system('clear')
print("Reading and plotting BNDO's data - Ilha Guaiba Terminal")

# location from tidal gauge
ilat = -23.00000000
ilon = -44.01916667

os.system('clear')
runs = 'run05 run06 run07 run08'.split(' ')
runs = ['run06','run08','run11','run10']

labels = {
	'run00': 'Original',
	'run05': r'$S_2$ - 5H',
	'run06': r'$M_2$ - 5H',
	'run10': r'$M_2$ - 12H',
	'run11': r'$M_2$ - 10H',
	'run07': r'$S_2$ - 6H',
	'run08': r'$M_2$ - 6H'
}


fig, axes = plt.subplots(nrows=len(runs),ncols=1, sharex=True)

i = 0 # contador dos eixos

for run in runs:

	observ = oceano.read_BNDO(DATA_DIR) 					# ler os dados BNDO
	ncfile = xr.open_dataset(SIMS_DIR+run+'/'+run+'.cdf')   # ler netcdf do modelo

	lon = ncfile['lon'].data 								# extrair lat e lon do modelo
	lat = ncfile['lat'].data

	lon[lon == 0.0] = np.nan
	lat[lon == 0.0] = np.nan

	iss, jss = oceano.find_nearest(lon,lat,ilon,ilat) 		# encontrar indices do ponto mais proximo

	time = ncfile['time'].values							# extrair tempo e elevacao
	elev = ncfile['elev'].data[:,iss,jss]
	elev = np.squeeze(elev)

	model = pd.DataFrame({'modeled':elev}, index=time) 		# criar dataframe

	# criar date_range
	if run=='run00' or run == 'run06' or run == 'run08' or run == 'run10' or run == 'run11':
		pass
	else:
		dtRange = pd.date_range(start=time[0],end=time[-1], freq='18360000000000ns')

	# resample dos dados BNDO para os mesmos instantes do modelo
	n = oceano.newTimerange(dtRange, observ['1997-01-06':'1997-02-06'].index, observ)

	# dataframe final
	df = pd.DataFrame({'sECOM': model.modeled.values[:-1],
					    'BNDO': n
						}, index = dtRange)

	# calculo da correlação e skill
	corr  = pearsonr(df.BNDO.values, df.sECOM.values)[0]
	skill = oceano.skill_willmott(df.BNDO.values, df.sECOM.values)

	corStr = 'Corr.: %0.2f'%df['1997-01-06':'1997-02-06'].corr()['sECOM'][0]
	skiStr = 'Skill: %0.2f'%skill

	title = labels[run]+': %s and %s' % (corStr, skiStr)

	axes[i].plot(df['1997-01-06':'1997-02-06'].index, df['1997-01-06':'1997-02-06'].values)
	axes[i].set_title(title, fontsize=15)

	axes[i].legend(['BNDO','sECOM'],loc='lower right')



	i = i + 1 # contador do eixo de plotagem

axes[-1].set_xlabel('Time in days')

plt.suptitle(u'Variação do nível do mar', fontsize=24, x=.515)
plt.show()

# criar novo date_range com frequencia de 18360000000000ns
# dtRange = pd.date_range(start=time[0], end=time[-1], freq='18360000000000ns')

