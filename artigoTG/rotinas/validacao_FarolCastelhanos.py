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
SAVE_DIR = BASE_DIR + 'artigoTG/figures/calibracao/'

os.system('clear')
print("Reading and plotting BNDO's data - Ilha Guaiba Termianl")
#print("Lendo e plotando dados do BNDO - Terminal da Ilha Guaíba")

# Tidal Gauge location:
location = [-23.16694444, -44.08500000]

# ler arquivos de 1997 fornecidos pelo BNDO: 50167002751204199713061997ALT.txt
lfiles = glob.glob(DATA_DIR+'50167002751204199713061997ALT.txt')

fname = pd.read_csv(lfiles[0], skiprows=11, delimiter=';', names=['level', 'x'])

fname.drop(fname.columns[len(fname.columns)-1], axis=1, inplace=True)

dtRange = pd.date_range(start='1997-04-12 00:00', end='1997-06-13 23:00', freq='H')
observ = pd.DataFrame({'level': fname['level'].values/100.}, index=dtRange)

# controle de qualidade
cond = observ['level'] > 4.
observ[cond] = np.nan

# removendo a média da série temporal
observ['level'] = observ['level'] - observ['level'].mean()

#observ.plot(title='Original Data from BNDO with Hourly Frequency- Farol Castelhanos')
#plt.show()


os.system('clear')
print("Reading and plotting sECOM's output")

try:
	df = pickle.load(open(SIMS_DIR+'run09/df.pickle','r'))
	wndDic = pickle.load(open(SIMS_DIR+'run09/wndDic.pickle','r'))
	wu = wndDic['wu']
	wv = wndDic['wv']
	time = wndDic['time']
except:
	print('Tre')
	# importar e tratar dados
	run = 'run09'

	ncfile = xr.open_dataset(SIMS_DIR+run+'/'+run+'.cdf')

	lon = ncfile['lon'].data
	lat = ncfile['lat'].data

	lon[lon == 0.0] = np.nan
	lat[lon == 0.0] = np.nan

	# localizacao do terminal da ilha guaiba
	ilat = location[0]
	ilon = location[1]

	iss, jss = oceano.find_nearest(lon,lat,ilon,ilat)

	# importing information from netCDF
	time = ncfile['time'].values # data from 1997-01-06 00h onwards
	elev = ncfile['elev'].data[:,iss,jss]
	elev = np.squeeze(elev)

	model = pd.DataFrame({'modeled':elev}, index=time)

	dtRange = pd.date_range(start=time[0],end=time[-1], freq='18359s')

	#model.plot()
	#plt.show()


	wu = ncfile['wu'].data[:,iss,jss]
	wv = ncfile['wv'].data[:,iss,jss]
	spd = np.sqrt(wu**2, wv**2)                            

	d = pd.DataFrame({'wndSpd':np.squeeze(spd)},index=time)
	d = d['1997-04-12':'1997-06-13']

	modelo = model['1997-04-12':'1997-06-13']

	df = pd.DataFrame({'sECOM': modelo.modeled.values, 'BNDO': observ.level.values, 'Wind': d.wndSpd.values}, index=modelo.index.values)



os.system('clear')
print('plotting data')

fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(15,10), sharex=True)

ax[0].plot(df.index.values, df.BNDO.values, label='BNDO observed data')
ax[0].plot(df.index.values, df.sECOM.values, label='sECOM modeled data')
ax[0].set_ylabel('Elevation [m]',fontsize=20)

ax[0].legend(loc='upper left')

ax[1].plot(df.index.values, df.Wind.values, label='Wind Speed from sECOM')
ax[1].set_xlabel('Time [days]',fontsize=20)
ax[1].set_ylabel(r'Intensity [$m s^{-1}$]',fontsize=20)

ax[1].legend(loc='upper left')

skill = oceano.skill_willmott(df.BNDO.values, df.sECOM.values)
corre = pearsonr(df.BNDO.values, df.sECOM.values)[0]

plt.suptitle('Farol Castelhanos [r = %0.2f, skill = %0.2f] \n 12/04/1997 to 13/06/1997' % (corre,skill), fontsize=24,x=0.52)

plt.show()



############

fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(15,10), sharex=True)

ax[0].plot(df['1997-05-12':].index.values, df['1997-05-12':].BNDO.values, label='BNDO observed data')
ax[0].plot(df['1997-05-12':].index.values, df['1997-05-12':].sECOM.values, label='sECOM modeled data')
ax[0].set_ylabel('Elevation [m]',fontsize=20)

ax[0].legend(loc='upper left')

ax[1].plot(df['1997-05-12':].index.values, df['1997-05-12':].Wind.values, label='Wind Speed from sECOM')
ax[1].set_xlabel('Time [days]',fontsize=20)
ax[1].set_ylabel(r'Intensity [$m s^{-1}$]',fontsize=20)

ax[1].legend(loc='upper left')

skill = oceano.skill_willmott(observ['1997-05-12':].level.values, modelo['1997-05-12':].modeled.values)
corre = pearsonr(observ.level.values, modelo.modeled.values)[0]

plt.suptitle('Farol Castelhanos [r = %0.2f, skill = %0.2f] \n 12/05/1997 to 13/06/1997' % (corre,skill), fontsize=24,x=0.52)

plt.show()

############

# detectando se é ou não shelfwaves


m=df['sECOM'].rolling(40, center=True).mean()
o=df['BNDO'].rolling(40, center=True).mean()

newdf = pd.DataFrame({'BNDO':o,'sECOM':m},index=df.index)

wind = pd.DataFrame({'wu':np.squeeze(wu),'wv':np.squeeze(wv)},index=time)

w = wind['1997-04-12':'1997-06-13']

fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(15,10))

ax[0].plot(newdf.index.values, newdf.BNDO.values, label='BNDO')
ax[0].plot(newdf.index.values, newdf.sECOM.values, label='sECOM')

ax[0].legend(loc='lower right')

ax[1].plot(w.index.values, w.wu.values, label='Zonal')
ax[1].plot(w.index.values, w.wv.values, label='Meridional')

ax[1].legend(loc='lower right')

ax[0].set_title(u'Elevação do Nível do mar (filtragem de 40h)')
ax[1].set_title(u'Componente Zonal [vermelho] e Meridional [azul] do vento')

plt.show()

##################

# # comparar ventos (componentes) com elevacao
# fig, ax = plt.subplots(nrows=3,ncols=1,figsize=(15,10))

# ax[0].plot(df.index.values, df.BNDO.values, label='BNDO observed data')
# ax[0].plot(df.index.values, df.sECOM.values, label='sECOM modeled data')

# wind = pd.DataFrame({'wu':np.squeeze(wu),'wv':np.squeeze(wv)},index=time)

# w = wind['1997-04-12':'1997-06-13']

# ax[1].plot(w.index.values, w.wu.values/2,label='wu')
# ax[2].plot(w.index.values, w.wv.values/2, label='wv')

# plt.show()