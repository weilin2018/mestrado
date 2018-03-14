#-*-coding: utf-8-*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import pickle

import sys
sys.path.append('../rotinas/artigoTGpack/')

import artigoTGpack as oceano

# criar BASE_DIR

BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR + 'dados_bndo/'

name = [u'Terminal Ilha Guaíba (50165)', 
        u'Farol Castelhanos (50167)', 
        u'Colegio Naval (501812)', 
        u'Parati (50193)']

lats = [-23.00000000, -23.16694444, -23.01694444, -23.21750000]
lons = [-44.01916667, -44.08500000, -44.31833333, -44.70138889]

# criar mapa

fig, ax = plt.subplots(figsize=(10,8))
m=pickle.load(open("/home/tparente/danilo/mestrado/github/artigoTG/rotinas/bigmap.p","rb"))
m.ax = ax 

# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.4) #linha de costa em alta resolução
m.drawmapboundary(fill_color='#e5f2ff') # fill_color colore o oceano
m.fillcontinents(color='#ffd480') # colorir o continente

# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

# plotar estações
for x,y,label in zip(lons,lats,name):
    m.plot(x,y,'s',label=label,latlon=True)

lg = plt.legend(loc='upper left',fontsize=12,numpoints=1,scatterpoints=1)
lg.get_frame().set_alpha(.4)

plt.title(u'Estações BNDO para validação')

plt.show()