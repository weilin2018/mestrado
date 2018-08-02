"""
funcao para plotar a dispersão da pluma d material nuclear no
cenário 6, com ventos variando no espaço e no tempo, plotando
também os vetores de direção do vento.
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
# from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib import dates
# import datetime
import cmocean as cmo

# pacotes para minimap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

os.system('clear')

## functions

def convertData2Bq(c,unit='mg/L'):

    # constante = 0.3e+22
    c[c<0.0] = np.nan # valores negativos são removidos
    if unit == 'mg/L':
        # converter de mg para g
        c = c*.1e-3
    else:
        # converter de kg para g
        c = c*.1e3

    c = c*35334440479135.016 # converter de g para Bq

    return c


## main
base = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/sims_dispersao/wcptec/'

conc = pickle.load(open(base+'concentracaoInt.pkl','rb'))
conf = pickle.load(open(base+'config.pkl','r'))
lon  = conf['lon']
lat  = conf['lat']
time = conf['time']

wind = pickle.load(open(base+'windComp.pkl','r'))
wu   = wind['wu']
wv   = wind['wv']


lon[lon == 0.] = np.nan
lat[lat == 0.] = np.nan

llon = -44.9;
ulon = -43.45;
llat = -23.6;
ulat = -22.75;

# skip vectors
winint = 15
curint =  5
skipcurr = (slice(None,None,curint),slice(None,None,curint))
skipwind = (slice(None,None,winint),slice(None,None,winint))


# converter o dado de Kg para Bq/l
k=convertData2Bq(conc)
# liberar espaço na memória
# del c
# normalizar os valores
k2=k/.1e+7

# linhas de contorno para usar no contourf
maxC =np.nanmax(k2[-1,:,:])+100
contour_levels=np.arange(0.,maxC,maxC/1e+3)

# plotting
os.system('clear')
for t in range(wu.shape[0]):
    print('Plotting for %s'%(str(time[t])))
    fig,ax = plt.subplots(figsize=(12,8))

    uw,vw = wu[t,:,:],wv[t,:,:]
    m = oceano.make_map(ax,llat=llat,ulat=ulat,llon=llon,ulon=ulon,resolution='f')

    pc = m.contourf(lon,lat,k2[t,:,:],contour_levels,latlon=True,cmap=cmo.cm.matter)
    qw = m.quiver(lon[skipwind],lat[skipwind],uw[skipwind],vw[skipwind],latlon=True,alpha=.3,scale=150,width=0.005,pivot='middle')

    # definir quantidade de marcadores da colorbar - muito util para valores com muitos zeros =)
    ticks_cb  = np.linspace(contour_levels[0],contour_levels[-1],6,endpoint=True)
    #ticks_cb  = [0.,1e10+3,1e10+4,1e10+5,1e10+6,1e10+7,1e10+8,1e10+9,1e10+10]

    cbaxes = fig.add_axes([0.2, 0.82, 0.1, 0.01])

    # plotar a colorbar, usando os ticks estabelecidos e o tamanho definido acima.
    cb=plt.colorbar(pc,orientation='horizontal',ticks=ticks_cb,format='%.1d',cax=cbaxes)
    # plotar a legenda da colorbar
    cb.set_label('Integrated Tritium Concentration \n (x$10^{7}Bq.m^{-3}$)',fontsize=10)
    cb.ax.tick_params(labelsize=10)

    plt.suptitle(str(time[t]),fontsize=24)

    # saving figure
    outName = str(t).zfill(5)+'.png'
    plt.savefig('/media/danilo/Danilo/mestrado/artigo_data/simulacoes/figures_sims/wcptec/'+outName)
    plt.close('all')
