#!/usr/bin/env python2
#-*-coding:utf-8-*-
'''
Rotina para plotar dispersão da pluma de trítio integrada na coluna em 4 instantes

PLOTA AS 4 IMAGENS EM UMA MESMA FIGURA

Cenário w50 e w210

'''

%reset -f
#################################################
#   MÓDULOS IMPORTAÇÃO/EXPORTAÇÃO ARQUIVOS      #
#################################################
import xray as xr # ler netcdf
from Tkinter import Tk
from tkFileDialog import askopenfilename,askdirectory
import pickle

#################################################
#           MÓDULOS GRÁFICOS                    #
#################################################
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import cmocean as cmo
import matplotlib.ticker as ticker
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib


#################################################
#           MÓDULOS CIENTÍFICOS                 #
#################################################
import scipy.io as sio
import numpy as np
import pandas as pd

#################################################
#           MÓDULOS MISCELANEOUS                #
#################################################
import os
import re
import sys

reload(sys)
sys.setdefaultencoding('utf8')

#### Criando o mapa da região
def make_map(location='bigmap.p'):
    # definindo coordenadas para plotar a região desejada
    # llon = -44.9;
    # ulon = -43.45;
    # llat = -23.6;
    # ulat = -22.75;
    fig, ax = plt.subplots(figsize=(12,8))

    # m=Basemap(projection='merc',llcrnrlat=llat,urcrnrlat=ulat,llcrnrlon=llon,urcrnrlon=ulon,resolution='f')
    # pra agilizar o processo, puxa-se um mapa já feito e salvo em pickle
    m=pickle.load(open("/media/danilo/Danilo/TG/from_tripoli/rotinas/"+location,"rb"))

    m.ax = ax
    return fig, m

def make_minimap():
    # create small region figure

    m2=pickle.load(open("/media/danilo/Danilo/TG/from_tripoli/rotinas/smallmap.p","rb"))

    # #### plotar batimetria
    # lon3,lat3=m2(lon,lat)
    # contorno = np.arange(-80,0,0.1)
    # m2.contourf(lon3,lat3,-depth,contorno,cmap=cmo.cm.bathy)

    # plotar outras coisas do mapa
    m2.drawcoastlines(linewidth=0.4) #linha de costa em alta resolução
    #m2.drawmapboundary(fill_color='#e5f2ff') # fill_color colore o oceano
    m2.fillcontinents(color='#D8D8D8',lake_color='#99ffff') # colorir o continente

    return m2


def extract_data():
    # entra no diretório em que o código foi rodado
    os.chdir(directory)
    # pede para usuário selecionar o arquivo que deseja carregar
    filename = askopenfilename()
    # carrega o arquivo para a variável f_xr
    f_xr = xr.open_dataset(filename)
    # extrai informações importantes (longitude, latitude e tempo)
    lon = f_xr['lon'].data
    lat = f_xr['lat'].data
    time = f_xr['time'].data
    ## ======== RETIRANDO LAT E LON IGUAL A 0 ====

    ind=np.where(lon==0)
    lon[ind]=np.nan
    ind=np.where(lat==0)
    lat[ind]=np.nan

    return lon,lat,time,f_xr


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

################################################################################################################################################################
#                                           PROGRAMA PRINCIPAL                                                                                                 #
################################################################################################################################################################

# determina qual o diretório atual (muito importante para não precisar rodar a rotina somente da pasta em que ela está)
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
directory = os.getcwd()
outputDir = askdirectory(initialdir='/media/danilo/Danilo/mestrado/github/artigoTG/figures/',title='Selecione o diretório para salvar as imagens:')

# selecionando o titulo fixo dos plots
# txt = 'Dispersão de $^{3}H$ forçado pelo:\n'+ r'Vento (de NE $5 m s^{-1}$), Maré (TPXO) e Descarga Fluvial - '
txt = ['Dispersão da pluma de $^{3}H$ na superfície - Após %s horas do vazamento']
# txt = selectTitle(directory)

#### extrair dados
lon,lat,time,f_xr = extract_data()

# extrair concentracao, mas tomando a soma da coluna e não apenas superficial
concInt = f_xr['conc'].data.sum(axis=1)

# remover f_xr pra liberar espaço na memória
del f_xr

# converter o dado de Kg para Bq/l
k=convertData2Bq(concInt)
# liberar espaço na memória
# del c
# normalizar os valores
k=k/.1e+7

# linhas de contorno para usar no contourf
maxC =np.nanmax(k[-1,:,:])+45
contour_levels=np.arange(0.,150,maxC/1e+3)

# contour_levels=np.arange(0.,0.01,0.00005)

max_img = time.shape[0] - 1
cont_time=np.arange(0,max_img+1,1)

# definir o vetor de dia: de 0 a 10 dias a cada 10./120 tempo (a cada 2 horas aproximadamente)
# timeD = np.arange(0,10,10./120)
t0    = pd.to_datetime(str(time[0])).strftime('%H')
t1    = pd.to_datetime(str(time[1])).strftime('%H')
dtImg = int(t1) - int(t0)
timeH = np.arange(dtImg,(time.shape[0])*dtImg+dtImg,dtImg)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################


def plotar_dado(fig,ax,k1,tmpDecor,horas,cbPlot=False,location='bigmap.p'):
    #### plotar dados
    #fig,m=make_map(location=location)
    m=pickle.load(open("/media/danilo/Danilo/mestrado/github/artigoTG/rotinas/"+location,"rb"))
    m.ax = ax

    lon2,lat2=m(lon,lat)

    m.drawcoastlines(linewidth=.4)
    m.fillcontinents(color='#D8D8D8',lake_color='#99ffff')

    # title_out = txt+datetime2title(time[i])
    # plt.title(txt,fontsize=15)

    parallels=np.arange(-23.6,-22.75,0.2)
    m.drawcoastlines(linewidth=0.4)
    m.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.1,fontsize=10)
    meridians=np.arange(-44.9,-33.4,0.5)
    m.drawmeridians(meridians,labels=[1,0,0,1],linewidth=0.1,fontsize=10)
    # plotando as informações de dispersão do material radioativo
    pc=m.contourf(lon2,lat2,k1[:,:],contour_levels,cmap=cmo.cm.matter,extend='both')

    if cbPlot:
        # definir quantidade de marcadores da colorbar - muito util para valores com muitos zeros =)
        # ticks_cb  = np.linspace(contour_levels[0],contour_levels[-1],7,endpoint=True)
        ticks_cb  = [0,30,60,90,120,150]

        #  ordem do add_axes [left, bottom, width, height]
        cbaxes = fig.add_axes([0.335, 0.05, 0.36, 0.02])

        # plotar a colorbar, usando os ticks estabelecidos e o tamanho definido acima.
        cb=plt.colorbar(pc,orientation='horizontal',ticks=ticks_cb,format='%.1d',cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Integrated Concentration (x$10^{7}Bq.m^{-3}$)',fontsize=10)
        cb.ax.tick_params(labelsize=8)

    #########################################################
    #              Plotar dados dentro do mapa              #
    # #######################################################
    # ### datetime
    # tempo = '%s horas após \n vazamento' % txt
    # plt.text(.18,.75,tempo,ha='center',va='center',transform=m.ax.transAxes,fontsize=20)
    if horas:
        tempo = 'After \n %s hours'%(tmpDecor)
        plt.text(0.9,0.85,tempo,ha='center',va='center',transform=m.ax.transAxes,fontsize=10)
    else:
        tempo = 'After \n %s days'%(tmpDecor)
        plt.text(0.9,0.85,tempo,ha='center',va='center',transform=m.ax.transAxes,fontsize=10)

    #########################################################
    #          Plotar mapa detalhado da região              #
    #########################################################
    axins = zoomed_inset_axes(m.ax, 2.5, loc=4)
    axins.set_xlim(-10,0)
    axins.set_ylim(3,10)

    plt.xticks(visible=False)
    plt.yticks(visible=False)

    # create small figure of detailed region
    m2=make_minimap()

    #### plotar dado no mapa pequeno
    lon2,lat2=m2(lon,lat)
    pc=m2.contourf(lon2,lat2,k1[:,:],contour_levels,cmap=cmo.cm.matter,extend='both')

    plt.subplots_adjust(hspace=0.10) # ajustar a distancia vertical entre os gráficos pra ficar mais compacto



os.system('clear')

# 2 horas, 12 horas, 3 dias, 7 dias
print('EVOLUÇÃO DA PLUMA ... Plotando imagens para 0, 30 e 60 dias')
# instantes = [0,39,119,159,199,239]
instantes = [0,120,-1]

fig,ax = plt.subplots(nrows=3,ncols=1,figsize=(10,20))

#def plotar_dado(fig,ax,k1,location='bigmap.p',cbPlot=False):
k1 = k[instantes[0],:,:]
plotar_dado(fig,ax[0],k1,tmpDecor=str(timeH[instantes[0]]),horas=True,cbPlot=True)

k1 = k[instantes[1],:,:]
plotar_dado(fig,ax[1],k1,tmpDecor=str(timeH[instantes[1]]/24),horas=False,cbPlot=False)

k1 = k[instantes[2],:,:]
plotar_dado(fig,ax[2],k1,tmpDecor=str(timeH[instantes[2]]/24),horas=False,cbPlot=False)

# k1 = k[instantes[3],:,:]
# plotar_dado(fig,ax[1][1],k1,tmpDecor=str(timeH[instantes[3]]/24),horas=False,cbPlot=False)

# plt.savefig(outputDir+'/fig5.png')
plt.show()
