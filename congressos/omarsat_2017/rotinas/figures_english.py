"""
Rotina para gerar imagens do resumo expandido para OMARSAT 2017

Figura 1:
  comparação da corrente média nos três cenários de forçantes individuais

Figura 2:
  comparação da dispersão do poluente no instante de 40 dias nos cenários de
  dispersão com vento constante (de SW e de NE)

Figura 3:
  arrumar colorbar e outros aspectos da imagem de evolução da pluma de material
  radioativo no domínio estudado

"""

# importar pacotes
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

# plt.style.use('ggplot')

#################################################
#           MÓDULOS CIENTÍFICOS                 #
#################################################
import scipy.io as sio
import numpy as np
import pandas as pd
import pickle

#################################################
#           MÓDULOS MISCELANEOUS                #
#################################################
import os
import re
import sys
reload(sys)
sys.setdefaultencoding('utf8')

os.system('clear')

def extract_data(f_xr):
    # extrai informações importantes (longitude, latitude e tempo)
    lon = f_xr['lon'].data
    lat = f_xr['lat'].data
    time = f_xr['time'].data
    ## ======== RETIRANDO LAT E LON IGUAL A 0 ====

    ind=np.where(lon==0)
    lon[ind]=np.nan
    ind=np.where(lat==0)
    lat[ind]=np.nan

    return lon,lat,time

def plotar_ventocorrente(fig,ax,item,u,v,spd,curint,location='bigmap.p',save=False,cbPlot=False):
    #########################################################
    # PROPRIEDADE DO MAPA (Linha de costa, estados e outros)#
    #########################################################
    #fig,m=make_map(location=location)
    m=pickle.load(open("/home/tparente/danilo/TG/rotinas/"+location,"rb"))
    m.ax = ax

    lon2,lat2=m(lon,lat)

    m.fillcontinents(color='#D8D8D8',lake_color='#99ffff')
    # title_out = txt+datetime2title(time[i])
    # plt.title(txt,fontsize=16)
    parallels=np.arange(-23.6,-22.75,0.2)
    m.drawcoastlines(linewidth=0.4)
    m.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.1,fontsize=10)
    meridians=np.arange(-44.9,-33.4,0.5)
    m.drawmeridians(meridians,labels=[1,0,0,1],linewidth=0.1,fontsize=10)

    #### criar o padrão de plots dos vetores
    skipcurr = (slice(None,None,curint),slice(None,None,curint))

    #########################################################
    #         PLOTAR CORRENTE  (QUIVER E CONTOURF)          #
    #########################################################
    cs1 = m.contourf(lon2,lat2,spd[:,:],contour_levels,cmap=cmo.cm.speed,extend='max')
    # ticks_cb = np.linspace(contour_levels[0],contour_levels[-1],7,endpoint=True)
    if cbPlot:
        ticks_cb = np.arange(0,0.8,0.1)
        # brackets[left, bottom, width, height]
        cbaxes = fig.add_axes([0.91, 0.365, 0.01, 0.29])
        cb = plt.colorbar(cs1, orientation='vertical',ticks=ticks_cb,format='%1.1f', cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Mean Current [m.s$^{-1}$]',fontsize=10,labelpad=-1)
        cb.ax.tick_params(labelsize=10)

    maxcurr = 'Max curr: \n %.2fm.s$^{-1}$'%(np.nanmax(spd))
    plt.text(0.9,0.90,maxcurr,ha='center',va='center',transform=m.ax.transAxes,fontsize=10)

    q=m.quiver(lon2[skipcurr],lat2[skipcurr],u[skipcurr],v[skipcurr],alpha=.7,scale=30,minshaft=2)

    # plot number subplot
    plt.text(0.9, 0.1, item, ha='center', va='center', transform=m.ax.transAxes, fontsize=10)

def tratarVento(wu,wv):
    indu = np.where(wu == wu.min())
    wu[indu] = np.nan
    wv[indu] = np.nan
    indv = np.where(wv == wv.min())
    wu[indv] = np.nan
    wv[indv] = np.nan
    wndspd = np.sqrt(wu*wu + wv*wv)
    wu     = wu/2
    wv     = wv/2
    wndspd = wndspd/2

    return wu,wv,wndspd

def tratarCorrente(u,v):
    ## substituir o menor valor da componente por np.nan, pois este é um flag de dado ruim
    ind = np.where(u == u.min())
    u[ind] = np.nan
    v[ind] = np.nan
    # media da corrente
    mu = u.mean(axis=0)
    mv = v.mean(axis=0)
    # gerar vetor de velocidade
    spd = np.sqrt(mu*mu + mv*mv)
    # normalizar vetores da corrente
    mun = mu/spd
    mvn = mv/spd

    del u,v,mu,mv,ind

    return mun,mvn, spd

def loadVelocity(BASE_DIR):
    """ load velocity data (components and speed) from pickle file

    important: normalized components by speed
    """
    data = pickle.load(open(BASE_DIR + 'currComp.pkl','rb'))
    mun = data['mun']
    mvn = data['mvn']
    spd = data['spd']

    U,V,SPD = mun[:,:],mvn[:,:],spd[:,:]

    return U,V,SPD

os.system('clear')

#########################################################
#                PROGRAMA PRINCIPAL                     #
#########################################################
print('Carregando configurações principais')

# diretorios onde estão armazenados os arquivos gcmplt.cdf de cada simulação
TIDE_DIR = '/home/tparente/danilo/mestrado/congressos/omarsat_2017/rotinas/arquivos/Tides/'
W_NE_DIR = '/home/tparente/danilo/mestrado/congressos/omarsat_2017/rotinas/arquivos/NEwind/'
W_SW_DIR = '/home/tparente/danilo/mestrado/congressos/omarsat_2017/rotinas/arquivos/SWwind/'
SAVE_DIR = '/home/tparente/danilo/mestrado/congressos/omarsat_2017/figuras/'

def criar_figura1(savefig=False):

    os.system('clear')
    #########################################################
    #                      FIGURA 1                         #
    #    comparação da corrente média nos três cenários de  #
    #              forçantes individuais                    #
    #########################################################
    print("#########################################################")
    print("#                      FIGURA 1                         #")
    print("#    comparação da corrente média nos três cenários de  #")
    print("#              forçantes individuais                    #")
    print("#########################################################")
    ###### creating fig for subplots
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,10))

    #############################################################################
    print('Select file for tide ...')
    tideFile = xr.open_dataset('/home/tparente/danilo/TG/simulacoes/rel_final/forcantesInd/Tides/gcmplt.cdf')

    # i = 210 # enchente de sizígia: mais intensa
    # uTide,vTide,spdTide = tratarCorrente(tideFile['u'].data[i,:,:,:],tideFile['v'].data[i,:,:,:])
    uTide, vTide, spdTide = loadVelocity(TIDE_DIR)

    ### primeira imagem - SW wind
    maxC = np.nanmax(spdTide)
    contour_levels=np.arange(0.,maxC,maxC/1e+3)
    #contour_levels = pickle.load(open(SAVE_DIR.replace('figuras/', 'rotinas/arquivos/contour.pkl'),'r'))

    print("Plotting tide simulation output ... ")
    # plotar experimento III - Tide
    plotar_ventocorrente(fig,axes[2],item='(c)',u=uTide[:,:],v=vTide[:,:],spd=spdTide[:,:],curint=15,save=False,cbPlot=True)
    axes[2].title.set_text(u'Experiment III - Flood Spring Tide')

    # remove variables
    del uTide, vTide, spdTide

    ############################################################################
    print('Select file for ne wind ... ')

    # load data from simulation
    uNE, vNE, spdNE = loadVelocity(W_NE_DIR)

    print("Plotting NE simulation output ... ")
    # plotar experimento I - w210
    plotar_ventocorrente(fig,axes[0],item='(a)',u=uNE[-1,:,:],v=vNE[-1,:,:],spd=spdNE[-1,:,:],curint=15,save=False,cbPlot=False)
    axes[0].title.set_text(u'Experiment I - Northeasterly Wind [5$ m s^{-1}$]')

    # remove variables
    del uNE, vNE, spdNE

    #############################################################################
    print('Select file for sw wind ... ')

    # load data from simulation
    uSW, vSW, spdSW = loadVelocity(W_SW_DIR)

    print("Plotting SW simulation output ... ")
    # plotar experimento II - w50
    plotar_ventocorrente(fig,axes[1],item='(b)',u=uSW[-1,:,:],v=vSW[-1,:,:],spd=spdSW[-1,:,:],curint=15,save=False,cbPlot=False)
    axes[1].title.set_text(u'Experiment II - Southwesterly Wind [5$ m s^{-1}$]')

    # remove variables
    del uSW, vSW, spdSW

    if savefig:
        # name's fig
        output1 = SAVE_DIR + 'figura1.png'
        # savefig
        plt.savefig(output1, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output1, output1))
    else:
        plt.show()

# lon and lat data - just one time, is the same in all simulations
dic = pickle.load(open(SAVE_DIR.replace("figuras/", "rotinas/arquivos/lonlat.pkl"), "rb"))
lon, lat, time = dic['lon'], dic['lat'], dic['time']

# criar_figura1(savefig=False)

####################################################################################################################
os.system('clear')

def extract_concentration(fname):
    """ extract concentration data from fname given """

    f_xr = xr.open_dataset(fname)

    c   = f_xr['conc'].data

    lon, lat, time = extract_data(f_xr)

    return lon, lat, time, c

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

def surfaceConcentration():
    """ import surface concentration from pickle file already extracted from .cdf """
    # set path and file from each simulation
    ne = '/home/tparente/danilo/TG/simulacoes/rel_final/w210/conc_sup.pkl'
    sw = '/home/tparente/danilo/TG/simulacoes/rel_final/w50/conc_sup.pkl'

    print("extracting concentration from NE simulation ... ")
    neConc = pickle.load(open(ne, 'rb'))[i,:,:]
    neConc = np.squeeze(neConc)

    print("extracting concentration from SW simulation ... ")
    swConc = pickle.load(open(sw, 'rb'))[i,:,:]
    swConc = np.squeeze(swConc)

    return neConc, swConc

def totalConcentration(i=159, pkl=None):
    """
    import all data of concentration and sum all sigma levels
     """

    if pickle:
        # importar pickle file
        ### IMPORTANT: data already converted to Bq and divide by 1e7
        totalConc = pickle.load(open(pkl,'rb'))
        neConc = totalConc['neConc']
        swConc = totalConc['swConc']

    else:
        print("extracting concentration from NE simulation and summing all sigma levels ... ")
        fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w210/vento_dobrado.cdf'
        f_xr = xr.open_dataset(fname)
        conc = np.squeeze(f_xr['conc'].data[i,:,:,:])

        cTotal = conc.sum(axis=0)

        neConc = cTotal

        print("extracting concentration from SW simulation and summing all sigma levels ... ")
        fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w50/vento_dobrado.cdf'
        f_xr = xr.open_dataset(fname)
        conc = np.squeeze(f_xr['conc'].data[i,:,:,:])

        cTotal = conc.sum(axis=0)
        swConc = cTotal

        # converting concentration to a friendly format
        neConc = convertData2Bq(neConc)/.1e7
        swConc = convertData2Bq(swConc)/.1e7

    return neConc, swConc

def plotar_dado2(fig,ax,k1,item,horas,cbPlot=False,location='bigmap.p'):
    #### plotar dados
    #fig,m=make_map(location=location)
    m=pickle.load(open("/home/tparente/danilo/TG/rotinas/"+location,"rb"))
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
        ticks_cb = [0, 40, 80, 120, 160, 200, 240 ]

        cbaxes = fig.add_axes([0.91, 0.29, 0.01, 0.42])

        # plotar a colorbar, usando os ticks estabelecidos e o tamanho definido acima.
        cb=plt.colorbar(pc,orientation='vertical',ticks=ticks_cb,format='%d',cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Interated Concentration of Tritium \n (x$10^{7}Bq.m^{-3}$)',fontsize=15)
        cb.ax.tick_params(labelsize=13)

    # plot number subplot
    plt.text(0.9, 0.1, item, ha='center', va='center', transform=m.ax.transAxes, fontsize=10)

    return pc

#########################################################
#                       FIGURA 2                        #
# comparação da dispersão do poluente no instante de    #
# 40 dias nos cenários de dispersão com vento constante #
#                  (de SW e de NE)                      #
#########################################################
def criar_figura2(savefig=False):
    os.system('clear')
    print("#########################################################")
    print("#                       FIGURA 2                        #")
    print("# comparação da dispersão do poluente no instante de    #")
    print("# 40 dias nos cenários de dispersão com vento constante #")
    print("#                  (de SW e de NE)                      #")
    print("#########################################################")

    # set timestep for 40 days of simulation
    i = [159] # 40 dias

    # extract data
    neConc, swConc = totalConcentration(i, pkl="/home/tparente/Dropbox/mestrado/artigoTG/data/totalConc.pkl")

    # linhas de contorno para usar no contourf
    maxC = np.nanmax(neConc[:,:])+55
    contour_levels=np.arange(0.,round(maxC),round(maxC)/1e+3)

    # mandar plotar os subplots

    fig,ax = plt.subplots(nrows=1,ncols=2,figsize=(12,8))

    # primeira imagem com vento de SW
    pc = plotar_dado2(fig,ax[0],neConc,item='(a)',horas=False,cbPlot=True)
    ax[0].set_title(u'Experimento IV - Vento de NE [5$ m s^{-1}$]')
    # segunda imagem com vento de NE
    pc = plotar_dado2(fig,ax[1],swConc,item='(b)',horas=False,cbPlot=False)
    ax[1].set_title(u'Experimento IV - Vento de SW [5$ m s^{-1}$]')

    plt.suptitle(u"Evolução da Pluma de Dispersão (Integrada) de Trítio após 40 dias de simulação", fontsize=25, y=0.78)

    if savefig:
        output = SAVE_DIR + "figura2_integrada.png"
        plt.savefig(output, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output, output))
    else:
        plt.show()

    del neConc, swConc, pc

# criar_figura2(savefig=False)

#########################################################
#                       FIGURA 3                        #
#     arrumar colorbar e outros aspectos da imagem de   #
#       evolução da pluma de materialradioativo no      #
#                   domínio estudado                    #
#########################################################
def make_minimap():
    # create small region figure

    m2=pickle.load(open("/home/tparente/danilo/TG/rotinas/smallmap.p","rb"))

    # #### plotar batimetria
    # lon3,lat3=m2(lon,lat)
    # contorno = np.arange(-80,0,0.1)
    # m2.contourf(lon3,lat3,-depth,contorno,cmap=cmo.cm.bathy)

    # plotar outras coisas do mapa
    m2.drawcoastlines(linewidth=0.4) #linha de costa em alta resolução
    #m2.drawmapboundary(fill_color='#e5f2ff') # fill_color colore o oceano
    m2.fillcontinents(color='#D8D8D8',lake_color='#99ffff') # colorir o continente

    return m2

def plotar_dado3(fig,ax,k1,tmpDecor,cbPlot=False,location='bigmap.p'):
    #### plotar dados
    #fig,m=make_map(location=location)
    m=pickle.load(open("/home/tparente/danilo/TG/rotinas/"+location,"rb"))
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
    pc=m.contourf(lon2,lat2,k1[:,:],contour_levels,cmap=cmo.cm.matter,extend='max')

    if cbPlot:
        # definir quantidade de marcadores da colorbar - muito util para valores com muitos zeros =)
        ticks_cb  = np.linspace(contour_levels[0],contour_levels[-1],7,endpoint=True)

        cbaxes = fig.add_axes([0.15, 0.87, 0.1, 0.01])

        # plotar a colorbar, usando os ticks estabelecidos e o tamanho definido acima.
        cb=plt.colorbar(pc,orientation='horizontal',ticks=ticks_cb,format='%.1d',cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Concentração de Trítio \n (x$10^{7}Bq.m^{-3}$)',fontsize=10)
        cb.ax.tick_params(labelsize=10)

    #########################################################
    #              Plotar dados dentro do mapa              #
    # #######################################################
    # ### datetime
    # tempo = '%s horas após \n vazamento' % txt
    # plt.text(.18,.75,tempo,ha='center',va='center',transform=m.ax.transAxes,fontsize=20)
    tempo = 'Após \n %s' % (tmpDecor)
    plt.text(0.9, 0.9, tempo, ha='center', va='center', transform=m.ax.transAxes, fontsize=10)

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

    return pc

def extractConc(i,data):
    # importar somente os dados do instante i
    conc = data['conc'].data[i,:,:,:]
    # somar na coluna
    conc = conc.sum(axis=0)
    # converter para bq
    k = convertData2Bq(conc)/.1e7

    return k

def criar_figura3(savefig=False):
    os.system('clear')

    print("#########################################################")
    print("#                       FIGURA 3                        #")
    print("#     arrumar colorbar e outros aspectos da imagem de   #")
    print("#       evolução da pluma de materialradioativo no      #")
    print("#                   domínio estudado                    #")
    print("#########################################################")

    fname     = '/home/tparente/danilo/TG/simulacoes/rel_final/wcptecJULAGO/gcmplt.cdf' # definir arquivo gcmplt.cdf a ser lido
    data = xr.open_dataset(fname)               # import netcdf file
    conc = data['conc'].data[120:,0,:,:]        # import only concentration data
                                                # cut period without pollution
    # conc = conc.sum(axis=1)                     # sum each levels data
    conc = convertData2Bq(conc)/.1e7            # convert to Bq

    lon, lat, time = extract_data(data)         # extract lon,lat,time info
    time = time[120:]                           # cut time without leaking

    instantes = [12, 120, 240, 360]              # referentes a 6 horas e 10,
                                                # 20 e 30 dias
    instTitle = ['01 dia', '10 dias', '20 dias', '30 dias']# titulo do instante plotado

    # # save all concentrations data into a pickle file
    # fpickle = SAVE_DIR.replace("figuras/", "rotinas/arquivos/concentracao.pkl")
    # pickle.dump(conc, open(fpickle, 'wb'))
    #
    # concentracoes = {}
    #
    # for i,title in zip(instantes, instTitle):
    #     c = conc[i,:,:]
    #     concentracoes.update({title: c})
    #
    # fpickle = SAVE_DIR.replace("figuras/", "rotinas/arquivos/plot_concentracao.pkl")
    # pickle.dump(concentracoes, open(fpickle, 'wb'))

    # criar contour_levels baseado no instante de maior concentração
    # maxC =np.nanmax(conc[-1,:,:])+45
    # contour_levels = np.asarray([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70])    # contour levels of contourf

    maxC = np.nanmax(conc[-1,:,:])+55
    contour_levels=np.arange(0.,70,70/1e+3)

    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(12,8))

    k = conc[instantes[0]+3,:,:]
    pc1 = plotar_dado3(fig,ax[0,0],k,tmpDecor=instTitle[0],cbPlot=False)

    k = conc[instantes[1]+3,:,:]
    pc2 = plotar_dado3(fig,ax[0,1],k,tmpDecor=instTitle[1],cbPlot=False)

    k = conc[instantes[2]+3,:,:]
    pc3 = plotar_dado3(fig,ax[1,0],k,tmpDecor=instTitle[2],cbPlot=False)

    k = conc[instantes[3]+3,:,:]
    pc4 = plotar_dado3(fig,ax[1,1],k,tmpDecor=instTitle[3],cbPlot=False)


    ticks_cb = np.asarray([0, 10, 20, 30, 40, 50, 60, 70])
    cb = plt.colorbar(pc4, ax=ax.ravel().tolist(), ticks=ticks_cb, extend='max') # inserir colorbar ocupado o eixo y completo
    cb.set_label('Concentração de Trítio (x$10^{7}Bq.m^{-3}$)',fontsize=20)

    plt.suptitle(u"Evolução (em superfície) da Pluma de Dispersão de Trítio liberado em 01/08/2016", fontsize=25, y=0.94)

    if savefig:
        output = SAVE_DIR + "figura3_superficie.png"
        plt.savefig(output, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output, output))
    else:
        plt.show()

# criar_figura3(savefig=False)

#########################################################
#                       FIGURA 4                        #
#      plotar campo de corrente médio gerado nos        #
#   cenários V e VI (maré, descarga fluvial e vento)    #
#########################################################
def load_data4(f_xr):
    u = f_xr['u'].values[:,0,:,:]
    v = f_xr['v'].values[:,0,:,:]

    mun,mvn,spd = tratarCorrente(u,v)

    return mun,mvn,spd


def load_pickle(nfile,u='u',v='v'):
    data = pickle.load(open(nfile,'r'))

    return data[u],data[v]

def criar_figura4(savefig=False):

    os.system('clear')
    #########################################################
    #                      FIGURA 1                         #
    #    comparação da corrente média nos três cenários de  #
    #              forçantes individuais                    #
    #########################################################
    print("#########################################################")
    print("#                      FIGURA 4                         #")
    print("#      plotar campo de corrente médio gerado nos        #")
    print("#   cenários V e VI (maré, descarga fluvial e vento)    #")
    print("#########################################################")


    ###### creating fig for subplots
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20,10))
    #############################################################################
    print('Select file for sw wind ... ')

    # load data from simulation
    fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w50/vento_dobrado.cdf'
    #f_xr = xr.open_dataset(fname)
    #uSW, vSW, spdSW = load_data4(f_xr)
    #uSW,vSW = load_pickle('/home/tparente/Dropbox/mestrado/artigoTG/data/expVI_current.pickle',u='uSW',v='vSW')
    #spdSW = np.sqrt(uSW**2 + vSW**2)
    maxC = np.nanmax(spdSW)
    contour_levels=np.arange(0.,maxC,0.01)

    print("Plotting SW simulation output ... ")
    # plotar experimento II - w50
    plotar_ventocorrente(fig,axes[1],item='(b)',u=uSW[:,:],v=vSW[:,:],spd=spdSW[:,:],curint=15,save=False,cbPlot=True)
    axes[1].title.set_text(u'Experiment VI - Southwesterly Wind [5$ m s^{-1}$]')

    ############################################################################
    print('Select file for ne wind ... ')
    # load data from simulation
    fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w210/vento_dobrado.cdf'
    #f_xr = xr.open_dataset(fname)
    #uNE, vNE, spdNE = load_data4(f_xr)
    #uNE,vNE = load_pickle('/home/tparente/Dropbox/mestrado/artigoTG/data/expV_current.pickle',u='uNE',v='vNE')
    #spdNE = np.sqrt(uNE**2 + vNE**2)
    
    print("Plotting NE simulation output ... ")
    # plotar experimento I - w210
    plotar_ventocorrente(fig,axes[0],item='(a)',u=uNE[:,:],v=vNE[:,:],spd=spdNE[:,:],curint=15,save=False,cbPlot=False)
    axes[0].title.set_text(u'Experiment V - Northeasterly Wind [5$ m s^{-1}$]')

    plt.suptitle('Surface Current Generated by Wind, Tide \n and Fluvial Discharge',fontsize=30)


    # remove variables
    del uNE, vNE, spdNE
    # remove variables
    del uSW, vSW, spdSW

    if savefig:
        # name's fig
        output1 = SAVE_DIR + 'figura1.png'
        # savefig
        plt.savefig(output1, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output1, output1))
    else:
        plt.show()

#########################################################
#                       FIGURA 5                        #
#      plotar campo de corrente médio gerado nos        #
#   cenários V e VI (maré, descarga fluvial e vento)    #
#   e a dispersão integrada nos mesmos cenários         #
#########################################################
def criar_figura5(savefig=False):

    # importar dados de corrente já tratados e em arquivo pickle
    return False
    # importar dados de concentração já tratados e em arquivo pickle

os.system('clear')

choose = input("Selecione qual figura voce deseja plotar: [1], [2] ou [3]? ")
savefig = input("Deseja salvar ou apenas visualizar a imagem: [1] - salvar, [2] visualizar ? ")

if savefig == 1:
    savefig = True
else:
    savefig = False

if choose == 1:
    criar_figura1(savefig=savefig)

elif choose == 2:
    criar_figura2(savefig=savefig)

elif choose == 3:
    criar_figura3(savefig=savefig)
elif choose == 4:
    criar_figura4(savefig=savefig)