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

def plotar_ventocorrente(fig,contour_levels,ax,item,u,v,spd,curint,location='bigmap.p',save=False,cbPlot=False,axesPos=[0.91, 0.365, 0.01, 0.29]):
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
        cbaxes = fig.add_axes(axesPos)
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

# def criar_figura3(savefig=False):
#     os.system('clear')

#     print("#########################################################")
#     print("#                       FIGURA 3                        #")
#     print("#     arrumar colorbar e outros aspectos da imagem de   #")
#     print("#       evolução da pluma de materialradioativo no      #")
#     print("#                   domínio estudado                    #")
#     print("#########################################################")

#     fname     = '/home/tparente/danilo/TG/simulacoes/rel_final/wcptecJULAGO/gcmplt.cdf' # definir arquivo gcmplt.cdf a ser lido
#     data = xr.open_dataset(fname)               # import netcdf file
#     conc = data['conc'].data[120:,0,:,:]        # import only concentration data
#                                                 # cut period without pollution
#     # conc = conc.sum(axis=1)                     # sum each levels data
#     conc = convertData2Bq(conc)/.1e7            # convert to Bq

#     lon, lat, time = extract_data(data)         # extract lon,lat,time info
#     time = time[120:]                           # cut time without leaking

#     instantes = [12, 120, 240, 360]              # referentes a 6 horas e 10,
#                                                 # 20 e 30 dias
#     instTitle = ['01 dia', '10 dias', '20 dias', '30 dias']# titulo do instante plotado

#     # # save all concentrations data into a pickle file
#     # fpickle = SAVE_DIR.replace("figuras/", "rotinas/arquivos/concentracao.pkl")
#     # pickle.dump(conc, open(fpickle, 'wb'))
#     #
#     # concentracoes = {}
#     #
#     # for i,title in zip(instantes, instTitle):
#     #     c = conc[i,:,:]
#     #     concentracoes.update({title: c})
#     #
#     # fpickle = SAVE_DIR.replace("figuras/", "rotinas/arquivos/plot_concentracao.pkl")
#     # pickle.dump(concentracoes, open(fpickle, 'wb'))

#     # criar contour_levels baseado no instante de maior concentração
#     # maxC =np.nanmax(conc[-1,:,:])+45
#     # contour_levels = np.asarray([0, 2, 4, 6, 8, 10, 20, 30, 40, 50, 60, 70])    # contour levels of contourf

#     maxC = np.nanmax(conc[-1,:,:])+55
#     contour_levels=np.arange(0.,70,70/1e+3)

#     fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(12,8))

#     k = conc[instantes[0]+3,:,:]
#     pc1 = plotar_dado3(fig,ax[0,0],k,tmpDecor=instTitle[0],cbPlot=False)

#     k = conc[instantes[1]+3,:,:]
#     pc2 = plotar_dado3(fig,ax[0,1],k,tmpDecor=instTitle[1],cbPlot=False)

#     k = conc[instantes[2]+3,:,:]
#     pc3 = plotar_dado3(fig,ax[1,0],k,tmpDecor=instTitle[2],cbPlot=False)

#     k = conc[instantes[3]+3,:,:]
#     pc4 = plotar_dado3(fig,ax[1,1],k,tmpDecor=instTitle[3],cbPlot=False)


#     ticks_cb = np.asarray([0, 10, 20, 30, 40, 50, 60, 70])
#     cb = plt.colorbar(pc4, ax=ax.ravel().tolist(), ticks=ticks_cb, extend='max') # inserir colorbar ocupado o eixo y completo
#     cb.set_label('Concentração de Trítio (x$10^{7}Bq.m^{-3}$)',fontsize=20)

#     plt.suptitle(u"Surface Evoluion of Tritium Released in Aug, 01", fontsize=25, y=0.94)
#     #plt.suptitle(u"Evolução (em superfície) da Pluma de Dispersão de Trítio liberado em 01/08/2016", fontsize=25, y=0.94)

#     if savefig:
#         output = SAVE_DIR + "figura3_superficie.png"
#         plt.savefig(output, dpi=150)
#         # let's clean this figure
#         os.system('convert -trim %s %s' % (output, output))
#     else:
#         plt.show()
#########################################################
#                       FIGURA 1                        #
#       Area de Estudo                                  #
#########################################################
# desenhar um circulo tracejado em um mapa
def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

# desenhar um circulo tracejado em um mapa
def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)

def create_figure1(savefig=False):

    # load grid
    data = xr.open_dataset('/home/tparente/danilo/TG/simulacoes/rel_final/forcantesInd/Tides/gcmplt.cdf')
    # extract variables
    depth = data['depth'].values
    lon = data['lon'].values
    lat = data['lat'].values
    lon[lon == 0.] = np.nan
    lat[lon == 0.] = np.nan

    # definindo coordenadas para plotar a região desejada
    llon = -44.9;
    ulon = -43.45;

    llat = -23.6;
    ulat = -22.75;

    TermGuaiba = [-44.01916667, -23.00000000, '^', 'black']
    ColNaval  =  [-44.31833333, -23.01694444, 'v', 'black']

    ######### DETALHES PARA PLOTAR
    CNAAAbase     = [-44.460789, -23.006711]
    ilhagrande    = [-44.613300, -23.141362]#-23.137504, -44.613300
    sepetiba      = [-43.864138, -23.003667]
    piraquaraf    = [-44.441523, -23.011983]
    itaorna       = [-44.472463, -23.006116]#-23.017116,
    porOeste      = [-44.508007, -23.194821]
    porLeste      = [-44.113577, -23.070955]
    canCentral    = [-44.333599, -23.083153]
    picinguaba    = [-44.830344, -23.378115]
    guaratiba     = [-43.561320, -23.064399]
    riomambuca    = [-44.521393, -23.026489]
    angradosreis  = [-44.316713, -22.986787]

    marambaia     = [-43.989471, -23.083472]
    mangaratiba   = [-44.056691, -22.944938]
    paraty        = [-44.757599, -23.232831]
    angradosreis  = [-44.358047, -23.020794]

    # tarituba      = [-44.598025, -23.039031]
    ## pontos de descarga/captação
    fozMambucaba  = [-44.521247, -23.026884]
    descargaCNAAA = [-44.445079, -23.012130]
    captacaoCNAAA = [-44.459445, -23.009400]

    ### quadrado
    lat0 = -22.980121
    lat1 = -23.054853
    lon0 = -44.374144
    lon1 = -44.558726

    fig = plt.figure(figsize=(20,15))
    ax = fig.add_subplot(111)

    #fig,m=make_map(location=location)
    m=pickle.load(open("/home/tparente/danilo/TG/rotinas/bigmap.p","rb"))
    m.ax = ax

    m.drawcoastlines(linewidth=0.2) #linha de costa em alta resolução
    m.drawmapboundary() # fill_color colore o oceano
    m.fillcontinents(color='#ffd480') # colorir o continente
    m.drawrivers()
    # definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(-44.9,-33.4,0.25)
    parallels=np.arange(-23.6,-22.75,0.1)
    # desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,
                                        fontweight='bold',color='gray',linewidth=0.5)
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,
                                        fontweight='bold',color='gray',linewidth=0.5)

    m.drawmapscale(-43.65, -23.55, -43.65, -23.55, 20, barstyle='fancy', yoffset=1000)
    # batimetria em forma de grade
    pc = m.pcolor(lon,lat,depth,latlon=True,cmap=cmo.cm.deep)

    # plotar detalhes em texto
    x,y=m(CNAAAbase[0], CNAAAbase[0])
    m.plot(x,y, 'bo', markersize=14)

    x,y=m(ilhagrande[0]+.100,ilhagrande[1]+.02)
    ax.text(x,y,u'ILHA GRANDE\nBAY',color='#000000', fontsize=15, ha='center',va='center')

    x,y=m(sepetiba[0]+.080,sepetiba[1])
    ax.text(x,y,u'SEPETIBA\nBAY',color='#000000', fontsize=15, ha='center',va='center')

    x,y=m(porOeste[0],porOeste[1])
    ax.text(x,y,u'WEST\nPORTION',color='#800000', fontsize=12, ha='center',va='center')

    x,y=m(porLeste[0],porLeste[1])
    ax.text(x,y,u'EAST\nPORTION',color='#800000', fontsize=12, ha='center',va='center')

    x,y=m(canCentral[0],canCentral[1])
    ax.text(x,y,u'CENTRAL CHANNEL',color='#800000', fontsize=12, ha='center',va='center')

    x,y=m(picinguaba[0]+0.05, picinguaba[1]+0.015)
    ax.text(x,y,u'Picinguaba',color='k',fontsize=8, ha='center',va='center')

    x,y=m(guaratiba[0]+0.035, guaratiba[1]+0.025)
    ax.text(x,y,u'Barra de\nGuaratiba',color='k',fontsize=8, ha='center',va='center')

    x,y=m(riomambuca[0], riomambuca[1]+0.01)
    ax.text(x,y,u'Mambucaba\nRiver',color='k',fontsize=10, ha='center',va='center')

    x,y=m(angradosreis[0]+0.01, angradosreis[1]+0.03)
    ax.text(x,y,u'Angra dos\nReis', color='k', fontsize=11, ha='center',va='center')

    x,y=m(paraty[0]+0.01, paraty[1]+0.03)
    ax.text(x,y,u'Paraty', color='k', fontsize=12, ha='center',va='center')

    x,y=m(-44.230458, -23.151578)
    ax.text(x,y,u'Ilha Grande',color='k',fontsize=12,ha='center', va='center')


    # -----------------------------------------------------------------------------

    loc_w210 = [[-23.010962, -44.440706], [-22.9741811,-44.323831],
                [-23.012479, -44.307350], [-23.059381, -44.241407],
                [-22.994208,-44.0397950], [-23.088036,-44.0980590],
                [-23.096075,-44.010663]]
    nam_w210 = ['Piraquara de Fora Inlet','Angra dos Reis (A)',
                'Angra dos Reis (B)', 'Canal Central',
                'Mangaratiba','Porção Leste (BIG)',
                'Marambaia']

    x,y = m(loc_w210[0][1],loc_w210[0][0])
    plt.scatter(x,y,marker='*',color='black',s=80,label=nam_w210[0],zorder=2)

    x,y = m(TermGuaiba[0], TermGuaiba[1])
    plt.scatter(x,y, marker=TermGuaiba[2], color=TermGuaiba[3], s=60, label='Guaiba Island Terminal', zorder=2)

    x,y = m(ColNaval[0], ColNaval[1])
    plt.scatter(x,y, marker=ColNaval[2], color=ColNaval[3], s=60, label='Naval College', zorder=2)


    kwargs = {'linestyle': '--', 'alpha':0.4, 'color': 'black'}
    equi(m, CNAAAbase[0], CNAAAbase[1], 15.,lw=1., **kwargs)
    # create legend with semi transparent background
    lg = plt.legend(loc='upper left',fontsize=12,numpoints=1,scatterpoints=1)
    lg.get_frame().set_alpha(.4)

    cbar = plt.colorbar(pc, orientation='horizontal', shrink=0.625, aspect=20, fraction=0.2,pad=0.04)
    cbar.set_label('Bathymetry [m]', fontsize=20)

    if savefig:
        # name's fig
        output1 = '/home/tparente/danilo/mestrado/github/congressos/iwmo2018/figs/fig1.png'
        # savefig
        plt.savefig(output1, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output1, output1))
    else:
        plt.show()

#########################################################
#                       FIGURA 3                        #
#     arrumar colorbar e outros aspectos da imagem de   #
#       evolução da pluma de materialradioativo no      #
#                   domínio estudado                    #
#########################################################
def create_figure3(savefig=False):

    os.system('clear')
    #########################################################
    #                      FIGURA 3                         #
    #    comparação da corrente média nos três cenários de  #
    #              forçantes individuais                    #
    #########################################################
    print("#########################################################")
    print("#                      FIGURA 3                         #")
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
    plotar_ventocorrente(fig,contour_levels,axes[2],item='(c)',u=uTide[:,:],v=vTide[:,:],spd=spdTide[:,:],curint=15,save=False,cbPlot=True)
    axes[2].title.set_text(u'Experiment III - Flood Spring Tide')

    # remove variables
    del uTide, vTide, spdTide

    ############################################################################
    print('Select file for ne wind ... ')

    # load data from simulation
    uNE, vNE, spdNE = loadVelocity(W_NE_DIR)

    print("Plotting NE simulation output ... ")
    # plotar experimento I - w210
    plotar_ventocorrente(fig,contour_levels,axes[0],item='(a)',u=uNE[-1,:,:],v=vNE[-1,:,:],spd=spdNE[-1,:,:],curint=15,save=False,cbPlot=False)
    axes[0].title.set_text(u'Experiment I - Northeasterly Wind [5$ m s^{-1}$]')

    # remove variables
    del uNE, vNE, spdNE

    #############################################################################
    print('Select file for sw wind ... ')

    # load data from simulation
    uSW, vSW, spdSW = loadVelocity(W_SW_DIR)

    print("Plotting SW simulation output ... ")
    # plotar experimento II - w50
    plotar_ventocorrente(fig,contour_levels,axes[1],item='(b)',u=uSW[-1,:,:],v=vSW[-1,:,:],spd=spdSW[-1,:,:],curint=15,save=False,cbPlot=False)
    axes[1].title.set_text(u'Experiment II - Southwesterly Wind [5$ m s^{-1}$]')

    plt.suptitle('Mean Current Generated by Wind (Experiments I and II) \n and Tidal Components (Experiment III)', fontsize=23, y=0.78)

    # remove variables
    del uSW, vSW, spdSW

    if savefig:
        # name's fig
        output1 = '/home/tparente/danilo/mestrado/github/congressos/iwmo2018/figs/fig3.png'
        # savefig
        plt.savefig(output1, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output1, output1))
    else:
        plt.show()

#########################################################
#                       FIGURA 4                        #
#      plotar campo de corrente médio gerado nos        #
#   cenários IV e V (maré, descarga fluvial e vento)    #
#    e concentração integrada nos mesmos experimentos   #
#########################################################
def load_data4(f_xr):
    u = f_xr['u'].values[:,0,:,:]
    v = f_xr['v'].values[:,0,:,:]

    mun,mvn,spd = tratarCorrente(u,v)

    return mun,mvn,spd

def load_pickle(nfile,u='u',v='v'):
    data = pickle.load(open(nfile,'r'))

    return data[u],data[v]

def plot_concData(fig,ax,contour_levels,k1,item,horas,cbPlot=False,location='bigmap.p',axesPos=[0.91, 0.09, 0.01, 0.42]):
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
        # ticks_cb  = np.linspace(contour_levels[0],contour_levels[-1],7,endpoint=True)
        ticks_cb = [0, 40, 80, 120, 160, 200, 240 ]

        cbaxes = fig.add_axes(axesPos)

        # plotar a colorbar, usando os ticks estabelecidos e o tamanho definido acima.
        cb=plt.colorbar(pc,orientation='vertical',ticks=ticks_cb,format='%d',cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Integrated Concentration of Tritium \n [x$10^{7}Bq.m^{-3}$]',fontsize=15)
        cb.ax.tick_params(labelsize=13)

    # plot number subplot
    plt.text(0.9, 0.1, item, ha='center', va='center', transform=m.ax.transAxes, fontsize=10)

    return pc

def plot_currData(fig,contour_levels,ax,item,u,v,spd,curint,location='bigmap.p',save=False,cbPlot=False,axesPos=[0.91, 0.365, 0.01, 0.29]):
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
        cbaxes = fig.add_axes(axesPos)
        cb = plt.colorbar(cs1, orientation='vertical',ticks=ticks_cb,format='%1.1f', cax=cbaxes)
        # plotar a legenda da colorbar
        cb.set_label('Surface Current [m.s$^{-1}$]',fontsize=15,labelpad=-1)
        cb.ax.tick_params(labelsize=10)

    maxcurr = 'Max curr: \n %.2fm.s$^{-1}$'%(np.nanmax(spd))
    plt.text(0.9,0.90,maxcurr,ha='center',va='center',transform=m.ax.transAxes,fontsize=10)

    q=m.quiver(lon2[skipcurr],lat2[skipcurr],u[skipcurr],v[skipcurr],alpha=.7,scale=30,minshaft=2)

    # plot number subplot
    plt.text(0.9, 0.1, item, ha='center', va='center', transform=m.ax.transAxes, fontsize=10)

def create_figure4(savefig=False):

    #########################################################
    #                      FIGURA 4                         #
    #    comparação da corrente média nos experimento IV    #
    #                V [corrente e concentracao]            #
    #########################################################
    print("#########################################################")
    print("#                      FIGURA 4                         #")
    print("#    comparação da corrente média nos experimento IV    #")
    print("#                V [corrente e concentracao]            #")
    print("#########################################################")

    # 1st: loading current data 
    # load data from experiment IV
    fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w50/vento_dobrado.cdf'
    f_xr = xr.open_dataset(fname)
    uSW, vSW, spdSW = load_data4(f_xr)
    # load data from experiment V
    fname = '/home/tparente/danilo/TG/simulacoes/rel_final/w210/vento_dobrado.cdf'
    f_xr = xr.open_dataset(fname)
    uNE, vNE, spdNE = load_data4(f_xr)

    debug = False # set as True if debugging this code
    if debug:
        del f_xr, fname

    # creating contour_levels for current, based on expV (with most intense currents)
    maxC = np.nanmax(spdSW)
    contour_levels_current=np.arange(0.,maxC,0.01)

    # 2nd: loading concentration data (already in pickle files)
    i = [159] # set timestep for 40 days of simulation

    # extract data
    swConc, neConc = totalConcentration(i, pkl="/home/tparente/Dropbox/mestrado/artigoTG/data/totalConc.pkl")

    # linhas de contorno para usar no contourf
    maxC = np.nanmax(neConc[:,:])+55
    contour_levels_concentration=np.arange(0.,round(maxC),round(maxC)/1e+3)

    # 3rd: creating axis to plot
    fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(20,11))

    # 4th: plotting current data
    print("Plotting experiment IV current data ... ")
    # plotar experimento I - w210
    plot_currData(fig,contour_levels_current,axes[0,0],item='(a)',u=uNE[:,:],v=vNE[:,:],spd=spdNE[:,:],curint=15,save=False,cbPlot=False)
    #axes[0,0].title.set_text(u'Experiment IV')

    print("Plotting experiment V current data ... ")
    # plotar experimento II - w50
    plot_currData(fig,contour_levels_current,axes[0,1],item='(b)',u=uSW[:,:],v=vSW[:,:],spd=spdSW[:,:],curint=15,save=False,cbPlot=True,axesPos=[0.89, 0.54, 0.01, 0.36])
    #axes[0,1].title.set_text(u'Experiment V')

    # 5th: plotting concentration data
    print("Plotting experiment IV concentration data ... ")
    pc = plot_concData(fig,axes[1,0],contour_levels_concentration,swConc,item='(c)',horas=False,cbPlot=False)
    #ax[1].set_title(u'Experiment VI - Southwesterly Wind [5$ m s^{-1}$]')

    print("Plotting experiment V concentration data ... ")
    pc = plot_concData(fig,axes[1,1],contour_levels_concentration,neConc,item='(d)',horas=False,cbPlot=True,axesPos=[0.89, 0.10, 0.01, 0.36])
    #ax[0].set_title(u'Experiment V - Northeasterly Wind [5$ m s^{-1}$]')

    title = 'Mean Surface Current (above) and Integrated Concentration of Tritium After \n 40 Days (below) in Experiments IV (left) and V (right)'
    plt.suptitle(title,fontsize=25)

    if debug:
        # remove variables
        del uNE, vNE, spdNE
        # remove variables
        del uSW, vSW, spdSW

    if savefig:
        # name's fig
        output1 = '/home/tparente/danilo/mestrado/github/congressos/iwmo2018/figs/fig4.png'
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
def plot_data5(fig,ax,contour_levels,k1,tmpDecor,cbPlot=False,location='bigmap.p'):
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
    tempo = '%s' % (tmpDecor)
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

def importConc():

    try:
        data = pickle.load(open('/home/tparente/danilo/TG/simulacoes/rel_final/wcptecJULAGO/poster_iwmo2018/conc_surface.pkl','r'))

        return data['c']
    except:
        print('Can\'t import pickle file with concentration data')

def create_figure5(savefig=False):
    os.system('clear')

    print("#########################################################")
    print("#                       FIGURA 5                        #")
    print("#         plotar evolucao da pluma de trítio no         #")
    print("#                   experimento VI                      #")
    print("#########################################################")

    fname     = '/home/tparente/danilo/TG/simulacoes/rel_final/wcptecJULAGO/gcmplt.cdf' # definir arquivo gcmplt.cdf a ser lido
    data = xr.open_dataset(fname)               # import netcdf file
    #conc = data['conc'].data[120:,0,:,:]       # import only concentration data
                                                # cut period without pollution
    # conc = conc.sum(axis=1)                   # sum each levels data
    #conc = convertData2Bq(conc)/.1e7           # convert to Bq

    conc = importConc()/.1e7
    
    lon, lat, time = extract_data(data)         # extract lon,lat,time info
    time = time[120:]                           # cut time without leaking

    instantes = [12, 120, 240, 360]              # referentes a 6 horas e 10,
                                                # 20 e 30 dias
    instTitle = ['Leakage Day', '10 Days After', '20 Days After', '30 Days After']# titulo do instante plotado

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

    fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(20,11))

    k = conc[instantes[0]+3,:,:]
    pc1 = plot_data5(fig,ax[0,0],contour_levels,k,tmpDecor=instTitle[0],cbPlot=False)

    k = conc[instantes[1]+3,:,:]
    pc2 = plot_data5(fig,ax[0,1],contour_levels,k,tmpDecor=instTitle[1],cbPlot=False)

    k = conc[instantes[2]+3,:,:]
    pc3 = plot_data5(fig,ax[1,0],contour_levels,k,tmpDecor=instTitle[2],cbPlot=False)

    k = conc[instantes[3]+3,:,:]
    pc4 = plot_data5(fig,ax[1,1],contour_levels,k,tmpDecor=instTitle[3],cbPlot=False)


    ticks_cb = np.asarray([0, 10, 20, 30, 40, 50, 60, 70])
    cb = plt.colorbar(pc4, ax=ax.ravel().tolist(), ticks=ticks_cb, extend='max') # inserir colorbar ocupado o eixo y completo
    cb.set_label('Surface Concentration of Tritium (x$10^{7}Bq.m^{-3}$)',fontsize=20)

    title = u'Final Surface Concentration of Tritium After 30 Days of a Nuclear Leakage'

    plt.suptitle(title, fontsize=25, y=0.94)

    if savefig:
        output = '/home/tparente/danilo/mestrado/github/congressos/iwmo2018/figs/fig5.png'
        plt.savefig(output, dpi=150)
        # let's clean this figure
        os.system('convert -trim %s %s' % (output, output))
    else:
        plt.show()

os.system('clear')

choose = input("Selecione qual figura voce deseja plotar: [1],[3],[4] or [5]? ")
savefig = input("Deseja salvar ou apenas visualizar a imagem: [1] - salvar, [2] visualizar ? ")

if savefig == 1:
    savefig = True
else:
    savefig = False

if choose == 1:
    create_figure1(savefig=savefig)
elif choose == 3:
    create_figure3(savefig=savefig)
elif choose == 4:
    create_figure4(savefig=savefig)
elif choose == 5:
    create_figure5(savefig=savefig)