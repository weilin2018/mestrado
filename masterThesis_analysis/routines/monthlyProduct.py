'''

'''

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

#### funcoes

##############################################################################
#                                                                            #
#                          FUNCOES GERAIS                                    #
#                                                                            #
##############################################################################

def read_month(date,DATA_DIR):
    '''
        Funcao que le os arquivos .nc baixados do NCAR/UCAR, calcula a média
        diária e mensal e retorna a média mensal pronta para plotar

        date = YEARMONTH (ex: 201411)

    '''
    nfiles = glob.glob(DATA_DIR + '*' + date + '*')
    nfiles.sort()

    # checar o shape dos dados para utilizar na criação da matriz de u,v
    ncdata = xr.open_dataset(nfiles[0])
    wu     = ncdata['U_GRD_L103'].values[0,:,:]
    shape  = np.array([ncdata.dims['lat'], ncdata.dims['lon']])

    matriz_u, matriz_v = np.zeros([len(nfiles),shape[0],shape[1]]), np.zeros([len(nfiles),shape[0],shape[1]])

    cont = 0        # contador para o dia

    for f in nfiles: # loop pelos arquivos do mes escolhido
        print('Reding file: %s \n' % (f))
        ncdata = xr.open_dataset(f)

        # extrair componentes do vento
        u = ncdata['U_GRD_L103'].values
        v = ncdata['V_GRD_L103'].values

        # tomar a media diaria
        umean = u.mean(axis=0)
        vmean = v.mean(axis=0)

        # armazenar a media diária no dia correspondete da matriz vazia
        matriz_u[cont,:,:] = umean[:,:]
        matriz_v[cont,:,:] = vmean[:,:]

        cont += 1

    # retorna a media mensal
    return matriz_u.mean(axis=0), matriz_v.mean(axis=0)

def getLonLat(DATA_DIR):
    ''' '''
    # extrair longitude e latitude
    nfiles = glob.glob(DATA_DIR+"*.nc")[0]
    ncdata = xr.open_dataset(nfiles)
    lon    = ncdata['lon'].values - 360
    lat    = ncdata['lat'].values

    lon,lat = np.meshgrid(lon,lat)

    return lon,lat

def make_map(ax, llat, ulat, llon, ulon):
    '''
        Plot a map using Basemap.

        Args:
            ax (matplotlib axes): axis to plot data
            llat,ulat (float): lower latitude, upper latitude
            llon,ulon (float): lower longitude, upper longitude
    '''
    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='l')
    # m = pickle.load(open("/media/danilo/Danilo/mestrado/ventopcse/rotinas/sudesteBR.pkl", "r"))
    m.ax = ax
    m.drawcoastlines(linewidth=.8)
    m.drawmapboundary()
    # definir meridianos e paralelos para plotar no mapa
    meridians=np.linspace(llon,ulon,4)
    parallels=np.linspace(llat,ulat,5)
    # desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

    return m

def plotar(saveData, llat, ulat, llon, ulon, dctTitles, DATA_DIR):
    '''
        receber os dados e plotar as
    '''
    # vetor com localizacao do subplots
    locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)]

    # obter grade griddada
    lon,lat = getLonLat(DATA_DIR)

    for mes,loc in zip(['nov', 'dec', 'jan', 'feb', 'mar'], locs):
        print('Plotting climatology for %s' % (mes))
        # dados:
        wu,wv = saveData[mes]['wu'], saveData[mes]['wv']
        # calcular velocidade
        spd = np.sqrt(wu**2 + wv**2)

        # realizar os plots
        ax = plt.subplot2grid(shape=(2,6), loc=loc, colspan=2)

        m = make_map(ax, llat, ulat, llon, ulon)

        contour_levels = np.arange(0,10.001,0.001)

        c = m.contourf(lon,lat,spd,contour_levels,latlon=True,extend='max')
        q = m.quiver(lon[::3,::3],lat[::3,::3],wu[::3,::3], wv[::3,::3], latlon=True,
                                        alpha=.7,scale=150,width=0.005,minshaft=2)
        m.ax.set_title(dctTitles[mes])

        cb = plt.colorbar(c,orientation='horizontal',ticks=[0,2,4,6,8,10],format='%d',fraction=.057,pad=.06)
        cb.set_label(r'Wind [$m.s^{-1}$]',fontsize=8, labelpad=-1)
        cb.ax.tick_params(labelsize=8)

    plt.show()

def regionParameters(region, BASE_DIR, DATABASE):
    '''

    Args:
        region (string): the study area. Could be PCSE or ATSW
        BASE_DIR (string): directory
        DATABASE (string): which database you're using. Could be CFSR or CFSv2

    Returns
        DATA_DIR    (string): directory where all data are stored
        llat,ulat   (float,float): lower latitude and upper latitude
        llon,ulon   (float,float): lower longitude and upper longitude

    '''
    if DATABASE == 'CFSR':
        if region == 'PCSE':
            DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSR/1992_2011/')
            # coordenadas para plot
            llat, ulat = -30, -20
            llon, ulon = -50, -40

        if region == 'ATSW':
            DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/atlanticoSW/wnd10m/')
            # coordenadas para plot
            llat, ulat = -45, -15
            llon, ulon = -60, -30

    if DATABASE == 'CFSv2':
        if region == 'PCSE':
            DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao')
            # coordenadas para plot
            llat, ulat = -30, -20
            llon, ulon = -50, -40

        if region == 'ATSW':
            DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/atlanticoSW/verao')
            # coordenadas para plot
            llat, ulat = -45, -15
            llon, ulon = -60, -30


    return DATA_DIR, llat, ulat, llon, ulon

##############################################################################
#                                                                            #
#                          FUNCOES CLIMATOLOGIA                              #
#                                                                            #
##############################################################################

def ler_climato(DATA_DIR):
    '''
        ler arquivos e calcular as médias para climatologia

        Args
            DATA_DIR (string): diretorio contendo os dados

        Returns
            saveData (dictionary): dicionário contendo a média mensal, onde
                cada mês, sendo um chave, é relacionado a uma matriz de dados médios.
    '''

    meses = ['nov', 'dec','jan', 'feb', 'mar']
    locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)] # vetor com localizacao do subplots
    dates = '11 12 01 02 03'.split(" ")
    # dicionario para armazenar os dados
    saveData = {}

    for mes,loc,date in zip(meses,locs,dates):
        print('Processing %s' % (mes))
        print('-----------------------------------------------')
        # ler arquivos referentes ao mes
        wu,wv = read_month(date,DATA_DIR)

        # calcular velocidade
        spd = np.sqrt(wu**2 + wv**2)

        data = { 'wu': wu, 'wv': wv }  	# organizando dados para armazenamento
        saveData[mes] = data   			# dicionario para salvar os dados em pickle

    return saveData

def climatologia(BASE_DIR=oceano.make_dir(),read=False,plot=True,region='PCSE'):
    '''
        Ler e tratar os dados baixados para gerar a climatologia

        Args
            BASE_DIR (string) = diretorio base
            read     (boolean) = True (ler arquivos e calcular as médias),
                                 False (ler arquivo .pickle)

            plot     (boolean) = True (plotar)
            region   (string)  = região a ser calculado:
                    'PCSE' [Plataforma Continental Sudeste] default ou
                    'ATSW' [Atlantico Sudoeste]

        Returns
            saveData (dictionary): médias mensais calculadas
    '''

    # definindo o diretorio contendo os dados da regiao de interesse
    DATA_DIR, llat, ulat, llon, ulon = regionParameters(region,BASE_DIR,'CFSR')

    # calculando ou nao a climatologia
    if read:
        saveData = ler_climato(DATA_DIR)
        # save data to a pickle file
        pickle.dump(saveData, open(BASE_DIR.replace('github/', '/ventopcse/data/pickles/climatology.pickle'), 'w'))
    else:
        saveData = pickle.load(open(BASE_DIR.replace('github/', '/ventopcse/data/pickles/climatology.pickle'), 'r'))

    dctTitles = {'nov': 'Nov/1992 - Nov/2010','dec': 'Dez/1992 - Dez/2010','jan': 'Jan/1993 - jan/2010','feb': 'Fev/1993 - Fev/2010','mar': 'Mar/1993 - Mar/2010'}

    if plot:
        fig = plt.figure(figsize=(16,8))

        plotar(saveData,llat,ulat,llon,ulon,dctTitles,DATA_DIR)

    return saveData

##############################################################################
#                                                                            #
#                          FUNCOES VERAO                                     #
#                                                                            #
##############################################################################

def ler_summer(DATA_DIR):
    '''

    '''
    meses = ['nov', 'dec','jan', 'feb', 'mar']
    locs  = [(0,0), (0,2), (0,4), (1,1), (1,3)] # vetor com localizacao do subplots
    dates = '11 12 01 02 03'.split(" ")
    # dicionario para armazenar os dados
    saveData = {}

    for mes,loc,date in zip(meses,locs,dates):

        # ler os dados do mes em loop
        wu,wv = read_month(date,DATA_DIR)
    	# calcular velocidade
        spd = np.sqrt(wu**2 + wv**2)

        data = { 'wu': wu, 'wv': wv }  	# organizando dados para armazenamento
        saveData[mes] = data   			# dicionario para salvar os dados em pickle

    return saveData

def summer(BASE_DIR=oceano.make_dir(),read=False,plot=True,region='PCSE',summer='2014'):
    '''
        Ler e tratar os dados baixados para gerar o campo médio dos verões

        Args
            BASE_DIR (string) = diretorio base
            read     (boolean) = True (ler arquivos e calcular as médias),
                                 False (ler arquivo .pickle)

            plot     (boolean) = True (plotar)
            region   (string)  = região a ser calculado:
                    'PCSE' [Plataforma Continental Sudeste] default ou
                    'ATSW' [Atlantico Sudoeste]
            summer   (string): which summer you want to plot [2014 or 2015]

        Returns
            saveData (dictionary): médias mensais calculadas
    '''

    # definindo o diretorio contendo os dados da regiao de interesse
    DATA_DIR, llat, ulat, llon, ulon = regionParameters(region,BASE_DIR,'CFSv2')

    # como temos dois cenarios de verao, precisamos completar o DATA_DIR com o verao em escolha
    DATA_DIR = DATA_DIR + str(summer) + '/'

    # calculando ou nao a climatologia
    if read:
        saveData = ler_summer(DATA_DIR)
        # save data in a pickle file
        pickleBase = BASE_DIR.replace('github/', '/ventopcse/data/pickles/')
        pickleName = pickleBase + 'summer'+str(summer)+'_'+str(region)+'.pickle'

        pickle.dump(saveData,open(pickleName, 'w'))
    else:
        pickleBase = BASE_DIR.replace('github/', '/ventopcse/data/pickles/')
        pickleName = pickleBase + 'summer'+str(summer)+'_'+str(region)+'.pickle'

        saveData = pickle.load(open(pickleName, 'r'))

    os.system('clear')
    print('#################################################################')
    print('#                        PLOTTING                               #')
    print('#################################################################')

    year = int(summer)

    dctTitles = {'nov': 'Nov/%s'%(str(year-1)),'dec': 'Dez/%s'%(str(year-1)),'jan': 'Jan/%s'%(summer),'feb': 'Fev/%s'%(summer),'mar': 'Mar/%s'%(summer)}

    if plot:
        fig = plt.figure(figsize=(16,8))

        plotar(saveData,llat,ulat,llon,ulon,dctTitles,DATA_DIR)

    return saveData

##############################################################################
#                                                                            #
#                          FUNCOES ANOMALY                                   #
#                                                                            #
##############################################################################


def anomaly(climatology,month):
    '''
        climatology: campo de vento médio tomado como estado básico
        monmth: campo de vento médio para o mês em questão
    '''
    # return climatology - month
    return month - climatology

def Anomaly(climatology,data,calculate=False,plot=True,summer='2014',region='PCSE',BASE_DIR=oceano.make_dir()):
    '''
        Function to calculate an anomaly for a specific periodo (data),
        using the climatology, taken from climatologia() function.

        Args
            climatology (dict): dictionary with monthly climatology, where
                                each key is a month
            data        (dict): dictionary with monthly mean, where
                                each key is a month
            summer    (string): which summer we're analysing [2014 or 2015]
            region    (string): which region we're working

        Returns
            anomaly     (dict): dictionary with monthly anomaly, where
                                each key is a month

        Example
            >>> climatology = climatologia(read=True,plot=False)
            >>> summer2014  = summer(read=True,plot=False,summer='2014')
            >>> anomaly2014 = anomaly(climatology, summer2014, plot=False)
    '''

    # importar as grades que serão utilizadas na interpolação ou plotagem
    coarsedFile = glob.glob(BASE_DIR.replace('github/', 'ventopcse/data/CFSR/1992_2011/*.nc'))[0]
    refinedFile = glob.glob(BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao2014/*.nc'))[0]

    coarsedGrid, refinedGrid = oceano.extrair_grades(coarsedFile, refinedFile)

    # definindo o diretorio contendo os dados da regiao de interesse (grade maior)
    DATA_DIR, llat, ulat, llon, ulon = regionParameters(region,BASE_DIR,'CFSR')

    if calculate:
        dataAnomaly = calcAnomaly(climatology,data,coarsedGrid,refinedGrid)

        # save data in a pickle file
        pickleBase = BASE_DIR.replace('github/', '/ventopcse/data/pickles/')
        pickleName = pickleBase + 'anomlay_summer'+str(summer)+'_'+str(region)+'.pickle'

        pickle.dump(dataAnomaly,open(pickleName, 'w'))
    else:
        pickleBase = BASE_DIR.replace('github/', '/ventopcse/data/pickles/')
        pickleName = pickleBase + 'anomlay_summer'+str(summer)+'_'+str(region)+'.pickle'

        saveData = pickle.load(open(pickleName, 'r'))

    if plot:
        fig = plt.figure(figsize=(16,8))

        year = int(summer)

        dctTitles = {
            'nov': 'November/%s - Climatology'%(str(year-1)),
            'dec': 'December/%s - Climatology'%(str(year-1)),
            'jan': 'January/%s - Climatology'%(summer),
            'feb': 'February/%s - Climatology'%(summer),
            'mar': 'March/%s - Climatology'%(summer)
        }

        plotar(saveData,llat,ulat,llon,ulon,dctTitles,DATA_DIR)

    return saveData

def calcAnomaly(climatology, summer,coarsedGrid,refinedGrid):
    ''' calcular '''

    meses = ['nov', 'dec','jan', 'feb', 'mar']

    # variable to store the anomaly data
    saveData = {}

    for mes in meses:
        # importar os dados do mes em analise
        mes_climato = climatology[mes]
        mes_summer  = summer[mes]

        # interpolar os dados para a mesma grade
        wui = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wu'])
        wvi = oceano.interpolar_grade(coarsedGrid,refinedGrid,mes_summer['wv'])

        # calculando as anomalias
        wu_anomaly = anomaly(mes_climato['wu'], wui)
        wv_anomaly = anomaly(mes_climato['wv'], wvi)

        # store data
        saveData[mes] = {'wu': wu_anomaly, 'wv': wv_anomaly}

        # calcular a velocidade anomala
        spd_anomaly = np.sqrt(wu_anomaly**2 + wv_anomaly**2)

    return saveData


##############################################################################
#                                                                            #
#                          FUNCAO MAIN                                       #
#                                                                            #
##############################################################################
