#!/usr/bin/env python2
#-*-coding:utf-8-*-

# arquivo contendo funcoes utilizadas no artigo
import numpy as np
import xray as xr
import pandas as pd
from scipy.spatial import cKDTree
from scipy import signal, fftpack
import scipy
import socket
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# gerar diretorio base
def make_dir():
    '''
        Funcao para gerar o BASE_DIR baseado no computador em que esta rodando
    '''

    hostname = socket.gethostname()

    if hostname == 'oceano': # estou no meu pc pessoal
        BASE_DIR = '/media/danilo/Danilo/mestrado/github/'

        return BASE_DIR
    if hostname == 'tripoli':
        BASE_DIR = '/home/tparente/danilo/mestrado/github/'

        return BASE_DIR

# funcao de Skill: validacao modelo
def skill_willmott(re,m):
    """
    Analise de Skill (Willmott, 1981) em Python
    Esta funcao esta definida como no trabalho de Miranda et al 2012 (BRAZILIAN JOURNAL OF OCEANOGRAPHY, 60(1):11-23, 201)
    CIRCULATION AND SALT INTRUSION IN THE PIACAGUERA CHANNEL, SANTOS (SP)
    Based on the MSE , a quantitative model skill was presented by Willmott (1981)
    The highest value, WS = 1, means perfect agreement between model and observation, while the lowest value,  WS = 0,
    indicates   complete     disagreement. Recently, this was used to evaluate ROMS in the simulation of multiple parameters in the
    Hudson River estuary [ Warner et al., 2005b] and on the southeast New England Shelf [Wilkin ,2006].
    The Willmott skill will be used to  quantify model performance in simulating different parameters from the best model run
    skill parameter (WILLMOTT, 1981)
    Parameters:
    re - real data
    m - model data
    skill - Skill parameter
    funcao traduzida por: Paula Birocchi
    """
    dif   = re - m
    soma  = np.nansum(abs(dif)**2)
    somam = m - np.nanmean(re)
    c     = re - np.nanmean(re)
    d     = np.nansum((abs(somam) + abs(c))**2)
    skill = 1 - (soma/d)
    return skill

# encontrar indices dos pontos mais proximo a uma coordenada
def find_nearest(lon,lat,ilon,ilat):
    '''
        lon,lat = lat e lon da grade
        ilon,ilat = ponto a ser encontrado
    '''

    lo = lon.ravel()
    la = lat.ravel()

    coords = []

    for i,j in zip(la,lo):
        coords.append([i,j])

    coords = np.array(coords)

    locations_name = ['Terminal Ilha Guaiba']
    locations_posi = [[ilat,ilon]]

    locs = np.asarray(locations_posi)

    tree = cKDTree(coords)
    # procura em tree os pontos mais próximos dos pontos definidos acima
    dists,indexes = tree.query(locs,k=1)

    pontos = []

    for index in indexes:
        pontos.append(coords[index])

    # converter de lista para array
    pontos = np.asarray(pontos)

    # findind indexes from lat and lon
    ind = []

    for p in pontos:
        ind.append(np.where(lon == p[1]))

    ind = np.asarray(ind)

    # vetores para separar i e j para facilitar resgatar os dados de concentração
    iss=[]
    jss=[]

    for i,j in ind:
        iss.append(int(i))
        jss.append(int(j))

    return iss,jss

# baixar dados do Climate Forecast System Version 2
def downloadCFSv2(year,month):
    ''' '''

    try:
        import wget
    except:
        print('Need to instasll wget package. Try pip install wget in the terminal.')

    HTTP_BASE = 'https://nomads.ncdc.noaa.gov/modeldata/cfsv2_analysis_timeseries/'

    date = year+month

    http = HTTP_BASE+year+'/'+date+'/'

    fname = http + 'wnd10m.gdas.'+date+'.grb2'

    print('Downloading %s\n'%(fname))
    wget.download(fname)

# baixar dados do Climate Forecast System Reanalysis
def downloadCFSR(start='1979',final='2011',MONTHS=np.arange(1,13,1)):
    '''
    '''

    try:
        import wget
    except:
        print('Need to instasll wget package. Try pip install wget in the terminal.')

    HTTP_BASE = 'https://nomads.ncdc.noaa.gov/data/cfsr/'

    # create list of year and months
    YEARS = np.arange(int(start), int(final)+1, 1)

    # create list of dates
    DATES = []
    for year in YEARS:
        for month in MONTHS:
            DATES.append(str(year)+str(month))

    for date in DATES:
        fname = HTTP_BASE+date+'wnd10m.gdas.'+date+'.grb2'
        print('Downloading %s\n' % (fname))

        wget.download(fname)

# Recortar uma area de arquivos grib2, converter em netCDF usando CDO
def cut_grb2File(DATA_DIR, SAVE_DIR, PROCESSED_DIR, GRIB_DIR, box=[-55., -35., -15., -35.]):

    '''
        funcao para recortar os arquivos grb2 para um box especifico.

        converte para netcdf

        dependencia: software CDO instalado

        args:
            DATA_DIR: diretorio com arquivos .grb2
            SAVE_DIR: diretorio para armazenar os dados finais
            PROCESSED_DIR: diretorio para armzenar arquiqvos processados
            GRIB_DIR: diretorio para separar os gribs recortados dos netcdf
            box: lista com coordenadas [llon ulon ulat llat]

        Estrutura de diretorios para organizar os dados:

            HOME_DIR/
                |
                ------/*.grb2
                |
                ------/PROCESSED_DIR/*.grb2
                |
                ------/SAVE_DIR/
                |           |
                |           -------/*.nc
                            -------/GRIB_DIR/*cuted_.grb2



    '''

    # extraindo coordenadas do box
    llon,llat = box[0], box[3]
    ulon,ulat = box[1], box[2]

    nfiles = glob.glob(DATA_DIR+'*.grib2') # ler arquivos grb2 no diretorio passado
    nfiles.sort()                         # ordenar alfabeticamente

    os.system('clear')

    # processar os arquivos para recortar
    for f in nfiles:
        print("Cutting: %s"%(f.split('/')[-1]))
        # criar nome do arquivo de saida com diretorio, baseado no proprio
        # arquivo de entrada
        outputFile = SAVE_DIR+f.split("/")[-1]
        # rodar CDO via terminal linux
        os.system('cdo -sellonlatbox,%s,%s,%s,%s %s %s' % (llon,ulon,ulat,llat,f,outputFile))
        # mover o arquivo processado para o diretório correspondente
        os.system('mv %s %s'%(f,PROCESSED_DIR))

    # listar arquivos recortados para converter
    nfiles = glob.glob(SAVE_DIR+'*.grib2')
    nfiles.sort()

    os.system('clear')
    os.system('Converting grib2 files to netCDF files')

    for f in nfiles:
        os.system('cdo -f nc copy %s %s' % (f,f.replace('.grib2', '.nc')))
        os.system('mv %s %s' % (f, GRIB_DIR))

# Quando baixar arquivos já recortados (Bob Dattore), usar a seguinte funcao
def read_month(date,DATA_DIR):
    '''
        Funcao que le os arquivos .nc baixados do NCAR/UCAR, calcula a média
        diária e mensal e retorna a média mensal pronta para plotar

        date = YEARMONTH (ex: 201411)

    '''
    nfiles = glob.glob(DATA_DIR + 'cdas1.' + date + '*')
    nfiles.sort()

    # checar o shape dos dados para utilizar na criação da matriz de u,v
    ncdata = xr.open_dataset(nfiles[0])
    wu     = ncdata['U_GRD_L103'].values[0,:,:]
    shape  = np.array([ncdata.dims['lon'], ncdata.dims['lat']])

    matriz_u, matriz_v = np.zeros([len(nfiles),shape[0],shape[1]]), np.zeros([len(nfiles),shape[0],shape[1]])

    cont = 0        # contador para o dia

    for f in nfiles: # loop pelos arquivos do mes escolhido
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

# interpolar grade mais refinada para uma grade menos refinada
def interpolar_grade(coarsedGrid,refinedGrid,variable):
    '''
        Funcao para interpolar da grade do CFSv2 (mais refinada) para a grade
        do CSFR (menos refinada).

        input:
            coarseGrid = DICIONARIO com lons e lats da grade menos refinada
            refineGrid = DICIONARIO com lons e lats da grade mais  refinada
            variable   = variavel a ser interpolada para nova grade

        output:
            interpolated_variable

        Importante: a coarseGrid deve ser maior que a refineGrid, ou seja,
        a grade menos refinada deve englobar totalmente a grade mais refinada.
        Assim a interpolação não irá gerar valores NaN.
    '''

    try:
        from scipy.interpolate import griddata
    except:
        print('Need to install scipy package to import griddata')

    # extrair as grades: coarsed (menos refinada) e refined (mais refinada)
    lon_coarsed = coarsedGrid['lon']
    lat_coarsed = coarsedGrid['lat']
    lon_refined = refinedGrid['lon']
    lat_refined = refinedGrid['lat']

    # devemos remover a primeira e ultima linha da grade refinada, pois elas estão
    # para fora da grade grosseira e isso pode gerar problemas na interpolação
    lon_refined = lon_refined[1:-1]
    lat_refined = lat_refined[1:-1]
    variable    = variable[1:-1,1:-1]

    # preparar os dados para interpolacao
    lon_r,lat_r = np.meshgrid(lon_refined, lat_refined)
    points = np.array([lon_r, lat_r])               # shape da variavel
    lon_c,lat_c = np.meshgrid(lon_coarsed, lat_coarsed)
    xi     = np.array([lon_c, lat_c])               # shape de interpolacao

    # interpolar
    interpolated_variable = griddata((lon_r.flatten(),lat_r.flatten()),variable.flatten(), (lon_c,lat_c), method='cubic')

    return interpolated_variable

# importar as grades refinadas e grosseiras
def extrair_grades(coarsedFile, refinedFile):
    '''
    args:
        coarsedFile = arquivos .nc contendo dados com lat e lon
        refinedFile = arquivos .nc contendo dados com lat e lon

    output:
        coarseGrid, refinedGrid = dicionários com lons e lats
    '''

    try:
        import xarray as xr
    except:
        print('Need to install xarray to use this function.')

    coarsed     = xr.open_dataset(coarsedFile)
    lon_coarsed = coarsed['lon'].values - 360
    lat_coarsed = coarsed['lat'].values

    refined = xr.open_dataset(refinedFile)
    lon_refined = refined['lon'].values - 360
    lat_refined = refined['lat'].values

    coarsedGrid = {'lon':lon_coarsed,'lat':lat_coarsed}
    refinedGrid = {'lon':lon_refined,'lat':lat_refined}

    return coarsedGrid, refinedGrid

# remover arquiqvos baixados do NCEP/NCAR
def removeFiles(years,months,prefix,posfix):
    '''
        Remove file downloaded from cfsr/cfsv2.

        Args:
            years (list): list with years to be removed
            months (list): list of months to be removed
            prefix (string): full path + some pattern
            posfix (string): basically the extension

        Returns:
            nothing, because this function will remove files
            from your computed.

            Be carefull ...


        Examples:

            >>> years  = np.arange(1993,2011,1)
            >>> months = np.arange(04,11,1)
            >>> prefix = '/media/danilo/Danilo/mestrado//ventopcse/data/CFSR/atlanticoSW/wnd10m/wnd10m.gdas.'
            >>> posfix = '.grb2.nc'

            >>> removeFiles(years,months,prefix,posfix)

    '''
    dates = []

    for y in years:
        for m in months:
            d = prefix + str(y) + str(m) + posfix

            dates.append(d)

    for f in dates:
        os.system('rm %s' % (f))

def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-40,resolution='l'):

	m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

	m.ax = ax

	m.drawcoastlines(linewidth=0.2)
	m.drawmapboundary(fill_color='#e5f2ff')
	m.fillcontinents(color='gray')

	m.drawcoastlines(linewidth=.8)
	m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
	meridians=np.arange(llon,ulon,3)
	parallels=np.arange(llat,ulat,2)
	# desenhar meridianos e paralelos conforme definido acima
	m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
	m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

	return m

def spdir2uv(spd, ang, deg=False):
    """
    Computes u, v components from speed and direction.
    Parameters
    ----------
    spd : array_like
          speed [m s :sup:`-1`]
    ang : array_like
          direction [deg]
    deg : bool
          option, True if data is in degrees. Default is False
    Returns
    -------
    u : array_like
        zonal wind velocity [m s :sup:`-1`]
    v : array_like
        meridional wind velocity [m s :sup:`-1`]
    """

    if deg:
        ang = np.deg2rad(ang)

    # Calculate U (E-W) and V (N-S) components
    u = spd * np.sin(ang)
    v = spd * np.cos(ang)

    return u, v
