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
import pickle

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
t=
def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-40,resolution='l'):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

    m.drawcoastlines(linewidth=0.1)
    m.drawmapboundary(fill_color='#e5f2ff')
    m.fillcontinents(color='#c0c0c0')

    m.drawcoastlines(linewidth=.1)
    m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
    meridians=np.arange(llon,ulon,3)
    parallels=np.arange(llat,ulat,2)
	# desenhar meridianos e paralelos conforme definido acima
    m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
    m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

    return m

def spdir2uv(spd, ang, deg=False, convention=None):
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

def rotaciona(u,v,angulo):
    """Rotacionar as componentes de velocidade (u e v) para paralelo e perpendicular.

    Routine based on the book Data Analysis in Physical Oceanography,
    page 425.

    Credits
    -------
    Created by Paula Birocchi (paula.birocchi@gmail.com)

    Parameters
    ----------
    u : array
        East-West component.
    v : array
        North-South component.
    angulo : float
        Angle (in radians) to rotate the components, related to the coast. Could
        be some constant, assuming a rectilinar coast,  or an array with the
        same size of u and v, in case of a curvilinear coast.

    Returns
    -------
    urotated : array
        Cross-shore velocity (perpendicular component to the coast).
    vrotated : array
        Along-shore velocity (parallel component to the coast)

    """
    urotated = u * np.cos(angulo) - v * np.sin(angulo)
    vrotated = u * np.sin(angulo) + v*np.cos(angulo)

    return urotated, vrotated

def rotateVectors(u,v,angle):
    """Rotate velocities components to along and cross shore, based
    on a given angle of rotation.

    This is based on the file available in:
     <link>

    where we develop the mathematical steps.

    Parameters
    ----------
    u : numpy.ndarray
        Eastward velocity component [m/s].
    v : numpy.ndarray
        Northward velocity component [m/s].
    angle : float
        Angle of rotation in degrees.

    Returns
    -------
    cross,along : numpy.ndarray
        Cross and Along shore components of velocity.

    """
    # convert from degrees to radians
    angle = np.deg2rad(angle)

    # compute cross and along shore components
    cross = u*np.cos(angle) - v*np.sin(angle)
    along = u*np.sin(angle) + v*np.cos(angle)

    return cross,along

def compass2uv(direction, speed, kind):
    """Converts speed and direction (of wind) from meteorological convention
    to oceanographic convention (u,v).

    This function was translated from compass2uv writed to MATLAB,
    founded in:
    https://marine.rutgers.edu/~codaradm/metstation/scripts/clean_met_data/compass2uv

    Converting directions to mathematical convention (increasing CCW from x-axis):

    We change angles to progress CCW from north (360-dir), then
    the angles are rotated by +90degrees ((360-dir)+90).

    After that, this function use pol2cart's function, to convert from polar
    coordinates to cartesian coordinates.

    Finally, rounding erros are removed from u and v components.

    Credits
    -------
    Translated by Danilo Augusto Silva <nilodna@gmail.com>

    Parameters
    ----------
    direction : array-like
        compass direction, in degrees (0 to 360).
    speed : array-like
        Vector length, in any units.

    Returns
    -------
    u   : array-like
        east vector component, in the same units as speed.
    v   : array-like
        north vector component, in the same units as speed
    """

    if kind=='meteo2uv':
        direction = 90 - (direction - 180)
        u = (-1) * speed * np.cos(direction)
        v = (-1) * speed * np.sin(direction)

    if kind=='ocean2uv':
        direction = 90 - direction
        u = speed * np.cos(direction)
        v = speed * np.sin(direction)


    # direction = (360 - np.asarray(direction)) + 90
    # direction = angle360(direction) # ensure angles range from 0 to 360
    #
    # # converting from polar to cartesian coordinates system
    # d = (direction*np.pi)/180
    # u,v = pol2cart(d, speed)

    # # removing any rounding error
    # roundoff = 1e-14
    # u = round(u/roundoff)*roundoff
    # v = round(v/roundoff)*roundoff

    return u,v

def angle360(ang,shift=0):
    """This function ensures angles (in degrees) range from 0 to 360. An angle shift can be added.

    This function was originally created to MatLab by Mike Whitney (link in the final
    of this documentation) and translated to python by Danilo Augusto Silva
    (nilodna@gmail.com).

    If you're using direction's wind and converting then to ocenographic convention,
    you need to ensure that the directions range between 0 to 360, so the
    function compass2uv can be used properly.

    Source code: https://marine.rutgers.edu/~codaradm/metstation/scripts/clean_met_data/angle360.m

    Parameters
    ----------
    ang : array-like
        Array containing all angles in degrees.
    shift : float
        A shift value to add to the angles.The default is zero.

    Returns
    -------
    ang
        Angles, in degrees, in the interval of 0 to 360.

    """
    # add an angle shift
    ang += shift

    ang[ang < 0] = ang[ang < 0]+360
    ang[ang >= 360] = ang[ang >= 360] - 360

    return ang

def pol2cart(rho, phi):
    """Convert from polar coordinates to cartesian coordinates.

    Parameters
    ----------
    rho : array-like
        Direction in polar coordinantes.
    phi : array-like
        Intensity/speed of vectors.

    Returns
    -------
    x,y : array-like
        Components in cartesian coordinates.

    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x,y

def dirmag2uv(direction,speed,decl_mag,ref_direcao):
    """ Transforma direcao (referenciada no LESTE - circulo trigonometrico ou Referenciada no NORTE) e magnitude/speed/velocidade em componentes U
    (leste-oeste) e V (norte-sul) - CHECK! Tudo correto aqui.
    ATENCAO!!! Se a sua direcao esta referenciada em NORTE, use: u = seno(dir)*speed e v=cos(dir)*speed
    Se a sua direcao esta referenciada em LESTE (circulo trignometrico) use: u = cosseno(dir)*speed e v=seno(dir)*speed.
    ref_direcao = 'trignometrico' ou 'norte'. 'trignometrico' se o dir=0 eh no LESTE. 'norte' se dir=0 eh no NORTE.
    """
    import math
    direction = direction + decl_mag
    # print("direction = " + str(direction))
    # direction = np.mod(direction, 360)
    # print("direction = " + str(direction))
    direction = direction * np.pi / 180
    if ref_direcao=='trigonometrico':
        u = np.cos(direction)*speed
        v = np.sin(direction)*speed
    if ref_direcao=='norte':
        u = (-1)*np.sin(direction)*np.abs(speed) # esta assim no BoB.
        v = (-1)*np.cos(direction)*np.abs(speed)
    return u,v

# baixar dados do Climate Forecast System Reanalysis
# def downloadCFSR(start='1979',final='2011',MONTHS=np.arange(1,13,1)):

# convert velocities components to direction and intensity
def uv2intdir(u,v):
    """Converting velocities components to intensity and direction.

    Parameters
    ----------
    u : numpy.ndarray
        East vector component [m.s$^{1}$].
    v : numpy.ndarray
        North vector component [m.s$^{1}$].

    Returns
    -------
    speed : array-like
        Vector length, in any units.
    direction : array-like
        Compass direction, in degrees (0 to 360).
    """

    from math import atan2

    # create list to store the data calculated
    ws,wd = [],[]

    for i,j in zip(u,v):
        intensity = np.sqrt(i**2 + j**2)
        direction = 270 - (180/np.pi)*atan2(j,i)

        if direction > 360:
            direction -= 360

        ws.append(intensity)
        wd.append(direction)


    ws = np.asarray(ws)
    wd = np.asarray(wd)

    return ws,wd

def uv2dirmag(u,v):
    """ Transforma componente U (leste-oeste) e V (norte-sul) em
    magnitude e direcao em relacao AO NORTE!!!.
    CUIDADO: Para o vento, a direcao que esta funcao retorna é pra onde ELE SOPRA e nao da direcao da onde ele vem !!!
    Exemplo: U = +1 e V = +1 -> gera direcao = 45 graus!!!
    Paula Birocchi
    magnitude - corrente resultante
    direcaorad - direcao em radianos
    direcaol - direcao em graus em relacao a direcao oeste-leste (0 graus em leste)
    direcaon - direcao em graus em relacao a direcao norte (0 graus em norte)! --> DIRECAO FINAL!!!!
    """
    magnitude     = np.sqrt(u**2+v**2)
    direcaorad = np.arctan2(v,u)
    direcaon = direcaorad/np.pi*180 # direcao em relacao a direcao
    direcaon = np.mod(90 - direcaon,360) #esta linha está na rotina do Bob, e só deve ser usada no caso dos dados do canal de São Sebastião e nos dados do NCEP!.
    return magnitude, direcaon

#### plotting directional histrogram
def polarPlot(ws,wd,title):
    """plot a polar plot with intensity and direction.

    Check the convention you're passing the data,
    because this may influence in your analyze. Remembering that
    meteorological convention (show the direction from which the wind is
    blowing) and oceanographic convention (show the direction towards which wind
    is blowing).

    Credits
    -------
    Function created by Danilo A. Silva <nilodna@gmail.com>,
    from Coastal Hydrodynamics Lab.

    Parameters
    ----------
    ws : numpy.ndarray
        Wind's intensity.
    wd : numpy.ndarray
        Wind's direction.
    """

    from windrose import WindroseAxes

    ax = WindroseAxes.from_ax()
    ax.bar(wd,ws, normed=True,opening=0.8)
    ax.set_legend()

    ax.set_title(title)

    plt.show()

######################################

def uv2compass(wve,wvn):
    """Short summary.

    Parameters
    ----------
    wve : array_like
        eastward component of velocity.
    wvn : array_like
        northwart component of velocity.

    Returns
    -------
    wmag
        magnitude of velocity vector in same units as wve and wvn.
    wdir
        compass direction of vector (0.deg north)

    """

    from math import atan2

    wmag = np.real(np.sqrt(wve**2 + wvn**2))

    ratio = wvn/wve

    whereLESS = np.where(ratio < 0.)
    ratio[whereLESS] *= -1.0

    d = 180/np.pi

    wdir = atan2(wvn,wve) * d

    wdir[np.where((wve >= 0.) & (wvn >= 0.))] =  90. - wdir[np.where((wve >= 0.) & (wvn >= 0.))]
    wdir[np.where((wve >= 0.) & (wvn >= 0.))] = 270. + wdir[np.where((wve <  0.) & (wvn >= 0.))]
    wdir[np.where((wve >= 0.) & (wvn >= 0.))] =  90. + wdir[np.where((wve >= 0.) & (wvn  < 0.))]
    wdir[np.where((wve >= 0.) & (wvn >= 0.))] = 270. - wdir[np.where((wve <  0.) & (wvn  < 0.))]

    return wmag,wdir

def save_data(data,directory):
    """save some variable in a pickle file located in directory pass as argument.

    Parameters
    ----------
    data : variable
        Some variable to store.
    directory : string
        Full path + filename.

    """

    import pickle

    pickle.dump(data,open(directory,'w'))

def load_data(directory):
    """Short summary.

    Parameters
    ----------
    directory : string
        Full path + filename.

    Returns
    -------
    data : variable
        Some saved variable.

    """

    import pickle

    data = pickle.load(open(directory,'r'))

    return data

def stickplot(df,ax):
    """Create a stickplot.

    With the u and v components given as argument in pd.DataFrame df,
    this function plot a stickplot, using MPL quiver.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe containing wu and wv components and datetimeindex.
    ax : matplotlib.axis
        axis to stickplot

    Example
    -------
    >>>
    >>>

    Credits
    -------
    Stephane Raynaud
    http://permalink.gmane.org/gmane.comp.python.matplotlib.general/24155
    """

    # creating the date axis
    # dates = pd.to_datetime(df.index)
    # extracting components from dataframe
    u = df['wu'].values
    v = df['wv'].values

    # calculating speed
    spd = np.sqrt(u**2 + v**2)
    maxSpd = np.nanmax(spd)

    # plotting
    # fig, ax = plt.subplots()

    qiv = ax.quiver(df.index, [[0]*len(df)], u, v, headlength=0, headwidth=0, headaxislength=0 )
    key = ax.quiverkey(qiv, 0.25, 0.75, maxSpd, "%0.2f $m^{2}s^{-1}$"%(maxSpd), labelpos='N', coordinates='axes' )

    # plot a horizontal line in y=0.0
    ax.axhline(y=0.0,xmin=df.index[0],xmax=df.index[-1],linewidth=1.,color='black')

    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis_date()

    # ax.set_xticks(['2012-01-07', '2012-01-21', '2012-02-04', '2012-02-18', '2012-03-03', '2012-03-17', '2012-03-31'])

    return ax

###### functions to store data into pickle
def createPickle(data,namePickle):
    """Short summary.

    Parameters
    ----------
    data : dict
        dicionary with all data to store into pickle.
    namePickle : str
        full path and pickle file name.
    """

    if not type(data) == dict:
        print('The argument \'data\' must be a dicionary.')
    else:
        pickle.dump(data,open(namePickle,'w'))

def removeLabelsFromDataset(ncdata):
    """ Function to remove some labels from xarray.dataset.

    Parameters
    ----------
    ncdata : xarray.Dataset
        Dataset
    Returns
    -------
    nc : xarray.Dataset
        Dataset with labels removed
    """

    labels = ['xpos','ypos','time','date','layer_bnds','x','y','h1','h2','depth','ang','FSM','DUM','DVM','lon_bnds','lat_bnds','cbc']

    return ncdata.drop(labels)
