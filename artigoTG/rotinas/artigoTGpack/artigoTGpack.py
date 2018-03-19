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

    # localizacao do terminal da ilha guaiba

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

# encontrar o lag de maxima correlacao entre duas series
def max_correlacao(modelo, observado):
    '''
    obtem a diferenca temporal (em horas) da maxima correlação entre
    duas séries temporais.
    modelo: série temporal saída do modelo, em pandas.Series
    observado: série temporal obtida in situ, em pandas.Series

    return
        diferença de fase em horas
    '''
    z = signal.fftconvolve(modelo, observado[::-1])
    lags = np.arange(z.size) - (observado.size - 1)

    return ( lags[np.argmax(np.abs(z))] )

# desenhar um quadrado colorido em uma imagem basemap
def draw_square(lats,lons,m):
    x,y=m(lons,lats)
    xy=zip(x,y)
    poly = Polygon(xy,facecolor='blue', alpha=0.4, edgecolor='black')
    plt.gca().add_patch(poly)

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

def read_BNDO(DATA_DIR):
    '''
        ler arquivos BNDO no diretório passado como argumento

        

        retorna os dados de 1997 como pandas.DataFrame
        
    '''
    lfiles = glob.glob(DATA_DIR+'1997/*')   
    lfiles.sort()

    # ler, inicialmente, os dois primeiros arquivos para ampliar uma série
    files = lfiles[:2]

    file1 = pd.read_csv(files[0], skiprows=11, delimiter=';', names=['nivel', 'x'])
    file1.drop(file1.columns[len(file1.columns)-1], axis=1, inplace=True)

    file2 = pd.read_csv(files[1], skiprows=11, delimiter=';', names=['nivel', 'x'])
    file2.drop(file2.columns[len(file2.columns)-1], axis=1, inplace=True)

    # criar os dataframes
    dtRange = pd.date_range(start=file1.index[0], end=file1.index[-1], freq='H')
    df1 = pd.DataFrame({'nivel': file1['nivel'].values/100.}, index=dtRange)

    dtRange = pd.date_range(start="1997-02-02 00:00", end="1997-03-05 23:00", freq='H')
    df2 = pd.DataFrame({'nivel': file2['nivel'].values/100.}, index=dtRange)

    dtRange = pd.date_range(start='1997-02-01 00:00', end='1997-02-01 23:00', freq='H')
    df3 = pd.DataFrame({'nivel': np.ones(dtRange.shape[0])*np.nan}, index=dtRange)

    # concatenar as séries
    observ = pd.concat([df1, df3, df2])

    # controle de qualidade
    cond = observ['nivel'] > 4.
    observ[cond] = np.nan

    # removendo a média da série temporal
    observ['nivel'] = observ['nivel'] - observ['nivel'].mean()

    return observ



def newTimerange(tm,to,observ):

    ''' pegar os dados mais próximos dos instantes passados do modelo '''

    modelTimestamp = []

    # converting to julian date the data from model
    for i in tm:
        modelTimestamp.append(pd.Timestamp(i).to_julian_date())

    observTimestamp = []

    # converting to julian date the data from BNDO
    for i in to:
        observTimestamp.append(pd.Timestamp(i).to_julian_date())
        
    modelTimestamp  = np.asarray(modelTimestamp) 
    observTimestamp = np.asarray(observTimestamp)

    # search for the closest values between model and observTimestamp
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx

    indices = []

    for i in modelTimestamp:
        indices.append(find_nearest(observTimestamp, i))

    indices = np.asarray(indices)

    ob=observ
    ob['ind'] = np.arange(0,len(observ))
    ob['datetime'] = observ.index.values
    ob.set_index('ind', inplace=True)

    # selected
    n = []
    for i in indices:
        value = ob[ob.index == i].nivel.values
        n.append(value)

    return np.squeeze(np.asarray(n))