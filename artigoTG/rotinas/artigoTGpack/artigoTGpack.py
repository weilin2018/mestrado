# arquivo contendo funcoes utilizadas no artigo

import numpy as np
import xray as xr
import pandas as pd
from scipy.spatial import cKDTree
from scipy import signal, fftpack
import scipy

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

def find_nearest(glon,glat,ilon,ilat):
    '''
        glon,glat = lat e lon da grade
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

def criar_diretorio_BASE_DIR(pwd):

    BASE_DIR = pwd.split('/')[:-3]

    BASE_DIR = '/'.join((BASE_DIR))
    return BASE_DIR

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
