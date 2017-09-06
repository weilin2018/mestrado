"""
    Lista 4 - mapeamento de funções de corrente a partir de dados sinóticos e a partir de dados observados

    Cruzeiro: WESTRAX 2

    1.A: anomalia do geopotencial
        . calcular anomalia do geopotencial dos pontos observados para 10dbar relativamente a 1200dbar
        . dividir por f0 calculado para latitude central de 5N
        . remover a média

    1.B: condições de contorno
        Condição de contorno no-slip
            . tomar lat/lon para isóbata de 200m (GEBCO ou ETOPO)
            . suavizar usando média móvel, janela móvel ou interpolação cúbica
    1.C:
        . selecionar latlon das isóbatas e atribuir valor ZERO na borda

    1.D:
        . combinar os vetores de velocidade com o vetor gerado pela isóbata

    1.E:
        . scaloa

    1.F:
        . plotar


    A tarefa referente ao mapeamento de função de corrente a partir dos dados sinóticos é a seguinte.
    O cruzeiro a ser utilizado é o WESTRAX 2 e tenham em mãos o meu velho artigo de 2000. Adotem para
    a AOE e AOV os parâmetros que forneço. Cuidado com as unidades das saídas! Sugiro a conversão da
    velocidade do perfilador Pegasus para graus/s quando entrarem na AOV. Depois da interpolação,
    reconvertam para m/s.

    1)  Mapeamento de Função de Corrente Geostrófica​

    ​a) ​ Calculem a função de corrente geostrófica nos pontos dos dados, ou seja, peguem os pontos de
    obs. e calculem a anomalia do geopotencial para 10 dbar relativamente a 1200 dbar. A posteriori,
    dividam por f0 com f avaliado a 5N (tal qual no artigo da minha tese). Removam a média do vetor
    formado pelos valores de psi_g.

    b) Precisamos satisfazer as condições de contorno. Comecemos com a mais fácil de implementar: no
    slip. Obtenham de arquivo a isóbata de 200 m, seja pelo ETOPO ou GEBCO. Suavizem-na com uma média
    corrida, janela móvel ou reinterpolem usando um esquema cúbico. Vcs escolhem.

    c) Selecionem lats e lons dessa isóbata e atribuam o valor de zero na borda.

    d) Componha o vetor juntando os dados de verdade com os dados falsos da c.c.

    e) Rode a rotina scaloa.m ou a versão python que o Iurizinho.

    ​f) Plote e verifique se o padrão é semelhante ao do artigo.

    g) Agora repita o mesmo caso, mas use imagens. Faça a aproximação linear da isóbata de 200 m e a
    use como eixo de simetria.​


    ​2) Mapeamento de Função de Corrente Observada

    a) o procedimento é semelhante ao da AOE, mas não removam a média dos vetores. Mapeiem a psi_obs
    em 10 dbar. A rotina tem um parâmetro especial que estabelece o nível médio de psi_obs=0.

    b) Usem tanto a no-slip como as imagens (free slip)​ como cc.

    c) Comparemos os padrões.

    d) se quiserem, podem fazer também a velocidade relativa a 1200 dbar para direta comparação entre
    psi_g e psi_obs. Para tanto, basta subtrairem o vetor de 1200 dbar do vetor de 10 m.

    ​NÃO EMPAQUEM. EM CASO DE DÚVIDA, PROCUREM-ME. ESTOU SEM AULAS PELA TARDE, A MENOS DAS 2 QUARTAS
    DE REPOSIÇÃO DE DFG2.​

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import seawater as sw
import gsw
import pandas as pd
import xarray as xr
import cmocean as cmo
import scipy.interpolate as scint
from mpl_toolkits.basemap import Basemap # módulo de mapas
from math import factorial

import matplotlib
matplotlib.style.use('ggplot')

###############################################################################
#                                                                             #
#                                 FUNÇÕES                                     #
#                                                                             #
###############################################################################
def create_VerticalStation(Z=4861):
    """
    create a vector that represent a vertical profile with a Z maximum depth

    INPUT:
        Z length

    OUTPUT:
        vector1D with (Z,) shape
    """

    vector1D = np.zeros((Z,))*np.nan

    return vector1D

def plotar_mapa(files,fEtopo='',savefig=""):
    """
    plotar mapa com a localização das estações de coleta do cruzeiro
    WESTRAX2
    """
    # extract data from etopo1

    etopo = xr.open_dataset(fEtopo)
    etopo_lat = etopo['lat']
    etopo_lon = etopo['lon']
    etopo_bat = etopo['Band1']

    # start by extract the position of each station from files list
    lons,lats = [],[]
    for fname in files:
        f = pd.read_csv(fname, delim_whitespace=True, names=['lat', 'lon'], nrows=1, usecols=[4,5])
        lons.append(float(f['lon']*(-1)))
        lats.append(float(f['lat']))

    lons = np.asarray(lons)
    lats = np.asarray(lats)

    longitude_min,longitude_max = -60.0,-40.0
    latitude_min,latitude_max   = -04.0, 16.0

    fig, ax = plt.subplots()

    m = Basemap(llcrnrlon=longitude_min, llcrnrlat=latitude_min,
                urcrnrlon=longitude_max, urcrnrlat=latitude_max,
                resolution='h',projection='merc')
    m.ax = ax

    m.drawcoastlines(linewidth=.5)
    m.fillcontinents(color='white')

    # plotar batimetria
    x,y = np.meshgrid(etopo_lon, etopo_lat)
    b_x, b_y = m(x,y)
    m.contourf(b_x, b_y, etopo_bat, cmap=cmo.cm.deep)
    m.contour(b_x, b_y, etopo_bat, [-200])
    lbs=[True, False, False, True]
    m.drawmeridians(range(-180,180,4),labels=lbs,linewidth=.3);
    m.drawparallels(range(-90,90,4),labels=lbs,linewidth=.3);

    # # plotar as estações
    # m.scatter(lons, lats, s=50, color='red', latlon=True, label='WESTRAX2')
    # plt.legend(loc='best', scatterpoints=1)

    names = list_StationsName(files)
    for i,j,n in zip(lons,lats,names):
        m.scatter(i,j,s=30,color='red', latlon=True)
        i,j = m(i,j)
        ax.text(i,j,n,ha='left')

    plt.title(u"37 Estações de Coleta - Cruzeiro WEXTRAX2")

    if savefig == "":
        plt.show()
    else:
        os.system('clear')
        print('saving fig in: %s' % SAVE_DIR)
        plt.savefig(SAVE_DIR+savefig, dpi=150)

    return lons, lats, etopo_lon, etopo_lat, etopo_bat, m

# criar matrizes para armazenar os dados
def createArrays(ndim=4861, mdim=37, dz=1):
    ### parte crucial: juntandos os dados em uma unica MATRIZ
    # criar matriz para receber dados
    ### importante o len(data) deve corresponder ao máximo de dados do maior perfil
    v = [] # dimensao
    for i in range(0, ndim,dz):
        v.append(np.nan)

    v = np.vstack(v)
    # matriz
    m =[]
    for i in range(0, mdim):
        m.append(v)
    # matriz (len(allFiles), len(data))
    newArray = np.squeeze((m), axis=2)

    return newArray

def list_StationsName(files):
    """ criar lista com os nomes das estações a partir da lista glob """
    names = []
    files.sort()

    for f in files:
        names.append(f.split('/')[-1][5:8])

    return names

BASE_DIR = "/home/tparente/danilo/mestrado/disciplinas/mestrado/IOC5817/lista3/"
DATA_DIR = BASE_DIR + "data/wx2"
SAVE_DIR = BASE_DIR + 'outputs/'

# read files
files = glob.glob(DATA_DIR+'/*.pro')
files.sort()
files = np.asarray(files)              # convert list in np.ndarray

# extract and store temp, salt and pres data into a dataframe for each
# parameter
stations = files.shape[0]              # how many stations
maxDepth = 5006                        # maximum depth registrated

# problem configuration
latitude = 5.                          # positive because is North Hemisphere
f0 = sw.f(latitude)                    # calculate coriolis parameter

# create matrixes
Prof_deCorte = 1200

indexes = []
# check the maximum depth of each station, only stations higher than 600m will keep
for station in range(0,stations-1):
    f = pd.read_csv(files[station], delim_whitespace=True, skiprows=1, names=['pres', 'temp', 'theta', 'salt', 'dens', 'long', 'lat'], usecols=[1,2,3,4,5,8,9])

    # check len() >= 600m
    if len(f['pres']) < 600:
        # files = np.delete(files, station)
        indexes.append(station)        # save index's station to remove later

# remover as estações muito rasas
files = np.delete(files, indexes)
stations = files.shape[0]                   # update stations variable

tempArray = createArrays(ndim=maxDepth, mdim=stations)
saltArray = createArrays(ndim=maxDepth, mdim=stations)
presArray = createArrays(ndim=maxDepth, mdim=stations)
nameArray = list_StationsName(files)        # extract name of each station

# extract and populate each array with real data instead of a lot of nans
for station in range(0,stations):
    # open file
    f = pd.read_csv(files[station], delim_whitespace=True, skiprows=1, names=['pres', 'temp', 'theta', 'salt', 'dens', 'long', 'lat'], usecols=[1,2,3,4,5,8,9])

    for z in range(0,len(f['pres'])):
        tempArray[station,z] = f['temp'][z]
        saltArray[station,z] = f['salt'][z]
        presArray[station,z] = f['pres'][z]

# cortar as profundidades no nível de referencia
temp = tempArray[:,:Prof_deCorte]
salt = saltArray[:,:Prof_deCorte]
pres = presArray[:,:Prof_deCorte]


# extrapolar as estações menores que 1200dbar usando expansão da série de fourier
def meanBlock(P,bloco,dZ):
    """
    calcular a média dentro de um bloco para uma propriedade P
    """
    v  = P[:,bloco-dZ:bloco]                    # valores dentro do bloco
    n  = np.where(~np.isnan(v))[0].shape[0]     # qntos valores são
    Pn = np.nansum(v)                           # média desses valores

    return Pn/n

def expansaoSerieTaylor(S,T,dZ,window=2):
    """
        (i) calculando o perfil médio da área de estudo, com um deltaZ de 10m.
        Importante lembrar que deve-se tratar das propriedades de forma escalar,
        ou seja, T, S e p independentemente.

        (ii) alisar o perfil médio com média corrida ou janela móvel,
        preservando a variância. Usa-se 'normalized root mean square' (NRMS) e
        então corrige-se o perfil médio alisado, através de
            Ps = (1-alpha)Ps, onde alpha é o resultado de NRMS

        (iii) calcular a primeira e segunda derivada parcial de Ps em z,
        atentando-se para não obter derivadas com valores erroneos

        (iv) série de taylor da forma:
            P*(1:n) = P*(1:n) - aonde tenho dado, mantenho o dado
            P*(n+1:m) = P*(n) + dPs/dz deltaZ + d²Ps/dz² (deltaZ/2!)

            sendo m o indice referente a 1200dbar
    """
    os.system('clear')

    # (i) obtendo o perfil médio da área de estudo para cada propriedade
    newSalt = []
    newTemp = []
    newPres = np.arange(0,Prof_deCorte,10)

    for bloco in range(10,Prof_deCorte+10,10): #deltaZ = 10m
        # somar todos os valores dentro do limite do bloco (somatória de Pn)
        newSalt.append(meanBlock(S,bloco,dZ))
        newTemp.append(meanBlock(T,bloco,dZ))

    newSalt = np.asarray(newSalt)
    newTemp = np.asarray(newTemp)

    # plotar o perfil vertical de temperatura e salinidade para verficar o
    # perfil médio obtido até agora
    print("Plotando perfil médio")
    fig, ax = plt.subplots(ncols=2,nrows=1,sharey=True)

    ax[0].plot(newTemp, -newPres, 'k')
    # ax[0].axhline(y=-1200, color='gray', linestyle='-')
    ax[0].set_xlabel(r'Temperature [$^{o}C$]')
    ax[0].set_ylabel(u'Depth [m]')

    ax[1].plot(newSalt, -newPres, 'k',label=u'Perfil Médio')
    ax[1].set_xlabel(r'Salinity')

    # plt.suptitle(u"Perfil Médio sem Alisar",fontsize=18)

    # (ii) alisar o perfil médio obtido em (i) mantendo a variância deste
    # para suavizar o perfil usando alguma ferramenta, convem transformar o
    # perfil em um pandas.dataframe:

    Tsmoothed = pd.rolling_mean(newTemp,window,center=True)
    Ssmoothed = pd.rolling_mean(newSalt,window,center=True)
    # fig, ax = plt.subplots(ncols=2,nrows=1,sharey=True)
    print("Plotando perfil médio suavizado")
    ax[0].plot(Tsmoothed, -newPres, 'r')
    # ax[0].axhline(y=-1200, color='gray', linestyle='-')
    ax[0].set_xlabel(r'Temperature [$^{o}C$]')
    ax[0].set_ylabel(u'Depth [m]')

    ax[1].plot(Ssmoothed, -newPres, 'r', label=u'Janela móvel de 2m')
    ax[1].set_xlabel(r'Salinity')

    # compute the normalized root mean square error (NRMS)
    alphaT = np.sqrt( (np.nanmean(newTemp - Tsmoothed))**2 )/ np.nanmean(newTemp)
    alphaS = np.sqrt( (np.nanmean(newSalt - Ssmoothed))**2 )/ np.nanmean(newSalt)

    # manter a variância no perfil suavizado
    Ts = (1-alphaT)*Tsmoothed
    Ss = (1-alphaS)*Ssmoothed

    print("Plotando perfil médio suavizado e corrigido")
    ax[0].plot(Ts, -newPres, 'g')
    # ax[0].axhline(y=-1200, color='gray', linestyle='-')
    ax[0].set_xlabel(r'Temperature [$^{o}C$]')
    ax[0].set_ylabel(u'Depth [m]')

    ax[1].plot(Ss, -newPres, 'g', label=r'Perfil Corrigido com $\alpha$')
    ax[1].set_xlabel(r'Salinity')

    plt.suptitle(u"Perfil Médio da área de estudo",fontsize=18)
    plt.legend(loc='lower right')
    plt.show()

    # calcular as derivadas parciais
    # dP/dz
    dsdz = np.diff(Ss)/dZ
    dTdz = np.diff(Ts)/dZ
    # d²P/dz²
    d2sdz = np.diff(dsdz)/dZ
    d2Tdz = np.diff(dTdz)/dZ

    # verificar o comportamento das derivadas
    print("Plotando as derivadas do perfil médio")
    fig, ax = plt.subplots(ncols=2, nrows=1, sharey=True)

    #plotar temperatura
    ax[0].plot(dTdz, -newPres[:-1], 'k', label=r'$\frac{\partial T}{\partial z}$')
    ax[0].plot(d2Tdz, -newPres[:-2], 'r', label=r'$\frac{\partial^{2} T}{\partial z^{2}}$')
    plt.legend(loc='lower right')

    ax[1].plot(dsdz, -newPres[:-1], 'k', label=r'$\frac{\partial S}{\partial z}$')
    ax[1].plot(d2sdz, -newPres[:-2], 'r', label=r'$\frac{\partial^{2} S}{\partial z^{2}}$')

    ax[0].set_xlabel(r'Temperature [$^{o}C$]')
    ax[0].set_ylabel(u'Depth [m]')
    ax[1].set_xlabel(r'Salinity')

    plt.suptitle(u"Derivadas Parciais do Perfil Médio da Temperatura (esquerda) e Salinidade (direita)",fontsize=18)
    plt.legend(loc='lower right')
    plt.show()

    return dsdz, d2sdz, dTdz, d2Tdz

# calcular as derivada dos perfis médios (expansão por série de Taylor)
dsdz, d2sdz, dTdz, d2Tdz = expansaoSerieTaylor(salt, temp, 10, window=2)

# iniciar o processo de iteração das estações e extrapolação destas
# a chave aqui é verificar os pontos nan, dada a forma como eu estabeleci o problema
def media(P):

    nP = []
    for bloco in range(10,Prof_deCorte+10,10):
        nP.append(np.nanmean(P[bloco-10:bloco]))

    return np.asarray(nP)

# regridar os dados para um valor a cada 10m - media em bloco de 10m
temp2 = createArrays(ndim=Prof_deCorte, mdim=stations, dz=10)
salt2 = createArrays(ndim=Prof_deCorte, mdim=stations, dz=10)
pres2 = np.arange(0,Prof_deCorte,10)

os.system('clear')
for station in range(0,stations):
    # calcular a media em bloco de 10metros de profundidade
    temp2[station,:] = media(temp[station,:])
    salt2[station,:] = media(salt[station,:])

    print('Plotting station #%i' %(station))
    fig, ax = plt.subplots(ncols=2, nrows=1, sharey=True)
    ax[0].plot(temp2[station,:], -pres2, 'k')
    ax[1].plot(salt2[station,:], -pres2, 'k')

    # extrapolar as estações incompletas
    # pegar os indices de nan
    indT = np.where(np.isnan(temp2[station,:]))[0]
    indS = np.where(np.isnan(salt2[station,:]))[0]
    if len(indT) > 0:
        # pegar o valor do ultimo dado registrado
        Tn = temp2[station, indT[0]-1]
        Sn = salt2[station, indS[0]-1]
        for n in indT[:-2]: # aqui vamos até o penultimo item, pois ao derivar a expansão, perde-se dados
            # completar os dados por expansão da série de Taylor pra ordem 2
            temp2[station,n] = Tn + dTdz[n]*10 + d2Tdz[n]*5
            salt2[station,n] = Sn + dsdz[n]*10 + d2sdz[n]*5

    # plotar a estação extrapolada em cima da incompleta
    ax[0].plot(temp2[station,:], -pres2, 'r--')
    ax[1].plot(salt2[station,:], -pres2, 'r--')

    plt.show()

# calcular o gpan
gpanArray = createArrays(ndim=Prof_deCorte/10, mdim=stations)

for station in range(0,stations):
    # calc geopotential anomaly
    s = salt2[station,:]
    t = temp2[station,:]
    p = pres2[:]

    # calculate geopotential anomaly for each station
    gpanArray[station,:] = sw.gpan(s,t,p)


# calculate geopotential anomaly with respect to 1200dbar as reference level
gpanRef = gpanArray[7, -1]              # reference level of 1200dbar

gpan = gpanArray[:, :] - gpanRef        # gpan with reference level
gpan = -gpan                            # corrigir sentido da integral

psi_g = gpan/f0                         # psi: geostrophic streamfunction

psi10dbar = psi_g[:, 1]                 # selecting 10dbar as level

psi = psi10dbar - np.nanmean(psi_g) # remove mean from vector

"""
DANILO DE SEGUNDA:

SÓ CONTINUAR A CALCULAR GPAN E AS OUTRAS PARAFERNALHAS. A PARTE DE
EXTRAPOLAÇÃO IS DONE!
"""

# obter dados da profundidade de 200m para implementação da condição de contorno
# dados obtidos em: https://maps.ngdc.noaa.gov/viewers/wcs-client/

# plotar mapa com a localização das estações do cruzeiro
fEtopo = BASE_DIR+'data/etopo1.nc'
lons, lats, blon,blat,bathy, m =  plotar_mapa(files,fEtopo=fEtopo,savefig='batimetria.png')

# extrair coordenadas dos pontos de 200m
ind = np.where(bathy == -200)[0]
lonBoundary = blon.data[ind]
latBoundary = blat.data[ind]

# gerar vetor com psi = 0 para contorno
psiBoundary = np.zeros(lonBoundary.shape)

"""
Comentário importante:

    lonBoundary, latBoundary dizem respeito as coordenadas da isolinha de
    200m de profundidade, que serão utilizadas como condição de contorno
    na resolução do problema.

    lons,lats são as coordenadas das estações realizadas


    Elaborar uma grade regular para inserir os dados, sabendo que:
        lonBoundary,latBoundary forma um vetor de dados zero [condição de contorno]
        lons,lats form outro vetor com dados de gpan/f0

        interpolar os dados dentro dessa grade para calcular a
        velocidade geostrófica
"""

# import Iury's grid
grid = np.load(BASE_DIR+'data/westrax_grd_IT.npy')

gLon,gLat = m(grid[0], grid[1])
m.plot(gLon,gLat,'k', alpha=.4)
m.plot(gLon.T, gLat.T,'k', alpha=.4)

plt.savefig(SAVE_DIR+'grade_iury.png')

# plotar a grade na imagem de coleta

# juntar os dados das estações com os dados da condição de contorno
# bem como as latitudes e longitude, tudo em 3 vetores (lon,lat,psi)
x,y,t = [],[],[]

for i,j,k in zip(lons,lats,psi):
    x.append(i)
    y.append(j)
    t.append(k)

for i,j,k in zip(lonBoundary,latBoundary,psiBoundary):
    x.append(i)
    y.append(j)
    t.append(k)

x,y,t = map(np.asarray, (x,y,t))    # convertendo as listas em ndarray

def usar_griddata(x,y,t,xi,yi):
    """
    x,y = coordenadas atuais x(n,), y(n,)
    t = vetor de dados 1D t(n)
    xi,yi = coordenadas para interpolação xi(n,m), yi(n,m)

    return data(n,m)
    """
    from matplotlib.mlab import griddata

    data = griddata(x,y,t,xi,yi,interp='linear')

    return data

# use scaloa function

def scaloa(xc, yc, x, y, t=None, corrlenx=None,corrleny=None, err=None, zc=None):
  """
      Scalar objective analysis. Interpolates t(x, y) into tp(xc, yc)
      Assumes spatial correlation function to be isotropic and Gaussian in the
      form of: C = (1 - err) * np.exp(-d**2 / corrlen**2) where:
      d : Radial distance from the observations.
      Parameters
      ----------
      corrlen : float
      Correlation length.
      err : float
      Random error variance (epsilon in the papers).
      Return
      ------
      tp : array
      Gridded observations.
      ep : array
      Normalized mean error.
      Examples
      --------
      See https://ocefpaf.github.io/python4oceanographers/blog/2014/10/27/OI/
      Notes
      -----
      The funcion `scaloa` assumes that the user knows `err` and `corrlen` or
      that these parameters where chosen arbitrary. The usual guess are the
      first baroclinic Rossby radius for `corrlen` and 0.1 e 0.2 to the sampling
      error.
  """
  corrlen = corrleny
  xc = xc*( corrleny*1./corrlenx)
  x = x*(corrleny*1./corrlenx)

  n = len(x)
  x = np.reshape(x, (1, n))
  y = np.reshape(y, (1, n))
  # Squared distance matrix between the observations.
  d2 = ((np.tile(x, (n, 1)).T - np.tile(x, (n, 1))) ** 2 +
  (np.tile(y, (n, 1)).T - np.tile(y, (n, 1))) ** 2)
  nv = len(xc)
  xc, yc = np.reshape(xc, (1, nv)), np.reshape(yc, (1, nv))
  # Squared distance between the observations and the grid points.
  dc2 = ((np.tile(xc, (n, 1)).T - np.tile(x, (nv, 1))) ** 2 +
  (np.tile(yc, (n, 1)).T - np.tile(y, (nv, 1))) ** 2)
  # Correlation matrix between stations (A) and cross correlation (stations
  # and grid points (C))
  A = (1 - err) * np.exp(-d2 / corrlen ** 2)
  C = (1 - err) * np.exp(-dc2 / corrlen ** 2)
  if 0: # NOTE: If the parameter zc is used (`scaloa2.m`)
    A = (1 - d2 / zc ** 2) * np.exp(-d2 / corrlen ** 2)
    C = (1 - dc2 / zc ** 2) * np.exp(-dc2 / corrlen ** 2)
  # Add the diagonal matrix associated with the sampling error. We use the
  # diagonal because the error is assumed to be random. This means it just
  # correlates with itself at the same place.
  A = A + err * np.eye(len(A))
  # Gauss-Markov to get the weights that minimize the variance (OI).
  tp = None
  ep = 1 - np.sum(C.T * np.linalg.solve(A, C.T), axis=0) / (1 - err)
  if t!=None:
    t = np.reshape(t, (n, 1))
    tp = np.dot(C, np.linalg.solve(A, t))
    #if 0: # NOTE: `scaloa2.m`
    #  mD = (np.sum(np.linalg.solve(A, t)) /
    #  np.sum(np.sum(np.linalg.inv(A))))
    #  t = t - mD
    #  tp = (C * (np.linalg.solve(A, t)))
    #  tp = tp + mD * np.ones(tp.shape)
    return tp, ep

  if t==None:
    print("Computing just the interpolation errors.")
    #Normalized mean error. Taking the squared root you can get the
    #interpolation error in percentage.
    return ep

# para facilitar a compreensão do código, vamos renomear as variáveis:
data = t
xi = grid[0]
yi = grid[1]

# usar um griddata estranho
data2 = usar_griddata(x,y,t,grid[0],grid[1])

# gridar x,y dos dados que temos (estações e contorno)
x,y = np.meshgrid(x,y)
# gridar os dados para o shape novo (80,80)
points = np.array([x,y])
points = points.T

data = scint.griddata(points, t, (grid[0], grid[1]), method='linear')

tp, ep = scaloa(xi,yi, x, y, data, corrlenx=2, corrleny=2, err=0.1)


u,v = psi2uv(lon,lat,psi10dbar)


#### griddata
x = np.asarray(lons)
y = np.asarray(lats)

xx, yy = np.meshgrid(x,y)
xx, yy = xx.T, yy.T

D = psi10dbar

DI = scint.griddata( (xx,yy), D.ravel(), (xx,yy), method='cubic')
