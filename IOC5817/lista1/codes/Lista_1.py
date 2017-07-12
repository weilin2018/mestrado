#-*-coding:utf-8-*-
"""

    Peço, que usando a radial do DEPROAS I2001, sigam os seguintes procedimentos.

    [OK] 1) Sigam o protocolo de tratamento básico de dados de CTD:  despike, loopedit, binning e filtering;

    [OK] 2) Plote perfis verticais de theta e S para cada estação [paineis duplos com suplots esquerda e direita]. Plotem em 
    pontilhado a curva bruta e linha sólida, a curva tratada;

    3) [OK] Plote o diagrama theta-S espalhado e arco-íris com as isopicnais como curvas paramétricas. 
       [OK] Plotem a curva theta-S espalhada da climatologia WOA para a região da radial.

    4) Apresentem seções verticais da radial para theta, salinidade e sigma_0.  Use como  topografia  a última profundidade 
    amostrada. Para criar o perfil topográfico, cria um vetor distância em km e usem uma interpolação cubica ou spline para 
    suavizar.  Usem o griddata para interpolar a seção. Escolha um método (eu gosto do v4, que é baseado no Akima).


    TODO::

    . deploy(): keep coding


    Comparar as seções vertical com a página 36, do pdf, ou 24 do arquivo em:
	ftp://ftp.io.usp.br/lado/bsctheses/pdf_version/BScDiogo.pdf

"""

# clear screen
import os
os.system('clear')

## --------------------------- ##
##    1. import packages       ##
## --------------------------- ##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import matplotlib
from glob import glob
import gsw 
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
import cmocean as cmo
import scipy.interpolate as scint
import seawater as sw

# function of pyoceans, created by Filipe Fernandes, available in: https://github.com/pyoceans/python-ctd
#from ctd import gen_topomask 

# import my personal package
import sys
sys.path.append('/home/danilo/Dropbox/pythocean/daniloToolbox/')
# function to create absolute path for some directories
#from misc import create_directories_path

# :)

matplotlib.style.use('ggplot')
# Yuri's pack
from OceanLab import CTD # misc of tools to use in CTD data

from matplotlib import rcParams

## --------------------------- ##
##    2. global variables      ##
## --------------------------- ##
DEPLOY = True

if DEPLOY:
    BASE_DIR                 = '/home/tparente/danilo/mestrado/disciplinas/sinoticos/lista1/'
    DATA_DIR                 = BASE_DIR + 'Deproas2/'
    SAVE_DIR                 = BASE_DIR + 'outputs/'
else:
    dirs = create_directories_path(SAVE_DIR=True, DATA_DIR=True) # retorna base_dir, save_dir e data_dir
    BASE_DIR                 = dirs['BASE_DIR']
    DATA_DIR                 = dirs['DATA_DIR']
    SAVE_DIR                 = dirs['SAVE_DIR']

VerticalProfile_basename = 'VertProfile_'
TSDiagram_name           = 'TSDiag_'
VerticalSection_name     = 'verticalSection_subplots.png'

os.chdir(BASE_DIR)

## --------------------------- ##
##    3. define functions      ##
## --------------------------- ##

""" funcoes para interpolacao da secao vertical """
def extrap1d(interpolator):
    """
    http://stackoverflow.com/questions/2745329/
    How to make scipy.interpolate return an extrapolated result beyond the
    input range.
    """
    xs, ys = interpolator.x, interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
        elif x > xs[-1]:
            return (ys[-1] + (x - xs[-1]) * (ys[-1] - ys[-2]) /
                    (xs[-1] - xs[-2]))
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike

def extrap_sec(data, dist, depth, w1=1., w2=0):
    """
    Extrapolates `data` to zones where the shallow stations are shadowed by
    the deep stations.  The shadow region usually cannot be extrapolates via
    linear interpolation.
    The extrapolation is applied using the gradients of the `data` at a certain
    level.
    Parameters
    ----------
    data : array_like
          Data to be extrapolated
    dist : array_like
           Stations distance
    fd : float
         Decay factor [0-1]
    Returns
    -------
    Sec_extrap : array_like
                 Extrapolated variable
    """
    from scipy.interpolate import interp1d

    # interpolação em cada nível de profundidade (x)
    new_data1 = [] # variavel output
    for row in data: # leitura de cada nível de profundidade
        mask = ~np.isnan(row) # mascara para dados que não são np.nan em cada nível
        if mask.any(): # se tiver algum ponto que não seja nan
            y = row[mask] # atribui esses valores a uma variável auxiliar y
            if y.size == 1: # se a quantidade de pontos for 1
                row = np.repeat(y, len(mask)) # repete o mesmo valor len(mask) vezes no nível
            else:  # se a quantidade de pontos for maior que 1
                x = dist[mask] # pega os pontos de distancia referentes aos pontos que não são nan
                f_i = interp1d(x, y) # interpola esses pontos
                f_x = extrap1d(f_i) # cria a funcao para extrapolar esses pontos com os pontos interpolados
                row = f_x(dist) # extrapola esses pontos para o tamanho de dist
        new_data1.append(row) # insere o novo nível na variavel output

    # interpolação em cada nível de distancia (y)
    new_data2 = []
    for col in data.T: # le a matriz de dados transposta (colocando o que era distancia no eixo vertical0)
        mask = ~np.isnan(col) # verifica dados not nan
        if mask.any(): # se existir pontos não nan
            y = col[mask] # atribui os valores a uma variavel auxliar
            if y.size == 1: # checa o tamanho da varivel auxliar: se tiver apenas 1 ponto
                col = np.repeat(y, len(mask)) # repete esse ponto ao longo de y
            else: # se tiver mais de 1 pontos
                z = depth[mask] # recupera os pontos de profundidade
                f_i = interp1d(z, y) # interpola essa coluna
                f_z = extrap1d(f_i) # cria funcao para extrapolar a coluna
                col = f_z(depth) # extrapola efetivamente a coluna
        new_data2.append(col) # insere a nova coluna na variavel output
 
    new_data = np.array(new_data1) * w1 + np.array(new_data2).T * w2 # tira a média ponderada das duas variaveis output
    # retorna dados
    return new_data

# ler o arquivo e armazenar em um pandas.DataFrame
def readCTD(fname, down_cast=True):
    """ 
    ########################################
    #                                      #
    #   function to read .cnv data (fname) #
    #                                      #
    #   input: string                      #
    #   output: pandas.dataframe           #
    #                                      #
    ########################################
    """
    # names = ['press', 'lat', 'lon', 'x', 'y', 'z', 'w', 'salt', 'temp', 'a', 'b']
    names = ['press', 'temp', 'salt'] # for cols 1, 8 and 9
    cast = pd.read_csv(fname, header=None, delim_whitespace=True, names=names, usecols=[1,2,8], skiprows=10)

    # remove negative pressure values
    ind = np.where(cast['press'] < 0)
    cast = cast.drop(cast.index[ind])    

    press_tmp = cast['press'].values

    # set pressure as a new index 
    cast.set_index('press', drop=True, inplace=True)
    # rename index to 'Pressure [db]'
    cast.index.name = 'Pressure [db]'

    cast['press'] = press_tmp

    # slice based on down and upcast
    dwn, up = CTD.abv_water(cast)

    if down_cast:
        return dwn
    else:
        return up

# realizar tratamento/processadmento dos dados, conforme apresentado em aula
def processCTD(data, looped=True, hann_f=False, hann_block=11, hann_times=2):
    """ 
    ########################################
    #                                      #
    #                                      #
    #      Apply some basic treatments:    #
    #           . cut down and up          #
    #           . remove spikes            #
    #           . bin averaging            #
    #           . hanning filter           #
    #                                      #
    #                                      #
    ########################################    
    """
    # remover voltas que o aparelho dá quando na água
    if looped:
        data = CTD.loopedit(data)

    # remover spikes usando um
        # removing spikes of data
    if (data.shape[0]<101)&(data.shape[0]>10): # se o tamanho do perfil for com menos de 101 medidas

        if (data.shape[0]/2)%2 == 0: # caso a metade dos dados seja par
            blk = (data.shape[0]/2)+1 # bloco = a metade +1
        else:
            blk = data.shape[0]/2 # se for impar o bloco e a metade

        # remove spikes dos perfis de temperatura e condutividade
        data = CTD.despike(data,propname='temp',block=blk,wgth=2)
        data = CTD.despike(data,propname='salt',block=blk,wgth=2)
    elif data.shape[0]>=101:
        # para perfis com mais de 101 medidas, utiliza-se blocos de 101
        data = CTD.despike(data,propname='temp',block=101,wgth=2)
        data = CTD.despike(data,propname='salt',block=101,wgth=2)
    else:
             print 'radial muito rasa'


    # bin averaging, with 1m of box - 1m pq a velocidade de descida do instrumento é de 1m/s
    data = CTD.binning(data,delta=1.)


    # hann filtering - alisando o perfil
    if hann_f:
        times=0
        while times<hann_times:
                data = CTD.hann_filter(data,'temp',hann_block)
                data = CTD.hann_filter(data,'salt',hann_block)
                times +=1


    # return data
    return data

# obter o ID da estacao a partir do nome do arquivo sendo processado
def get_profileID(fname):
    """ 
        return only the name file (without full path and extension)
    """

    return fname.split('/')[-1][3:-4]

# obter a coordenada da estacao a partir do ID, pesquisando no arquivo posicoes.txt
def get_latlon(fname):

    posicoesFiles = DATA_DIR + '/posicoes.txt'

    with open(posicoesFiles) as f:
        lines = f.readlines()
        for line in lines:
            elementos = line.split(' ')
            # confere se o ID consta na lista de estacoes que estamos trabalhando
            if elementos[-1][:-1] == get_profileID(fname):
                # resgata as coordenadas e já transforma de degº min secs em coordenadas decimais
                lat_estacao = float(elementos[0]) + float(elementos[1])/60
                lon_estacao = float(elementos[2]) + float(elementos[3])/60
                # insere no vetor. O (-1) * por estarmos no HS
                return lon_estacao*(-1), lat_estacao*(-1)

# criar perfil vertical com os dados tratados e originais
def verticalProfiles(original, filtered, fname,savefig=False):
    """ 
    Plot in subplots temperature and salinity, original (dashed line) and filtered 
    (solid line)
    """
    # create figure
    fig = plt.figure()

    # create subplots
    axTemp = plt.subplot2grid((1,2), (0,0))
    axSalt = plt.subplot2grid((1,2), (0,1))

    ### plot temperature
    axTemp.plot(original['temp'], original.index, 'k--');
    axTemp.plot(filtered['temp'], filtered.index, 'k');
    # subplot settings
    axTemp.set_ylim(axTemp.get_ylim()[::-1])
    axTemp.set_xlim([0., 30.])
    # axTemp.set_xlim(xmin=0., xmax=30.)
    axTemp.set_ylabel('Pressure [db]', size=15)
    axTemp.set_xlabel(u'Temperature [$^oC$]')
    axTemp.xaxis.set_label_position('top')
    axTemp.xaxis.set_ticks_position('top')
    axTemp.legend(['Original', 'Filtered'], loc='best')

    ### plot salinity
    axSalt.plot(original['salt'], original.index, 'k--');
    axSalt.plot(filtered['salt'], filtered.index, 'k');
    # subplot settings
    axSalt.set_ylim(axSalt.get_ylim()[::-1])
    axSalt.set_xlim([24., 40.])
    axSalt.set_xlabel('Salinity [psu]')
    axSalt.xaxis.set_label_position('top')
    axSalt.xaxis.set_ticks_position('top')
    axSalt.legend(['Original', 'Filtered'], loc='best')

    # global title
    plt.suptitle("Station " + get_profileID(fname), size=20)

    # savefig
    if savefig:
        plt.savefig(SAVE_DIR + VerticalProfile_basename + get_profileID(fname) +'.png', dpi=150)

# calcular as isopicnais paramétricas
def parametricIsopicnal(data=None):

    """ """

    # creating a grid to calculate density

    # Figure out boudaries (mins and maxs) - automatically
    # smin = data['salt'].values.min() - (0.01 * data['salt'].values.min())
    # smax = data['salt'].values.max() + (0.01 * data['salt'].values.max())
    # tmin = data['temp'].values.min() - (0.1 * data['temp'].values.max())
    # tmax = data['temp'].values.max() + (0.1 * data['temp'].values.max())

    try:
        # load pickle file with si, ti and sigma0 already calculated
        dic = pickle.load(open('/home/danilo/Documents/mestrado/disciplinas/sinoticos/lista1/parametricIsopicnals.pkl', 'r'))
        si = dic['si']
        ti = dic['ti']
        sigma0 = dic['sigma0']

    except:
        # calculate

        # Figure out boundaries, mins and maxs, manually
        smin, smax = 23., 39.
        tmin, tmax = -10., 28.

        # Calculate how many gridcells we need in the x and y dimensions
        xdim = round((smax-smin)/0.1+1,0)
        ydim = round((tmax-tmin)+1,0)

        # Create temp and salt vectors of appropiate dimensions
        ti = np.linspace(1,ydim-1,ydim)+tmin
        si = np.linspace(1,xdim-1,xdim)*0.1+smin

        # Create empty grid of zeros
        sigma0 = np.zeros((int(ydim),int(xdim)))
        # Loop to fill in grid with densities
        for j in range(0,int(ydim)):
            for i in range(0, int(xdim)):
                # for z in range(0, zdim):
                #     rho[j ,i] = gsw.rho_t_exact(si[i], ti[j], z)
                sigma0[j,i]=gsw.sigma0(si[i],ti[j])

        # save into pickle file
        dic = {"si": si, "ti": ti, "sigma0": sigma0}
        pickle.dump(dic, open("/home/danilo/Documents/mestrado/disciplinas/sinoticos/lista1/parametricIsopicnals.pkl", "w"))

    return si, ti, sigma0

# criar o diagrama TS espalhado: um por estação
def diagramaTS_espalhado(ax, data, fname):
    """ 
        inputs:
            ax      = eixo do plot
            data    = pandas.dataframe com dados
            fname   = nome do arquivo sendo processado (para pegar ID da estação)

        output:
            o output é inserido automaticamente no axis da figura criada

        tip:: após rodar a função, apenas usar plt.show() ou plt.savefig()
    """
    si, ti, sigma0 = parametricIsopicnal(data)

    # 2. plotar contorno com a densidade
    # plot isopicnals
    CS = ax.contour(si, ti, sigma0, linestyles='dashed', colors='k', alpha=.5)
    # CS = ax.contour(si, ti, rho, linestyles='dashed', colors='k')
    ax.clabel(CS, fontsize=12, inline=1, fmt='%1.0f')

    # plotar os registros da estação
    cs = ax.scatter(data['salt'], data['temp'], c='k', linewidth=0.1, edgecolor='k')

    ax.grid()
    ax.set_xlim([33, 38])
    ax.set_ylim([0, 28])
    ax.set_xlabel(u'Salinity [psu]')
    ax.set_ylabel(u'Temperature [$^\circ$C]')
    # plt.title("Station " + get_profileID(fname))

# criar o diagrama TS arco íris: por enquanto um por estação
def diagramaTS_arcoiris(ax, data, fname):
    """
        inputs:
            ax      = eixo do plot
            data    = pandas.dataframe com dados
            fname   = nome do arquivo sendo processado (para pegar ID da estação)

        output:
            o output é inserido automaticamente no axis da figura criada

        tip:: após rodar a função, apenas usar plt.show() ou plt.savefig()

    """
    ################################################################################################
    #                                                                                              #
    #   plotar isopicnais como background, inserindo label para cada contorno em pontos discretos  #
    #                                                                                              #
    ################################################################################################
    si, ti, sigma0 = parametricIsopicnal(data)

    manual_locations = [(33.9, 26.5), (35.3, 26.5), (36.7, 26.5), 
                        (37.2, 24.4), (37.3, 21.3), (40.0, 20.0),
                        (40.0, 20.0), (40.0, 20.0), (40.0, 20.0),
                        (37.6, 14.5), (37.5, 8.80), (37.7, 2.11)]

    levels = np.arange(22, 31, 1)
    CS = ax.contour(si, ti, sigma0, levels, linestyles='dashed', colors='k', alpha=.5)
    ax.clabel(CS,  fontsize=12, inline=1, inline_spacing=1, fmt='%1.0f', manual=manual_locations)

    ################################################################################################
    #                                                                                              #
    #   criar informacoes das massas d'agua (nome, localizacao da isopicnal, densidade, cor)       #
    #   calcula sigma0 para ter um grid de sigma0 para plotar                                      #
    #   plota somente as linhas de massas_dagua, colorindo os intervalos com massas_cores          #
    #                                                                                              #
    ################################################################################################

    # preparando o background com as cores de cada massa dagua
    """
        reference: ftp://ftp.io.usp.br/lado/papers/bub_brown_1996.pdf
        South Atlantic Subtropical Tropospheric - SAST (sigma_theta médio de 25.96)
        South Atlantic Central Waters - SACW
        North Atlantic Deep Water - NADW
        Antartic Intermediate Waters - AAIW

	segundo Silveira et al (2000), as massas d'agua que ocupam os primeiros 3000m são:
		AT, ACAS, AIA, ACS, APAN

    """
    massas_nomes = ['AT', 'ACAS', 'AIA', 'ACS', 'APAN']
    massas_locat = [(37.5, 25.8), (37.4, 20.), (37.4, 17.4), (37.4, 16.0), (37.3, 14.6)]
    massas_dagua = [24.7, 25.7, 26.9, 27.32, 27.57, 27.88]
    massas_cores = (    'r', 'darkorange',   'y',   'g',   'b'    )

    # preparando o grid de sigma0 para plotar
    Smassas, Tmassas = np.arange(32, 40., 0.5), np.arange(0, 28.9, 0.01)

    S, T = np.meshgrid(Smassas, Tmassas)
    densmassas = gsw.sigma0(S,T)
    # coordenadas x,y para inserir a label da isopicnal
    manual_locations = [(33.5, 9.7), (34., 4.), (34.47, 3.13), (34.7, 2.6), (35., 2.)]

    # plotar somente densmassas iguais as que eu quero em massas_dagua
    CC = ax.tricontour(S.ravel(), T.ravel(), densmassas.ravel(), massas_dagua[1:], linestyles='solid', colors='k')
    ax.clabel(CC, fontsize=12, inline=1, inline_spacing=1, fmt='%1.2f', manual=manual_locations)
    # pintar os intervalos
    CS = ax.contourf(S, T, densmassas, massas_dagua, colors=massas_cores, alpha=.5)

    # plotar nome das massas d'agua TODO::
    # for loc,name in zip(massas_locat, massas_nomes):
    #     x,y = loc[0], loc[1]
    #     ax.text(x,y, name, rotation=20)

    ax.grid()
    ax.set_xlim([33, 38])
    ax.set_ylim([0, 28])
    ax.set_xlabel(u'Salinity [psu]')
    # ylabel comentado para o caso de subplot
    # plt.ylabel(u'Temperature [$^\circ$C]')

    # plotar os registros da estação
    # cs = ax.scatter(data['salt'], data['temp'], c='k', s=7, linewidth=0.1, edgecolor='w')
    cs = ax.plot(data['salt'], data['temp'], 'k', linewidth=3.)

# plotar diagrama TS espalhado com a climatologia WOA
def diagramaTS_woa(savefig=False):
    """ 
        input: 
            savefig: boolean para salvar ou não a imagem
    """
    # load WOA climatology from Martim
    woa = pickle.load(open('/home/danilo/Documents/mestrado/disciplinas/sinoticos/lista1/perfil_woa_prox7048.pkl', 'r'))

    tclim = woa['Temp']
    sclim = woa['Salt']

    si, ti, sigma0 = parametricIsopicnal()

    plt.figure(figsize=(18, 10))
    # plot isopicnals
    levels = np.arange(22, 31, 1)
    CS = plt.contour(si, ti, sigma0, levels, linestyles='dashed', colors='k', alpha=.5)

    manual_locations = [(33.9, 26.5), (35.3, 26.5), (36.7, 26.5), 
                        (37.2, 24.4), (37.3, 21.3), (37.4, 18.0),
                        (40.0, 20.0), (40.0, 20.0), (40.0, 20.0),
                        (37.6, 14.5), (37.5, 8.80), (37.7, 2.11)]

    plt.clabel(CS, fontsize=12, inline=1, inline_spacing=1,fmt='%1.0f', manual=manual_locations)

    # plot scatter with TS pairs
    plt.scatter(sclim, tclim, c='k', linewidth=0.1, edgecolor='w')

    plt.grid()
    plt.xlim([33, 38])
    plt.ylim([0, 28])
    plt.xlabel('Salinity [psu]', fontsize=18)
    plt.ylabel(r'Temperature [$^\circ$C]', fontsize=18)
    plt.title("WOA Climatology for [-23.5, -41.5]", fontsize=20)
    # savefig
    if savefig:
        plt.savefig(SAVE_DIR + TSDiagram_name + 'WOAclimatology.png', dpi=150, bbox_inches='tight')

# calcular distancia entre estações
def distance(origin, destination):
    #Source one: http://www.platoscave.net/blog/2009/oct/5/calculate-distance-latitude-longitude-python/
    #Stoled from: https://oceanpython.org/2013/03/21/otn-gliders-more-advanced-plotting/
    import math
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 # km
 
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
 
    return d

# resgatar as coordenadas de cada estação
def coordenadasEstacoes(lista_estacoes):
    """ resgatar as coordenadas das estacoes no arquivo posicoes.txt """

    coords_estacoes = []

    with open(DATA_DIR + '/posicoes.txt') as f:
        lines = f.readlines()
        for line in lines:
            elementos = line.split(' ')
            # confere se o ID consta na lista de estacoes que estamos trabalhando
            if elementos[-1][:-1] in lista_estacoes:
                # resgata as coordenadas e já transforma de degº min secs em coordenadas decimais
                lat_estacao = float(elementos[0]) + float(elementos[1])/60
                lon_estacao = float(elementos[2]) + float(elementos[3])/60
                # insere no vetor. O (-1) * por estarmos no HS
                coords_estacoes.append([ (-1) * lat_estacao, (-1) * lon_estacao])

    coords_estacoes = np.asarray(coords_estacoes)

    return coords_estacoes

# calcular a distancia entre as estações e criar vetor de topografia
def distanciaEstacoes(coords_estacoes):
    """ 
    . calcula a distância entre cada estação
    """
    dist = [0. ]

    for i in np.arange(1,len(coords_estacoes)):
        dist.append( distance(coords_estacoes[0], coords_estacoes[i]) )

    return dist

# criar secao vertical
def createVerticalSection(xm, hm, xplot, yplot, props, lista_estacoes, name='', suptitle=None, interpolate=True, topColor='k', savefig=False, showfig=True, alpha=.6, comparar=False):
    """ 
    plotar seção vertical, baseado nos dados enviados:
        xm, hm = x,y para topografia
        x,y = distancia do transecto, profundidade para plotar o contourf
        props = dicionario com propriedades a serem plotadas (temp, salt, dens)

        x,y,props = mesma dimensão

    """

    if comparar:
        topDepth = -1000
        yticks_markers = np.arange(topDepth, 0, 100)
        cmapTemp = matplotlib.cm.jet
        cmapSalt = matplotlib.cm.jet
        cmapDens = matplotlib.cm.jet
    else:
        topDepth = -2280
        yticks_markers = np.arange(topDepth+280, 0, 500)
        cmapTemp = cmo.cm.thermal
        cmapSalt = cmo.cm.haline
        cmapDens = cmo.cm.dense

    # extract data from dictionary
    tempArray = props['tempArray']
    saltArray = props['saltArray']
    densArray = props['densArray']

    # create figure with 3 rows in subplots, sharing x axis
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)
    
    # plot temperature
    levels = np.arange(0, 28, .5)    
    temperaturPlot = ax1.contourf(xplot, -yplot, tempArray.T, levels, cmap=cmapTemp)
    topografiaPlot = ax1.fill_between(xm, topDepth, -hm, color=topColor, interpolate=interpolate, alpha=alpha)
    cbar = plt.colorbar(temperaturPlot, ax=ax1, ticks=np.arange(0, 30, 5.))
    cbar.set_label(u'Temperature [$^o$C]')
    ax1.set_xlim([0, 174])
    ax1.set_ylim([topDepth, 0])
    ax1.yaxis.set_ticks(yticks_markers)
    ax1.set_ylabel('Pressure [dbar]')

    ## salinidade
    levels = np.arange(33., 38., 0.1)
    salinityPlot   = ax2.contourf(xplot, -yplot, saltArray.T, levels, cmap=cmapSalt)
    topografiaPlot = ax2.fill_between(xm, topDepth, -hm, color=topColor, interpolate=interpolate, alpha=alpha)
    cbar = plt.colorbar(salinityPlot, ax=ax2, ticks=np.arange(33, 38.1, .5))
    cbar.set_label(u'Salinity [psu]')
    ax2.yaxis.set_ticks(yticks_markers)
    ax2.set_xlim([0, 174])
    ax2.set_ylim([topDepth, 0])
    ax2.set_ylabel('Pressure [dbar]')

    ## densidade
    levels = np.arange(23, 28.5, 0.1)
    densityPlot    = ax3.contourf(xplot, -yplot, densArray.T, levels, cmap=cmapDens)
    topografiaPlot = ax3.fill_between(xm, topDepth, -hm, color=topColor, interpolate=interpolate, alpha=alpha)
    cbar = plt.colorbar(densityPlot, ax=ax3, ticks=[23, 24, 25, 26, 27, 28, 29, 30])
    cbar.set_label(u'$\sigma_{\Theta}$ [kg m$^{-3}$]')
    ax3.yaxis.set_ticks(yticks_markers)
    ax3.set_xlim([0, 174])
    ax3.set_ylim([topDepth, 0])
    ax3.set_xlabel('Distance along transect [km]')
    ax3.set_ylabel('Pressure [dbar]')

    if suptitle == None:
        subtitle = "Seção Vertical - Deproas I2001"

    plt.suptitle(suptitle, fontsize=20)

    # insert small map with stations
    m2 = insertTransectLocation(ax3, lista_estacoes, bathy=True)

    if savefig:
        plt.savefig(SAVE_DIR + VerticalSection_name[:-4] + name + '.png', dpi=150)

    if showfig:
        plt.show()

# inserir um pequeno mapa no subplot 'ax', com as estações que compoem o transecto plotado
def insertTransectLocation(ax, lista_estacoes, bathy=False):

    """ criar minimapa com a localização do transecto sendo plotado """

    inset = inset_axes(ax, width="40%", height=1.4, loc=3)

    m2 = pickle.load(open(BASE_DIR + '/cabofrio.pkl', 'r'))
    m2.ax = inset

    m2.fillcontinents(color='coral')

    if bathy:
        import xarray as xr
        etopo = xr.open_dataset('/home/tparente/danilo/mestrado/cabofrio/dados/etopo5.nc')
        tlat = etopo['topo_lat']
        tlon = etopo['topo_lon']
        tdep = etopo['topo']

        # gridar coordenadas
        x,y = np.meshgrid(tlon, tlat)

        levels = [-2000, -1000, -200, -100]
        CS = m2.contour(x, y, tdep, levels, colors='k', latlon=True, linestyle='dashed', linewidth=.5)
        ax.clabel(CS, fontsize=8, inline=1, fmt='%1.0f')

    coords_estacoes = coordenadasEstacoes(lista_estacoes)

    for coord,IDstation in zip(coords_estacoes, lista_estacoes):
        lat_estacao = coord[0]
        lon_estacao = coord[1]
        m2.scatter(lon_estacao, lat_estacao, latlon=True, label=IDstation, color='k')

    return m2

# interpolar seção vertical, usando np.griddata, method='nearest'
def gridarSecao(props, lastDepth, dist, method='nearest'):
    """ 

        reference: http://christopherbull.com.au/python/scipy-interpolate-griddata/
    """
    temp = props['tempArray'].T
    salt = props['saltArray'].T
    dens = props['densArray'].T

    # old grid dimension
    loni = np.asarray(dist)
    depi = np.arange(0, int(lastDepth.max()))

    # new grid dimension
    lon  = np.arange(loni.min(), loni.max(), 5)
    dep  = np.arange(0, int(lastDepth.max()))

    # create grid
    X, Y   = np.meshgrid(loni, depi) # old grid
    XI, YI = np.meshgrid(lon, dep) # new grid

    # grid data
    px, py = X.flatten(), Y.flatten()
    temp = scint.griddata((px, py), temp.flatten(), (XI, YI), method=method)
    salt = scint.griddata((px, py), salt.flatten(), (XI, YI), method=method)
    dens = scint.griddata((px, py), dens.flatten(), (XI, YI), method=method)

    # organize output in a dictionary
    dicOutput = {
                "XY": [XI, YI],
        "tempArray": temp.T,
        "saltArray": salt.T,
        "densArray": dens.T
    }

    return dicOutput

# plotar topografia interpolada e não interpolada
def plotarTopografia(xorig, yorig, xinterp, yinterp, savefig=False):
    """
        observar a melhora da topografia com a interpolação
    """
    # plotar topografia
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    # plot crude topography
    ax1.fill_between(xorig, -yorig.max(), -yorig, color='k', interpolate=True, alpha=.6)
    ax1.set_xlim([0, int(np.nanmax(xorig))+1])
    ax1.set_ylim([-int(yorig.max()), 0])
    ax1.set_title("Crude Topography")
    ax1.set_ylabel("Depth [m]")

    # plot interpolated topography
    ax2.fill_between(xinterp, -2270, -yinterp, color='k', alpha=0.6)
    ax2.set_xlim([0, int(np.nanmax(xorig))+1])
    ax2.set_ylim([-int(yorig.max()), 0])
    ax2.set_title("Interpolated Topography [cubic]")
    ax2.set_ylabel("Depth [m]")
    ax2.set_xlabel("Distance along transect [km]")

    if savefig:
        plt.savefig(SAVE_DIR + '/topography_created.png', dpi=150)
        plt.close()
    else:
        plt.show()

# realizar todos os processos - equivalente a um main()
def realize_process(fname, gerarImagens=False):

    print("Processando a estacao %s" % (get_profileID(fname)))
    # leitura do arquivo, retorna somente downcast
    cast = readCTD(fname)
    # copia do dataframe
    data = cast.copy()
    # processamento dos dados
    data = processCTD(data, hann_f=True)
    # resgatar lon,lat da estacao sendo lida
    lon, lat = get_latlon(fname)

    # criar novas variáveis
    # calculate absolute salinity
    data['SA'] = gsw.SA_from_SP(data['salt'].values, data['press'].values, lon, lat)
    # calculate conservative temperature
    data['CT'] = gsw.CT_from_t(data['SA'].values, data['temp'].values, data['press'].values)

    # condição para gerar ou nao imagens, útil durante os testes das secoes verticais
    if gerarImagens:
        print("Gerando perfil vertical")
        # criar o perfil vertical
        verticalProfiles(cast, data, fname, savefig=True)

        print("Gerando Diagrama TS espalhado e arco íris")
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(18,10))
        # criar o diagrama TS espalhado
        diagramaTS_espalhado(ax1, data, fname)
        # criar o diagrama TS arco iris
        diagramaTS_arcoiris(ax2, data, fname)

        plt.suptitle("Station " + get_profileID(fname), fontsize=20)

        savefig = False
        if savefig:
            plt.savefig(SAVE_DIR + TSDiagram_name + get_profileID(fname) +'.png', dpi=150)

        plt.show()

    return data

# realize a simple test, only to check the consistency of the code
def deploy(allFiles):

    """ only for tests """

    fname = allFiles[-1]

    cast = readCTD(fname)
    data = cast.copy()
    data = processCTD(data, hann_f=True)
    lon, lat = get_latlon(fname)
    data['SA'] = gsw.SA_from_SP(data['salt'].values, data['press'].values, lon, lat)
    data['CT'] = gsw.CT_from_t(data['SA'].values, data['temp'].values, data['press'].values)

    return data, lon, lat

# como estou  com preguiça ...
def plotar_todas_estacoes_emTS_arcoiris(allFiles):
    fig, ax = plt.subplots()

    si, ti, sigma0 = parametricIsopicnal()

    manual_locations = [(33.9, 26.5), (35.3, 26.5), (36.7, 26.5), 
                        (37.2, 24.4), (37.3, 21.3), (40.0, 20.0),
                        (40.0, 20.0), (40.0, 20.0), (40.0, 20.0),
                        (37.6, 14.5), (37.5, 8.80), (37.7, 2.11)]

    levels = np.arange(22, 31, 1)
    CS = ax.contour(si, ti, sigma0, levels, linestyles='dashed', colors='k', alpha=.5)
    ax.clabel(CS,  fontsize=12, inline=1, inline_spacing=1, fmt='%1.0f', manual=manual_locations)

    ################################################################################################
    #                                                                                              #
    #   criar informacoes das massas d'agua (nome, localizacao da isopicnal, densidade, cor)       #
    #   calcula sigma0 para ter um grid de sigma0 para plotar                                      #
    #   plota somente as linhas de massas_dagua, colorindo os intervalos com massas_cores          #
    #                                                                                              #
    ################################################################################################

    # preparando o background com as cores de cada massa dagua
    massas_nomes = ['TW', 'SACW', 'AAIW', 'NADW', 'AAIW']
    massas_locat = [(37.5, 25.8), (37.4, 20.), (37.4, 17.4), (37.4, 16.0), (37.3, 14.6)]
    massas_dagua = [24.7, 25.7, 26.9, 27.32, 27.57, 27.88]
    massas_cores = (    'r', 'darkorange',   'y',   'g',   'b'    )

    # preparando o grid de sigma0 para plotar
    Smassas, Tmassas = np.arange(32, 40., 0.5), np.arange(0, 28.9, 0.01)

    S, T = np.meshgrid(Smassas, Tmassas)
    densmassas = gsw.sigma0(S,T)
    # coordenadas x,y para inserir a label da isopicnal
    manual_locations = [(33.5, 9.7), (34., 4.), (34.47, 3.13), (34.7, 2.6), (35., 2.)]

    # plotar somente densmassas iguais as que eu quero em massas_dagua
    CC = ax.tricontour(S.ravel(), T.ravel(), densmassas.ravel(), massas_dagua[1:], linestyles='solid', colors='k')
    ax.clabel(CC, fontsize=12, inline=1, inline_spacing=1, fmt='%1.2f', manual=manual_locations)
    # pintar os intervalos
    CS = ax.contourf(S, T, densmassas, massas_dagua, colors=massas_cores, alpha=.5)

    # plotar nome das massas d'agua TODO::
    # for loc,name in zip(massas_locat, massas_nomes):
    #     x,y = loc[0], loc[1]
    #     ax.text(x,y, name, rotation=20)

    ax.grid()
    ax.set_xlim([33, 38])
    ax.set_ylim([0, 28])
    ax.set_xlabel(u'Salinity [psu]')

    cores_estacoes = [ 'black', 'maroon', 'sienna', 'darkgreen', 'darkblue', 'darkmagenta' ]

    for fname,c, IDstation in zip(allFiles, cores_estacoes, lista_estacoes):
        data = realize_process(fname)

        ax.plot(data['salt'], data['temp'], linewidth=4., c=c, label=IDstation)

    plt.legend(loc=2)

    plt.show()

## --------------------------- ##
##    4. run this program      ##
## --------------------------- ##

# list all files you want to plot -> TODO:: operacionalizar essa listagem
allFiles = glob(DATA_DIR + '/d2_*')
allFiles.sort()

# set deploy (below) to True only if you want to test this code
deploy = False
if deploy:
    data, lon, lat = deploy(allFiles)

# criar o diagrama TS para a climatologia
# diagramaTS_woa(savefig=True)

# variaveis para armazenar dados de cada estacao conforme processamento
temp = []
salt = []
pres = []

# variavel para criar lista de IDs das estacoes
lista_estacoes = []

####################################################
#                     PARTE I                      #
#           i) tratar dados do CTD                 #
#                                                  #
#           ii) plotar seção  vertical             #
#                                                  #
#          iii) plotar diagramas TS                #
#                                                  #
####################################################
for fname in allFiles:
    data = realize_process(fname)

    # armazenar dados em variáveis auxiliares durante o processamento de cada
    # estação, para evitar ter que rodar tudo novamente
    temp.append(np.asarray(data['temp'].values))
    salt.append(np.asarray(data['salt'].values))
    pres.append(np.asarray(data['press'].values))

    outName = SAVE_DIR + get_profileID(fname)+'_data.pkl'
    pickle.dump(data, open(outName, 'w'))

    lista_estacoes.append(get_profileID(fname))

####################################################
#                     PARTE II                     #
#       i) calcular distancia do transecto         #
#                                                  #
#       i) calcular profundidades e gerar          #
#      topografia (interpolada cubic e linear)     #
#                                                  #
#     iii) juntar as propriedades para plotar      #
#                                                  #
####################################################
# coordenadas das estações do arquivo posicoes.txt
coords_estacoes = coordenadasEstacoes(lista_estacoes)

# calcular a distancia do transecto - criar eixo X [km]
dist = distanciaEstacoes(coords_estacoes)
# separa para facilitar a vida mais pra frente
lons, lats = coords_estacoes[:,1], coords_estacoes[:,0]

# calcular as profundidades máximas amostradas - criar eixo Y [m]
lastDepth = []
for stationPress,lat in zip(pres, lats):
    # pegar somente a maior pressão
    lastDepth.append(stationPress.max())

lastDepth = np.asarray(lastDepth)

# True: interpolará cubicamentel, False: mantem
interpTopography = False

if interpTopography:
    dic = pickle.load(open(BASE_DIR + 'topografia.pkl', 'r'))
    xm, hm = dic['xm'], dic['hm']
else:
    xm, hm = dist, lastDepth

# plotar topografia para visualizar a suavização da interpolação
plotarTopografia(dist, lastDepth, xm, hm, savefig=False)

def createArrays():

    ### parte crucial: juntandos os dados em uma unica MATRIZ
    # criar matriz para receber dados
    ### importante o len(data) deve corresponder ao máximo de dados do maior perfil
    v = [] # dimensao
    for i in range(0, int(lastDepth.max())):
        v.append(np.nan)

    v = np.vstack(v)
    # matriz
    m =[]
    for i in range(0, len(allFiles)):
        m.append(v)
    # matriz (len(allFiles), len(data))
    tempArray = np.squeeze((m), axis=2)
    saltArray = np.squeeze((m), axis=2)
    densArray = np.squeeze((m), axis=2)

    return tempArray, saltArray, densArray

#
    ################################################################################################
    #                                                                                              #
    #                         PLOTAR DADOS ORIGINAIS, SEM INTERPOLAR E COM GAPS                    #
    #                                                                                              #
    ################################################################################################

tempArray, saltArray, densArray = createArrays()

## leitura e armazenamento dos perfis
# variavel para controlar qual perfil esta sendo lido em questao de coluna da matriz de dados
perfil = 0 

for fname in allFiles:
    data = realize_process(fname)

    data['sigma0'] = gsw.sigma0(data['SA'].values, data['CT'].values)

    for i in range(0, len(data)):
        tempArray[perfil, i] = data.temp.values[i]
        saltArray[perfil, i] = data.salt.values[i]
        densArray[perfil, i] = data.sigma0.values[i]

    perfil += 1

## plotar as seções
depth = np.arange(0, int(lastDepth.max()))
xplot, yplot = np.meshgrid(dist, depth)
xplot, yplot = xplot.T, yplot.T

props_original = {
    'tempArray': tempArray.T,
    'saltArray': saltArray.T,
    'densArray': densArray.T
}

createVerticalSection(xm, hm, xplot, yplot, props_original, lista_estacoes, 
                      name='original', # para diferencia cada seção quando salvar a figura
                      suptitle='Vertical Section - Deproas I2001', # titulo do plot
                      interpolate=True, topColor='k', savefig=False, #interpolar (sim), cor da topografia (preto), savefig(nao salvar)
                      showfig=True, alpha=1, comparar=True)

    ################################################################################################
    #                                                                                              #
    #                    PLOTAR DADOS INTERPOLADOS, com rotinas FF, com gaps                       #
    #                                                                                              #
    ################################################################################################
print("interpolando seção vertical para preencher buracos")

distancia = np.asarray(dist)
lastDepth = np.asarray(lastDepth)

#contourf levels
levels = np.arange(0, 25, 0.1)


M = tempArray.T
dist = np.arange(distancia[0], distancia[-1], 10) #np.asarray(distancia)
depth = np.arange(0, int(lastDepth.max()),1)

ntemp = extrap_sec(tempArray.T, dist, depth)
nsalt = extrap_sec(saltArray.T, dist, depth)
ndens = extrap_sec(densArray.T, dist, depth)

props_interpolated = {
    'tempArray': ntemp,
    'saltArray': nsalt,
    'densArray': ndens
}

createVerticalSection(xm, hm, xplot, yplot, props_interpolated, lista_estacoes, 
                      name='nearest', # para diferencia cada seção quando salvar a figura
                      suptitle='Vertical Section - Deproas I2001', # titulo do plot
                      interpolate=True, topColor='k', savefig=False, #interpolar (sim), cor da topografia (preto), savefig(nao salvar)
                      showfig=True, alpha=1, comparar=True)


    ################################################################################################
    #                                                                                              #
    #               PLOTAR DADOS INTERPOLADOS, com rotinas FF sem gaps(martim)                     #
    #                   preenchendo gaps com o gradiente do nível conhecido                        #
    #                     somado ao valor do nível mais próximo conhecido                          #                                                                                 #
    ################################################################################################

ntempArray = tempArray.copy()
nsaltArray = saltArray.copy()
ndensArray = densArray.copy()

def completeColumn(prop):

    """ """
    # localizar os indices de nao nans
    ind = np.where(~np.isnan(prop))
    # localizar o ultimo valor registrado
    lastValue = prop[ind][-1]
    # reproduzir ultimo valor em todos os nan's
    ind = np.where(np.isnan(prop))
    prop[ind] = lastValue

    return prop

for i in np.arange(0, 6): # a cada iteração avanço uma estação/coluna da matriz
    # propriedade sendo lida no momento
    ntempArray[i,:] = completeColumn(tempArray[i,:])
    nsaltArray[i,:] = completeColumn(saltArray[i,:])
    ndensArray[i,:] = completeColumn(densArray[i,:])


ntemp_ff = extrap_sec(ntempArray.T, dist, depth)
nsalt_ff = extrap_sec(nsaltArray.T, dist, depth)
ndens_ff = extrap_sec(ndensArray.T, dist, depth)

props_interpolated_v2 = {
    'tempArray': ntempArray.T,
    'saltArray': nsaltArray.T,
    'densArray': ndensArray.T
}

createVerticalSection(xm, hm, xplot, yplot, props_interpolated_v2, lista_estacoes, 
                      name='gradient', # para diferencia cada seção quando salvar a figura
                      suptitle='Vertical Section - Deproas I2001', # titulo do plot
                      interpolate=True, topColor='k', savefig=False, #interpolar (sim), cor da topografia (preto), savefig(nao salvar)
                      showfig=True, alpha=1, comparar=True)

    ################################################################################################
    #                                                                                              #
    #           PLOTAR DADOS INTERPOLADOS, com rotinas scint.griddata, sem gaps(martim)            #
    #                                                                                              #
    ################################################################################################


#### testar por griddata
# kx, ky - knowed X and Y
# xi, yi - interpolated X and Y
# values - property to grid
# method - nearest

# important: kx, ky and values must have same dimension (2D)
kx = np.asarray(dist)
ky = np.arange(0, int(lastDepth.max()))

kx, ky = xplot.flatten(), yplot.flatten()

XI, YI = np.meshgrid( np.arange(0, 174, 10), np.arange(0, int(lastDepth.max())))
XI, YI = XI.T, YI.T

ntemp_v2 = scint.griddata((kx, ky), ntempArray.flatten(), (XI, YI), method='nearest')

    ################################################################################################
    #                                                                                              #
    #        PLOTAR DADOS INTERPOLADOS, preenchendo gaps com o gradiente do nível conhecido        #
    #       somado ao valor do nível mais próximo conhecido: Tb = Ta + gradT                       #
    #                                                                                              #
    ################################################################################################

"""
Algoritmo
    por estação e por nível vertical (profundidade), testar se o valor é NaN, caso seja devo calcular a
    variação da propriedade média no nível analisado e somar o valor mais proximo conhecido a essa variação

    mesmo esquema para se utilizar no calculo do método dinâmico, caso o nível de referência do movimento
    nulo seja topografia

"""
