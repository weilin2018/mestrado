#-*-coding:utf-8-*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import cmocean as cmo
import glob
import gsw
import seawater as sw
import scipy.interpolate as scint
import matplotlib
import matplotlib.patches as mpatches

matplotlib.style.use('ggplot')

from OceanLab import CTD


# importar os dados
# tratar os dados: loopedit, despike, binning and hanning filter 
BASE_DIR = '/home/tparente/danilo/mestrado/disciplinas/sinoticos/lista1/'


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

# obter o ID da estacao a partir do nome do arquivo sendo processado
def get_profileID(fname):
    """ 
        return only the name file (without full path and extension)
    """

    return fname.split('/')[-1][3:-4]

# obter a coordenada da estacao a partir do ID, pesquisando no arquivo posicoes.txt
def get_latlon(fname):

    posicoesFiles = BASE_DIR + 'Deproas2/posicoes.txt'

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

# resgatar as coordenadas de cada estação
def coordenadasEstacoes(lista_estacoes):
    """ resgatar as coordenadas das estacoes no arquivo posicoes.txt """

    coords_estacoes = []

    with open(BASE_DIR + 'Deproas2/posicoes.txt') as f:
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

# criar matrizes para armazenar os dados
def createArrays(ndim=2270, mdim=6):

    ### parte crucial: juntandos os dados em uma unica MATRIZ
    # criar matriz para receber dados
    ### importante o len(data) deve corresponder ao máximo de dados do maior perfil
    v = [] # dimensao
    for i in range(0, ndim):
        v.append(np.nan)

    v = np.vstack(v)
    # matriz
    m =[]
    for i in range(0, mdim):
        m.append(v)
    # matriz (len(allFiles), len(data))
    tempArray = np.squeeze((m), axis=2)
    saltArray = np.squeeze((m), axis=2)
    densArray = np.squeeze((m), axis=2)

    return tempArray, saltArray, densArray

# interpola uma propriedade para novas dimensoes, usando interp1d
def interpolate2newdimension(prop, dist, depth, new_dist, new_depth):
	""" """
	new_prop = createArrays(ndim=len(new_depth), mdim=len(new_dist))[0]

	for z in depth:
		f = scint.interp1d(dist, prop[:,z], kind='linear')
		new_prop[:,z] = f(new_dist)

	return new_prop

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
        dic = pickle.load(open(BASE_DIR + "pickleData/parametricIsopicnals.pkl", 'r'))
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
        pickle.dump(dic, open(BASE_DIR + "pickleData/parametricIsopicnals.pkl", "w"))

    return si, ti, sigma0

# gera imagem do diagrama TS rainbow com todos os dados 
def diagramaTS_arcoiris(allFiles):
	""" 
	recebe lista de arquivos pickle com os dados a serem plotados
	"""

	fig, ax = plt.subplots()

	# isopicnais que serão utilizadas como linhas paramétricas
	si, ti, sigma0 = parametricIsopicnal()



	manual_locations= [(33.9, 26.5), (35.3, 26.5), (36.7, 26.5), 
                        (37.2, 24.4), (37.3, 21.3), (40.0, 20.0),
                        (40.0, 20.0), (40.0, 20.0), (40.0, 20.0),
                        (37.6, 14.5), (37.5, 8.80), (37.7, 2.11)]

	massas_nomes = ['TW', 'SACW', 'AAIW', 'UCDW', 'NADW']

	levels = np.arange(22, 31, 1)
	CS = ax.contour(si, ti, sigma0, levels, linestyles='dashed', colors='k', alpha=.5)
	ax.clabel(CS,  fontsize=12, inline=1, inline_spacing=1, fmt='%1.0f', manual=manual_locations)

	TW_patch = mpatches.Patch(color='r', label='TW', alpha=.5)
	SACW_patch = mpatches.Patch(color='darkorange', label='SACW', alpha=.5)
	AAIW_patch = mpatches.Patch(color='y', label='AAIW', alpha=.5)
	UCDW_patch = mpatches.Patch(color='g', label='UCDW', alpha=.5)
	NADW_patch = mpatches.Patch(color='b', label='NADW', alpha=.5)

	leg_massas = plt.legend(handles=[TW_patch, SACW_patch, AAIW_patch, UCDW_patch, NADW_patch], loc=2)

    ################################################################################################
    #                                                                                              #
    #   criar informacoes das massas d'agua (nome, localizacao da isopicnal, densidade, cor)       #
    #   calcula sigma0 para ter um grid de sigma0 para plotar                                      #
    #   plota somente as linhas de massas_dagua, colorindo os intervalos com massas_cores          #
    #                                                                                              #
    ################################################################################################

	# preparando o background com as cores de cada massa dagua
	massas_nomes = ['TW', 'SACW', 'AAIW', 'NADW', 'AAIW']
	massas_locat = [(37.5, 25.8), (37.4, 20.), (37.5, 17.4), (37.4, 17.0), (37.3, 15.6)]
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
	ax.set_ylabel(u'Potential Temperature [$^oC$]')

	cores_estacoes = [ 'black', 'maroon', 'sienna', 'darkgreen', 'darkblue', 'darkmagenta' ]

	# importar salinidade e temperatura da climatologia WOA
	woa = pickle.load(open(BASE_DIR + 'pickleData/perfil_woa_prox7048.pkl', 'r'))

	sclim, tclim = woa['Salt'], woa['Temp']

	ax.scatter(sclim, tclim, c='gray', label='WOA', linewidth=0.)

	for fname,c, IDstation in zip(allFiles, cores_estacoes, lista_estacoes):
		data = pickle.load(open(fname, 'r'))

		# ax.plot(data['salt'], data['temp'], linewidth=4., c='k', label=IDstation)
		ax.scatter(data['salt'], data['CT'], c=c, label=IDstation, linewidth=0.)

	plt.legend(loc=4, scatterpoints=1)

	plt.gca().add_artist(leg_massas)

	plt.title('Diagrama TS')

	plt.show()

#
    ################################################################################################
    #                                                                                              #
    #                      PARTE I - TRATAR E PLOTAR OS PERFIS VERTICAIS E TS                      #
    #                                                                                              #
    ################################################################################################

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
    axTemp.set_xlim([original['temp'].min(), original['temp'].max()])
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
    axSalt.set_xlim([original['salt'].min(), original['salt'].max()])
    axSalt.set_xlabel('Salinity [psu]')
    axSalt.xaxis.set_label_position('top')
    axSalt.xaxis.set_ticks_position('top')
    axSalt.legend(['Original', 'Filtered'], loc='best')

    # global title
    plt.suptitle("Station " + get_profileID(fname), size=12)

    # savefig
    if savefig:
        plt.savefig(BASE_DIR + "outputs/perfis/Perfil_" + get_profileID(fname) +'.png')

    plt.close("ALL")

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

        # print("Gerando Diagrama TS espalhado e arco íris")
        # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(18,10))
        # # criar o diagrama TS espalhado
        # diagramaTS_espalhado(ax1, data, fname)
        # # criar o diagrama TS arco iris
        # diagramaTS_arcoiris(ax2, data, fname)

        # plt.suptitle("Station " + get_profileID(fname), fontsize=20)

        # savefig = False
        # if savefig:
        #     plt.savefig(SAVE_DIR + TSDiagram_name + get_profileID(fname) +'.png', dpi=150)

        # plt.show()

    return data

def save2pickle(data, fname):
	""" store pd.DataFrame into pickles files """

	outputFile = fname.split('/')[-1][:-4] + '.pickle'

	pickle.dump(data, open(BASE_DIR + 'pickleData/' + outputFile, 'w'))

def readFrompickle(fname):

	""" read a pickle file with dataframe information and return """

	inputFile = fname.split('/')[-1][:-7] + '.pickle'

	return pickle.load(open(BASE_DIR + 'pickleData/' + inputFile, 'r'))

###### gerar perfis verticais e salvar os dados tratados em pickles

# list all files you want to plot
allFiles = glob.glob(BASE_DIR + 'Deproas2/d2_*')
allFiles.sort()

lista_estacoes = []
pres = []

for fname in allFiles:
	data = realize_process(fname, gerarImagens=False)
	pres.append(np.asarray(data['press'].values))

	# store data into pickle file 
	save2pickle(data, fname)

	lista_estacoes.append(get_profileID(fname))

####################################################
# codigos intermediários

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

xm, hm = dist, lastDepth

####################################################

    ################################################################################################
    #                                                                                              #
    #                      				PARTE III - SEÇÕES VERTICAIS                               #
    #                                                                                              #
    ################################################################################################

# listar arquivos pickles com dados ja tratados
allFiles = glob.glob('/home/tparente/danilo/mestrado/disciplinas/sinoticos/lista1/pickleData/data/*.pickle')
allFiles.sort()

# gerar diagrama TS rainbow:
diagramaTS_arcoiris(allFiles)

# gerar matriz para os dados [6, 2270]
temp, salt, dens = createArrays()

# importar os dados e inserir nas matrizes
perfil = 0

for fname in allFiles:
	#data = readFrompickle(fname) # temp, salt, press, SA, CT
	data = pickle.load(open(fname, 'r'))

	# calcular densidade sigma0
	data['sigma0'] = gsw.sigma0(data['SA'].values, data['CT'].values)

	# armazenar dados ate a profundidade em que ha informacoes
	for i in range(0, len(data)):
		temp[perfil, i] = data.temp.values[i]
		salt[perfil, i] = data.salt.values[i]
		dens[perfil, i] = data.sigma0.values[i]

	perfil += 1

# Salvar informacoes em um dicionário para plotar os dados originais posteriormente
depth = np.arange(0, int(lastDepth.max()))
xplot, yplot = np.meshgrid(dist, depth)
xplot, yplot = xplot.T, yplot.T

origData = {
	"depth": depth,
	"xplot": xplot,
	"yplot": yplot,
	"tempO": temp,
	"saltO": salt,
	"densO": dens
}

save2pickle(origData, "secao_vertical_original.pkl")

# interpolar a matriz inteira, aumentando a quantidade do espaçamento horizontal (de 6 para 31, por exemplo)
ndist = np.arange(0, 170, 5) # ndim: 34
ndept = np.arange(0, 2270, 1) # mdim: 2270

# criar nova matriz para armazenar dados
new_temp = interpolate2newdimension(temp, dist, depth, ndist, ndept)
new_salt = interpolate2newdimension(salt, dist, depth, ndist, ndept)
new_dens = interpolate2newdimension(dens, dist, depth, ndist, ndept)

# plotar para conferir
xplot, yplot = np.meshgrid(ndist, ndept)
xplot, yplot = xplot.T, yplot.T

# salvar informacoes em um dicionario para plotar os dados interpolados posteriormente
interpData = {
	"ndept": ndept,
	"ndist": ndist,
	"xplot": xplot,
	"yplot": yplot,
	"tempI": new_temp,
	"saltI": new_salt,
	"densI": new_dens
}

save2pickle(interpData, "secao_vertical_interpolada.pkl")

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(8,11))

tempPlot = ax1.contourf(xplot, -yplot, new_temp, cmap=cmo.cm.thermal)
ax1.fill_between(xm, -1000, -hm, color='k')
ax1.set_title(u'Potential Temperature [$^o$C]')
ax1.set_ylim([-800, 0])
ax1.set_xlim([0, 165])
ax1.set_ylabel('Pressure [dbar]')

cbar = plt.colorbar(tempPlot, ax=ax1, ticks=np.arange(0., 28., 3))
cbar.set_label(u'Temperature [$^o$C]')

saltPlot = ax2.contourf(xplot, -yplot, new_salt, cmap=cmo.cm.haline)
ax2.fill_between(xm, -1000, -hm, color='k')
ax2.set_ylim([-800, 0])
ax2.set_xlim([0, 165])
ax2.set_title('Salinity')
ax2.set_ylabel('Pressure [dbar]')

cbar = plt.colorbar(saltPlot, ax=ax2, ticks=np.arange(34., 39., 0.5))
cbar.set_label('Salinity [psu]')

densPlot = ax3.contourf(xplot, -yplot, new_dens, cmap=cmo.cm.dense)
# l = np.asarray([26.9])
# cs = ax3.contour(xx, -yy, DI, l, colors='k')
# plt.clabel(cs, fontsize=9, inline=1, fmt='%1.2f')

ax3.fill_between(xm, -1000, -hm, color='k')
ax3.set_ylim([-800, 0])
ax3.set_xlim([0, 165])
ax3.set_title('sigma_0')
cbar = plt.colorbar(densPlot, ax=ax3, ticks=np.arange(24., 28., .4))
cbar.set_label(u'$\sigma_{\Theta}$ [kg m$^{-3}$]')

ax3.set_xlabel('Distance along transect [km]')
ax3.set_ylabel('Pressure [dbar]')

plt.suptitle("Interpolation from [6, 2270] to [34, 2270]", size=20)

plt.savefig(BASE_DIR + 'outputs/secao_vertical/secao_vertical_interpolacaohorizontal.png', dpi=150)
plt.show()


#############################################
#				GRIDDATA 					#
#############################################

# pontos conhecidos
x = np.asarray(dist)
y = np.asarray(depth)

xx, yy = np.meshgrid(x,y)
xx, yy = xx.T, yy.T

# localizar indices de não nan
ind = np.where(~np.isnan(temp))

x1 = xx[ind]
y1 = yy[ind]

T = temp[ind]
S = salt[ind]
D = dens[ind]

TI = scint.griddata( (x1, y1), T.ravel(), (xx, yy), method='nearest')
SI = scint.griddata( (x1, y1), S.ravel(), (xx, yy), method='nearest')
DI = scint.griddata( (x1, y1), D.ravel(), (xx, yy), method='nearest')

# salvar informacoes em um dicionario para plotar os dados interpolados posteriormente
griddataData = {
	"depth": depth,
	"dist": dist,
	"xplot": xx,
	"yplot": yy,
	"tempG": TI,
	"saltG": SI,
	"densG": DI
}

save2pickle(griddataData, "secao_vertical_griddata.pkl")

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(8,11))

tempPlot = ax1.contourf(xx, -yy, TI, cmap=cmo.cm.thermal)
ax1.fill_between(xm, -1000, -hm, color='k')
ax1.set_title(u'Potential Temperature [$^o$C]')
ax1.set_ylim([-800, 0])
ax1.set_xlim([0, 165])
ax1.set_ylabel('Pressure [dbar]')

cbar = plt.colorbar(tempPlot, ax=ax1, ticks=np.arange(0., 28., 3))
cbar.set_label(u'Temperature [$^o$C]')

saltPlot = ax2.contourf(xx, -yy, SI, cmap=cmo.cm.haline)
ax2.fill_between(xm, -1000, -hm, color='k')
ax2.set_ylim([-800, 0])
ax2.set_xlim([0, 165])
ax2.set_title('Salinity')
ax2.set_ylabel('Pressure [dbar]')

cbar = plt.colorbar(saltPlot, ax=ax2, ticks=np.arange(34., 39., 0.5))
cbar.set_label('Salinity [psu]')

densPlot = ax3.contourf(xx, -yy, DI, cmap=cmo.cm.dense)
l = np.asarray([26.9])
cs = ax3.contour(xx, -yy, DI, l, colors='k')
plt.clabel(cs, fontsize=9, inline=1, fmt='%1.2f')

ax3.fill_between(xm, -1000, -hm, color='k')


ax3.set_ylim([-800, 0])
ax3.set_xlim([0, 165])
ax3.set_title('sigma_0')
cbar = plt.colorbar(densPlot, ax=ax3, ticks=np.arange(23., 28.4, .4))
cbar.set_label(u'$\sigma_{\Theta}$ [kg m$^{-3}$]')

ax3.set_xlabel('Distance along transect [km]')
ax3.set_ylabel('Pressure [dbar]')

plt.suptitle("Interpolation using griddata [nearest]", size=20)

#plt.savefig(BASE_DIR + 'outputs/secao_vertical/secao_vertical_griddata_nearest.png', dpi=150)
plt.show()


#
    ################################################################################################
    #                                                                                              #
    #                      		   PARTE IV - VELOCIDADE GEOSTRÓFICA                               #
    #                                                                                              #
    ################################################################################################


def createPressArray(densArray, pres):

	for i in range(0, 6):
		densArray[i, :] = pres[-1]

	return densArray.T

# calcular geopotencial baseado no nível de refência
presArray, gaArray, gvelArray = createPressArray(dens, pres).T, dens*np.nan, dens*np.nan

""" RODAR CADA ESTACAO, RODAR CADA NIVEL VERTICAL  """
for perfil in range(0,6): # importante destacar que só vamos até penultima estação
	# z variando de 0 ao nível vertical de 560dbar
	for z in range(0, len(TI[perfil, :])): 
		if ~np.isnan(TI[perfil, z]): # testa se o valor é nan
			# calcular geopotential anomaly
			gaArray[perfil, z] = sw.gpan( SI[perfil, z], TI[perfil, z], pres[perfil, z] )

# calcular velocidade geostrófica usando sw.gvel
gvelArray = sw.gvel(gaArray.T, lats, lons)


fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)

cs1 = ax1.contourf(xx, -yy, gaArray, cmap=cmo.cm.amp)
ax1.fill_between(xm, -2270, -hm, color='k', interpolate=True, alpha=1.)
c = plt.colorbar(cs1, ax=ax1)
c.set_label('Geopotential anomaly')
# ax1.set_title('Geopotential Anomaly')
ax1.set_ylim([-560,0])
ax1.set_xlim([0, 140])


xplot, yplot = np.meshgrid(dist[:5], depth)
xplot, yplot = xplot.T, yplot.T
levels = np.arange(np.nanmin(gvelArray), np.nanmax(gvelArray), 0.001)

cs2 = ax2.contourf(xplot[:,:], -yplot[:,:], gvelArray.T[:,:], levels, cmap='coolwarm')

l = np.asarray([0.])
cs = ax2.contour(xplot, -yplot, gvelArray.T, l, colors='k')
plt.clabel(cs, fontsize=9, inline=1, fmt='%1.0f')

ax2.fill_between(xm, -560, -hm, color='k', interpolate=True, alpha=1.)

cbar = plt.colorbar(cs2, ax=ax2)
cbar.set_label(r'Geostrophic Velocity [m s$^(-1)$]')

ax2.set_ylim([-560,0])
ax2.set_xlim([0, 140])

plt.suptitle('Geopotential Anomaly and Geostrophic Velocity - without RL')



### calcular usando a referencia de 560dbar

### preciso remover o valor da velocidade geostrófica com referencia a 560dbar
# calcular geopotencial baseado no nível de refência
ga560, gvel560 = dens*np.nan, dens*np.nan
for perfil in range(0,6): # importante destacar que só vamos até penultima estação
	# z variando de 0 ao nível vertical de 560dbar
	for z in range(0, 560): 
		if ~np.isnan(TI[perfil, z]): # testa se o valor é nan
			# calcular geopotential anomaly
			ga560[perfil, z] = sw.gpan( SI[perfil, z], TI[perfil, z], pres[perfil, z] )



gvel560 = sw.gvel(ga560.T, lats, lons)

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)

cs1 = ax1.contourf(xx, -yy, ga560, cmap=cmo.cm.amp)
ax1.fill_between(xm, -2270, -hm, color='k', interpolate=True, alpha=1.)
c = plt.colorbar(cs1, ax=ax1)
c.set_label('Geopotential anomaly')
# ax1.set_title('Geopotential Anomaly')
ax1.set_ylim([-560,0])
ax1.set_xlim([0, 140])


xplot, yplot = np.meshgrid(dist[:5], depth)
xplot, yplot = xplot.T, yplot.T
levels = np.arange(np.nanmin(gvel560), np.nanmax(gvel560), 0.001)

cs2 = ax2.contourf(xplot[:,:], -yplot[:,:], gvel560.T[:,:], levels, cmap='coolwarm')

l = np.asarray([0.])
cs = ax2.contour(xplot, -yplot, gvel560.T, l, colors='k')
plt.clabel(cs, fontsize=9, inline=1, fmt='%1.0f')

ax2.fill_between(xm, -560, -hm, color='k', interpolate=True, alpha=1.)

cbar = plt.colorbar(cs2, ax=ax2)
cbar.set_label(r'Geostrophic Velocity [m s$^(-1)$]')

ax2.set_ylim([-560,0])
ax2.set_xlim([0, 140])

plt.suptitle('Geopotential Anomaly and Geostrophic Velocity - with RL=560dbar')

plt.show()



###########3 calcular o geopotencial relativo ao nível de referencia

""" 
em sw.gpan ele calcula relativo a superficie ao fundo (0 a pfundo)

preciso calcular relativo ao meu nivel de referencia, então ao inves de 0 será
560dbar

"""
# recortar dados de 0 a 560dbar
t_cut = TI[:,:559]
s_cut = SI[:,:559]
p_cut = pres[:, :559]

# inverter para ficar de 560 a 0dbar
t_inv = np.fliplr(t_cut)
s_inv = np.fliplr(s_cut)
p_inv = np.fliplr(p_cut)

# calcular gpan
gpan = sw.gpan(s_inv, t_inv, p_inv)

# calcular gvel
gvel = sw.gvel(gpan.T,  lats, lons)

def gpan(s, t, p):
    """
    Geopotential Anomaly calculated as the integral of svan from the
    the sea surface to the bottom. THUS RELATIVE TO SEA SURFACE.
    Adapted method from Pond and Pickard (p76) to calculate gpan relative to
    sea surface whereas Pond and Pickard calculated relative to the deepest
    common depth.  Note that older literature may use units of "dynamic
    decimeter" for above.
    Parameters
    ----------
    s(p) : array_like
           salinity [psu (PSS-78)]
    t(p) : array_like
           temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    Returns
    -------
    gpan : array_like
           geopotential anomaly
           [m :sup:`3` kg :sup:`-1`
           Pa = m :sup:`2` s :sup:`-2` = J kg :sup:`-1`]
    Examples
    --------
    >>> # Data from Unesco Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> s = [[0, 1, 2], [15, 16, 17], [30, 31, 32], [35,35,35]]
    >>> t = [[15]*3]*4
    >>> p = [[0], [250], [500], [1000]]
    >>> sw.gpan(s, t, p)
    array([[   0.        ,    0.        ,    0.        ],
           [  56.35465209,   54.45399428,   52.55961152],
           [  84.67266947,   80.92724333,   77.19028933],
           [ 104.95799186,   99.38799979,   93.82834339]])
    References
    ----------
    .. [1] S. Pond & G.Pickard 2nd Edition 1986 Introductory Dynamical
       Oceanography Pergamon Press Sydney. ISBN 0-08-028728-X
    """

    s, t, p = list(map(np.asanyarray, (s, t, p)))
    s, t, p = np.broadcast_arrays(s, t, p)
    s, t, p = list(map(atleast_2d, (s, t, p)))

    svn = svan(s, t, p)

    # NOTE: Assumes that pressure is the first dimension!
    mean_svan = (svn[1:, ...] + svn[0:-1, ...]) / 2.
    top = svn[0, ...] * p[0, ...] * db2Pascal
    bottom = (mean_svan * np.diff(p, axis=0)) * db2Pascal
    ga = np.concatenate((top[None, ...], bottom), axis=0)
    return np.cumsum(ga, axis=0).squeeze()
