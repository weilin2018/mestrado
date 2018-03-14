"""
Conjunto de funçoes para extrapolar estações correntográficas através
da expansão da série de Taylor.

funções:
    meanBlock: calcula a média ponderada pela quantidade de estações no bloco,
    sendo P uma propriedade escalar e dZ o intervalo vertical/tamanho do bloco

    expansaoSerieTaylor: realiza uma sequência de calculos:
        . calculo do perfil médio da propriedade P, utilizando um intervalo
        de profundidade dZ

        . alisa o perfil médio utilizando pandas.rolling_mean

        . preserva a variância do perfil ao calcular o NRMS e
            recalculando o perfil médio como: P = (1-alpha)Psmoothed,
            sendo alpha o resultado da NRMS

        . calcula a primeira e segunda derivada do perfil médio suavizado (Psmoothed)

    createArrays: cria uma matriz na direção N sendo a quantidade de estações
     e na direção M a profundidade máxima amostrada pela estação mais funda.
     A matriz retornada é preenchida com np.nan

"""

import numpy as np
import matplotlib.pyplot as plt

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


# extrapolar as estações menores que 1200dbar usando expansão da série de fourier
def meanBlock(P,bloco,dZ):
    """
    calcular a média, ponderada pela quantidade de estações realizadas em
    cada bloco, dentro de um bloco para uma propriedade P
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

Prof_deCorte = 1200
stations = 35

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
