"""
    Lista 3 de DFGII
        Exercício 2 - solução analítica para a camada de Ekman bêntica
        Exercício 3 - solução analítica e numérica, considerando Av(z), para a camada de 
                        Ekman de superfície

    Elaborado por Danilo A. Silva <nilodna@gmail.com>, com algumas funções
    obtidas em:
        Filipe Fernandes: https://ocefpaf.github.io/python4oceanographers/blog/2015/02/09/compass/
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.gridspec as gridspec
import seawater as sw
import matplotlib
matplotlib.style.use('ggplot')

import sys
reload(sys)
sys.setdefaultencoding('utf8')

""" converter de coordenadas cartesianadas para polares """
def cart2pol(x, y):
    """Convert from Cartesian to polar coordinates.

    Example
    -------
    >>> theta, radius = pol2cart(x, y)
    """
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)

    return theta, radius

""" calcular velocidade na solução analítica da camada de ekman bêntica """
def calcular_UeV(Av=1e-2, ug=0.1, fo=-1e-4):
	he=np.sqrt(Av/abs(fo))
	z=np.linspace(0,6*he,80)
	u=ug*(1-np.exp(-z/he)*np.cos(z/he))
	s=fo/abs(fo)
	v=s*ug*np.exp(-z/he)*np.sin(z/he)
	return u,v,z

""" extraída de Filipe Fernandes, simular a função compass do matlab, gerando um hodógrafo """
def compass(ax, u, v, color='k', arrowprops=None):
    """
    Compass draws a graph that displays the vectors with
    components `u` and `v` as arrows from the origin.

    Examples
    --------
    >>> import numpy as np
    >>> u = [+0, +0.5, -0.50, -0.90]
    >>> v = [+1, +0.5, -0.45, +0.85]
    >>> compass(u, v)
    """

    angles, radii = cart2pol(u, v)

    kw = dict(arrowstyle="->", color=color)
    if arrowprops:
        kw.update(arrowprops)
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),
                 arrowprops=kw) for
     angle, radius in zip(angles, radii)]

    ax.set_ylim(0, np.max(radii))

    return ax

""" rodar o exercício 2 """
def exercicio2(Av=4*1e-2, theta=-np.pi/6, l=int(5e4), y=0):
    """
        plotar hodógrafo do vetor velocidade em função da profundidade

        Av: coeficiente de viscosidade vertical
        ug: velocidade do interior geostrófico
        fo: parâmetro de coriolis
    """
    ue,ve,z = calcular_UeV() # calculo da velocidade na camada de ekman bentica

    #####################
    #  plotar hodógrafo #
    #####################
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2,2)

    ax1 = fig.add_subplot(gs[:,0], projection='polar') # cobrir toda primeira linha
    # se quiser plotar uma linha e nao vetores, descomente as próximas linhas e
    # comente a linha ax1=compass()

    #angles, radii = cart2pol(ue,ve)
    # ax1.plot(angles, radii, 'k', label=u'Hemisfério Sul [s < 0]')
    # ax1.set_ylim(0, 0.13)
    ax1 = compass(ax1, ue, ve)
    ax1.set_title(u"Ex1.b - Hodógrafo da Espiral de Ekman na \n camada de Ekman Bêntica", fontsize=15)
    # ax1.set_title(u'Ex1.b - Hodógrafo do vetor velocidade \n em função da profundidade', fontsize=13)

    ##################
    #  plotar perfil #
    ##################
    fo = 2*7.2921*(1e-5)*np.sin(theta)
    he = np.sqrt(Av/abs(fo))
    ug = np.exp(-y**2/l**2)
    z = np.linspace(0,6*he,80)

    ue,ve,ze = calcular_UeV(Av=Av, ug=ug, fo=fo)

    ax2 = fig.add_subplot(gs[0,1])
    ax2.grid(True)
    ax2.plot(ue,z/he,'k--', label='U')
    ax2.plot(ve,z/he, 'k-.', label='V')

    ax2.axvline(0, color='black', alpha=.5)
    ax2.set_xlabel(r'Velocidade normalizada ($\frac{u}{U_0}$, $\frac{v}{U_0}$)', fontsize=17)
    ax2.set_ylabel(r'Profundidade normalizada ($\frac{z}{he}$)', fontsize=17)

    ax2.set_title(u'Ex1.c - Perfil vertical das correntes', fontsize=15)
    plt.legend(loc='best')

    ##################
    #  plotar we     #
    ##################
    ys = np.arange(-100*1e3,100*1e3,1)
    U  = np.exp(-ys**2/l**2)
    s = fo/abs(fo)
    we = he*s*ys*np.exp(-ys**2/l**2)

    ax3 = fig.add_subplot(gs[1,1])
    ax3.plot(we*1e-5, ys/1e3, 'k', label=u'$we$')
    ax3.plot(U,ys/1e3, 'k--', label=u'$u_g$')
    ax3.set_xlabel(u'Velocidade do Bombeamento de Ekman ($we$) [x$10^{-5} m s^{-1}$] e \n da Corrente Geostrófica ($u_g$) [$ m s^{-1}$]', fontsize=17)
    ax3.set_ylabel(u'Distância Meridional [km]', fontsize=17)
    ax3.set_title(u'Ex1.e - Bombeamento de Ekman no topo da camada de Ekman Bêntica', fontsize=15)
    plt.legend(loc='best')
    ######################
    #  plotar corrente   #
    ######################
    # ax4 = fig.add_subplot(gs[1,1])
    # ax4.plot(U,ys/1e3, 'k')
    # ax4.set_xlabel(u'Velocidade da Corrente Geotrófica [$m s^{-1}$]', fontsize=20)
    # ax4.set_ylabel(u'Distância Meridional [m]', fontsize=20)
    # ax4.set_title(u'Perfil meridional da corrente geostrófica', fontsize=13)

    ax1.text(np.pi * 1.749, 0.101, '(A)')
    ax2.text(1.09104, 0.478974, '(B)')
    ax3.text(5.01755, -78.6419, '(C)')


    plt.suptitle(r'Exercício 1 - Modelo Clássico da Camada de Ekman Bêntica para o Hemisfério Sul ($\theta_0 = -30^o$)', fontsize=20)
    #plt.savefig('../outputs/exercicio2.png', dpi=250)
    plt.show()

""" calcular a velocidade na camada de Ekman segundo solução analítica clássica """
def calcular_ekman(rho,fo,he,z,taux,tauy):
    tau_u = (taux*np.cos((z/he)-np.pi/4) - tauy*np.sin((z/he)-np.pi/4))
    tau_v = (taux*np.sin((z/he)-np.pi/4) + tauy*np.cos((z/he)-np.pi/4))
    ue=(fo/np.abs(fo))*(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_u
    ve=(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_v
    return ue,ve

""" criar matriz A e vetor solução S para resolver problema de álgebra linear na solução numérica de Ekman na superfície """
def create_Amatrix(Av, JJ, dz, z, fo, Tau, rho):
    """
    Criar a matriz A e o vetor solução S, baseados nas equações de (15) a (20) da lista

    ::TODO colocar aqui as equaçoes e esquema da matrizs

    """
    # definição da variação de Av
    dAv = np.diff(Av)/np.diff(JJ)

    # recriar o número hemisférico
    s = fo/abs(fo)

    # Definição das diagonais principais baseado nas equações (16), (17) e (18)
    G1 = -1+((dAv/Av[:-1])*(dz/2))
    G2 = (np.tile(2,z.shape[0]))+1j*((dz**2)*s*np.abs(fo))/Av
    G3 = (-1-((dAv/Av[:-1])*(dz/2)))

    # POG pra dar certo
    if np.isinf(G1[0]) and G1[0] > 0:
        G1[0]=100000000
    elif np.isinf(G1[0]) and G1[0] < 0:
        G1[0]=-100000000

    G3[0]=0
    G2[0]=1

    # Definição do vetor solução S
    S = np.append(np.zeros(z.shape[0]-1), (2*Tau*dz)/(rho*(Av[-1]+Av[-2])))

    # Definição da matriz A
    A = np.diag(G2,0) + np.diag(G1, -1) + np.diag(G3, 1)
    # Adicionando os contornos nas últimas linhas, conforme (15)
    A[0,0],A[-1,-1],A[0,1],A[-1,-2] = 1,1,0,-1

    return A, S

""" calcular a velocidade complexa na camada de Ekman de superfície """
def calcular_Vcomplexa(A,S):
    """ 
    Calculo do vetor das incógnitas \nu (velocidade complexa), usando o 
    problema linear dado por (14)

    calcula-se a matriz inversa de A para que possamos realizar o produto
    vetorial entre o vetor solução S e a matriz inversa, obtendo \hat{\nu}.

    """
    # Cálculo das velocidades, invertendo a matriz A e multiplicando,
    # matricialmente, pelo vetor solução S
    A2 = np.linalg.inv(A)
    Vcomplexo  = np.dot(A2,S) # resultados da velocidade (complexa)

    # extrair componentes da velocidade (u,v) de V (velocidade complexa)
    # velocidade complexa: V = u + iv
    return np.real(Vcomplexo)[::-1], np.imag(Vcomplexo)[::-1]

""" plotar diversos tipos de vetores no mesmo hodógrafo """
def multiple_compass(ax, us, vs, colors, arrowprops=None):
    """
        Criar hodógrafo utilizando compass() extraído de Filipe Fernandes.

        Parameters
        ----------
            ax: eixo de plotagem
            us: os dois tipos de velocidade zonal a serem plotados
            vs: os dois tipos de velocidade meridional a serem plotados
            colors: cores para cada tipo de dado

        Return
        ------
            ax: eixo com dados plotados

    """

    U,V,C    = us[0], vs[0], colors[0]                          # extrair solução numérica
    ue,ve,ce = us[1], vs[1], colors[1]                          # extrair solução analítica

    # plotar primeiro analítico
    angles, radii = cart2pol(ue,ve)                             # converter u,v para polar
    kw = dict(arrowstyle="->", color=ce)                        # keyword com informações das setas
    if arrowprops:
        kw.update(arrowprops)
    
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),         # plotar as setas
                 arrowprops=kw) for
     angle, radius in zip(angles, radii)]

    ax.set_ylim(0, np.max(radii))                               # setar limite do raio do plot

    # plotar numérico
    angles, radii = cart2pol(U,V)
    kw = dict(arrowstyle="->", color=C)
    if arrowprops:
        kw.update(arrowprops)
    [ax.annotate("", xy=(angle, radius), xytext=(0, 0),
                 arrowprops=kw) for
     angle, radius in zip(angles, radii)]

                                                                # plotar vetor do vento
    kw = dict(arrowstyle="->", edgecolor='black', facecolor='darkgreen', lw=2, alpha=0.6)
    ax.annotate("", xy=(0., 0.2), xytext=(0,0), arrowprops=kw)  # características do vetor
    
    ax.text(1.917*np.pi, 0.105, r'$\tau_0 = 0.2Pa$', ha='center', fontsize=13)# expressão do vento

    return ax

""" encontrar o índice nos dados do valor mais próximo ao requisitado """
def near(dat,val,how_many=1):
    """
        find the closest value from datas
    """
    dif = np.abs(dat-val)           # calcula a diferença entre os pontos
    idx = np.argsort(dif)           # organiza as diferenças da menor a maior
    return dat[idx][:how_many][0]      # retorna  'how_many' dados de menor diferença

""" calcular a prof efetiva na solução numérica da camada de ekman de superfície """
def calcular_profundidadeEfetiva(SPDe,SPD,z,hE):

    """ 
        Calcula a profundidade efetiva de uma solução numérica da camada de Ekman, 
        utilizando como referência a solução clássica.

        Parameters
        ----------
            SPDe    - velocidade na solução analítica clássica 
            SPD     - velocidade na solução numérica
            z       - domínio vertical
            hE      - espessura característica, dada por $hE = \sqrt{\frac{2A_v}{|f_o|}}$

        Returns
        -------
            prof_camada - indice de z da profundidade efetiva calculada

        Algoritmo
        ---------

        1) calcular profundidade efetiva na camada de Ekman na solução clássica
        2) estimar uma razão entre a velocidade na superfície e na profundidade
            efetiva
        3) encontrar uma profundidade na velocidade numérica que tenha a mesma razão

    """
    pE = np.pi * hE                                             # profundidade efetiva da solução clássica

    rat=SPDe[0]/SPDe[np.argwhere(z==near(z,-pE,1))][0]          # razão entre velocidade na superfície e 
                                                                # na profundidade efetiva clássica
    vel_teoric=SPD[0]/rat                                       # velocidade teorica encontrada com relação a velocidade
                                                                # de Ekman para cada um dos Av(z)
    
    prof_camada = np.argwhere(SPD==near(SPD,vel_teoric[0],1))   # indice da profundidade da camada efetiva para cada
                                                                # experimento

    He = z[prof_camada][0][0]                                   # Profundidade efetiva calculada com razao entre velocidades
    
    return He                                                   # retornar somente o índice

""" rodar exercício 3 """
def exercicio3():

    """ """
    ###############################################
    #   Definição das constantes para o exercicio #
    ###############################################

    theta = -25.                                # latitude central do fenômeno
    fo = sw.f(theta)                            # calculo do parâmetro de Coriolis
    s = fo/abs(fo)                              # número hemisférico: fo/|fo|
    dz = 0.2                                    # variação de z em metros
    uar=10                                      # velocidade do vento a 10m de altura da superfície
    Tau=1.225*0.0015*uar**2                     # tensão de cisalhamento do vento
    rho=1027                                    # densidade do oceano: homogêneo
    hE = np.sqrt((2*0.5e-2)/np.abs(fo))         # Espessura característica hE = sqrt(2Av/|fo|)
    pE= np.pi*hE                                # profundidade efetiva teórica
    z = np.arange(0,-2*np.pi*hE,-dz)            # domínio vertical (z) dado pelo exercício: -2piHe < z < 0
    JJ = np.arange(z.shape[0])+1                # índices baseado no domínio vertical

    ue,ve = calcular_ekman(rho,fo,hE,z,Tau,0)   # cálculo da velocidade na camada de Ekman teórica
                                                # como ela independe das variações de Av, será constante nos hodógrafos

    fator = 10e-4                               # fator de multiplicação para os valores de Av

    # definição dos tipos de variação de Av(z)
    experimento_real = {
        "Av": np.array([0.05, 0.04, 0.03, 0.02, 0.01]),
        "zr": np.array([10,20,40,50,80]),
    }
    experimentos = {
        u'Ex 2.C - Av(z) constante ($5x10^{-2} m^2 s^{-1}$)': np.repeat(0.5e-2, z.shape[0]),
        u'Ex 2.D - Variação Linear de Av(z) = $5.0 + 0.02z$': (5.+0.02*z)*fator,
        u'Ex 2.E - Variação Linear de Av(z) = $-0.0625z$ - Madsen(JPO, 1970)': (-0.0625*z)*fator,
        r'Ex 2.F - Variação Exponencial de Av(z) = $5e^{\frac{z}{d}}$, para d=30m': (5*np.exp(z/30))*fator,
        # u'Variação de Av(z) baseado em Yu and O\'Brien (JPO, 1991)': True 
    }

    for key in experimentos.keys():                                 # loop para trabalhar com cada tipo de Av(z)
        print("Working with: %s"%key)
        # checar se é hora de utilizar valores reais de Av
        # if experimentos[key] == True:                               # usar o segundo dicionário
        #     Av = experimento_real['Av']                             # extrair valores reais de Av
        #     zr = experimento_real['zr']                             # extrair as profundidades associadas

        #     z = np.arange(0,-80,-0.2)                               # calcular novo domínio vertical
        #     Av = np.interp(z,zr,Av)                                 # interpolar os valores de Av obtidos
        #     Av *= 10e-4                                             # fator 
        #     z *= -1                                                 # domínio vertical negativo
        #     JJ = np.arange(z.shape[0]) +1                           # novo vetor com indices
        #     title = key                                             # título do plot

        # else:
        Av    = experimentos[key]                               # definindo Av(z) baseado no experimento
        title = key                                             # definindo o título para o plot

        A, S = create_Amatrix(Av, JJ, dz, z, fo, Tau, rho)          # criação da matriz tridiagonal A e vetor solução S 
                                                                    # para resolver sistema linear

        u,v = calcular_Vcomplexa(A,S)                               # calculo da velocidade complexa utilizando A e S no sistema linear
                                                                    # retornando somente a componente u (real) e v (imaginária)

        splt = 0                                                    # se desejar limitar os dados a um nível z

        if splt != 0:
            U,V,Z,ue,ve,Av = u[:splt], v[:splt], z[:splt], ue[:splt], ve[:splt], Av[:splt]
        else:
            U,V,Z,ue,ve,Av = u[:], v[:], z[:], ue[:], ve[:], Av[:]

        spd, spde = np.sqrt(U**2 + V**2), np.sqrt(ue**2 + ve**2)    # cálculo da velocidade na solução numérica (spd) e analítica (spde)

        prof_camada = calcular_profundidadeEfetiva(spde, spd, Z, hE)# calculo da profundidade efetiva com o perfil de Av do experimento
        
        plotar_exec3(U,V,Z,ue,ve,title,profEfetiva=prof_camada,savefig=key[:6].replace(' ', '').replace('.', '') + '.png')

""" plotar exercicio 3"""
def plotar_exec3(U,V,Z,ue,ve,title,profEfetiva,savefig=''):
    """
        Plotar hodógrafo e perfil vertical das componentes de velocidade, tanto da 
        solução analítica clássica, como da solução numérica para Av(z).

        Parameters
        ----------
            U,V         - numerical solution from velocity
            ue,ve       - analytial solution from velocity
            Z           - vertical domain or ... depth
            title       - plot title
            profEfetiva - index of effective depth from influence of the wind stress

    """
    ##################################################################
    ###                   configurações do gráfico                 ###
    ##################################################################

    figsize      = (15,10)                          # tamanho da figura
    gradeamento  = [1,3]                            # gradeamento dos subplots
    spaceStep    = 20                               # step para plotar os vetores (evitar imagem densa)

    limite_eixoX = [-0.15, 0.3]                     # xlim para os perfis verticais
    limite_eixoY = [-80.0, 0.0]                     # xlim para os perfis verticais
    he_plot      = profEfetiva                      # linha horizontal para profundidade efetiva, em metros e negativo

    xy_a         = [(270*np.pi)/180, 0.25]          # posição do '(A)'
    xy_b         = [0.2725, -75.0]                  # posição do '(B)'
    xy_c         = [0.2725, -75.0]                  # posição do '(C)'

    qualid_output= 180                              # qualidade da imagem de saída, em dpi
    save_dir     = '../outputs/'                    # diretório para salvar imagens
    
    # início da plotagem
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows=gradeamento[0], ncols=gradeamento[1])

    ##################################################################
    ### plotar hodógrafo confrontando solução analítica e numérica ###
    ##################################################################
    
    us = (U[::spaceStep],ue[::spaceStep])                           # juntando o U e ue com um step
    vs = (V[::spaceStep],ve[::spaceStep])                           # juntando o V e ve com um step

    ax1 = fig.add_subplot(gs[0,0], projection='polar')              
    ax1 = multiple_compass(ax1, us, vs, colors=['k', 'r'])          # plotar os vetores de velocidade
    
    ax1.text(xy_a[0], xy_a[1], '(A)', ha='center')                  # identificação do subplot

    ##################################################################
    ###                   plotar perfil vertical                   ###
    ##################################################################

    ax2 = fig.add_subplot(gs[0,1])
    ax2.plot(ue,Z,'r-.',label=u'U analítico')                       # plotar U analítico
    ax2.plot(U,Z,'k--',label=u'U numérico')                         # plotar U numérico
    ax2.axhline(y=he_plot, color='black', alpha=.5)                 # He numérico - linha
    ax2.text(0.15, he_plot-3.,                                      # He numérico - valor
        r'$H_e = %2.2f$m'%(he_plot*(-1)), ha='center')            

    ax2.set_title(u'Perfil Vertical da Componente U', fontsize=13)  # titutlo subplot

    ax2.set_ylabel(u'Profundidade [m]', fontsize=13)                # label eixo Y
    ax2.set_xlabel(u'Velocidade Zonal [$m s^{-1}$]', fontsize=13)   # label eixo X

    ax2.set_xlim(limite_eixoX)                                      # limite do eixo X
    ax2.set_ylim(limite_eixoY)                                      # limite do eixo Y

    ax2.text(xy_b[0], xy_b[1], '(B)', ha='center')                  # identificação do plot

    plt.legend(loc='best')                                          # geração da legenda

    ax3 = fig.add_subplot(gs[0,2])
    ax3.plot(ve,Z,'r-.', label=u'V analítico')                      # plotar V analítico
    ax3.plot(V,Z,'k--', label=u'V numérico')                        # plotar V numérico
    ax3.axhline(y=he_plot, color='black', alpha=.5)                 # He numérico - linha
    ax3.text(0.15, he_plot-3., 
        r'$H_e = %2.2f$m'%(he_plot*(-1)), ha='center')              # He numérico - valor

    ax3.set_title(u'Perfil Vertical da Componente V', fontsize=13)  # titulo do subplot

    ax3.tick_params(axis='y', labelleft='off', )                    # esconder ticks do eixo Y
    ax3.set_xlabel(u'Velocidade Meridional [$m s^{-1}$]', fontsize=13) # label eixo x

    ax3.set_xlim(limite_eixoX)                                      # limite do eixo X
    ax3.set_ylim(limite_eixoY)                                      # limite do eixo Y
    ax3.text(xy_c[0], xy_c[1], '(C)', ha='center')                  # identificação do subplot

    plt.legend(loc='best')                                          # geração da legenda

    plt.suptitle(title, fontsize=18)                                # super título

    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])                   # tight_subplot incorporando o
                                                                    # suptitle da imagem

    if savefig == '':
        plt.show()
    else:
        plt.savefig(save_dir+savefig, dpi=qualid_output)

    plt.close("ALL")

# exercicio2()

exercicio3()


import glob
import os

lstFile = glob.glob('../outputs/Ex**')
for fname in lstFile:
    os.system('convert -trim %s %s' % (fname, fname))
