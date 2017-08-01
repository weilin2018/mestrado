"""
Funções:


    def exec2() - plotar hodógrafo
    def exec3() - plotar perfil
    def exec6() - perfil de distribuição meridional
        do bombeamento de ekman bentico

    compass() based on https://ocefpaf.github.io/python4oceanographers/blog/2015/02/09/compass/

    How to use:

        type:
            $ exercicio2() <enter>



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

def cart2pol(x, y):
    """Convert from Cartesian to polar coordinates.

    Example
    -------
    >>> theta, radius = pol2cart(x, y)
    """
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)

    return theta, radius

def calcular_UeV(Av=1e-2, ug=0.1, fo=-1e-4):
	he=np.sqrt(Av/abs(fo))
	z=np.linspace(0,6*he,80)
	u=ug*(1-np.exp(-z/he)*np.cos(z/he))
	s=fo/abs(fo)
	v=s*ug*np.exp(-z/he)*np.sin(z/he)
	return u,v,z

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

# exercicio2()

def calcular_ekman(rho,fo,he,z,taux,tauy):
    tau_u = (taux*np.cos((z/he)-np.pi/4) - tauy*np.sin((z/he)-np.pi/4))
    tau_v = (taux*np.sin((z/he)-np.pi/4) + tauy*np.cos((z/he)-np.pi/4))
    ue=(fo/np.abs(fo))*(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_u
    ve=(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_v
    return ue,ve

def exercicio3():

    """ """
    ###############################################
    #   Definição das constantes para o exercicio #
    ###############################################

    theta = -25. # latitude central do fenômeno
    fo = sw.f(theta) # calculo do parâmetro de Coriolis
    s = fo/abs(fo) # número hemisférico
    dz = 0.2 # metros
    uar=10
    Tau=1.225*0.0015*uar**2
    rho=1027
    hE = np.sqrt((2*0.5e-2)/np.abs(fo))
    pE= np.pi*hE

    # domínio vertical: -2piHe < z < 0
    z = np.arange(0,-2*np.pi*hE,-dz)
    JJ = np.arange(z.shape[0])+1

    # definição dos títulos dos experimentos e variações de Av
    experimentos = {
        u'Ex 3.C - Av(z) constante ($5x10^{-2} m^2 s^{-1}$)': np.repeat(0.5e-2, z.shape[0]),
        # u'Ex 3.D - Av variando linearmente ($5.0 + 0.00625z$)': 5.+0.0625*z,
        # u'Ex 3.E - Variação Linear de Av - Madsen(JPO, 1970)': -0.0625*z,
        # u'Ex 3.F - Escala de Decaimento de Av': 5*np.exp(z/30),
        # u'Variação de Av baseado em Yu&O\'Brien (JPO, 1991)': 10
    }

    for key in experimentos.keys():
        Av    = experimentos[key]
        title = key

        # definição da variação de Av
        dAv = np.diff(Av)/np.diff(JJ)

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

        # Cálculo das velocidades, invertendo a matriz A e multiplicando,
        # matricialmente, pelo vetor solução S
        A2 = np.linalg.inv(A)
        Vcomplexo  = np.dot(A2,S) # resultados da velocidade (complexa)

        # extrair componentes da velocidade (u,v) de V (velocidade complexa)
        # velocidade complexa: V = u + iv
        u,v = np.real(Vcomplexo)[::-1], np.imag(Vcomplexo)[::-1]

        # calcular velocidade no modelo clássico de Ekman
        ue,ve = calcular_ekman(rho,fo,hE,z,Tau,0)

        splt = 0 # se quiser recortar o perfil em algum ponto

        if splt != 0:
            U,V,Z,ue,ve,Av = u[:splt], v[:splt], z[:splt], ue[:splt], ve[:splt], Av[:splt]
        else:
            U,V,Z,ue,ve,Av = u[:], v[:], z[:], ue[:], ve[:], Av[:]


        # calculo da velocidade na solução numérica (spd) e analitica (spde)
        spd, spde = np.sqrt(U**2 + V**2), np.sqrt(ue**2 + ve**2)

        # calcular a profundidade efetiva
        # He = calcular_profundidadeEfetiva(spde, spd, z, pE)
        
        plotar_exec3(U,V,Z,ue,ve,title,savefig='')

def plotar_exec3(U,V,Z,ue,ve,title,savefig=''):

    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(1,3)

    ### plotar hodógrafo confrontando solução analítica e numérica
    
    ax1 = fig.add_subplot(gs[0,0], projection='polar') # cobrir toda primeira linha
    ax1 = compass(ax1, ue[::20], ve[::20], color='r')
    ax1 = compass(ax1, U[::20], V[::20], color='k')
    
    # ax1.set_title(u"Hodógrafo", fontsize=15)    

    ### plotar perfil vertical
    ax2 = fig.add_subplot(gs[0,1])
    ax2.plot(ue,Z,'r-.',label=u'U analítico')
    ax2.plot(U,Z,'k--',label=u'U numérico')
    ax2.set_title(u'Perfil Vertical da Componente U', fontsize=13)

    ax2.set_ylabel(u'Profundidade [m]', fontsize=13)
    ax2.set_xlabel(u'Velocidade Zonal [$m s^{-1}$]', fontsize=13)

    plt.legend(loc='best')

    ax3 = fig.add_subplot(gs[0,2])
    ax3.plot(ve,Z,'r-.', label=u'V analítico')
    ax3.plot(V,Z,'k--', label=u'V numérico')    
    ax3.set_title(u'Perfil Vertical da Componente V', fontsize=13)

    ax3.tick_params(axis='y', labelleft='off', )
    ax3.set_xlabel(u'Velocidade Meridional [$m s^{-1}$]', fontsize=13)

    plt.legend(loc='best')

    plt.suptitle(title, fontsize=18)

    # tight layout using GridSpec and inserting suptitle in the box
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95]) 

    if savefig == '':
        plt.show()
    else:
        plt.savefig('../outputs/'+savefig, dpi=250)

def near(dat,val,how_many=1):
    """ to find the closest value to datas 

    return index """
    dif = np.abs(dat-val) # calcula a diferença entre os pontos
    idx = np.argsort(dif) # organiza as diferenças da menor a maior
    return idx # retorna as how_many menores diferenças

# def calcular_profundidadeEfetiva(spde,spd,z,hE):

#     """ 
#         Algoritmo:

#         1) calcular profundidade efetiva na camada de Ekman na solução clássica
#         2) estimar uma razão entre a velocidade na superfície e na profundidade
#             efetiva
#         3) encontrar uma profundidade na velocidade numérica que tenha a mesma razão

#     """
#     # calcular profundidade efetiva da solução clássica
#     pE = np.pi * hE

#     # calcular a razão entre a velocidade na solução clássica na superfície e em pE
#     ind = near(spde,-pE,1)[:1] 
#     razao = spde[0] / spde[ind]

#     # encontrar a profundidade em que a razao é respeitada

#     # 
#         plt.show()
