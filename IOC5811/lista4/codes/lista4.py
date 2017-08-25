"""

    Lista 4 de DFGII
        Exercício 2 - modelo de Stommel na presença de topografia

    Elaborado por Danilo A. Silva <nilodna@gmail.com>, traduzido do código
    desenvolvido por Ilson da Silveira em Matlab


    Este programa calcula numericamente a solução de Stommel para
    campos de rotacional do vento e de topografia de fundo teóricos
    especificados.

    A resolução se dá por esquemas iteraivos em forma não-dimensional
    (que é melhor descrita no notebook que acompanha este código),
    para uma bacia oceânica quadrada de dimensão de 100km.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.gridspec as gridspec
import seawater as sw
import cmocean as cmo
# import matplotlib
# matplotlib.style.use('ggplot')

import sys
reload(sys)
sys.setdefaultencoding('utf8')

import os
os.system('clear')

def TODO_create_gridSubplot(rows, cols, outerSpace, figsize):
    """
    Function based on: https://stackoverflow.com/questions/34933905/matplotlib-adding-subplots-to-a-subplot

    rows, cols = quantidade de subplots [int]
    figsize = tamanho da figura [tupla]
    outerSpace = espaçamento em width e height entre os subplots [tupla]



    """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=figsize)
    outer = gridspec.GridSpec(rows, cols, wspace=outerSpace[0], hspace=outerSpace[1])

    return fig, outer

BASE_DIR = '/home/tparente/danilo/mestrado/disciplinas/mestrado/IOC5811/lista4/'
SAVE_DIR = '/home/tparente/danilo/mestrado/disciplinas/mestrado/IOC5811/lista4/outputs/'

    ###############################################
    #   Definição das constantes para o exercicio #
    ###############################################

tol = 1e-3                      # tolerancia do erro medio durante iteração

# parâmetros do modelo
L  = 1e6                        # largura da bacia
H  = 5e3                        # profundidade
f0 = 7.3e-5                     # parâmetro de coriolis
beta = 2e-11                    #
Av = 1e-2                       # Av = 1e-2 [m^2/s^2]
hE = np.sqrt((2*Av)/f0)         # espessura da camada de Ekman teórica
r = (f0*hE)/H                   # parâmetro de fricção
tau0 = 1e-4                     # tensão de cisalhamento
y0 = 0.5                        # latitude central

# normalizações
r  = r/(beta*L)                 # parâmetro de fricção normalizado
ff = f0/(beta*L)                # parâmetro de coriolios normalizado

# construção da grade da bacia oceânica quadrada
nmax = 101
dx   = 0.01                     # deltax = passo de espaço

x = np.arange(0,nmax) * dx      # construção da dimensão X
y = x.T                         # construção da dimensão Y
[xg,yg] = np.meshgrid(x,y)        # construindo a grade

# xg=np.flipud(xg)
# yg=np.flipud(yg)

yL = max(y)

# criando o campo de rotacional do vento não-dimensional
curl = -(np.pi/yL) * np.sin((np.pi * yg)/yL) # giro único - subtropical

def create_VPambiente(y0, yg, fac, dx):
    """ calcula a vorticidade potencial ambiente """
    B = fac*(yg-y0)+ff*hb               #função VP ambiente
    [Bx,By] = np.gradient(B,dx,dx)      # calculo do gradiente de VP ambiente
    By = (-1)*By

    return B, Bx, By

def stommelSolution(Bx, By, curl, r, nmax, tol):
    """
    to calculate stommel's solution for a subtropical gyre
    """

    psi = np.zeros(xg.shape)                    # create field of psi

    # start the iteration process
    for i in range(2000):

        psi_old = psi                           # save the old psi to calculate a new one

        # define and calculate the 'av' parameter
        av=psi[0:nmax-2,1:nmax-1] + psi[2:nmax,1:nmax-1] +psi[1:nmax-1,0:nmax-2] + psi[1:nmax-1,2:nmax]
        # define and calculate the 'F' parameter
        F=1/r*(curl[1:nmax-1,1:nmax-1]*dx*dx - psi[1:nmax-1,2:nmax]*(By[1:nmax-1,1:nmax-1])*dx +psi[2:nmax,1:nmax-1]*(Bx[1:nmax-1,1:nmax-1])*dx)
        # define and calculate psi[i,j]
        psi[1:nmax-1,1:nmax-1]=(av-F)/(4+(dx/r)*(By[1:nmax-1,1:nmax-1]-Bx[1:nmax-1,1:nmax-1]))
        # check the size of the error
        crit=(np.abs(psi-psi_old));
        crit=crit.max()
        # if the maximum error is <= the tolerance value, contiune the iteration
        if crit <= tol:
            continue
        # else stop the iteration
        else:
            break

    # calculate the relative vorticity
    zeta=(psi[0:nmax-2,1:nmax-1]+psi[2:nmax,1:nmax-1]+psi[1:nmax-1,0:nmax-2]+psi[1:nmax-1,2:nmax]-4*psi[1:nmax-1,1:nmax-1])/dx/dx

    return psi, zeta

def plotar_campoVPAmb(ax, fig, B):
    """
        Plotar campo de vorticidade potencial ambiente

        B = função da VP ambiente
        ax = eixo do subplot
    """

    b1=0.1*np.floor(B.min()*10);
    b2=0.1*np.ceil(B.max()*10);

    lb=np.arange(b1,b2+0.05,0.05)

    # testar se lb é maior que 1 para plotar os contour_levels
    if lb.shape[0] > 1:
        cs = ax.contourf(x,y,np.flipud(B),levels=lb, cmap=cmo.cm.curl)
    else:
        cs = ax.contourf(x,y,np.flipud(B), cmap=cmo.cm.curl)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.35)

    cbar = plt.colorbar(cs,cax=cax)
    cbar.set_label(u'$q_a$ ', fontsize=20)

    # ax.axis('square')
    ax.set(adjustable='box-forced', aspect='equal')
    ax.set_xlabel(u'x [$10^6$ m]')
    ax.set_ylabel(u'y [$10^6$ m]')
    ax.set_title(u'Vorticidade Potencial Ambiente [$q_a$]')

    return ax,fig

def plotar_topografia2D(ax,hb):
    """
    plotar topografia dinâmica
    """

    lt=np.arange(hb.min(),hb.max()+0.005,0.005)
    hbPlot = np.flipud(hb)

    cs = ax.contourf(x,y,hbPlot,levels=lt,cmap=plt.get_cmap('coolwarm'))

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.35)

    cbar = plt.colorbar(cs,cax=cax)
    cbar.set_label('Profundidade [m]', fontsize=17)

    ax.set_xlabel(u'x [$10^6$ m]')
    ax.set_ylabel(u'y [$10^6$ m]')

    ax.set(adjustable='box-forced', aspect='equal')

    return ax
#
# def plotar_topografia3D(ax,hb):
#     """
#     plotar topografia 3D
#     """
#
#     xx,yy = np.meshgrid(x,y)
#
#

def plotar_solucaoStommel(ax, psi):
    """ """
    # determinar níveis de plotagem
    p1=0.1*np.floor(psi.min()*10)
    p2=0.1*np.ceil(psi.max()*10)
    lp = np.arange(p1, p2+0.3, 0.1)
    lp2 = np.arange(p1,p2, .4)
    psiPlot = np.flipud(psi)

    cs = ax.contourf(x,y,psiPlot,levels=lp,cmap=cmo.cm.amp)
    cb = ax.contour(x,y,psiPlot,levels=lp2, colors='k')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.35)

    cbar = plt.colorbar(cs,cax=cax)
    cbar.set_label(u'$\psi$ ', fontsize=20)

    # plt.axis('square')
    ax.set(adjustable='box-forced', aspect='equal')
    cbar.set_label(u'$\psi$', fontsize=20)

    ax.set_xlabel(u'x [$10^{6}$ m]')
    ax.set_ylabel(u'y [$10^{6}$ m]')

    return ax

def fundoPlano_planoBeta():
    """
    Fundo plano:
    plano Beta:  fac =1
    """
    # campo de topografia de fundo não-dimensional para H = 5e3m
    hb = np.zeros(xg.shape)
    # campo de gradiente de vorticidade ambiente
    fac = 1.0
    B=fac*(yg-y0)+ff*hb                     #função VP ambiente
    [By,Bx]=np.gradient(B,dx,dx)            # calcular dB/dx,dB/dy
    Bx = (-1) * Bx                          # adequar a equação discretizada

    # calcular a solução de stommel para este cenário
    psi, zeta = stommelSolution(Bx, By, curl, r, nmax, tol)

    # plotar
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(24,7))

    axes[0],fig = plotar_campoVPAmb(axes[0],fig, B)

    axes[1] = plotar_solucaoStommel(axes[1], psi)
    axes[1].set_title(u'Solução de Stommel')

    plt.suptitle(r'Ex2-C: Solução de Stommel para fundo plano e plano $\beta$', fontsize=20)

    plt.show()

def fundoPlano_planoF():
    """
    Fundo plano:
    plano f: fac=0
    """

    # campo de topografia de fundo não-dimensional para H = 5e3m
    hb = np.zeros(xg.shape)
    # campo de gradiente de vorticidade ambiente
    fac = .0
    B=fac*(yg-y0)+ff*hb                     #função VP ambiente
    [By,Bx]=np.gradient(B,dx,dx)            # calcular dB/dx,dB/dy
    Bx = (-1) * Bx                          # adequar a equação discretizada

    # calcular a solução de stommel para este cenário
    psi, zeta = stommelSolution(Bx, By, curl, r, nmax, tol)

    # plotar
    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(24,7))

    axes[0].set_title(u'VP Ambiente')
    axes[0],fig = plotar_campoVPAmb(axes[0],fig, B)

    axes[1].set_title(u'Solução de Stommel')
    axes[1] = plotar_solucaoStommel(axes[1], psi)

    plt.suptitle(r'Ex2-D: Solução de Stommel para fundo plano e plano $f$', fontsize=20)

    plt.show()

def fundoLinearMeridional_planoF():
    """
    topografia linear na direção y: b = b(y)
    beta_topografico = beta
    """

    # campo de topografia de fundo não-dimensional para H = 5e3m
    hb = -0.275 * yg
    # campo de gradiente de vorticidade ambiente
    fac = .0
    B=fac*(yg-y0)+ff*hb                     #função VP ambiente
    [By,Bx]=np.gradient(B,dx,dx)            # calcular dB/dx,dB/dy
    Bx = (-1) * Bx                          # adequar a equação discretizada

    # calcular a solução de stommel para este cenário
    psi, zeta = stommelSolution(Bx, By, curl, r, nmax, tol)

    # plotar
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(20,6))
    plt.subplots_adjust(wspace = .5)

    axes[0].set_title(u'Topografia [2D]')
    axes[0] = plotar_topografia2D(axes[0], hb)

    axes[1].set_title(u'VP Ambiente')
    axes[1],fig = plotar_campoVPAmb(axes[1],fig, B)

    axes[2].set_title(u'Solução de Stommel')
    axes[2] = plotar_solucaoStommel(axes[2], psi)

    plt.suptitle(r'Ex2-E: Solução de Stommel para fundo com variação linear meridional e $\beta_T$ = $\beta$', fontsize=20)

    plt.show()

def fundoLinearZonal_planoBeta():
    """
    fundo linear em x: b = b(x) = alpha x = beta/4 x

    """

    # campo de topografia de fundo não-dimensional para H = 5e3m
    hb = -3 * 0.275 * xg
    # campo de gradiente de vorticidade ambiente
    fac = 1.0
    B=fac*(yg-y0)+ff*hb                     #função VP ambiente
    [By,Bx]=np.gradient(B,dx,dx)            # calcular dB/dx,dB/dy
    # Bx = (-1) * Bx                          # adequar a equação discretizada

    # calcular a solução de stommel para este cenário
    psi, zeta = stommelSolution(Bx, By, curl, r, nmax, tol)

    # plotar
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(20,6))
    plt.subplots_adjust(wspace = .5)

    axes[0].set_title(u'Topografia [2D]')
    axes[0] = plotar_topografia2D(axes[0], hb)

    axes[1].set_title(u'VP Ambiente')
    axes[1],fig = plotar_campoVPAmb(axes[1],fig, B)

    axes[2].set_title(u'Solução de Stommel')
    axes[2] = plotar_solucaoStommel(axes[2], psi)

    plt.suptitle(r'Ex2-F: Solução de Stommel para fundo com variação linear zonal e $\beta$ = $\alpha x$ = $\frac{\beta}{4}x$', fontsize=20)

    plt.show()
