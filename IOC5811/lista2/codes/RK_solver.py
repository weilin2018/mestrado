#-*-coding:utf-8-*-
"""
Converting RK_soltar_lista2.m from matlab to python.

In the process, answering exercise 2.E from the second homework of DFGII

Some useful sources and links:

. Eigenvalue and eigenvector problem:
    mathematical theory:
        https://www.respondeai.com.br/resumos/33/capitulos/1
    numerical solution using Jacobi's Method:
        http://www.southampton.ac.uk/~feeg6002/lecturenotes/feeg6002_numerical_methods08.pdf


. https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html

"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use("ggplot")


def importMatfile(fname=''):

    import scipy.io as sio

    fname = '../outputs/mats/'+fname

    data = sio.loadmat(fname)

    cr  = np.squeeze(data['cr'])
    sig = np.squeeze(data['sig'])
    P   = data['P']

    return cr, sig, P

def plotar_tudo(k,y,cr,sig,Pphasemax, Pmax, figname, title):
    """
        plotar 4 subplots with:
            cr x k*1000
            sig x k*1000

            Pmax x 1e-3*y
            Pphasemax x 1e-3*y
    """

    kplot = k*1000
    yplot = y*1e-3

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(13,12))

    # plot cr x k
    axes[0, 0].plot(kplot, cr)
    axes[0, 0].set_title(u'Velocidade de Fase - $c_r$')
    axes[0, 0].set_xlabel(u'Número de onda [$km^{-1}$]')
    axes[0, 0].set_ylabel(u'Velocidade de Fase [$m s^{-1}$]')
    axes[0, 0].set_xlim([0.003, 0.04])
    axes[0, 0].set_ylim([0., 0.05])
    axes[0, 0].text(0.0375717, 0.00453675, '(A)', ha='center')

    axes[0, 1].plot(kplot, sig)
    axes[0, 1].set_title(u'Taxa de Crescimento - $\sigma = k c_i$')
    axes[0, 1].set_xlabel(u'Número de onda [$km^{-1}$]')
    axes[0, 1].set_ylabel(u'Taxa de Crescimento [$dias^{-1}$]')
    axes[0, 1].set_ylim([0, 0.02])
    axes[0, 1].set_xlim([0.003, 0.04])
    axes[0, 1].text(0.0374733, 0.00266874, '(B)', ha='center')

    # plotar o k onde sigma é maximo:
    if title[-1] != "G":
        imax   = np.where(sig == np.nanmax(sig))[0]
        axes[0, 1].scatter(kplot[imax], sig[imax], c=80, marker='o')
        axes[0, 1].text(kplot[imax] + 0.000280, sig[imax] + 0.0003, u'$k_{max} = %1.3f km^{-1}$' % (kplot[imax]))

    axes[1, 0].plot(Pmax, yplot)
    axes[1, 0].set_title(u'Amplitude do Modo mais Instável')
    axes[1, 0].set_xlabel(u'Amplitude do Modo [$\phi$]')
    axes[1, 0].set_ylabel(u'Distância Meridional [km]')
    axes[1, 0].set_ylim([-L*1e-3, L*1e-3])
    axes[1, 0].set_xlim([0., 0.18])

    axes[1, 0].text(0.16984, -176.77, '(C)', ha='center')

    axes[1, 1].plot(Pphasemax, yplot)
    axes[1, 1].set_title('Fase')
    axes[1, 1].set_xlabel(u'Fase [$^o$]')
    axes[1, 1].set_ylabel(u'Distância Meridional [km]')
    axes[1, 1].set_ylim([-L*1e-3, L*1e-3])
    axes[1, 1].set_xlim([-70, 0])
    axes[1, 1].text(-4.99454, -176.77, '(D)', ha='center')

    plt.suptitle(u'Solução da Equação de Rayleigh-Kuo - ' + title, fontsize=20)

    if figname!='': # então rola salvar
        plt.savefig('../outputs/'+figname, dpi=100)
    else:
        plt.show()


### SET GLOBAL VARIABLES
lat = 30
# set beta = 0 for f-plane
beta = (2 * 7.29e-5 * math.cos(lat * np.pi/180))/6371e3

L  = 200e3 # 1/2 domain
dy = 2e3
y  = np.arange(-L, L+dy, dy)
ny = len(y)

# jet profile
L0 = 50e3
U0 = 0.085
U  = np.asarray([U0 * 1/math.cosh(Y/L0)**2 for Y in y])
Uy = np.gradient(U, dy)
Uyy= np.gradient(Uy, dy)

Qy = beta - Uyy

# # plotar
# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

# ax1.plot(U, 1e-3*y, 'b')
# ax1.set_xlabel('jet velocity in m/s')
# ax1.set_ylabel('meridional distance in m')

# ax2.plot(Qy, 1e-3*y, 'g')
# ax2.axvline(0, color='black')
# ax2.set_xlabel('potential vorticity gradiente in 1/s')
# ax2.set_ylabel('meridional distance in m')
# plt.show()

#### THE EIGENVALUE PROBLEM


"""
Entendendo o algoritmo do Prof Ilson

Problema:
    (A-Bc)P = 0
        A - matriz de dados
        B - matriz identidade
        c - autovalor associado
        P - autovetor de A $P \neq 0$

    A e B são matrizes tridiagonais, ou seja, matrizes em que os elementos
    não nulos são da diagonal principal e 1 diagonal acima e 1 abaixo da diagonal
    principal:

        | a1 b1           |
        | c1 a2 b2        |
    A = |    c2 a3 b3     |
        |      c3 a4 b4   |
        |        c4 a5 b5 |

    Definição de algumas variáveis:
        k   - número de ondas em 1/km
        nk  - comprimento de k

    criar matrizes que armazenarão os dados de:
        cr  - velocidade de fase (real)
        sig - taxa de crescimento, dada por $\sigma = k ci$, sendo ci a parte imaginária de c
        P   - modos instáveis ou o autovetor da matriz de dados




"""

# wavenumbers in 1/km
k = np.asarray(np.arange(0.04, 0.002, -0.001)) * 1e-3
nk = len(k)

# initialize the phase speed (cr), growth rate (sig from sigma) and unstable modes (P) matrices
cr = np.zeros(nk)
sig = cr
P = np.zeros([ny,nk])

# for n in range(0,nk,1): # wavenumber loop

#     # setting up the matrix coefficients
#     L1 = 1/dy/dy
#     L2 = 2/dy/dy + k[n]*k[n]
#     L3 = L1

#     # building tridiagonal matrices A and B
#     A = np.diag(Qy - L2 * U) + np.diag( L3 * U[1:ny],1) + np.diag(L1*U[1:ny] ,-1)
#     B = np.diag(-L2 * np.squeeze(np.ones([1,ny]))) + np.diag(L3*np.squeeze(np.ones([1,ny-1])), 1) + np.diag(L1*np.squeeze(np.ones([1,ny-1])), -1)
#     # obtaining the eigenvalue matrix
#     C = np.linalg.lstsq(B, A)

#     # calculate eigenvalues and eigenvector for k[n]
#     # return evals/lambda, evecs/F
#     lamb, F  = np.linalg.eig(C[0])
#     # extract the max imaginary part
#     ci = np.nanmax(np.diag(np.imag(lamb)))

#     if ci == 0:
#         sig[n] = 0
#         cr[n]  = np.nan
#         P[:,n] = np.nan*np.ones([ny])
#     else:
#         # find the index of ci
#         jmax = np.where(np.imag(lamb) == ci)
#         sig[n] = k[n]*ci*86400          # growth rate in 1/days
#         cr[n]  = np.real(lamb[jmax,jmax])  # phase speed in m/s
#         P[:,n] = F[:, jmax]             # meridional mode structure

exercs = ['exF.mat', 'exG.mat', 'exH.mat']

for ex in exercs:
    print("Plotando dados do %s" % (ex))
    # import data
    cr, sig, P = importMatfile(fname=ex)

    Pamp   = np.abs(P)
    Pphase = 180/np.pi * np.angle(P) # np.angle retorna o angulo do argumento complexo

    # finding and normalizing most unstable mode
    imax   = np.where(sig == np.nanmax(sig))[0]
    sigmax = sig[imax]

    Pmax      = np.squeeze(Pamp[:, imax])

    # orthonomalization
    Pmax      = Pmax / (np.linalg.norm(Pmax))
    Pphasemax = Pphase[:, imax]

    plotar_tudo(k, y, cr, sig, Pphasemax, Pmax, figname=ex.replace('.mat', '.png'), title='Ex %s'%(ex[2:3]))


import glob
import os

lstFile = glob.glob('../outputs/ex*')
for fname in lstFile:
    os.system('convert -trim %s %s' % (fname, fname))
