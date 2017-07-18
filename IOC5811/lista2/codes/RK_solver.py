#-*-coding:utf-8-*-
""" 
Converting RK_soltar_lista2.m from matlab to python.

In the process, answering exercise 2.E from the second homework of DFGII
"""

import math
import numpy as np
import matplotlib.pyplot as plt



### SET GLOBAL VARIABLES
lat = 30
# set beta = 0 for f-plane
beta = (2 * 7.29e-5 * math.cos(lat * np.pi/180))/6371e3

L  = 200e3 # 1/2 domain
dy = 2e3
y  = np.arange(-L, L, dy)
ny = len(y)

# jet profile
L0 = 50e3
U0 = 0.085
U  = np.asarray([U0 * 1/math.cosh(Y/L0)**2 for Y in y])
Uy = np.gradient(U, dy)
Uyy= np.gradient(Uy, dy)

Qy = beta - Uyy

# plotar
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True)

ax1.plot(U, 1e-3*y, 'b')
ax1.set_xlabel('jet velocity in m/s')
ax1.set_ylabel('meridional distance in m')

ax2.plot(Qy, 1e-3*y, 'g')
ax2.axvline(0, color='black')
ax2.set_xlabel('potential vorticity gradiente in 1/s')
ax2.set_ylabel('meridional distance in m')
plt.show()

#### THE EIGENVALUE PROBLEM

# wavenumbers in 1/km
k = np.asarray(np.arange(0.04, 0.002, -0.001)) * 1e-3
nk = len(k)

# initialize the phase speed (cr), growth rate (sig from sigma) and unstable modes (P) matrices
cr = np.zeros(nk)
sig = cr
P = np.zeros([ny,nk])

for n in range(0,nk,1): # wavenumber loop

    # setting up the matrix coefficients
    L1 = 1/dy/dy
    L2 = 2/dy/dy + k[n]*k[n]
    L3 = L1

    # building tridiagonal matrices A and B
    A = np.diag(Qy - L2 * U) + np.diag( L3 * U[1:ny],1) + np.diag(L1*U[1:ny] ,-1)
    B = np.diag(-L2 * np.squeeze(np.ones([1,ny]))) + np.diag(L3*np.squeeze(np.ones([1,ny-1])), 1) + np.diag(L1*np.squeeze(np.ones([1,ny-1])), -1)
    # obtaining the eigenvalue matrix
    C = np.linalg.solve(B, A)

    # calculate eigenvalues and eigenvector for k[n]
    F,lamb  = np.linalg.eig(C)
    w = np.diag(lamb) # create matrix [ny,ny]
    ci = np.nanmax(np.diag(np.imag(w)))

    if ci != 0:
        print(ci)

    # if ci == 0:
    #     sig[n] = 0
    #     cr[n]  = np.nan
    #     P[:,n] = np.nan*np.ones([ny])
    # else:
    #     sig[n] = k[n]*ci*86400          # growth rate in 1/days
    #     cr[n]  = np.real(w[jmax,jmax])  # phase speed in m/s
    #     P[:,n] = F[:, jmax]             # meridional mode structure

