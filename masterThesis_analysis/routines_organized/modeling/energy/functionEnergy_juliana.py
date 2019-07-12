# -*- coding: utf-8 -*-

#script para calcular a energia cinética e potencial do modelo


import numpy as np
import matplotlib.pyplot as plt
import subprocess
import seaWaterDensity as sw

# Leitura de pacotes:
import xarray as xr

def imprimir(frase, n,j,i):
    print(frase)
    print(N)
    print(J)
    print(I)

def energia_no_tempo(time,depth,salt,temp,sigma,lat,elev,u,v,kin_s, kin_b, pot_s):

    for J in range(1,132):

        pos = np.where(~np.isnan(depth[J,:])) #pegando todos valores de prof que nao sao NaN

        if len(pos[0]>0):
            ISTART = np.min(pos)
            IEND   = np.max(pos)

    for I in range(ISTART,IEND):

        if ~np.isnan(depth[J,I]):
            ### surface ###
            N  = 0 #superficie
            Sa = salt[N,J,I] #nao preciso mais do T poruqe estou definindo ele no for
            Te = temp[N,J,I]
            if Sa==0:
                imprimir(u'Sa é zero', N, J, I)

            if Te==0:
                imprimir(u'Te é zero', N, J, I)

            if np.isnan(depth[J,I])==True:
                J=J+1

            prof = np.abs(depth[J,I]*sigma[N]) #sta pegando o nivel primeiro nivel sigma

            if np.isnan(prof)==True:
                J=J+1
                imprimir(u'Prof é nan', N, J, I)

            la = lat[J,I]
            if np.isnan(la)==True:
                imprimir(u'la é nan', N, J, I)

            rho = sw.seaWaterDensity(Sa, Te, prof, la)
            if np.isnan(rho)==True:
                imprimir(u'rho é nan', N, J, I)

            Z = np.abs(elev[J,I]) #np.abs retorna todos valores de forma positiva (+)

            # kin_s += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 ) #energia cinetica
            kin_s = kin_s + ((rho/2)*Z*(np.abs(u[N,J,I])**2) + (np.abs(v[N,J+1,I])**2) )
            pot_s += rho*9.8*(Z**2) #energia potencial

            ### bottom ###

            N = len(sigma)-2 #uma camada a cima do fundo
            Sa = salt[N,J,I]
            Te = temp[N,J,I]
            if np.isnan(Sa)==True:
                imprimir(u'Sa é nan', N, J, I)


            if np.isnan(Te)==True:
                imprimir(u'Te é nan', N, J, I)


            prof = np.abs(depth[J,I]*sigma[N]) #esta pegando o nivel primeiro nivel sigma
            if np.isnan(prof)==True:
                imprimir(u'Prof é nan', N, J, I)

            la = lat[J,I]
            if np.isnan(la)==True:
                imprimir(u'la é nan', N, J, I)

            rho = sw.seaWaterDensity(Sa, Te, prof, la)
            if np.isnan(rho)==True:
                imprimir(u'rho é nan', N, J, I)


            Z = np.abs(elev[J,I]) #np.abs retorna todos valores de forma positiva (+)

            # kin_b += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 ) #energia cinetica
            kin_b = kin_b + (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 )


    return kin_s, pot_s, kin_b

# a = xr.open_dataset('/data0/juliana/inverno/run6/gcmplt.cdf')

a = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/output/exp06.cdf')

mask = a['FSM'].data.copy() #mascara de pontos em terra
depth = a['depth'].data.copy()*mask # H
depth[depth==-0] = np.nan
time=a['time'].data.copy()


n = time.shape[0]
kin_s = 0
pot_s = 0
kin_b = 0

KIN_s = []
POT_s = []
KIN_b = []

import os
os.system('clear')

for i in range(n):
    print(i)
    time = a.time[i]
    salt = a.salt[i]
    temp = a.temp[i]
    sigma = a.sigma
    lat = a.lat
    elev = a.elev[i]
    u = a.u[i]
    v = a.v[i]

    kin_s, pot_s, kin_b = energia_no_tempo(time,depth,salt,temp,sigma,lat,elev,u,v,kin_s,kin_b, pot_s)

    KIN_s.append(kin_s.values)
    POT_s.append(pot_s.values)
    KIN_b.append(kin_b.values)

    if i % 100 == 0: #salva a cada 100 passos de tempo
                     #quando 1 = 0, entao o resto da divisao é zero.
                     #nesse momento salva.
        np.save('kin_s_wint2013.npy',KIN_s) #salva a energia cinetica de superficie
        np.save('kin_b_wint2013.npy',KIN_b) #salva a energia cinetica de fundo
        np.save('pot_s_wint2013.npy',POT_s) #salva a energia potencial de superficie
