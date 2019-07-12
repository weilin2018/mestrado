#!/usr/bin/env python2.7
#coding: utf8
import sys
sys.path.append('/media/paula/Paulinha/@OCEANO/mestrado/calculate_energy_model/')
reload(sys)

sys.setdefaultencoding('utf8')


import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')  # or whatever other backend that you want
import matplotlib.pyplot as plt
import subprocess
import seaWaterDensity as sw

# Leitura de pacotes:
import xarray as xr


# Leitura de dados:
a=xr.open_dataset('/media/paula/Paulinha/@OCEANO/mestrado/BACKUP_SHELFWAVES_PAULA_MAIO_2017/rodada_teste74/gcmplt74.cdf')

# Armazenamento de variáveis:

run='rodada_teste74'
mask = a['FSM'].data.copy()
depth = a['depth'].data.copy()*mask # H
depth[depth==-0] = np.nan
time=a['time'].data.copy()

def plot(run,time,depth,salt,temp,sigma,lat,elev,u,v):

	print("plotting energy graphs")

	# initiate kin and pot variables = 0 (surface and bottom)
	kin_s = 0
	pot_s = 0
	kin_b = 0

	# create variables to store KIN and POT energies over time, for surface and bottom layers
	KIN_s = np.empty(len(time))
	POT_s = np.empty(len(time))
	KIN_b = np.empty(len(time))
	KIN_s[:] = np.NAN
	POT_s[:] = np.NAN
	KIN_b[:] = np.NAN

	for T in np.arange(0,len(time)):

		print("\ntime = " + str(time[T]) + " days")
		#posj = np.where(~np.isnan(depth[:,ISTART:IEND]))
		#if len(posj[0] > 0):
		#	JSTART = np.min(posj)
		#	JEND = np.max(posj)
		# CUIDADO! APARECE J+1
		for J in np.arange(1,132):

			pos = np.where(~np.isnan(depth[J,:]))

			if len(pos[0] > 0):
				ISTART = np.min(pos)
				IEND = np.max(pos)

			for I in np.arange(ISTART,IEND):

				if ~np.isnan(depth[J,I]):

					### --- surface --- ###
					N = 0
					Sa = salt[T,N,J,I]
					Te = temp[T,N,J,I]
					# if Sa==0:
					# 	print 'Sa é zero'
					# 	print T
					# 	print N
					# 	print J
					# 	print I
					# if Te==0:
					# 	print 'Temp é zero'
					# 	print T
					# 	print N
					# 	print J
					# 	print I
					# if np.isnan(depth[J,I])==True:
					# 	J=J+1

					prof = np.abs(depth[J,I]*sigma[N])

					# if np.isnan(prof)==True:
					# 	J=J+1
					# 	print 'Prof é nan'
					# 	print T
					# 	print N
					# 	print J
					# 	print I

					la = lat[J,I]
					# if np.isnan(la)==True:
					# 	print 'la é nan'
					# 	print T
					# 	print N
					# 	print J
					# 	print I

					# I wrote this function from MatLab seawater toolbox
					# check again!
					rho = sw.seaWaterDensity(Sa, Te, prof, la)
					# if np.isnan(rho)==True:
					# 	print 'rho é nan'
					# 	print T
					# 	print N
					# 	print J
					# 	print I

					# print(rho)

					# I changed local depth for the depth of the layer considered; I think it makes sense
					# check again!
					# kin_s and pot_s: Z or HT????
					# check again!
					# HT = depth[J,I] + np.abs(elev[T,J,I])
					# HT = prof + np.abs(elev[T,J,I])
					Z = np.abs(elev[T,J,I])
					kin_s = kin_s + (rho/2)*Z*( (np.nanmean([np.abs(u[T,N,J,I]), np.abs(u[T,N,J,I+1])]))**2 + (np.nanmean([np.abs(v[T,N,J,I]), np.abs(v[T,N,J+1,I])]))**2 )
					# kin_s = kin_s + (rho/2)*HT*( (np.nanmean([np.abs(u[T,N,J,I]), np.abs(u[T,N,J,I+1])]))**2 + (np.nanmean([np.abs(v[T,N,J,I]), np.abs(v[T,N,J+1,I])]))**2 )
					pot_s = pot_s + rho*9.8*(Z**2)

					### --- bottom --- ###
					# 1 layer above bottom
					N = len(sigma)-2

					Sa = salt[T,N,J,I]
					Te = temp[T,N,J,I]
					if np.isnan(Sa)==True:
						print 'Sa é nan'
						print T
						print N
						print J
						print I
					if np.isnan(Te)==True:
						print 'Temp é nan'
						print T
						print N
						print J
						print I


					prof = np.abs(depth[J,I]*sigma[N])
					if np.isnan(prof)==True:
						print 'Prof é nan'
						print T
						print N
						print J
						print I
					la = lat[J,I]
					if np.isnan(la)==True:
						print 'la é nan'
						print T
						print N
						print J
						print I

					rho = sw.seaWaterDensity(Sa, Te, prof, la)
					if np.isnan(rho)==True:
						print 'rho é nan'
						print T
						print N
						print J
						print I

					# I wrote this function from MatLab seawater toolbox
					# check again!

					# print(rho)

					# HT = depth[J,I] + np.abs(elev[T,J,I])
					Z = np.abs(elev[T,J,I])
					kin_b = kin_b + (rho/2)*Z*( (np.nanmean([np.abs(u[T,N,J,I]), np.abs(u[T,N,J,I+1])]))**2 + (np.nanmean([np.abs(v[T,N,J,I]), np.abs(v[T,N,J+1,I])]))**2 )
					# kin = kin + (rho/2)*HT*( (np.nanmean([np.abs(u[T,N,J,I]), np.abs(u[T,N,J,I+1])]))**2 + (np.nanmean([np.abs(v[T,N,J,I]), np.abs(v[T,N,J+1,I])]))**2 )
					if np.isnan(kin_b)==True:
						kin_b = 0.0

		KIN_s[T] = kin_s
		POT_s[T] = pot_s
		KIN_b[T] = kin_b
		print("surface KIN = " + '{:.5f}'.format(KIN_s[T]))
		print("bottom KIN = " + '{:.5f}'.format(KIN_b[T]))
	return KIN_s, POT_s ,KIN_b


KIN_s, POT_S, KIN_b=plot(run,time,depth,a['salt'].data,a['temp'].data,a['sigma'].data,a['lat'].data,a['elev'].data,a['u'].data,a['v'].data)

# figures
plt.subplots_adjust(wspace=0, hspace=0)
plt.subplot(121);
plt.semilogy(time,KIN_s);
plt.ylim(ymax=10**9)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

plt.semilogy(time,KIN_b);
# plt.plot(time,KIN_s);
# plt.plot(time,KIN_b);
plt.legend(['surface','bottom'], loc='best');
plt.ylabel(r'Kinetic energy/area (kg/s$^2$)',fontsize=14);
plt.xlabel('Days',fontsize=14);
plt.xticks(rotation=25)
plt.subplot(122);
plt.semilogy(time,POT_S, color='red');
# plt.plot(time,POT_s, color='red');
plt.ylabel(r'Potential energy/area (kg/s$^2$)',fontsize=14);
plt.xlabel('Days',fontsize=14);
plt.xticks(rotation=25)
plt.ylim(ymax=10**9)
plt.tight_layout(pad=.5, w_pad=0.3, h_pad=1.0);

plt.savefig('/media/paula/Paulinha/@OCEANO/mestrado/calculate_energy_model/figures/' + run +'/semilogyenergy_paula8_74' + run +'.pdf', bbox_inches='tight', dpi = 130)

# Mudar o estilo de plot:
matplotlib.style.use('ggplot')

import pandas as pd
dt = pd.date_range('12-26-2004 03:00','01-29-2005 21:00',freq='6H')
df = pd.DataFrame({'KIN_s':KIN_s,'POT_s':POT_s,'KIN_b':KIN_b},index=dt)
df.plot(subplots=True)
