
"""

Adapted by Danilo A. Silva <nilodna@gmail.com>, based on Carine and Paula's version

"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle


import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')
import masterThesisPack as oceano

sys.path.append('modelling/')
import seaWaterDensity as sw

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def calcular(T,depth,salt,temp,sigma,lat,elev,u,v,kin_s,kin_b,pot_s):
    ''' '''

    print("\ntime = " + str(T) + " days")
    for J in np.arange(1,132):
        pos = np.where(~np.isnan(depth[J,:]))

        if len(pos[0] > 0):
            ISTART = np.min(pos)
            IEND = np.max(pos)

        for I in np.arange(ISTART,IEND):
            if ~np.isnan(depth[J,I]):

                ### calculando na superficie
                N     = 0
                Sa    = salt[N,J,I]
                Te    = temp[N,J,I]

                prof  = np.abs(depth[J,I]*sigma[N])
                la    = lat[J,I]

                rho   = sw.seaWaterDensity(Sa,Te,prof,la)

                Z     = np.abs(elev[J,I])
                kin_s += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 )
                pot_s += rho*9.8*(Z**2)

                ### calculando no fundo
                N      = len(sigma)-2
                Sa    = salt[N,J,I]
                Te    = temp[N,J,I]

                prof  = np.abs(depth[J,I]*sigma[N])
                la    = lat[J,I]

                rho   = sw.seaWaterDensity(Sa,Te,prof,la)

                Z     = np.abs(elev[J,I])

                kin_b += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 )
                if np.isnan(kin_b)==True:
                    kin_b = 0.0

    return kin_s,pot_s,kin_b

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

# selecting file
BASE   = '/media/danilo/Danilo/mestrado/ventopcse/output/'  # data directory
fname  = BASE + 'exp06.cdf'                                 # gcmplt name

ncdata = xr.open_dataset(fname)                             # opening file

mask   = ncdata['FSM'].data.copy()                          # copy FSM data
# extract depth data, using the mask matrix for land and ocean cells
depth  = ncdata['depth'].data.copy()*mask
time   = ncdata['time'].data.copy()
sigma  = ncdata['sigma'].values
lat    = ncdata['lat'].values

# depth with 0. replaced by nan values
depth[depth == 0.] = np.nan

# initializing variables
n      = time.shape[0]                             # timesteps

# defining vectors
KIN_s  = np.zeros(len(time))                       # kinetic energy at surface
KIN_b  = np.zeros(len(time))                       # kinetic energy at bottom
POT_s  = np.zeros(len(time))                       # potential energy at surface

kin_s,kin_b,pot_s = 0,0,0

os.system('clear')                                 # clear screen
print('Calculating energy ...\n')

def calculateRho(N,J,I,salt,temp,depth,sigma,lat):

    Sa = salt
    Te = temp
    la = lat

    prof = np.abs(depth*sigma)

    rho = sw.seaWaterDensity(Sa,Te,prof,la)

    return rho


for i in np.arange(0,n):

    # import data from ncdata
    T     = time[i]
    salt  = ncdata['salt'][i,:,:,:].values
    temp  = ncdata['temp'][i,:,:,:].values
    u     = ncdata['u'][i,:,:,:].values
    v     = ncdata['v'][i,:,:,:].values
    elev  = ncdata['elev'][i,:,:].values



    print("\ntime = " + str(T) + " days")
    for J in np.arange(1,132):
        pos = np.where(~np.isnan(depth[J,:]))

        if len(pos[0] > 0):
            ISTART = np.min(pos)
            IEND = np.max(pos)

        for I in np.arange(ISTART,IEND):
            if ~np.isnan(depth[J,I]):

                ### calculando na superficie
                N     = 0

                # rho = calculateRho(N,J,I,salt[N,J,I],temp[N,J,I],depth[J,I],sigma[N],lat[J,I])

                Sa    = salt[N,J,I]
                Te    = temp[N,J,I]

                prof  = np.abs(depth[J,I]*sigma[N])
                la    = lat[J,I]

                rho   = sw.seaWaterDensity(Sa,Te,prof,la)

                Z     = np.abs(elev[J,I])
                kin_s += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 )
                pot_s += rho*9.8*(Z**2)

                ### calculando no fundo
                N      = len(sigma)-2
                Sa    = salt[N,J,I]
                Te    = temp[N,J,I]

                prof  = np.abs(depth[J,I]*sigma[N])
                la    = lat[J,I]

                rho   = sw.seaWaterDensity(Sa,Te,prof,la)

                Z     = np.abs(elev[J,I])

                kin_b += (rho/2)*Z*( (np.nanmean([np.abs(u[N,J,I]), np.abs(u[N,J,I+1])]))**2 + (np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])]))**2 )
                if np.isnan(kin_b)==True:
                    kin_b = 0.0

    KIN_s[i] = kin_s
    POT_s[i] = pot_s
    KIN_b[i] = kin_b


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
plt.semilogy(time,POT_s, color='red');
# plt.plot(time,POT_s, color='red');
plt.ylabel(r'Potential energy/area (kg/s$^2$)',fontsize=14);
plt.xlabel('Days',fontsize=14);
plt.xticks(rotation=25)
plt.ylim(ymax=10**9)
plt.tight_layout(pad=.5, w_pad=0.3, h_pad=1.0);
