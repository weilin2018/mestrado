# add some description here

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

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

# selecting file
BASE   = '/media/danilo/Danilo/mestrado//artigo_data/simulacoes/2ndPhase/run01_validation/'  # data directory
fname  = BASE + 'run01_validation.cdf'                      # gcmplt name

ncdata = xr.open_dataset(fname)                             # opening file

mask   = ncdata['FSM'].data.copy()                          # copy FSM data
# extract depth data, using the mask matrix for land and ocean cells
depth  = ncdata['depth'].data.copy()*mask
time   = ncdata['time'].data.copy()

# depth with 0. replaced by nan values
depth[depth == 0.] = np.nan

# initializing variables
n      = time.shape[0]                             # timesteps

# defining vectors
KIN_s  = np.zeros(len(time))                       # kinetic energy at surface
KIN_b  = np.zeros(len(time))                       # kinetic energy at bottom
POT_s  = np.zeros(len(time))                       # potential energy at surface

# defining auxiliar variables
kin_s = 0

os.system('clear')                                 # clear screen
print('Calculating energy ...\n')


# running each timestep and calculating energy
# os.system('clear')
for i in range(n):
    print('Timestep: %i' % (i))

    # extracting data from ncdata for each timestep
    t     = time[i]
    salt  = ncdata.salt[i].values
    temp  = ncdata.temp[i].values
    elev  = ncdata.elev[i].values
    u     = ncdata.u[i].values
    v     = ncdata.v[i].values
    sigma = ncdata.sigma.values
    lat   = ncdata.lat.values


    kin   = 0
    pot   = 0
    kin_b = 0
    for J in np.arange(2,135,1):

        ISTART = np.min(np.where(~np.isnan(depth[J,:])))
        IEND   = np.max(np.where(depth[J,:] <= 200))

        # if you want to calculate the energy in the entire range of I,
        # replace IEND by 108 or uncomment the next line
        # IEND = 108

        for I in np.arange(ISTART,IEND):
            print('working on J=%i,I=%i'%(J,I))

            if ~np.isnan(depth[J,I]):

                S = salt[0,J,I]
                T = temp[0,J,I]
                p = depth[J,I]
                la= lat[J,I]

                rho = sw.seaWaterDensity(S,T,p,la)

                HT = depth[J,I] + np.abs(elev[J,I])
                Z  = np.abs(elev[J,I])

                kin += (rho/2) * Z * ( np.nanmean([np.abs(u[0,J,I]),np.abs(u[0,J,I+1])])**2 + np.nanmean([np.abs(v[0,J,I]), np.abs(v[0,J+1,I])])**2 )
                pot += rho*9.8*(Z**2)

                ### BOTTOM
                N = len(sigma)-2

                S = salt[N,J,I]
                T = temp[N,J,I]
                p = depth[J,I]
                la= lat[J,I]

                rho = sw.seaWaterDensity(S,T,p,la)

                HT = depth[J,I] + np.abs(elev[J,I])
                Z  = np.abs(elev[J,I])

                kin_b += (rho/2) * Z * ( np.nanmean([np.abs(u[N,J,I]),np.abs(u[N,J,I+1])])**2 + np.nanmean([np.abs(v[N,J,I]), np.abs(v[N,J+1,I])])**2 )

    KIN_s[i] = kin
    POT_s[i] = pot
    KIN_b[i] = kin_b


df_kin = pd.DataFrame({'kin_s':KIN_s,'kin_b':KIN_b},index=pd.DatetimeIndex(time))
df_pot_Z = pd.DataFrame({'pot_s':POT_s},index=pd.DatetimeIndex(time))


# figures
plt.subplot(121);
plt.semilogy(df_kin.resample('12H').mean());
# plt.semilogy(time,KIN_b);
# plt.plot(time,KIN_s);
# plt.plot(time,KIN_b);
plt.legend(['surface','bottom'], loc='best');
plt.ylabel('Kinetic energy/area (kg/s^2)');
plt.xlabel('Days');
plt.subplot(122);
plt.semilogy(df_pot_Z, color='red');
# plt.plot(time,POT_s, color='red');
plt.ylabel('Potential energy/area (kg/s^2)');
plt.xlabel('Days');

plt.tight_layout(pad=.5, w_pad=0.5, h_pad=1.0);
