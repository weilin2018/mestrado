# add some description here

import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import dates
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def domain_area(xe,xw,yn,ys,R):
    # area of domain, expressed in m^2
    A = R*R*(xw-xe)*((2*np.pi)/360.)*(np.sin(ys)-np.sin(yn))

    return A

def load_Surface_HeatFlux(nfile):
    cdf_hf = xr.open_dataset(nfile)
    sfc_hflux = cdf_hf['THFLX_L1_Avg_1'].values
    x_hf = cdf_hf['lon'].values - 360
    y_hf = cdf_hf['lat'].values

    dA_hf = R*R*.1*((2*np.pi)/360)*(np.sin(lat[1,0]) - np.sin(lat[0,0]))

    mask_hf = np.ones([len(y_hf),len(x_hf)])

    return sfc_hflux, dA_hf, mask_hf


def calculate_averageVelocity_centerGridCell(data,begin,end):

    out = np.zeros(data.shape)
    for i in np.arange(begin,end):
        out[:,:,i] = np.mean(data[:,:,i:i+1],axis=2)

    return out

def set_timeIndexes(times,startDate,finalDate):
    try:
        from dateutil import parser
    except:
        print("Please install python-dateutils package to use this function")
        pass

    refDate_start = parser.parse(str(startDate))
    refDate_final = parser.parse(str(finalDate))
    df = pd.DataFrame(np.arange(0,len(times)),index=times)

    # performing search
    startIndex = df.iloc[df.index.get_loc(refDate_start,method='nearest')]
    finalIndex = df.iloc[df.index.get_loc(refDate_final,method='nearest')]

    return startIndex,finalIndex


def temperatureBudget(year,DATA_DIR,MERCATOR):

    # defining some constants
    R = 6371*1000  # earth's radius
    izh = 25       # represeting z = -120m as my mixed layer depth
    h = 100        # mixed layer depth
    cp      = 3990 # J degC^-1 kg^-1: specific heat at constant pressure
    rho_ref = 1035 # kg m^-3, reference density

    ncin = xr.open_dataset(DATA_DIR + MERCATOR)

    # selecting index for the period (Jan and Feb) from the year defined in #81
    startIndex,finalIndex = set_timeIndexes(ncin.time.values,'%s-01-01'%(str(year)),'%s-03-05'%(str(year)))
    startIndex,finalIndex = startIndex[0],finalIndex[0]

    # extracting variables
    x    = ncin.longitude.values
    y    = ncin.latitude.values
    z    = ncin.depth.values
    time = ncin.time[startIndex:finalIndex].values
    temp = ncin.temperature[startIndex:finalIndex,:izh,:,:].values

    # depth of each cell
    DZ = -np.diff(z[:izh+1])

    nx = len(x)
    ny = len(y)
    nz = len(z)
    nt = len(time)

    # gridding lon,lat
    lon,lat = np.meshgrid(x,y)

    # defining indexes for boundaries
    iyn = 0        # index for north section
    iys = len(x)-1 # index for south section
    ixw = 0        # index for west section
    ixe = len(y) -1 # index for east section

    # calculating the area of the domain
    A = domain_area(x[ixe],x[ixw],y[iyn],y[iys],R)

    # calculating the area of each cell.
    # Considering a regular grid, the area will be the same
    dA =  R*R*.1*((2*np.pi)/360)*(np.sin(lat[1,0]) - np.sin(lat[0,0]))

    # calculating the width of each face along edges of domain
    # .25 is the spacing in deg
    dy_w = R*((.25*2*np.pi)/360)                    # m, meridional lenght
    dy_e = R*((.25*2*np.pi)/360)                    # m, meridional lenght
    dx_n = R*np.cos(y[iyn])*((.25*2*np.pi)/360)     # m, zonal lenght
    dx_s = R*np.cos(y[iys])*((.25*2*np.pi)/360)     # m, zonal lenght

    # area of each face along edge of domain
    int_dy_w = np.nansum(ny*dy_w)
    int_dy_e = np.nansum(ny*dy_e)
    int_dx_n = np.nansum(nx*dx_n)
    int_dx_s = np.nansum(nx*dx_s)

    #### LOADING HEAT FLUX DATA
    sfc_hflux, dA_hf, mask_hf = load_Surface_HeatFlux('/home/danilo/Dropbox/mestrado/data/data2model/JF2014/hflx/heatflux_JF2014.nc')
    mask_hf[np.isnan(np.squeeze(np.mean(sfc_hflux,axis=0)))] = 0

    ##################################################################
    ##### calculating fluxes on the faces of the domain
    ##################################################################
    # extracting u,v components for each face of the domain
    ue1 = ncin.u[startIndex:finalIndex,:izh,iyn:iys,ixe].values # east face
    uw1 = ncin.u[startIndex:finalIndex,:izh,iyn:iys,ixw].values # west face
    vn1 = ncin.v[startIndex:finalIndex,:izh,iyn,ixw:ixe].values # north face
    vs1 = ncin.v[startIndex:finalIndex,:izh,iys,ixw:ixe].values # south face

    # calculating average velocity onto the center of the grid cell
    ue = calculate_averageVelocity_centerGridCell(ue1,iyn,iys)
    uw = calculate_averageVelocity_centerGridCell(uw1,iyn,iys)
    vn = calculate_averageVelocity_centerGridCell(ue1,ixw,ixe)
    vs = calculate_averageVelocity_centerGridCell(uw1,ixw,ixe)

    # average temperature onto the grid cell
    tw = np.squeeze(np.mean(temp[:,:izh,iyn:iys,ixw:ixw+1],axis=3))
    te = np.squeeze(np.mean(temp[:,:izh,iyn:iys,ixe:ixe+1],axis=3))
    ts = np.squeeze(np.mean(temp[:,:izh,iys:iys+1,ixw:ixe],axis=2))
    tn = np.squeeze(np.mean(temp[:,:izh,iyn:iyn+1,ixw:ixe],axis=2))

    # Integrating over z
    def integrate_z(compponentVelocity,temperatureSide,nt,DZ):
        u = compponentVelocity
        t = temperatureSide
        int_data = np.zeros([nt,u.shape[2]])

        for i in np.arange(0,u.shape[1]):
            utw = u[:,:,i]*t[:,:,i]
            wet = ~np.isnan(utw[1,:])
            int_data[:,i] = np.nansum(utw * np.tile(np.transpose(DZ),[nt,1]),axis=1)/np.sum(DZ[wet])

        return int_data

    def integrate_horizontally(data,widthFace):
        # zonal or meridionally
        int_horizontally = np.nansum(data*widthFace,axis=1)

        return int_horizontally


    # integrating over z
    int_utw_d = integrate_z(uw,tw,nt,DZ)
    int_ute_d = integrate_z(ue,te,nt,DZ)
    int_vtn_d = integrate_z(vn,tn,nt,DZ)
    int_vts_d = integrate_z(vs,ts,nt,DZ)

    # integrating zonally or meridionally
    int_utw = integrate_horizontally(int_utw_d,dy_w)
    int_ute = integrate_horizontally(int_ute_d,dy_e)
    int_vtn = integrate_horizontally(int_vtn_d,dx_n)
    int_vts = integrate_horizontally(int_vts_d,dx_s)

    # calculating the contribution to the volume-averaged temperature tendency equation
    # on the right hand side. Note the minus sign
    int_adv_d = -(1/A)*(int_utw - int_ute + int_vtn - int_vts)

    ##########################################
    # volume-average tmeperature within domain and area-average surface heat flux
    DDZ           = np.zeros([nz-(nz-izh),ny,nx])
    for i in np.arange(0,nx):
        for j in np.arange(0,ny):
            DDZ[:,j,i] = DZ

    sfc_hflux_avg = np.zeros([nt,1])
    temp_avg      = np.zeros([nt,1])

    for i in np.arange(0,nt):
        sfc_hflux_avg[i,0] = np.nansum(np.nansum(np.squeeze(sfc_hflux[i,:,:])*dA_hf))/np.sum(np.sum(dA_hf*mask_hf))

        wet = np.squeeze(~np.isnan(temp[i,:,:,:]))
        temp_d = np.squeeze(np.nansum(np.squeeze(temp[i,:,:,:])*DDZ,axis=0)) / np.squeeze(np.sum(DDZ*wet,axis=0))
        mask_temp = np.ones(temp_d.shape)
        mask_temp[np.isnan(temp_d)]=0
        # print(np.sum(np.sum(dA*mask_temp)))
        temp_avg[i,:] = np.nansum(np.nansum(temp_d*dA))/np.sum(np.sum(dA*mask_temp))

    sfc_hflux_avg = np.squeeze(sfc_hflux_avg)
    temp_avg      = np.squeeze(temp_avg)

    # calculating Qv component
    Qsnet =  sfc_hflux_avg/(h*cp*rho_ref) # deg C s^-1

    # storing all info into a dictionary
    f = {}

    f['time'] = time
    f['int_adv_d'] = int_adv_d
    f['temp_avg'] = np.squeeze(temp_avg)
    f['Qsnet'] = np.squeeze(Qsnet)

    # individual contribution of each convection's component
    f['int_adv_de'] = -(1/A)*(int_ute)
    f['int_adv_dw'] = -(1/A)*(int_utw)
    f['int_adv_dn'] = -(1/A)*(int_vtn)
    f['int_adv_ds'] = -(1/A)*(int_vts)

    # calculating values for plot
    ## time-rate of change in volume-averaged temperature
    f['dTdt'] = f['temp_avg']- f['temp_avg'][1]
    ## Th = horizontal advection
    f['Th'] = np.matlib.cumsum(f['int_adv_d'])*(24*60*60)
    ## Tq = air-sea heat flux
    f['Tq'] = np.matlib.cumsum(f['Qsnet'])*(24*60*60)
    # estimating residual
    f['Residual'] = f['dTdt'] - f['Th'] - f['Tq']

    # plotting
    df = pd.DataFrame(f,index=f['time'])

    df['dTdt'].plot(color='k',label=r'Total (T$_{v}$)')
    df['Th'].plot(color='b',label=r'Advection (T$_{h}$)')
    df['Tq'].plot(color='r',label=r'Surface Heat Flux (T$_{q}$)')
    df['Residual'].plot(color='g',label=r'Residual')

    return df



##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code
DATA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/mercator/'
MERCATOR = 'data.nc'

# define which year to plot
year = 2014


df_2010 = temperatureBudget(2010,DATA_DIR,MERCATOR)
df_2014 = temperatureBudget(2014,DATA_DIR,MERCATOR)
#################
fig,ax = plt.subplots(nrows=2)

df_2010['dTdt'].plot(ax=ax[0],color='k',label=r'Total (T$_{v}$)')
df_2010['Th'].plot(ax=ax[0],color='b',label=r'Advection (T$_{h}$)')
df_2010['Tq'].plot(ax=ax[0],color='r',label=r'Surface Heat Flux (T$_{q}$)')
df_2010['Residual'].plot(ax=ax[0],color='g',label=r'Residual')

df_2014['dTdt'].plot(ax=ax[1],color='k',label=r'Total (T$_{v}$)')
df_2014['Th'].plot(ax=ax[1],color='b',label=r'Advection (T$_{h}$)')
df_2014['Tq'].plot(ax=ax[1],color='r',label=r'Surface Heat Flux (T$_{q}$)')
df_2014['Residual'].plot(ax=ax[1],color='g',label=r'Residual')

ax[1].legend()
ax[0].set_title('Temperature Budget Analysis over South Brazil Bight for 2010 (above) and 2014 (below)',fontsize=25)
