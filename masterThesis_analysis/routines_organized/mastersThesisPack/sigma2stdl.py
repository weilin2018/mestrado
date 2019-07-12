import numpy as np
import xarray as xr
from scipy import interpolate

fname = '/media/danilo/Danilo/mestrado/ventopcse/output/EC1.cdf'
ncin = xr.open_dataset(fname)

var = ncin.temp[0,:,:,:].values.copy()
sig = ncin.sigma.values.copy()
dep = ncin.depth.values.copy()
lon = ncin.lon.values.copy(); lon[lon == 0.]=np.nan
lat = ncin.lat.values.copy(); lat[lat == 0.]=np.nan

# only for a 3D variable, to maintain organization
nstdl           = 23         # new standard levels to interpolate
Ns,Mp,Lp        = var.shape  # old dimensions (sigma,JJ,II)
newvar          = np.zeros((nstdl,Mp,Lp))*np.nan # new array for data

# new Z dimension
stdl = [0, 10, 25, 40, 50, 60, 70, 80, 100, 150, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200, 1500, 1800, 2000]

for j in range(Mp): # running all lines
    for i in range (Lp): # running all columns
        localdep        = dep[j,i]

        # if we are in a valid depth
        if ~np.isnan(localdep):
            depsigma = -localdep*sig

    		# include surface with same value of first sigma level m to interpolate
            D = list(depsigma)
            D.insert(0, 0)

    		# select profile and include surface
            profile     = np.zeros(Ns+1)
            profile[1:] = var[:,j,i]
            profile[0]  = profile[1]

    		# watercolumn positions only
            watercolumn = stdl <= localdep
            stdl2interp = np.array(stdl)[watercolumn]

    		# interpolate to the same standard levels
            fsigma2stdl = interpolate.interp1d(D, profile)
            profileI = fsigma2stdl(stdl2interp)

    		# stores at vectorized variable
            newvar[watercolumn,j,i] = profileI
            newvar[~watercolumn,j,i] = np.nan

# corners to NaN
newvar[:, 1, -2] = np.nan
newvar[:, -2, -2] = np.nan
