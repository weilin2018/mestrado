
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import netCDF4

import octant
import octant.roms

ncfilename = '/media/danilo/Danilo/mestrado/ventopcse/output/EC1.cdf'



nc = netCDF4.Dataset(ncfilename)

Nskip = 2

eke = []
mke = []
tke = []

ekep = []
mkep = []
tkep = []

dx = nc.variables['h1'][:]
dy = nc.variables['h2'][:]

salt0 = nc.variables['salt'][0]
ds = salt0.ptp()
s0 = salt0.max()
s_crit = s0 - 0.05*ds

tidx = 0

salt = nc.variables['salt'][tidx]
salt_bar = salt[:, :, 1:-1].mean(axis=-1)[:, :, None]

u = nc.variables['u'][tidx]
u_bar = u[:, :, :-1].mean(axis=-1)[:, :, None]
up = u - u_bar
vp = nc.variables['v'][tidx]
up = octant.tools.shrink(up, (21, 137, 110))
vp = octant.tools.shrink(vp, (21, 137, 110))

u = octant.tools.shrink(u, (21, 137, 110))

# zr = octant.roms.nc_depths(nc, grid='rho')[tidx]
# zw = octant.roms.nc_depths(nc, grid='w')[tidx]
zw = nc.variables['depth']
dz = np.diff(zw, axis=0)
# dV = (dx*dy*dz)[:, 1:-1, 1:-1]
dV = (dx[:-1,:]*dy[:-1,:]*dz)
dV_bar = dV[:, :, :-1].mean(axis=-1)[:, :, None]

eke.append( ((up**2 + vp**2)*dV).sum() / dV.sum() )
tke.append( ((u**2 + vp**2)*dV).sum() / dV.sum() )
mke.append( (u_bar[:, 1:-1, :]**2 *dV_bar).sum() / dV_bar.sum() )

idx = salt[:, 1:-1, 1:-1] < s_crit
ekep.append( ((up**2 + vp**2)*dV)[idx].sum() / dV[idx].sum() )
tkep.append( ((u**2 + vp**2)*dV)[idx].sum() / dV[idx].sum() )

idx = salt_bar[:, 1:-1, :] < s_crit
mkep.append( (u_bar[:, 1:-1, :]**2 *dV_bar)[idx].sum() / dV_bar[idx].sum() )
