clear all; close all; clc;

# entering in the working directory
cd /media/danilo/Danilo/mestrado/ventopcse/data/mercator/

addpath(genpath('/home/danilo/Downloads/delft/delft3d_repository/src/third_party_open/netcdf/matlab/mexnc/'))
rehash


xt = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','longitude');
yt = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','latitude');
zt = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','depth');
t  = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','time');

temp = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','temperature');
u    = ncread('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/data.nc','u');

[nt, nzt, nyt, nxt] = size(temp);
[ntu, nzu, nyu, nxu] = size(u);

load('/media/danilo/Danilo/mestrado/ventopcse/data/mercator/MERC_grid_bounds');
st_bnds_u2 = flipud(st_bnds);
for iz = 1:nzu-1; DZ(iz) = diff(st_bnds_u2(iz,:)); end;


%%Create land mask
temp2   = squeeze(temp(1,:,:,:));
u2      = squeeze(u(1,:,:,:));
mask_ts = ones(nzt,nyt,nxt);
mask_uv = ones(nzu,nyu,nxu);
%Note that u,v,temp are NaN over land.
for iz = 1:nzu; for iy = 1:nyu; for ix = 1:nxu;
if isnan(temp2(iz,iy,ix)) == 1; mask_ts(iz,iy,ix) = NaN; end;
if isnan(u2(iz,iy,ix)) == 1; mask_u(iz,iy,ix) = NaN; end;
end; end; end;


%save grid here....
MERC_grid.z_bnds = st_bnds_u2;
MERC_grid.xt = xt;
MERC_grid.yt = yt;
MERC_grid.zt = zt;
MERC_grid.xu = xu;
MERC_grid.yu = yu;
MERC_grid.zu = zu;
MERC_grid.mask_ts = mask_ts; %1 over water, NaN over land
MERC_grid.mask_uv = mask_uv;

save('MERC_grid_data_SWBrazil.mat', 'MERC_GRID');
