import xray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap

plt.ion() 


def make_map(ax,llat=-30,ulat=-20,llon=-50,ulon=-39,resolution='l',nmeridians=3,nparallels=2):

    m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution=resolution)

    m.ax = ax

#    m.drawcoastlines(linewidth=0.1)
#    m.drawmapboundary()
#    m.fillcontinents(color='#c0c0c0')

    m.drawcoastlines(linewidth=.5)

    return m


def transect(ilat,ilon,flat,flon,data):
    d = 

DATA_DIR = '/media/danilo/Danilo/mestrado/ventopcse/data/HeatBudget/'
ncin = xr.open_dataset(DATA_DIR+'mercator_2006.nc')

time= 0

lon = ncin.longitude.values
lat = ncin.latitude.values
dep = ncin.depth.values # selecionando ate a profundidade 11 (que deve ser algo em torno de 20m)
temp = ncin.temperature[time,:,:,:]

x,y = np.meshgrid(lon,lat)

# ponto do comeco do transecto
ilat = 22
ilon = 11
# ponto do final do transecto
flat = ilat - 10
flon = ilon + 10

# recortando
newT = temp[:,flat:ilat,ilon:flon]
newlo= x[flat:ilat,ilon:flon]
newla= y[flat:ilat,ilon:flon]

# visualizando o novo campo e aonde sera feito o transecto
fig,ax = plt.subplots()
m = make_map(ax)
m.contourf(x[flat:ilat,ilon:flon],y[flat:ilat,ilon:flon],temp[0,flat:ilat,ilon:flon],latlon=True)
m.plot(x[flat:ilat,ilon:flon],y[flat:ilat,ilon:flon],'k',alpha=.3,latlon=True)
m.plot(x[flat:ilat,ilon:flon].T,y[flat:ilat,ilon:flon].T,'k',alpha=.3,latlon=True)

jindexes = np.arange(9,-1,-1)
iindexes = np.arange(0,10,1)
for j,i in zip(jindexes,iindexes):
    m.scatter(newlo[j,i],newla[j,i],s=30,c='r',latlon=True)


# criando a martriz do transecto
transect = np.zeros(newT.shape[0:2])*np.nan

for k in range(newT.shape[0]):
    for j,i in zip(jindexes,iindexes):
        transect[k,i] = newT[k,j,i]

# visualizando a secao vertical
fig,ax = plt.subplots()
ax.contourf(transect)
ax.invert_yaxis()

# linha do transecto no mapa
fig,ax = plt.subplots()
m = make_map(ax)
m.plot(newlo[jindexes,iindexes],newla[jindexes,iindexes],'r',latlon=True)

