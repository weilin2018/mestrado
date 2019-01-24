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

# importando arquivo
DATA_DIR = '/home/danilo/mestrado/arian/'
fname = DATA_DIR + 'CARS_2009_WORLD_MONTHLY_01.nc'

ncin = xr.open_dataset(fname,decode_times=False)

# extraindo uma variavel qualquer para trabalhar
lon = ncin.LONGITUDE.values
lon = ((lon - 180)%360) - 180 # convertendo de 0/360 para -180/180
lat = ncin.LATITUDE.values
dep = ncin.DEPTH.values
data = ncin.TEMP[0,:,:,:] # tempo 0, todas as profs, lons e lats

# selecionando o ponto de inicio do transecto: eu sei pq descobri olhando pros
# dados de lon e lat e encontrando os pontos que eu queria (np.where(), por exemplo)
ilat,flat = 103, 92
ilon,flon = 626, 637

# recortando: aqui tem um pulo do gato
ndata = data[:,flat:ilat,ilon:flon]
nlon1D= lon[ilon:flon]
nlat1D= lat[flat:ilat]

# gridando a matriz de coordenadas
nlon,nlat = np.meshgrid(nlon1D,nlat1D)

# visualizando o novo campo e aonde sera feito o transecto
fig,axes = plt.subplots(ncols=2)
ax=axes[0]
m = make_map(ax)
m.contourf(xplot,yplot,data[0,:,:],latlon=True)
plt.title('Grade Completa')

ax=axes[1]
m = make_map(ax)
m.contourf(xplot[flat:ilat,ilon:flon],yplot[flat:ilat,ilon:flon],data[0,flat:ilat,ilon:flon],latlon=True)
m.plot(xplot[flat:ilat,ilon:flon],yplot[flat:ilat,ilon:flon],'k',alpha=.3,latlon=True)
m.plot(xplot[flat:ilat,ilon:flon].T,yplot[flat:ilat,ilon:flon].T,'k',alpha=.3,latlon=True)
plt.title('Grade Subamostrada')

jindexes = np.arange(10,-1,-1)
iindexes = np.arange(0,11,1)
for j,i in zip(jindexes,iindexes):
    m.scatter(nlon[j,i],nlat[j,i],s=30,c='r',latlon=True)

# criando a martriz do transecto
transect = np.zeros(ndata.shape[0:2])*np.nan

for k in range(ndata.shape[0]):
    for j,i in zip(jindexes,iindexes):
        transect[k,i] = ndata[k,j,i]

# visualizando a secao vertical
fig,ax = plt.subplots()
ax.contourf(transect)
ax.set_ylim([0,3])
ax.invert_yaxis()

# linha do transecto no mapa
fig,ax = plt.subplots()
m = make_map(ax)
m.plot(nlon[jindexes,iindexes],nlat[jindexes,iindexes],'r',latlon=True)


# para visualizar as profundidades no eixo y, precisamos criar uma nota matriz pra bater com o
# shape da matriz transect que criamos acima
ndep = np.tile(dep[:3],(11,1))
l    = np.tile(nlat1D,(3,1))
# ainda, pra ficar mais legal, queremos colocar no eixo X a distancia do transecto desde o ponto
# de referencia inicial (ou desde a costa)
dist = [0]
for lo,la in zip(nlon1D,nlat1D[::-1]):
    dist.append(
inicio = [nlat1D[-1],nlon1D[0]]
