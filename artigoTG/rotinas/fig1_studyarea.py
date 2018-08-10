#!/usr/bin/env python2
#-*-coding:utf-8-*-

import numpy as np # atribuir as arrays do .cdf para numpy's arrays
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import xray as xr

import cmocean as cmo

import pickle

from artigoTGpack import artigoTGpack as oceano

#
def draw_screen_poly( lats, lons, m,**kwargs):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), **kwargs)
    plt.gca().add_patch(poly)

# # diretorio base
BASE_DIR = oceano.make_dir()

SIMS_DIR = BASE_DIR.replace('github/', 'artigo_data/simulacoes/')
SIMS_DIR = '/media/danilo/Danilo/mestrado/artigo_data/simulacoes/2ndPhase/'

############### FIGURE PARAMETERS
figWdith  = 17.4/2.54
figHeight = 10./2.54


# importar batimetria direto de uma saida do modelo
#data = xr.open_dataset('/media/danilo/Danilo/mestrado/artigo_data/simulacoes/val3/gcmplt.cdf')
data = xr.open_dataset(SIMS_DIR+'reRUN00.cdf')
ncin = xr.open_dataset(SIMS_DIR+'depth_rerun00.cdf')
depth = ncin['depth'].values
lon = data['lon'].values
lat = data['lat'].values
lon[lon == 0.] = np.nan
lat[lon == 0.] = np.nan

# definindo coordenadas para plotar a região desejada
llon = -44.9;
ulon = -43.45;

llat = -23.6;
ulat = -22.75;

######### DETALHES PARA PLOTAR
squareWP      = [[-23.25,-22.89,-22.89,-23.25],[-44.72,-44.72,-44.29,-44.29]]
squareEP      = [[-23.15,-22.94,-22.94,-23.15],[-44.15,-44.15,-43.99,-43.99]]
porOeste      = [-44.72,-22.89]#[-44.508007, -23.194821]
porLeste      = [-43.99,-22.94]#[-44.113577, -23.070955]

## scale position
scaleLon      = [-44.65]
scaleLat      = [-23.55]

RibeiraArrow  = [-44.4,-23.,10000,13000,15000,15000]
JacuecangaArro= [-44.26,-23.02,4000,4000,5000,6000]
MambucabaArrow= [-44.525,-23.03, -2000,2000,-5000,5000]
MarambaiaArrow= [-43.95,-23.07,5000,-5000,7000,-8000]
GuaratibaArrow= [-43.561320, -23.064399, -1000, 1500,0,5000]
PicinguabaArro= [-44.827, -23.376, 2000,2800,2000,5000]
MamanguaArrow = [-44.644, -23.285, -5000,0,-9000,0]

CNAAAbase     = [-44.460789, -23.006711]
ilhagrande    = [-44.613300, -23.141362]#-23.137504, -44.613300
sepetiba      = [-43.864138, -23.003667]
piraquaraf    = [-44.440706, -23.010962]
paraty        = [-44.757599, -23.232831]
TermGuaiba = [-44.01916667, -23.00000000, '^', 'black']

### quadrado
lat0 = -22.980121
lat1 = -23.054853
lon0 = -44.374144
lon1 = -44.558726

fig = plt.figure(figsize=(figWdith,figHeight))
ax = fig.add_subplot(111)

# -----------------------------------------------------------------------------
# ---------------------    PLOTANDO A BASE DO MAPA   --------------------------
# -----------------------------------------------------------------------------
m = Basemap(projection='merc', llcrnrlat=llat, urcrnrlat=ulat, llcrnrlon=llon, urcrnrlon=ulon, resolution='f')
m.ax = ax
# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.1) #linha de costa em alta resolução
m.drawmapboundary(fill_color='#f2f2f2') # fill_color colore o oceano
# m.fillcontinents(color='#ffd480') # colorir o continente
m.fillcontinents(color='0.85')
m.drawrivers()
# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=8,
                                    color='gray',linewidth=0.1)
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=8,
                                    color='gray',linewidth=0.1)

m.drawmapscale(scaleLon[0], scaleLat[0], scaleLon[0], scaleLat[0], 20, barstyle='fancy', yoffset=1000)
# batimetria em forma de grade
pc = m.contourf(lon,lat,depth,latlon=True,cmap=cmo.cm.deep)

# -----------------------------------------------------------------------------
# ------------------------RETANGULOS DA DIVISAO--------------------------------
# -----------------------------------------------------------------------------
kwargs = {'linestyle': '--', 'alpha':0.1, 'color': '#800000','facecolor':'none','fill':"#800000",'linewidth':0.2}
draw_screen_poly(squareWP[0],squareWP[1],m,**kwargs)
kwargs = {'linestyle': '--', 'alpha':0.1, 'color': '#800000','facecolor':'none','fill':"#800000",'linewidth':0.2}
draw_screen_poly(squareEP[0],squareEP[1],m,**kwargs)

# -----------------------------------------------------------------------------
# --------------------------LOCALIZACOES EM TEXTO -----------------------------
# -----------------------------------------------------------------------------

########### physiographic classification
x,y=m(ilhagrande[0]+.100,ilhagrande[1]+.02)
ax.text(x,y,u'Ilha Grande\nBay',color='#000000', fontsize=9, ha='center',va='center')

x,y=m(sepetiba[0]+.080,sepetiba[1])
ax.text(x,y,u'Sepetiba\nBay',color='#000000', fontsize=9, ha='center',va='center')

x,y=m(porOeste[0],porOeste[1])
ax.text(x,y,u'Western Portion',color='#800000', fontsize=8, ha='left',va='baseline')

x,y=m(porLeste[0],porLeste[1])
ax.text(x,y,u'Eastern Portion',color='#800000', fontsize=8, ha='right',va='baseline')

x,y=m(-44.230458, -23.151578)
ax.text(x,y,u'Ilha Grande',color='k',fontsize=10,ha='center', va='center')

########### Important Locations
x,y=m(CNAAAbase[0], CNAAAbase[0])
m.plot(x,y, 'bo', markersize=8)


#### PRINCIPAL LOCATIONS
x,y=m(RibeiraArrow[0],RibeiraArrow[1])
ax.arrow(x,y,RibeiraArrow[2],RibeiraArrow[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+RibeiraArrow[4],y+RibeiraArrow[5],u'Ribeira Bay',fontsize=8,ha='center',va='center')

x,y=m(JacuecangaArro[0],JacuecangaArro[1])
ax.arrow(x,y,JacuecangaArro[2],JacuecangaArro[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+JacuecangaArro[4],y+JacuecangaArro[5],u'Jacuecanga\nBay',fontsize=8,ha='center',va='center')

x,y=m(MambucabaArrow[0],MambucabaArrow[1])
ax.arrow(x,y,MambucabaArrow[2],MambucabaArrow[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+MambucabaArrow[4],y+MambucabaArrow[5],u'Mambucaba\nRiver Mouth',fontsize=8,ha='center',va='center')

x,y=m(MarambaiaArrow[0],MarambaiaArrow[1])
ax.arrow(x,y,MarambaiaArrow[2],MarambaiaArrow[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+MarambaiaArrow[4],y+MarambaiaArrow[5],u'Marambaia\nBarrier',fontsize=8,ha='center',va='center')

x,y=m(GuaratibaArrow[0],GuaratibaArrow[1])
ax.arrow(x,y,GuaratibaArrow[2],GuaratibaArrow[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+GuaratibaArrow[4],y+GuaratibaArrow[5],u'Barra de\nGuaratiba',fontsize=8,ha='center',va='center')

x,y=m(PicinguabaArro[0],PicinguabaArro[1])
ax.arrow(x,y,PicinguabaArro[2],PicinguabaArro[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+PicinguabaArro[4],y+PicinguabaArro[5],u'Picinguaba',fontsize=8,ha='center',va='center')

x,y=m(MamanguaArrow[0],MamanguaArrow[1])
ax.arrow(x,y,MamanguaArrow[2],MamanguaArrow[3],head_width=1.5,head_length=0.99,fc='k',ec='k',lw=.5)
ax.text(x+MamanguaArrow[4],y+MamanguaArrow[5],u'Mamangua\nSac',fontsize=8,ha='center',va='center')

x,y=m(paraty[0]+0.01, paraty[1]+0.03)
ax.text(x,y,u'Paraty', color='k', fontsize=8, ha='center',va='center')

# plt.scatter(x,y,marker='o',color='k',alpha=.2,label='Direct Impact Area (20km)')

# -----------------------------------------------------------------------------
# ---------------------    MARKERS PARA LEGENDA         -----------------------
# -----------------------------------------------------------------------------
x,y = m(piraquaraf[0],piraquaraf[1])
plt.scatter(x,y,marker='*',color='black',s=20,label='Piraquara de Fora Inlet',zorder=2)

x,y = m(TermGuaiba[0], TermGuaiba[1])
plt.scatter(x,y, marker=TermGuaiba[2], color=TermGuaiba[3], s=15, label='Guaiba Island Terminal', zorder=2)

# -----------------------------------------------------------------------------
# -----------------------        ID RADIUS            -------------------------
# -----------------------------------------------------------------------------
# create legend with semi transparent background
lg = plt.legend(loc='upper right',fontsize=8,numpoints=1,scatterpoints=1)
lg.get_frame().set_alpha(.4)

cbar = plt.colorbar(pc, orientation='horizontal', shrink=0.625, aspect=20, fraction=0.2,pad=0.04)
cbar.set_label('Bathymetry (m)', fontsize=8)

# -----------------------------------------------------------------------------
# -----------------------    MINI GLOBE LOCATION      -------------------------
# -----------------------------------------------------------------------------
axin = inset_axes(m.ax, width='30%', height='30%', loc=4)
# inmap = Basemap(projection='ortho', lon_0=-44.5,lat_0=-25.5,ax=axin)
inmap = Basemap(projection='merc', llcrnrlat=-35,urcrnrlat=7,llcrnrlon=-77,urcrnrlon=-32,ax=axin)
inmap.drawmapboundary(fill_color='#f2f2f2')
inmap.drawcoastlines(linewidth=.1)
inmap.drawcountries(linewidth=.1,color='white')
inmap.fillcontinents(color='gray')

bx,by = inmap(m.boundarylons,m.boundarylats)
xy = list(zip(bx,by))
mapboundary = Polygon(xy,edgecolor='b',linewidth=1,fill=True)
inmap.ax.add_patch(mapboundary)

plt.tight_layout()
plt.subplots_adjust(top=.99,bottom=-0.045,right=0.993,left=0.0,wspace=0.)

plt.savefig('/media/danilo/Danilo/mestrado/github/artigoTG/figures/Fig1.png',dpi=600)

# plt.show()
#
# os.system('evince /media/danilo/Danilo/mestrado/github/artigoTG/figures/Fig1.eps')
