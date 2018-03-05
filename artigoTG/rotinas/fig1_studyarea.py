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

def draw_square(lats,lons,m):
    x,y=m(lons,lats)
    xy=zip(x,y)
    poly = Polygon(xy,facecolor='blue', alpha=0.4, edgecolor='black')
    plt.gca().add_patch(poly)

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])

    #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    X,Y = m(X,Y)
    plt.plot(X,Y,**kwargs)

#
# # diretorio base
# pwd = os.getcwd() # get current directory
# BASE_DIR = pwd.split('/')[:-3]
#
# BASE_DIR = '/'.join((BASE_DIR))
BASE_DIR = oceano.make_dir()

SIMS_DIR = BASE_DIR.replace('github/', 'artigo_data/simulacoes/')


# importar batimetria direto de uma saida do modelo
#data = xr.open_dataset('/media/danilo/Danilo/mestrado/artigo_data/simulacoes/val3/gcmplt.cdf')
data = xr.open_dataset(SIMS_DIR+'val3/gcmplt.cdf')
depth = data['depth'].values
lon = data['lon'].values
lat = data['lat'].values
lon[lon == 0.] = np.nan
lat[lon == 0.] = np.nan

# definindo coordenadas para plotar a região desejada
llon = -44.9;
ulon = -43.45;

llat = -23.6;
ulat = -22.75;

TermGuaiba = [-44.01916667, -23.00000000, '^', 'black']
ColNaval  =  [-44.31833333, -23.01694444, 'v', 'black']

######### DETALHES PARA PLOTAR
CNAAAbase     = [-44.460789, -23.006711]
ilhagrande    = [-44.613300, -23.141362]#-23.137504, -44.613300
sepetiba      = [-43.864138, -23.003667]
piraquaraf    = [-44.441523, -23.011983]
itaorna       = [-44.472463, -23.006116]#-23.017116,
porOeste      = [-44.508007, -23.194821]
porLeste      = [-44.113577, -23.070955]
canCentral    = [-44.333599, -23.083153]
picinguaba    = [-44.830344, -23.378115]
guaratiba     = [-43.561320, -23.064399]
riomambuca    = [-44.521393, -23.026489]
angradosreis  = [-44.316713, -22.986787]

marambaia     = [-43.989471, -23.083472]
mangaratiba   = [-44.056691, -22.944938]
paraty        = [-44.757599, -23.232831]
angradosreis  = [-44.358047, -23.020794]

# tarituba      = [-44.598025, -23.039031]
## pontos de descarga/captação
fozMambucaba  = [-44.521247, -23.026884]
descargaCNAAA = [-44.445079, -23.012130]
captacaoCNAAA = [-44.459445, -23.009400]

### quadrado
lat0 = -22.980121
lat1 = -23.054853
lon0 = -44.374144
lon1 = -44.558726

fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(111)

# Plotando a base do mapa
# m=Basemap(projection='merc',llcrnrlat=llat,urcrnrlat=ulat,llcrnrlon=llon,urcrnrlon=ulon,resolution='f')

# pra agilizar o processo, puxa-se um mapa já feito e salvo em pickle
#m=pickle.load(open("/home/danilo/Dropbox/0TG_DANILO/rotinas/rotinas/bigmap.p","rb"))
m=pickle.load(open(BASE_DIR+"artigoTG/rotinas/bigmap.p","rb"))
# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.2) #linha de costa em alta resolução
m.drawmapboundary() # fill_color colore o oceano
m.fillcontinents(color='#ffd480') # colorir o continente
m.drawrivers()
# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,
                                    fontweight='bold',color='gray',linewidth=0.5)
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,
                                    fontweight='bold',color='gray',linewidth=0.5)

m.drawmapscale(-43.65, -23.55, -43.65, -23.55, 20, barstyle='fancy', yoffset=1000)
# batimetria em forma de grade
m.pcolor(lon,lat,depth,latlon=True,cmap=cmo.cm.deep)

# plotar detalhes em texto
x,y=m(CNAAAbase[0], CNAAAbase[0])
m.plot(x,y, 'bo', markersize=14)

x,y=m(ilhagrande[0]+.100,ilhagrande[1]+.02)
ax.text(x,y,u'ILHA GRANDE\nBAY',color='#000000', fontsize=15, ha='center',va='center')

x,y=m(sepetiba[0]+.080,sepetiba[1])
ax.text(x,y,u'SEPETIBA\nBAY',color='#000000', fontsize=15, ha='center',va='center')

x,y=m(porOeste[0],porOeste[1])
ax.text(x,y,u'WEST\nPORTION',color='#800000', fontsize=12, ha='center',va='center')

x,y=m(porLeste[0],porLeste[1])
ax.text(x,y,u'EAST\nPORTION',color='#800000', fontsize=12, ha='center',va='center')

x,y=m(canCentral[0],canCentral[1])
ax.text(x,y,u'CENTRAL CHANNEL',color='#800000', fontsize=12, ha='center',va='center')

x,y=m(picinguaba[0]+0.05, picinguaba[1]+0.015)
ax.text(x,y,u'Picinguaba',color='k',fontsize=8, ha='center',va='center')

x,y=m(guaratiba[0]+0.035, guaratiba[1]+0.025)
ax.text(x,y,u'Barra de\nGuaratiba',color='k',fontsize=8, ha='center',va='center')

x,y=m(riomambuca[0], riomambuca[1]+0.01)
ax.text(x,y,u'Mambucaba\nRiver',color='k',fontsize=10, ha='center',va='center')

x,y=m(angradosreis[0]+0.01, angradosreis[1]+0.03)
ax.text(x,y,u'Angra dos\nReis', color='k', fontsize=11, ha='center',va='center')

x,y=m(paraty[0]+0.01, paraty[1]+0.03)
ax.text(x,y,u'Paraty', color='k', fontsize=12, ha='center',va='center')

x,y=m(-44.230458, -23.151578)
ax.text(x,y,u'Ilha Grande',color='k',fontsize=12,ha='center', va='center')

# plt.scatter(x,y,marker='o',color='k',alpha=.2,label='Direct Impact Area (20km)')

# -----------------------------------------------------------------------------

loc_w210 = [[-23.010962, -44.440706], [-22.9741811,-44.323831],
            [-23.012479, -44.307350], [-23.059381, -44.241407],
            [-22.994208,-44.0397950], [-23.088036,-44.0980590],
            [-23.096075,-44.010663]]
nam_w210 = ['Piraquara de Fora Inlet','Angra dos Reis (A)',
            'Angra dos Reis (B)', 'Canal Central',
            'Mangaratiba','Porção Leste (BIG)',
            'Marambaia']

x,y = m(loc_w210[0][1],loc_w210[0][0])
plt.scatter(x,y,marker='*',color='black',s=80,label=nam_w210[0],zorder=2)

x,y = m(TermGuaiba[0], TermGuaiba[1])
plt.scatter(x,y, marker=TermGuaiba[2], color=TermGuaiba[3], s=60, label='Guaiba Island Terminal', zorder=2)

x,y = m(ColNaval[0], ColNaval[1])
plt.scatter(x,y, marker=ColNaval[2], color=ColNaval[3], s=60, label='Naval College', zorder=2)


kwargs = {'linestyle': '--', 'alpha':0.4, 'color': 'black'}
equi(m, CNAAAbase[0], CNAAAbase[1], 15.,lw=1., **kwargs)
# create legend with semi transparent background
lg = plt.legend(loc='upper left',fontsize=12,numpoints=1,scatterpoints=1)
lg.get_frame().set_alpha(.4)


plt.show()
