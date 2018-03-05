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

# importar batimetria direto de uma saida do modelo
#data = xr.open_dataset('/media/danilo/Danilo/mestrado/artigo_data/simulacoes/val3/gcmplt.cdf')
data = xr.open_dataset('/home/tparente/danilo/mestrado/artigo_data/simulacoes/val3/gcmplt.cdf')
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
CNAAAbase  = [-44.460789, -23.006711]
ilhagrande = [-44.613300, -23.141362]#-23.137504, -44.613300
sepetiba   = [-43.864138, -23.003667]
piraquaraf = [-44.441523, -23.011983]
itaorna    = [-44.472463, -23.006116]#-23.017116,
porOeste   = [-44.508007, -23.194821]
porLeste   = [-44.113577, -23.070955]
canCentral = [-44.333599, -23.083153]

marambaia     = [-43.989471, -23.083472]
mangaratiba   = [-44.056691, -22.944938]
paraty        = [-44.757599, -23.232831]
angradosreis  = [-44.358047, -23.020794]

# tarituba      = [-44.598025, -23.039031]
## pontos de descarga/captação
fozMambucaba  = [-44.521247, -23.026884]
descargaCNAAA = [-44.445079, -23.012130]
captacaoCNAAA = [-44.459445, -23.009400]

# séries temporais

loc_w210 = [[-23.010962, -44.440706], [-22.9741811,-44.323831], [-23.012479, -44.307350], [-23.059381, -44.241407], [-22.994208,-44.039795],[-23.088036,-44.098059],[-23.096075,-44.010663]]
nam_w210 = ['Piraquara de Fora','Angra dos Reis (A)', 'Angra dos Reis (B)', 'Canal Central', 'Mangaratiba','Porção Leste (BIG)', 'Marambaia']

loc_w50 = [[-23.010962, -44.440706], [-22.960510, -44.394883],[-22.9741811,-44.323831],[-23.034976, -44.525694],[-23.213933, -44.704401]]
nam_w50 = ['Piraquara de Fora','Ilha Comprida', 'Angra dos Reis (A)', 'Rio Mambucada', 'Paraty']

loc_wcptec = [[-23.010962, -44.440706], [-22.960510, -44.394883],[-22.9741811,-44.323831],[-23.012479, -44.307350], [-23.059381, -44.241407], [-22.994208,-44.039795],[-23.088036,-44.098059],[-23.096075,-44.010663]]
nam_wcptec = ['Piraquara de Fora', 'Ilha Comprida', 'Angra dos Reis (A)', 'Angra dos Reis (B)', 'Canal Central', 'Porção Leste', 'Mangaratiba', 'Marambaia']


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
m=pickle.load(open("/home/tparente/danilo/mestrado/github/artigoTG/rotinas/bigmap.p","rb"))
# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.4) #linha de costa em alta resolução
m.drawmapboundary() # fill_color colore o oceano
m.fillcontinents(color='#ffd480') # colorir o continente
m.drawrivers()
# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.5)
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray',linewidth=0.5)

# batimetria em forma de grade
m.pcolor(lon,lat,depth,latlon=True,cmap=cmo.cm.deep)

# plotar detalhes
x,y=m(CNAAAbase[0], CNAAAbase[0])
m.plot(x,y, 'bo', markersize=14)

x,y=m(ilhagrande[0],ilhagrande[1])
ax.text(x,y,u'ILHA GRANDE\nBAY',color='#000000', fontsize=15)

x,y=m(sepetiba[0],sepetiba[1])
ax.text(x,y,u'SEPETIBA\nBAY',color='#000000', fontsize=15)

x,y=m(porOeste[0],porOeste[1])
ax.text(x,y,u'WEST\nPORTION',color='#800000', fontsize=12)

x,y=m(porLeste[0],porLeste[1])
ax.text(x,y,u'EAST\nPORTION',color='#800000', fontsize=12,horizontalalignment='center',verticalalignment='center')

x,y=m(canCentral[0],canCentral[1])
ax.text(x,y,u'CENTRAL CHANNEL',color='#800000', fontsize=12)


#### desenhar quadrado
# lats = [lat0,lat1,lat1,lat0]
# lons = [lon0,lon1,lon1,lon0]
# draw_square(lats,lons,m)

loc_w210 = [[-23.010962, -44.440706], [-22.9741811,-44.323831], [-23.012479, -44.307350], [-23.059381, -44.241407], [-22.994208,-44.039795],[-23.088036,-44.098059],[-23.096075,-44.010663]]
nam_w210 = ['Piraquara de Fora','Angra dos Reis (A)', 'Angra dos Reis (B)', 'Canal Central', 'Mangaratiba','Porção Leste (BIG)', 'Marambaia']

loc_w50 = [[-22.960510, -44.394883],[-23.026489, -44.521393],[-23.213933, -44.704401]]
nam_w50 = ['Ilha Comprida', 'Rio Mambucada', 'Paraty']

x,y = m(loc_w210[0][1],loc_w210[0][0])
plt.scatter(x,y,marker='*',color='black',s=100,label=nam_w210[0])

x,y = m(loc_w210[1][1],loc_w210[1][0])
plt.scatter(x,y,marker='*',color='fuchsia',s=100,label=nam_w210[1])

# x,y = m(loc_w210[2][1],loc_w210[2][0])
# plt.scatter(x,y,marker='*',color='sienna',s=100,label=nam_w210[2])

# x,y = m(loc_w210[3][1],loc_w210[3][0])
# plt.scatter(x,y,marker='*',color='indigo',s=100,label=nam_w210[3])

# x,y = m(loc_w210[4][1],loc_w210[4][0])
# plt.scatter(x,y,marker='*',color='black',s=100,label=nam_w210[4])

x,y = m(loc_w210[5][1],loc_w210[5][0])
plt.scatter(x,y,marker='*',color='tan',s=100,label='Porcao Leste (BIG)')

x,y = m(loc_w210[6][1],loc_w210[6][0])
plt.scatter(x,y,marker='*',color='darkmagenta',s=100,label=nam_w210[6])


x,y = m(loc_w50[0][1],loc_w50[0][0])
plt.scatter(x,y,marker='*',color='lightcoral',s=100,label=nam_w50[0])

x,y = m(loc_w50[1][1],loc_w50[1][0])
plt.scatter(x,y,marker='.',color='darkorange',s=100,label=nam_w50[1])

x,y = m(loc_w50[2][1],loc_w50[2][0])
plt.scatter(x,y,marker='*',color='seagreen',s=100,label=nam_w50[2])

x,y = m(TermGuaiba[0], TermGuaiba[1])
plt.scatter(x,y, marker=TermGuaiba[2], color=TermGuaiba[3], s=40, label='Guaiba Island')

x,y = m(ColNaval[0], ColNaval[1])
plt.scatter(x,y, marker=ColNaval[2], color=ColNaval[3], s=40, label='Naval College')


kwargs = {'linestyle': '--', 'alpha':0.4, 'color': 'black'}
equi(m, CNAAAbase[0], CNAAAbase[1], 15.,lw=1., **kwargs)
# create legend with semi transparent background
lg = plt.legend(loc='upper left',fontsize=12,numpoints=1,scatterpoints=1)
lg.get_frame().set_alpha(.4)


plt.show()
