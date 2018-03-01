#!/usr/bin/env python2
#-*-coding:utf-8-*-

import numpy as np # atribuir as arrays do .cdf para numpy's arrays
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset, inset_axes
import cmocean as cmo

import pickle

def draw_square(lats,lons,m):
    x,y=m(lons,lats)
    xy=zip(x,y)
    poly = Polygon(xy,facecolor='blue', alpha=0.4, edgecolor='black')
    plt.gca().add_patch(poly)

# definindo coordenadas para plotar a região desejada
llon = -44.9;
ulon = -43.45;

llat = -23.6;
ulat = -22.75;

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
m=pickle.load(open("/home/danilo/Dropbox/0TG_DANILO/rotinas/rotinas/bigmap.p","rb"))

# plotar outras coisas do mapa
m.drawcoastlines(linewidth=0.4) #linha de costa em alta resolução
m.drawmapboundary(fill_color='#e5f2ff') # fill_color colore o oceano
m.fillcontinents(color='#ffd480') # colorir o continente

# definir meridianos e paralelos para plotar no mapa
meridians=np.arange(-44.9,-33.4,0.25)
parallels=np.arange(-23.6,-22.75,0.1)
# desenhar meridianos e paralelos conforme definido acima
m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

# plotar detalhes
x,y=m(CNAAAbase[0], CNAAAbase[0])
m.plot(x,y, 'bo', markersize=14)

x,y=m(ilhagrande[0],ilhagrande[1])
ax.text(x,y,u'Ilha Grande\n Bay',color='#000000', fontsize=15)

x,y=m(sepetiba[0],sepetiba[1])
ax.text(x,y,u'Sepetiba\n Bay',color='#000000', fontsize=15)

x,y=m(porOeste[0],porOeste[1])
ax.text(x,y,u'West \n Portion',color='#800000', fontsize=12)

x,y=m(porLeste[0],porLeste[1])
ax.text(x,y,u'East \n Portion',color='#800000', fontsize=12)

x,y=m(canCentral[0],canCentral[1])
ax.text(x,y,u'Central Channel',color='#800000', fontsize=12)


#### desenhar quadrado
# lats = [lat0,lat1,lat1,lat0]
# lons = [lon0,lon1,lon1,lon0]
# draw_square(lats,lons,m)

loc_w210 = [[-23.010962, -44.440706], [-22.9741811,-44.323831], [-23.012479, -44.307350], [-23.059381, -44.241407], [-22.994208,-44.039795],[-23.088036,-44.098059],[-23.096075,-44.010663]]
nam_w210 = ['Piraquara de Fora','Angra dos Reis (A)', 'Angra dos Reis (B)', 'Canal Central', 'Mangaratiba','Porção Leste (BIG)', 'Marambaia']

loc_w50 = [[-22.960510, -44.394883],[-23.034976, -44.525694],[-23.213933, -44.704401]]
nam_w50 = ['Ilha Comprida', 'Rio Mambucada', 'Paraty']

x,y = m(loc_w210[0][1],loc_w210[0][0])
plt.scatter(x,y,marker='*',color='black',s=100,label=nam_w210[0])

x,y = m(loc_w210[1][1],loc_w210[1][0])
plt.scatter(x,y,marker='*',color='fuchsia',s=100,label=nam_w210[1])

x,y = m(loc_w210[2][1],loc_w210[2][0])
plt.scatter(x,y,marker='*',color='sienna',s=100,label=nam_w210[2])

x,y = m(loc_w210[3][1],loc_w210[3][0])
plt.scatter(x,y,marker='*',color='indigo',s=100,label=nam_w210[3])

x,y = m(loc_w210[4][1],loc_w210[4][0])
plt.scatter(x,y,marker='*',color='black',s=100,label=nam_w210[4])

x,y = m(loc_w210[5][1],loc_w210[5][0])
plt.scatter(x,y,marker='*',color='tan',s=100,label='Porcao Leste (BIG)')

x,y = m(loc_w210[6][1],loc_w210[6][0])
plt.scatter(x,y,marker='*',color='darkmagenta',s=100,label=nam_w210[6])


x,y = m(loc_w50[0][1],loc_w50[0][0])
plt.scatter(x,y,marker='*',color='lightcoral',s=100,label=nam_w50[0])

x,y = m(loc_w50[1][1],loc_w50[1][0])
plt.scatter(x,y,marker='*',color='darkorange',s=100,label=nam_w50[1])

x,y = m(loc_w50[2][1],loc_w50[2][0])
plt.scatter(x,y,marker='*',color='seagreen',s=100,label=nam_w50[2])

# create legend with semi transparent background
# lg = plt.legend(loc='upper left',fontsize=12,numpoints=1,scatterpoints=1)
# lg.get_frame().set_alpha(.4)


plt.show()
