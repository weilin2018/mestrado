from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

class mapa(object):
    """Those examples of objects meridional and zonal limits"""
    def __init__(self):
        self.lon0 = None
        self.lat0 = None
        self.lonF = None
        self.latF = None

    def Mapa(self,ax):
        self.mapa =  Basemap(projection='mill',lat_ts=10,llcrnrlon=self.lon0, urcrnrlon=self.lon1, llcrnrlat=self.lat0, urcrnrlat=self.lat1, resolution=self.res)
        self.mapa.ax = ax

        self.mapa.drawparallels(self.parallel,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
        self.mapa.drawmeridians(self.meridian,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')

        self.mapa.drawcoastlines(linewidth=0.1)
        self.mapa.drawmapboundary()
        self.mapa.fillcontinents(color='#c0c0c0')

        self.mapa.drawcoastlines(linewidth=.1)
        self.mapa.drawmapboundary()

    def pcse(self):
        self.lat0 = -30
        self.lat1 = -20
        self.lon0 = -50
        self.lon1 = -40
        self.meridian = np.arange(self.lon0,self.lon1,3)
        self.parallel = np.arange(self.lat0,self.lat1,2)
        self.res  = 'i'

    def Araca(self):
        self.lat0 = -23.9
        self.lat1 = -23.7
        self.lon0 = -45.5
        self.lon1 = -45.3
        self.res  = 'f'
        self.parallel = [-23.75,-23.85]
        self.meridian = [-45.45,-45.35]
        # self.Mapa()

    def Canal_ssb(self):
        self.lat0 = -23.9
        self.lat1 = -23.7
        self.lon0 = -45.5
        self.lon1 = -45.3
        self.res  = 'i'
        self.parallel = [-23.75,-23.85]
        self.meridian = [-45.45,-45.35]
        # self.Mapa()

    def Embaiamento_sp(self):
        #boundaries - mill
        self.lat0 = -30
        self.lat1 = -21
        self.lon0 = -50
        self.lon1 = -38
        #basemap
        self.res  = 'i'
        self.parallel = [-27,-24]
        self.meridian = [-47,-43]
        # self.Mapa()

    def Sp_coast(self):
        #boundaries - mill
        self.lat0 = -26.5
        self.lat1 = -22.5
        self.lon0 = -49
        self.lon1 = -43.5
        #basemap
        self.res  = 'i'
        self.parallel = [-25,-23]
        self.meridian = [-48,-46,-44]
        # self.Mapa()

    def Tropical_atlantic(self):
        #boundaries - mill
        self.lat0 = -20
        self.lat1 = 25
        self.lon0 = -60
        self.lon1 = 15
        #basemap
        self.res  = 'l'
        self.parallel = [-20,0,20]
        self.meridian = [-55,-30,15]
        # self.Mapa()

    def Eq_south_atlantic(self):
        self.lat0 = -30
        self.lat1 = 15
        self.lon0 = -60
        self.lon1 = 15
        #basemap
        self.res  = 'l'
        self.parallel = [-40,-20,0,20]
        self.meridian = [-45,-15,15]
        # self.Mapa()

class location(object):
    # this class define coordinates location to plot graphics in
    # many places over the PCSE model_grid

    def __init__(self):
        self.ilon = None
        self.ilat = None

    def SBC(self):
        # sao sebastiao channel
        self.ilat = 55
        self.ilon = 7

    def cananeia(self):
        self.ilat = 19
        self.ilon = 70

    def santos(self):
        self.ilat = 28
        self.ilon = 70

    def ubatuba(self):
        self.ilat = 99
        self.ilon = 70

    def cabofrio(self):
        self.ilat = 126
        self.ilon = 70
