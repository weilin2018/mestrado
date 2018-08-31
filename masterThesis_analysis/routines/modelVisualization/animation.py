import matplotlib.pyplot as plt
import numpy as np
import cmocean as cmo

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


class Animation:

    def __init__(self,var='elev'):
        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart:timeEnd,:,:,:].values
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart:self.timeEnd,:,:].values

    def select_colormap(self):
        if self.varname == 'elev':
            cmap = 'RdBu_r'
        elif self.varname == 'temp':
            cmap = cmo.cm.thermal
        elif self.varname == 'salt':
            cmap = cmo.cm.haline

        return cmap

    def field(self,**kwargs):
        plt.ion()

        # cmap = self.select_colormap()

        fig,ax = plt.subplots()
        for t in np.arange(self.timeStart,self.timeEnd,1):
            ax.clear()

            data = self.var[t]
            m = oceano.make_map(ax)

            m.contourf(self.lon,self.lat,data,**kwargs)
            plt.pause(0.5)
