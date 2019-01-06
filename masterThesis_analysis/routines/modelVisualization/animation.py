import matplotlib.pyplot as plt
import numpy as np
import cmocean as cmo
import decomp

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano


def rotate_velocityField(u,v,ang):

    ur = np.zeros(u.shape)*np.nan
    vr = np.zeros(v.shape)*np.nan

    for j in range(u.shape[0]):
        U,V = u[j,:].values,v[j,:].values
        angle = ang[j,:]

        INT,DIR = decomp.uv2intdir(U,V,0,angle)
        uro,vro = decomp.intdir2uv(INT,DIR,0,angle)
        ur[j,:] = uro
        vr[j,:] = vro

    return ur,vr

def formatGrid_plot(grid,fname):
    import numpy as np
    ij=np.load(fname)
    # for a 2D array (lon,lat)
    if len(grid.shape)==2:
        grid=grid[ij[1], :]
        grid=grid[:, ij[0]]
    # if grid is a 3D array (temp,salt,speed)
    if len(grid.shape)==3:
        grid=grid[:,ij[1], ij[0]]
    return grid

class Animation:

    def __init__(self,var='elev',sigma=0):
        plt.ion()

        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,:,:]
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),:,:]

    def select_colormap(self):
        if self.varname == 'elev':
            cmap = 'RdBu_r'
        elif self.varname == 'temp':
            cmap = cmo.cm.thermal
        elif self.varname == 'salt':
            cmap = cmo.cm.haline

        return cmap

    def velocidade(self,index_file,intervalo=.2,sigma=0,wind=True):

        # checar se existe u e v no objeto
        if ~hasattr(self,'u'):
            self.u = self.ncin['u'][self.timeStart.item():self.timeEnd.item(),sigma,:,:]
            self.v = self.ncin['v'][self.timeStart.item():self.timeEnd.item(),sigma,:,:]

        # preparar as variaveis para plot
        xplot = oceano.formatGrid_plot(self.lon,index_file)
        yplot = oceano.formatGrid_plot(self.lat,index_file)
        contour_levels = np.arange(0,2,0.2)
        # animacao
        fig,ax = plt.subplots()
        for t in np.arange(self.timeStart.item(),self.timeEnd.item(),2):
            ax.clear()
            m = oceano.make_map(ax)
            plt.suptitle('Velocity at: '+ str(self.ncin.time[t].values))

            # rotate vectors based on cell's angle
            ur,vr = rotate_velocityField(self.u[t,:,:],self.v[t,:,:],self.ncin.ang.values)
            s = np.sqrt(ur**2 + vr**2)

            un = ur/s
            vn = vr/s

            uplot = oceano.formatGrid_plot(un[:,:],index_file)
            vplot = oceano.formatGrid_plot(vn[:,:],index_file)
            wuplot = formatGrid_plot(self.ncin.wu[t,:,:],index_file)
            wvplot = formatGrid_plot(self.ncin.wv[t,:,:],index_file)

            m.contourf(self.lon,self.lat,s[:,:],contour_levels,latlon=True,
                        cmap=cmo.cm.speed)
            m.quiver(xplot[:,::3],yplot[:,::3],uplot[:,::3],vplot[:,::3],
                     scale=80,width=0.001,headwidth=4,headlength=4,latlon=True,
                     pivot='middle')
            if wind:
                m.quiver(xplot[::2,::4],yplot[::2,::4],wuplot[::2,::4],
                         wvplot[::2,::4],latlon=True,color='k',alpha=.4,
                         pivot='middle',headwidth=4,headlength=4,minshaft=2)

            plt.pause(intervalo)

    def field(self,title='',**kwargs):

        # in case user don't send a cmap, we use the variable name
        # to select a proper colormap
        if 'cmap' not in kwargs.keys():
            kwargs['cmap'] = self.select_colormap()

        # we create contour_levels outside the for loop, so we can
        # use the extreme values in the entire time domain
        # dif = np.abs(np.nanmin(self.var.values)-np.nanmax(self.var.values))
        # contour_levels = np.arange(np.nanmin(self.var.values),np.nanmax(self.var.values),round(dif/10.))

        if (self.var.ndim > 2):
            # if self.var is a 3D array, plot an animation
            plt.ion()

            fig,ax = plt.subplots()
            for t in np.arange(self.timeStart.item(),self.timeEnd.item(),1):
                ax.clear()

                data = self.var[t]
                # m = oceano.make_map(ax)
                self.Mapa(ax)

                self.mapa.contourf(self.lon,self.lat,data,**kwargs)
                plt.title(title + '\n' +str(self.ncin.time[t].values))
                plt.pause(0.5)
        else:
            # if self.var is a 2D array, plot a static map
            plt.ion()

            fig,self.ax = plt.subplots()
            self.Mapa(self.ax)
            self.csf = self.mapa.contourf(self.lon,self.lat,self.var,**kwargs)
            plt.colorbar(self.csf)
            plt.title(title,fontsize=12)


    def cross_section(self,**kwargs):
        pass


class Visualization:

    def __init__(self,var='elev',sigma=0):
        self.varname = var

        if len(self.ncin[var].shape) == 4:
            # is a 3D var (temp,salt,others) with sigma level in axis 1
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),sigma,self.ilat,self.ilon]
        elif len(self.ncin[var].shape) == 3:
            # is a 2D var with time in axis 0
            self.var = self.ncin[var][self.timeStart.item():self.timeEnd.item(),self.ilat,self.ilon]

        self.time = self.ncin.time[self.timeStart.item():self.timeEnd.item()]

    def graph(self,title='',**kwargs):
        fig,ax = plt.subplots()
        ax.plot(self.time,self.var,'k')
        ax.set_ylabel(self.var.attrs['long_name'] + " [" +self.var.attrs['units'] + "]",fontsize=8)
        ax.set_xlabel(self.time.attrs['long_name'],fontsize=8)
        ax.set_title(title,fontsize=12)
        ax.margins(0)
