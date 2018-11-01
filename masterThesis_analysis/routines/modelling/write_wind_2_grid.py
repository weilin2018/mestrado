#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Modificado por Danilo:
    adicionado bloco de comandos para repetir o ultimo instante de dados para um tempo 9999.00

"""

#escreve os vento para a grade

#IMPORTANTE LEMBRAR DE ABRIR O PYTHON E COLOCAR NA PASTA CERTA,
#PARA QUE OS ARQUIVOS SEJAM SALVOS NO LOCAL CORRETO
#obs.: pode ser antes de abrir o python


import numpy as np
import xarray as xr
from scipy import interpolate



class Secom_model_grid2(object):
    #primeira linha reservada para o autor do model grid
    #segunda linha descreve os sigma levels
    #terceira linha apresenta os n sigma levels introduzidos
    #3+n+1 linha apresenta o 'Horizontal Segmentations'
    #3+n+2 linha apresenta as dimensões os tamanho de I e J (conferir se respectivamente)
    def __init__(self,filename):
        with open(filename) as f:
            model_grid = f.readlines()
        self.model_grid2_list(model_grid)
        self.sigma_level()
        self.IJ_size()
        self.model_grid_data()

        lonlat      = lambda x : np.array(self.mg_data)[:,x]
        self.i   = lonlat(0)
        self.j   = lonlat(1)
        self.lon = lonlat(7)
        self.lat = lonlat(6)
        self.dep = lonlat(4)

    def reshape(self):
        reshape     = lambda x: x.reshape(self.ni-2,self.nj-2)
        self.lon = reshape(self.lon)
        self.lat = reshape(self.lat)
        self.i   = reshape(self.i)
        self.j   = reshape(self.j)
        self.dep = reshape(self.dep)


    def model_grid2_list(self,model_grid):
        l =[]
        for t in model_grid:
            t1 = t.split()
            try:
                l.append([float(i) for i in t1 ])
            except ValueError:
                l.append([t])
            self.model_grid_list = l

    def sigma_level(self,n=2):
        self.nsigma_level = int(self.model_grid_list[n][0])

    def IJ_size(self,n=4):
        self.ni  = int(self.model_grid_list[self.nsigma_level+n][0])
        self.nj  = int(self.model_grid_list[self.nsigma_level+n][1])

    def model_grid_data(self,n=5):
        self.mg_data = self.model_grid_list[self.nsigma_level+n:]



class Wind_data_input(object):
    # def __init__(self):
    #     self.wnd = Object()

    def wind_nc(self,nc,x='x',y='y',t='t',u='u',v='v',offset=0):
        """
        Reads wavewatch reanalysis results
        nc: netcdf file
        hs: wave height ...
        hs: wave period ...
        hs: wave direction ...
        """
        self.f_wnd = xr.open_dataset(nc)
        var            = lambda x : self.f_wnd[x].data
        self.x     = var(x)-offset
        self.y     = var(y)
        self.u     = var(u)
        self.v     = var(v)
        self.t     = var(t)
        self.xm,self.ym = np.meshgrid(self.x,self.y)

class Interpolation(object):
    """
    """
    def __init__(self):
        self.griddata = interpolate.griddata

    def instancia_dados(self,dados):
        self.dados   = dados

    def instancia_modelgrid(self,model_grid):
        self.mg     = model_grid


    def interpolate_flag(self,var,flag=-999):
        """
        Retorna array com True e False de acordo com o flag.
        Nao funciona com nan.
        """
        var[var == flag] =-999
        var_list         = np.array([i.ravel() for i in var])
        var_not_flagged  = var_list[:]>-999
        return var_not_flagged

    def coarser2finer_ww3_to_secom(self,values,intp_flag,method = 'linear'):
        """
        Metodo para interpolar os dados do WaveWatchIII para a grade do sECOM.
        Ambos devem ser intanciados antes, o WW3 com instancia_dados e sECOM
        com instancia_modelgrid.
        Por algum motivo pode retornar grade com nan's.
        """
        x = self.dados.xm.ravel()
        y = self.dados.ym.ravel()
        xi= self.mg.lon
        yi= self.mg.lat

        try:
            interpolated = self.griddata((x[intp_flag],y[intp_flag]),
                              values.ravel()[intp_flag],
                             (xi,yi),
                              method=method)
        except:
            interpolated = self.griddata((x,y),values.ravel(),(xi,yi), method=method)

        return interpolated

    def fill_gaps_secom_grid(self,values,intp_flag):
        """
        Metodo criado para substituir valores nan na grade do sECOM.
        Ele procura substitui os flags por valores do ponto mais proximo da sua propria grade.
        """
        x = self.mg.lon
        y = self.mg.lat

        fill_gaps = self.griddata((x[intp_flag],y[intp_flag]),
                            values.ravel()[intp_flag],
                            (x,y),
                            method='nearest')
        return fill_gaps


    def to_secom(self,t, intp_flag, xwind_multiplier=1, ywind_multiplier=1,flag=-999):
        """
        t: numero de time steps do modelo utilizado
        wind_multiplier : multiplica o vento por um valor arbitrário
        flag: valore inválidos no modelo
        """
        u=[]
        v=[]
        size = (self.mg.nj-2, self.mg.ni-2)
        coarser2finer_reshape = lambda x : np.reshape(interp.coarser2finer_ww3_to_secom(x,intp_flag),(size))
        for i in xrange(t):
            u.append(coarser2finer_reshape(self.dados.u[i]))
            v.append(coarser2finer_reshape(self.dados.v[i]))
        p = np.full_like(u,1024) # pressure: pressao atmosferica igual a 1024 milibares, para alta pressao, comum na regiao.

        u=np.array(u)*xwind_multiplier
        v=np.array(v)*ywind_multiplier

        #os passos a seguir atribuem o valor -999 como flag e os usbstituem pelo valor mais proximo
        u[np.isnan(u)]= -999
        v[np.isnan(v)]= -999

        fill_gaps_reshape = lambda x,y : np.reshape(interp.fill_gaps_secom_grid(x,y),(size))
        for i in xrange(t):
            intp_flag = interp.interpolate_flag(u[i]).ravel()
            u[i]      = fill_gaps_reshape(u[i],intp_flag)
            intp_flag = interp.interpolate_flag(v[i]).ravel()
            v[i]      = fill_gaps_reshape(v[i],intp_flag)

        U = self.completa(u)
        V = self.completa(v)
        P = self.completa(p)

        return U,V,P

    def completa(self,matriz):
        """
        Funcao copia as primeiras e as ultimas linhas e colunas
         de cada intervalo de tempo das variaveis
        """
        matriz = np.insert(matriz,0,matriz[:,0,:],1)
        matriz = np.insert(matriz,matriz.shape[1],matriz[:,-1,:],1)
        matriz = np.insert(matriz,0,matriz[:,:,0],2)
        matriz = np.insert(matriz,matriz.shape[2],matriz[:,:,-1],2)
        return matriz

if __name__=='__main__':


    def wind_to_ascii(file_name,U,V,P,ano,mes,dt=3,t0=0,formato="w+"):
        """
        Funcao que salva as matrizes U,V e P no formato ascii
        com intervalo de tempo dt (tempo em horas entre dois passos).
        Os valores de ano e mes sao para nomear o arquivo.
        INPUT:
        U - velocidade zonal do vento (m/s)
        V - velocidade meridional do vento (m/s)
        P - Pressao atmosferica (mBar)

        """
        vento = open(file_name, formato)

        for t in range(U.shape[0]):
            vento.write("%10f" % (t*dt+t0))
            vento.write("\n")
            for i in range(U.shape[1]):
                for j in range(U.shape[2]):
                    vento.write("%5d%5d%10.3f%10.3f%10.3f" % (j+1,i+1,U[t,i,j],V[t,i,j],P[t,i,j]))
                    vento.write("\n")

        # Danilo: repetir os ultimos dados de vento para um tempo 9999.00
 #       vento.write("%10f"%(9999))
 #       vento.write("\n")
 #       for i in range(U.shape[1]):
 #           for j in range(U.shape[2]):
 #               vento.write("%5d%5d%10.3f%10.3f%10.3f" % (j+1,i+1,U[t,i,j],V[t,i,j],P[t,i,j]))
 #               vento.write("\n")


        vento.close()


    import matplotlib.pyplot as plt
    import glob

    #essa primeira parte e de leitura de dados
    arquivo  = '/media/danilo/Danilo/mestrado/dados_potranca/' #onde esta o dado de vento .nc
    nc    = '*.nc'
    f_mg  = '/home/danilo/Dropbox/mestrado/grade/model_grid_com_pontos_em_terra' #model grid com pontos em Terra.
                                             # só funcina com esse model grid, mas pra rodar o modelo é sem os pontos em terra depois

    ncfiles = glob.glob(arquivo+nc)
    ncfiles.sort()

    timestep        = 6 #horas (dado original)
    t0              = 0
    ano             = 2014
    mes             = 01
    wind_multiplier = 1.6 # atualizar ali embaixo na formula!
    file_name       = 'vento'#+str(ano)+str(mes)

    for file_i,nc in enumerate(ncfiles):

        if file_i==0:
            write_format    = "w+"
        else:
            write_format    = "a"


        smgrd = Secom_model_grid2(f_mg)
        wndnc = Wind_data_input()
        wndnc.wind_nc(nc,t='time',x='lon',y='lat',u='U_GRD_L103',v='V_GRD_L103',offset=360)

        interp = Interpolation()
        interp.instancia_dados(wndnc)
        interp.instancia_modelgrid(smgrd)
        intp_flag = interp.interpolate_flag(wndnc.u[0]).ravel()
        t = wndnc.u.shape[0]
        U,V,P = interp.to_secom(t,intp_flag,xwind_multiplier=1.6,ywind_multiplier=1.6)
        wind_to_ascii(file_name,U,V,P,ano,mes,dt=6,t0=t0,formato=write_format)
        t0 = timestep*t+t0
        number_time_steps = t0/6 # numero de timesteps usado no write_wind.py
