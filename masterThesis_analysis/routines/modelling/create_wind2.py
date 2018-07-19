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

class HeatFlux_data_input(object):

    def hflx_nc(self,nc,x='x',y='y',t='t',hflx='hflx',offset=0):
        """ read wind data from netcdf file """

        self.f_wind = xr.open_dataset(nc)
        var         = lambda x : self.f_wind[x].data
        self.x     = var(x)-offset
        self.y     = var(y)
        self.hflx     = var(hflx)
        self.t     = var(t)
        self.xm,self.ym = np.meshgrid(self.x,self.y)

        # changing sign, because Total Downward Heat Flux must to be Total Upward: from the ocean to the atmosphere.
        # So, following the suggestion in www.ccpo.odu.edu/POMWEB/FAQ.txt, we need to multiply these data by (-1).
        # for more information, check in these website, searching for the Question. "I have heat flux data in W/m2, how to I convert it to WTSURF":
        """
            Q:
               I have heat flux data in W/m2, how do I convert it to WTSURF?

            A:
               The model surface heat flux forcing WTSURF (as well as the short wave
            radiation SWRAD) are in units of Km/s; so you have to divide the fields by
            the factor 4.1876E6. Note also that WTSURF is the flux from the ocean upward
            into the atmosphere, so that WTSURF>0 -> cooling; some atmospheric data sets
            (e.g., Oberhuber's COADS) give downward fluxes from the atmosphere to the
            ocean, so you will need to change the sign too.

        """
        self.hflx *= -1

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

    def to_secom_hflx(self,t,intp_flag,flag=-999):

        hf = []
        size = (self.mg.nj-2,self.mg.ni-2)
        coarser2finer_reshape = lambda x : np.reshape(self.coarser2finer_ww3_to_secom(x,intp_flag),(size))
        for i in xrange(t):
            hf.append(coarser2finer_reshape(self.dados.hflx[i]))

        hf = np.array(hf)
        hf[np.isnan(hf)] = -999.

        fill_gaps_reshape = lambda x,y : np.reshape(interp.fill_gaps_secom_grid(x,y),(size))
        for i in xrange(t):
            intp_flag = interp.interpolate_flag(hf[i]).ravel()
            hf[i]      = fill_gaps_reshape(hf[i],intp_flag)

        HF = interp.completa(hf)

        return HF

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

        vento.close()

    def hflx_to_ascii(file_name,HF,ano,mes,dt=3,t0=0,formato="w+"):
        """
        Funcao que salva a matriz HF no formato ascii
        com intervalo de tempo dt (tempo em horas entre dois passos).
        Os valores de ano e mes sao para nomear o arquivo.
        INPUT:
        HF - heat flux [W/m^2]

        """
        fluxo = open(file_name, formato)

        for t in range(HF.shape[0]):
            fluxo.write("%10f" % (t*dt+t0))
            fluxo.write("\n")
            for i in range(U.shape[1]):
                for j in range(U.shape[2]):
                    fluxo.write("%5d%5d%10.3f" % (j+1,i+1,HF[t,i,j]))
                    fluxo.write("\n")

        # Danilo: repetir os ultimos dados de vento para um tempo 9999.00
 #       vento.write("%10f"%(9999))
 #       vento.write("\n")
 #       for i in range(U.shape[1]):
 #           for j in range(U.shape[2]):
 #               vento.write("%5d%5d%10.3f%10.3f%10.3f" % (j+1,i+1,U[t,i,j],V[t,i,j],P[t,i,j]))
 #               vento.write("\n")


        fluxo.close()


    import matplotlib.pyplot as plt
    import glob

    #essa primeira parte e de leitura de dados
    arquivo  = '/home/danilo/Dropbox/mestrado/data/data2model/JF2014/tuv/' #onde esta o dado de vento .nc
    nc    = '*.nc'
    f_mg  = '/home/danilo/Dropbox/mestrado/grade/model_grid_com_pontos_em_terra' #model grid com pontos em Terra.
                                             # só funcina com esse model grid, mas pra rodar o modelo é sem os pontos em terra depois

    ncfiles_wind = glob.glob(arquivo+nc)
    ncfiles_wind.sort()

    arquivo  = '/home/danilo/Dropbox/mestrado/data/data2model/JF2014/hflx/' #onde esta o dado de vento .nc
    ncfiles_hflx = glob.glob(arquivo+nc)
    ncfiles_hflx.sort()

    timestep        = 6 #horas (dado original)
    t0              = 0 # não começa em zero pq estou rodando hot start
    ano             = 2014
    mes             = 01
    wind_multiplier = 1.6 # atualizar ali embaixo na formula!
    file_name_wind       = 'vento'#+str(ano)+str(mes)
    file_name_fluxo      = 'calor'

    for file_i in np.arange(0,len(ncfiles_hflx),1):
        nc_hflx = ncfiles_hflx[file_i]
        nc_wind = ncfiles_wind[file_i]

        if file_i == 0:
            write_format    = "w+"
        else:
            write_format    = "a"

        # load model_grid with land cells
        smgrd = Secom_model_grid2(f_mg)

        # read wind file
        wndnc = Wind_data_input()
        wndnc.wind_nc(nc_wind,t='time',x='lon',y='lat',u='U_GRD_L103',v='V_GRD_L103',offset=360)

        # read heat flux file
        hflxnc = HeatFlux_data_input()
        hflxnc.hflx_nc(nc_hflx,t='time',x='lon',y='lat',hflx='THFLX_L1_Avg_1',offset=360)

        # interpolate reanalysis data into model_grid
        interp = Interpolation()
        interp.instancia_dados(wndnc)
        interp.instancia_modelgrid(smgrd)
        intp_flag = interp.interpolate_flag(wndnc.u[0]).ravel()
        t = wndnc.u.shape[0]
        U,V,P = interp.to_secom(t,intp_flag,xwind_multiplier=1.6,ywind_multiplier=1.6)

        # write interpolated wind data
        wind_to_ascii(file_name_wind,U,V,P,ano,mes,dt=6,t0=t0,formato=write_format)

        # interpolate heat flux data
        interp = Interpolation()
        interp.instancia_dados(hflxnc)
        interp.instancia_modelgrid(smgrd)
        intp_flag = interp.interpolate_flag(hflxnc.hflx[0]).ravel()
        t = hflxnc.hflx.shape[0]
        HF = interp.to_secom_hflx(t,intp_flag)

        # write interpolatead heat flux data
        hflx_to_ascii(file_name_fluxo,HF,ano,mes,dt=6,t0=t0,formato=write_format)

        # informations for the next routine
        t0 = timestep*t+t0
        number_time_steps = t0/6 # numero de timesteps usado no write_wind.py
