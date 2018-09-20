#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
coded by danilo, to use the meteorological option "2DMETRND" in run_data,
by forcing the simulation with:
	wind
	short wave radiation
	air temperature [2m]
	relative humidity
	barotropic pressure in mbar
	cloud cover fraction
	extinction coefficient
	amount of precipitation in m/year
	amount of evaporation in m/year
		but both are converted to mm/sec in the code

"""

#escreve os vento para a grade

#IMPORTANTE LEMBRAR DE ABRIR O PYTHON E COLOCAR NA PASTA CERTA,
#PARA QUE OS ARQUIVOS SEJAM SALVOS NO LOCAL CORRETO
#obs.: pode ser antes de abrir o python


import numpy as np
import xarray as xr
from scipy import interpolate
import glob

def convert_uv2intdir(wu,wv):
	from math import atan2

	ws,wd = np.zeros(wu.shape),np.zeros(wu.shape)

	# loop in time
	ws = np.sqrt(wu**2 + wv**2)

	for t in range(wu.shape[0]):
		u,v = wu[t,:,:],wv[t,:,:]

		for i in range(u.shape[0]):
			for j in range(u.shape[1]):
				direction = 270 - (180/np.pi)*atan2(v[i,j],u[i,j])
				if direction > 360:
					direction -= 360

				wd[t,i,j] = direction

	return ws,wd

# class to instanciate model_grid (FULL, with land nodes too)
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

class MeteorData_input(object):

	def wind_nc(self,nc,x='x',y='y',t='t',u='u',v='v',kind=None,offset=0):
		self.f_wnd = xr.open_dataset(nc)
		var        = lambda x : self.f_wnd[x].data
		self.x     = var(x) - offset
		self.y     = var(y)
		self.u     = var(u)
		self.v     = var(v)
		self.t     = var(t)
		self.xm,self.ym = np.meshgrid(self.x,self.y)

	def data_nc(self,nc,x='x',y='y',t='t',data='u',v='v',kind=None,offset=0):
		self.ncin = xr.open_dataset(nc)
		var            = lambda x : self.ncin[x].data
		self.x     = var(x)-offset
		self.y     = var(y)
		self.data  = var(data)
		if kind == 'wind':
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


    def wind_to_secom(self,t, intp_flag, xwind_multiplier=1, ywind_multiplier=1,flag=-999):
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

        fill_gaps_reshape = lambda x,y : np.reshape(self.fill_gaps_secom_grid(x,y),(size))
        for i in xrange(t):
            intp_flag = interp.interpolate_flag(u[i]).ravel()
            u[i]      = fill_gaps_reshape(u[i],intp_flag)
            intp_flag = interp.interpolate_flag(v[i]).ravel()
            v[i]      = fill_gaps_reshape(v[i],intp_flag)

        U = self.completa(u)
        V = self.completa(v)
        P = self.completa(p)

        return U,V,P

    def meteo_to_secom(self,t,intp_flag,flag=-999):

        pass

    def to_secom_data(self,t,intp_flag,flag=-999):

        raw_data = []
        size = (self.mg.nj-2,self.mg.ni-2)
        coarser2finer_reshape = lambda x : np.reshape(self.coarser2finer_ww3_to_secom(x,intp_flag),(size))
        for i in xrange(t):
            raw_data.append(coarser2finer_reshape(self.dados.data[i]))

        raw_data = np.array(raw_data)
        raw_data[np.isnan(raw_data)] = -999.

        fill_gaps_reshape = lambda x,y : np.reshape(interp.fill_gaps_secom_grid(x,y),(size))
        for i in xrange(t):
            intp_flag = interp.interpolate_flag(raw_data[i]).ravel()
            raw_data[i]      = fill_gaps_reshape(raw_data[i],intp_flag)

        interpolated_data = interp.completa(raw_data)

        return interpolated_data

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

def data_to_ascii(file_name,U,V,SW,AT,RH,BP,CC,QP,QE,ano,mes,dt=6,t0=0,formato='w+'):
    metdata = open(file_name, formato)

    for t in range(U.shape[0]):
        metdata.write("%10f" % (t*dt+t0))
        metdata.write("\n")
        for i in range(U.shape[1]):
            for j in range(U.shape[2]):
                # metdata.write("%5d%5d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f" % (j+1,i+1,U[t,i,j],V[t,i,j],SW[t,i,j],AT[t,i,j],RH[t,i,j],BP[t,i,j],CC[t,i,j],0.000,QP[t,i,j],QE[t,i,j]))
				metdata.write("%5d%5d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f" % (j+1,i+1,U[t,i,j],V[t,i,j],SW[t,i,j],AT[t,i,j]-273.15,RH[t,i,j],BP[t,i,j],CC[t,i,j],0.042))
				# using Extinction Coefficient = 0.042, like Optical Type Water described by Jerlov (1976)
				metdata.write("\n")

    metdata.close()

import matplotlib.pyplot as plt
import glob
import os

os.system('clear')

BASE_DIR = '/home/danilo/Dropbox/mestrado/data/data2model/JF2014/'

f_mg  = '/home/danilo/Dropbox/mestrado/grade/model_grid_com_pontos_em_terra' #model grid com pontos em Terra.

# define path and read files for each parameter
nfiles_wind = glob.glob(BASE_DIR + 'tuv/*.nc')  # wind files
nfiles_wind.sort()
nfiles_swra = glob.glob(BASE_DIR + 'swobs/*.nc') # short wave radiation
nfiles_swra.sort()
nfiles_airt = glob.glob(BASE_DIR + 't2m/*.nc')   # air temperature
nfiles_airt.sort()
nfiles_relh = glob.glob(BASE_DIR + 'relhum/*.nc')# relative humidity
nfiles_relh.sort()
nfiles_pmsl = glob.glob(BASE_DIR + 'prmsl/*.nc') # pressure mean sea level
nfiles_pmsl.sort()
nfiles_cldf = glob.glob(BASE_DIR + 'cld/*.nc')   # cloud fraction
nfiles_cldf.sort()
nfiles_qpre = glob.glob(BASE_DIR + 'qprec/*.nc') # amount of precipitation
nfiles_qpre.sort()
nfiles_evap = glob.glob(BASE_DIR + 'qevap/*.nc') # amount of evaporation
nfiles_evap.sort()

# define global parameters
timestep        = 6 #horas (dado original)
t0              = 0 # não começa em zero pq estou rodando hot start
ano             = 2014
mes             = 01
wind_multiplier = 1.6 # atualizar ali embaixo na formula!
file_name_data  = 'metData_2DMETLNP'

for file_i in np.arange(0,63,1):

	# check for create ou just open metData
	if file_i == 0:
		write_format   = 'w+'
	else:
		write_format   = 'a'

	# instanciate model_grid
	smgrd = Secom_model_grid2(f_mg)

	#-------------------- wind
	data = MeteorData_input()
	data.wind_nc(nfiles_wind[file_i],t='time',x='lon',y='lat',u='U_GRD_L103',v='V_GRD_L103',kind='wind',offset=360)

	# interpolate reanalysis data into model_grid

	interp = Interpolation()
	interp.instancia_dados(data)
	interp.instancia_modelgrid(smgrd)
	intp_flag = interp.interpolate_flag(data.u[0]).ravel()
	t = data.u.shape[0]
	# ao interpolar para a grade refinada do modelo, ja convertemos de componentes para direcao e magnitude
	U,V,P = interp.wind_to_secom(t,intp_flag,xwind_multiplier=1.6,ywind_multiplier=1.6)

	#-------------------- SW radiation
	data = MeteorData_input()
	data.data_nc(nfiles_swra[file_i],t='time',x='lon',y='lat',data='DSWRF_L1',offset=360)
	interp = Interpolation()
	interp.instancia_dados(data)
	interp.instancia_modelgrid(smgrd)
	intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	t = data.data.shape[0]
	SW = interp.to_secom_data(t,intp_flag)

	#-------------------- Air Temperature (AT)
	data = MeteorData_input()
	data.data_nc(nfiles_airt[file_i],t='time',x='lon',y='lat',data='TMP_L103',offset=360)
	interp = Interpolation()
	interp.instancia_dados(data)
	interp.instancia_modelgrid(smgrd)
	intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	t = data.data.shape[0]
	AT = interp.to_secom_data(t,intp_flag)

	#-------------------- Relative Humidity (RH)
	data = MeteorData_input()
	data.data_nc(nfiles_relh[file_i],t='time',x='lon',y='lat',data='R_H_L103',offset=360)
	interp = Interpolation()
	interp.instancia_dados(data)
	interp.instancia_modelgrid(smgrd)
	intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	t = data.data.shape[0]
	RH = interp.to_secom_data(t,intp_flag)

	#-------------------- Batropic Pressure (BP)
	data = MeteorData_input()
	data.data_nc(nfiles_pmsl[file_i],t='time',x='lon',y='lat',data='PRMSL_L101',offset=360)
	interp = Interpolation()
	interp.instancia_dados(data)
	interp.instancia_modelgrid(smgrd)
	intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	t = data.data.shape[0]
	BP = interp.to_secom_data(t,intp_flag)
	BP *= 0.01

	#-------------------- Cloud Cover (CC)
	# data = MeteorData_input()
	# data.data_nc(nfiles_pmsl[file_i],t='time',x='lon',y='lat',data='T_CDC_L244',offset=360)
	# interp = Interpolation()
	# interp.instancia_dados(data)
	# interp.instancia_modelgrid(smgrd)
	# intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	# t = data.data.shape[0]
	CC = BP* 0. #interp.to_secom_data(t,intp_flag)

	#-------------------- Amount of Precipitation (QP)
	# data = MeteorData_input()
	# data.data_nc(nfiles_pmsl[file_i],t='time',x='lon',y='lat',data='A_PCP_L1_Accum_1',offset=360)
	# interp = Interpolation()
	# interp.instancia_dados(data)
	# interp.instancia_modelgrid(smgrd)
	# intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	# t = data.data.shape[0]
	QP =  BP* 0. #interp.to_secom_data(t,intp_flag) * 0.
	# precisar converter de kg/m² para m/year

	#-------------------- Amount of Evaporation (QE)
	# data = MeteorData_input()
	# data.data_nc(nfiles_pmsl[file_i],t='time',x='lon',y='lat',data='PEVPR_L1_Avg_1',offset=360)
	# interp = Interpolation()
	# interp.instancia_dados(data)
	# interp.instancia_modelgrid(smgrd)
	# intp_flag = interp.interpolate_flag(data.data[0]).ravel()
	# t = data.data.shape[0]
	QE = BP* 0. # interp.to_secom_data(t,intp_flag) * 0.
	# precisar converter de W/m² para m/year

	data_to_ascii(file_name_data,U,V,SW,AT,RH,P,CC,QP,QE,ano,mes,dt=6,t0=t0,formato=write_format)

	# informations for the next routine
	t0 = timestep*t+t0
	number_time_steps = t0/6 # numero de timesteps usado no write_wind.py
