#-*-coding:utf-8-*-
# lista 2 de DFG II
import numpy as np
import matplotlib.pyplot as plt
from math import tanh, cosh, sin
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

import matplotlib
matplotlib.style.use('ggplot')

# exercicio 2.A - jato de bickley dado por u(y) =

def calcular_dqdy(y, Lo, U, beta):
	from math import tanh, cosh, sin

	argume = y/Lo # argumento das funcoes trigonometricas

	termo1 = (4/Lo**2) * ( tanh(argume)**2 * 1/cosh(argume)**2 )
	termo2 = (2/Lo**2) * ( 1/cosh(argume)**4 )
	
	return beta - U * ( (termo1) - (termo2) )

def plotar_ex2a(U, Lo, L, beta, figname=''):
	""" plotar velocidade zonal e gradiente da vorticidade potencial basica """

	print("Calculando o potencial de vorticidade potencial basico ...")
	y 	 = np.arange(-L, L, 100) # escala meridional // eixo y
	dqdy = [] # eixo x

	for i in y:
		dqdy.append(calcular_dqdy(i, Lo, U, beta))

	print("Calculando o perfil de velocidade zonal do jato de bickley ...")
	u = []
	for i in y:
		u.append( U * 1/cosh(i/Lo)**2 )


	u = np.asarray(u)
	dqdy = np.asarray(dqdy)

	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,14))

	ax1.plot(u,y)
	#ax1.set_ylim(200000, -200000)
	ax1.set_title('Ex.2-a) Velocidade Zonal (u)')
	#ax1.yaxis.set_major_formatter(FormatStrFormatter)
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y)
	ax2.set_title('Ex.2-a) Gradiente de VP Basica')
	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(u"Gradiente da Vorticidade Potencial Básica [$m s^{-1}$]")
	ax2.set_xlim([dqdy.min(), dqdy.max()])

	plt.suptitle(r'Jato de Bickley [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]', fontsize=20)

	if figname!='': # então rola salvar
		plt.savefig(figname)

	plt.show()

	return u, dqdy, y

# CONSTANTES FORNECIDAS NO EXERCÍCIO
U  = 0.05 # m/s # velocidade maxima do jato
Lo = 5e4 # m # escala de decaimento do jato
L  = 4e5  # m -> Largura do canal idealizado no exercício

latitude = 30

# beta = 2*omega*sin(latitude central) : estamos considerando a aproximação do plano f 
beta   = 2 * (2*np.pi / 86400) * sin(latitude)

# calcular e plotar ex 2a
u, dqdy, y = plotar_ex2a(U, Lo, L, beta, figname='../outputs/ex2_a.png')

# calcular e plotar ex 2b
