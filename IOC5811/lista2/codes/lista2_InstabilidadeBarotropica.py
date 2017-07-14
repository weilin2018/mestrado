#-*-coding:utf-8-*-
# lista 2 de DFG II
import numpy as np
import matplotlib.pyplot as plt
from math import tanh, cosh, sin, cos
import seawater as sw
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

import matplotlib
matplotlib.style.use('ggplot')

# exercicio 2.A - jato de bickley dado por u(y) =
def calcular_dqdy(y, Lo, U, beta):
	"""
		Calcula o gradiente de vorticidade potencial básica de um Jato
		de Bickley no plano f, onde:

		o gradiente de vorticidade potencial básica por:
			.. math::
				\frac{d\hat{q}}{dy} = \beta - \frac{d^2\hat{u}}{dy^2}

		a estrutura de velocidade do Jato de Bickley é dada por:
			.. math::
				\bar{u}(y) = \hat{U} sech^{2}(\left(\frac{y}{L_o} \right))

		Após derivação, considerando:
			.. math::
				\beta = 0

		Obtemos:
			.. math::
				\frac{d\hat{q}}{dy} = - \frac{2 * \hat{U}}{L_o^{2} cosh^{4}(\left(\frac{y}{L_o} \right))}
	"""

	from math import sinh, cosh

	# criar vetor para armazenar os valores de dqdy no intervalo de y
	dqdy = np.zeros(len(y)) * np.nan

	# calcular o termo 1: 2*U / (Lo^2 cosh^4(y/Lo))
	Lo_2  = Lo**2
	cosh4 = np.asarray([cosh(i/Lo)**4 for i in y])
	termo1 = 2*U / (Lo_2 * cosh4)

	# calcular o termo 2 : [ 2sinh^2(y/Lo) - 1 ]
	sinh2 = np.asarray([sinh(i/Lo)**2 for i in y])
	termo2 = 2 * sinh2 - 1

	dqdy = beta - termo1*termo2

	return dqdy

def plotar_ex2a(U, Lo, L, beta=0., figname=''):
	""" plotar velocidade zonal e gradiente da vorticidade potencial basica """

	print("Calculando o potencial de vorticidade potencial basico ...")
	y 	 = np.arange(-L, L, 100) # escala meridional // eixo y
	dqdy = calcular_dqdy(y, Lo, U, beta)  # eixo x

	print("Calculando o perfil de velocidade zonal do jato de bickley ...")
	u = np.asarray([ U*1/cosh(i/Lo)**2 for i in y ])

	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,14))

	ax1.plot(u,y, 'k')

	#ax1.set_ylim(200000, -200000)
	#ax1.set_title('Ex.2-a) Velocidade Zonal (u)')
	#ax1.yaxis.set_major_formatter(FormatStrFormatter)
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y, 'k')

	p_inflexao = [(0., 0.66*Lo), (0., -0.66*Lo)]
	for p in p_inflexao:
		ax2.scatter(p[0], p[1], color='k', marker='o')

	# ax2.text(2.66791e-15, 37228.7, r"$y = 0.66L_o$")
	# ax2.text(2.66791e-15, -60228.7, r"$y = 0.66L_o$")
	#ax2.set_title('Ex.2-a) Gradiente de Vorticidade Potencial Basica')
	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(r"$\frac{d\bar{q}}{dy}$ [$m s^{-1}$]")
	ax2.set_xlim([dqdy.min(), dqdy.max()])
	ax2.set_ylim([y.min(), y.max()])

	plt.suptitle(r" Ex 2a: Jato de Bickley [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]" + "\n" + r" no plano f e $\hat{U} = 0.05 m s^{-1}$", fontsize=20)

	ax1.text(0.0478912, -354741, '(A)')
	ax2.text(3.76764e-11, -353546, '(B)')

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	return u, dqdy, y

def calcular_fjortoft(u_barra, u_zero, dqdy):
	# calcular termo1: (ubarra - uzero)
	termo1 = u_barra - u_zero

	# calcular termo2: dqdy * termo1
	termo2 = dqdy * termo1

	return termo2

def plotar_ex2b(y, dqdy, u_barra, u_0, Lo, figname=""):

	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,14))

	# plotando o critério de Rayleigh-Kuo
	p_inflexao = [(0., 0.66*Lo), (0., -0.66*Lo)]

	ax1.plot(dqdy, y, 'k')
	for p in p_inflexao:
		ax1.scatter(p[0], p[1], color='k', marker='o')

	# ax1.text(2.66791e-15, 37228.7, r"$y = 0.66L_o$")
	# ax1.text(2.66791e-15, -60228.7, r"$y = 0.66L_o$")

	# ax1.set_title(u'Ex.2-b) Critério de Rayleigh-Kuo')
	#ax1.yaxis.set_major_formatter(FormatStrFormatter)
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Critério de Rayleigh-Kuo -" + r" $\frac{d\bar{q}}{dy}$ [$m s^{-1}$]")
	ax1.set_xlim([dqdy.min(), dqdy.max()])
	ax1.set_ylim(y.min(), y.max())

	# plotando o critério de Fjortoft
	fjor = calcular_fjortoft(u_barra, u_0, dqdy)
	ax2.plot(fjor, y, 'k')

	# for p in p_inflexao:
	# 	ax2.scatter(p[0], p[1], color='k', marker='o')
	#
	# ax2.set_title(u'Ex.2-b) Critério de Fjortoft')
	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(u"Critério de Fjortoft - " + r"$(\bar{u} - u_o) \frac{dq}{dy}$ [$m s^{-1}$]")
	ax2.set_xlim([fjor.min(), fjor.max()])
	ax2.set_ylim(y.min(), y.max())

	ax1.text(3.73795e-11, -368410, '(A)')
	ax2.text(7.68869e-13, -368410, '(B)')

	plt.suptitle(u" Ex 2b: Avaliação dos critérios de Rayleigh-Kuo e Fjortoft \n para o Jato de Bickley "
					r"do item 2a: [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]", fontsize=20)

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	#plt.show()

def plotar_ex2c(U, Lo, L, beta, dqdy_min, figname=""):
	""" plotar: velocidade zonal do jato no plano beta e dq/dy no plano beta """

	# calcular o gradiente de VP basico no plano beta
	y 	 = np.arange(-L, L, 100) # escala meridional // eixo y
	dqdy = calcular_dqdy(y, Lo, U, beta)  # eixo x
	# calcular a velocidade zonal no plano beta
	u = np.asarray([ U*1/cosh(i/Lo)**2 for i in y ])

	# plotar no grafico superior o perfil de velocidade zonal e no inferior
	# o gradiente de VP basico
	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,14))

	# plotando o critério de Rayleigh-Kuo
	ax1.plot(u, y, 'k')
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y, 'k')

	#ax2.set_title('Ex.2-a) Gradiente de Vorticidade Potencial Basica')
	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(u"Critério de Rayleigh-Kuo -" + r" $\frac{d\bar{q}}{dy}$ [$m s^{-1}$]")
	ax2.set_xlim([dqdy_min, dqdy.max()])
	ax2.set_ylim([y.min(), y.max()])

	ax1.text(0.0478912, -354741, '(A)')
	ax2.text(5.64334e-11, -353546, '(B)')

	plt.suptitle(r" Ex 2c: Jato de Bickley [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]" + "\n" +
				 r" no plano $\beta$ e $\hat{U} = 0.05 m s^{-1}$", fontsize=20)

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	return u, dqdy, y

def plotar_ex2d(U, Lo, L, beta, figname=''):
	""" plotar
		U = 0.085 m/s
		plano Beta

		perfil de velocidade zonal
		gradiente de VP básico
	"""
	# calcular dqdy
	y 	 = np.arange(-L, L, 100) # escala meridional // eixo y
	dqdy = calcular_dqdy(y, Lo, U, beta)  # eixo x
	# calcular velocidade zonal
	u = np.asarray([ U*1/cosh(i/Lo)**2 for i in y ])

	fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,14))

	ax1.plot(u,y,'k')

	p_inflexao = [(0., 0.95*Lo), (0., 1.39*Lo), (0., -0.95*Lo), (0., -1.39*Lo)]
	for p in p_inflexao:
		ax2.scatter(p[0], p[1], color='k', marker='o')

	# ax2.text(2.66791e-15, 37228.7, r"$y = 0.66L_o$")
	# ax2.text(2.66791e-15, -60228.7, r"$y = 0.66L_o$")

	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y, 'k')

	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(u"Critério de Rayleigh-Kuo -" + r" $\frac{d\bar{q}}{dy}$ [$m s^{-1}$]")
	ax2.set_xlim([dqdy.min(), dqdy.max()])
	ax2.set_ylim([y.min(), y.max()])

	plt.suptitle(r" Ex 2d: Jato de Bickley [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]" + "\n" +
			     r" no plano $\beta$ e $\hat{U} = 0.085 m s^{-1}$", fontsize=20)

	ax1.text(0.0817097, -352464, '(A)')
	ax2.text(8.3876e-11, -351269, '(B)')

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	return u, dqdy, y

# localizar o indice da array do ponto mais próximo do value passado
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def calcular_u0(dqdy, U, Lo, y):
	""" baseado no vetor de gradiente de VP basico,
	procurar onde dq/dy = 0 e utilizar este y para calcular u0

	"""
	# localizar o indice do ponto mais proximo a zero
	idy = find_nearest(dqdy, 0.0)
	# extrair somente os y de interesse
	ys = np.asarray(y[idy])

	if ys.size == 1:
		u = U * (1/cosh(ys/Lo))**2
	else:
		u = np.asarray([ U*1/cosh(i/Lo)**2 for i in ys ])

	return u




# CONSTANTES FORNECIDAS NO EXERCÍCIO
U  = 0.05 # m/s # velocidade maxima do jato
Lo = 5e4 # m # escala de decaimento do jato
L  = 4e5  # m -> Largura do canal idealizado no exercício

latitude = 30

# beta = 2*omega*sin(latitude central) : estamos considerando a aproximação do
# plano f
f   = 2 * (2*np.pi / 86400) * sin(latitude)

#
u, dqdy, y = plotar_ex2a(U, Lo, L, beta=0., figname='../outputs/Fig2_1.png')

dqdy_min = dqdy.min()

# calcular e plotar ex 2b
plotar_ex2b(y, dqdy, u, 0.03, Lo, figname="../outputs/Fig2_2.png")

# calcular e plotar ex 3c
lat_radians = np.deg2rad(latitude) # converter a latitude para radianos
beta = 2 * sw.constants.OMEGA*cos(lat_radians) / sw.constants.earth_radius
u, dqdy, y = plotar_ex2c(U, Lo, L, beta, dqdy_min, figname="../outputs/Fig2_3.png")

# calcular e plotar ex 3d - ../outputs/Fig2_4.png
U = 0.085
u, dqdy, y = plotar_ex2d(U, Lo, L, beta=beta, figname="")
