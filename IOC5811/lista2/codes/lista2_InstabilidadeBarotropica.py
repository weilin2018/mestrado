#-*-coding:utf-8-*-
# lista 2 de DFG II
import numpy as np
import matplotlib.pyplot as plt
from math import tanh, cosh, sin, cos
import seawater as sw
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches

import matplotlib
matplotlib.style.use('ggplot')


# set to True if you're plotting the final graphs, otherwise keep as False
# and no image will be save, only showed to you.
PLOTandSAVE = True

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

	ax1.plot(u,y)

	#ax1.set_ylim(200000, -200000)
	#ax1.set_title('Ex.2-a) Velocidade Zonal (u)')
	#ax1.yaxis.set_major_formatter(FormatStrFormatter)
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y)

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

	ax1.plot(dqdy, y)
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
	ax2.plot(fjor, y)

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

	return fjor

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
	ax1.plot(u, y)
	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y)

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

	ax1.plot(u,y)

	ax1.set_ylabel(u"Distância meridional [m]")
	ax1.set_xlabel(u"Velocidade zonal do jato [$m s^{-1}$]")
	ax1.set_xlim([u.min(), u.max()])

	ax2.plot(dqdy, y)

	# calcular os pontos de inflexão e a velocidade nestes pontos
	u0, ys, idy = calcular_u0(dqdy, U, Lo, y)

	ax2.scatter(dqdy[idy[:2]], ys[:2], color='k', s=25, marker='o', label=r'$u_0 = %1.3f m s^{-1}$'%(u0[0]))
	ax2.scatter(dqdy[idy[2:]], ys[2:], color='k', s=25, marker='*', label=r'$u_0 = %1.3f m s^{-1}$'%(u0[-1]))

	ax2.set_ylabel(u"Distância meridional [m]")
	ax2.set_xlabel(u"Critério de Rayleigh-Kuo -" + r" $\frac{d\bar{q}}{dy}$ [$m s^{-1}$]")
	ax2.set_xlim([dqdy.min(), dqdy.max()])
	ax2.set_ylim([y.min(), y.max()])

	plt.suptitle(r" Ex 2d: Jato de Bickley [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]" + "\n" +
			     r" no plano $\beta$ e $\hat{U} = 0.085 m s^{-1}$", fontsize=20)

	ax1.text(0.0817097, -352464, '(A)')
	ax2.text(8.3876e-11, -351269, '(B)')

	# colocar box com os valores de u_0
	ax2.legend(loc='best', scatterpoints=1)

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	return u, dqdy, y, u0

def plotar_ex2d_parteII(y, dqdy, u_barra, u_0, Lo, figname=""):

	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,14))

	style = ['k', 'k--', 'k:', 'k-.']

	cont = 0

	for u0 in u_0:
		fjor = calcular_fjortoft(u_barra, u0, dqdy)
		ax.plot(fjor, y, style[cont],label=r'$u_0 = %1.3f m s^{-1}$'%(u0))

		cont += 1

	ax.set_ylabel(u"Distância meridional [m]")
	ax.set_xlabel(u"Critério de Fjortoft - " + r"$(\bar{u} - u_o) \frac{dq}{dy}$ [$m s^{-1}$]")
	ax.set_xlim([fjor.min(), fjor.max()])
	ax.set_ylim(y.min(), y.max())

	# vertical line in zero
	plt.axvline(0, color='black')

	# create legend box
	plt.legend(loc='best')

	plt.suptitle(u" Ex 2d: Avaliação do critério de Fjortoft \n para o Jato de Bickley "
					r"do item 2a: [$\bar{u}(y) = \hat{U} sech^2(\frac{y}{L_o})$]", fontsize=20)

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

	return fjor

def plotar_ex2d_parteIII(jatoA, jatoD, y, figname=""):

	jatoA = fjor_jatoA
	jatoD = fjor_jatoD

	# calculate area
	areaA = np.trapz(jatoA, y)
	areaD = np.trapz(jatoD, y)

	labelA = u'Área: {:0.3e}'.format(areaA)
	labelD = u'Área: {:0.3e}'.format(areaD)

	# plot

	fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,8), sharey=True)
	# plot jatoA
	ax1.plot(jatoA, y)
	ax1.fill_between(jatoA, -100, y, alpha=.5, label=labelA)
	ax1.set_title(u'Critério de Fjortoft - Jato A')
	ax1.set_xlabel(u"Critério de Fjortoft - " + r"$(\bar{u} - u_o) \frac{dq}{dy}$ [$m s^{-1}$]")
	ax1.set_ylabel(u'Distância Meridional [m]')
	ax1.axvline(0., color='gray')
	ax1.legend()

	ax2.plot(jatoD, y)
	ax2.fill_between(jatoD, -100, y, alpha=.5, label=labelD)
	ax2.set_title(u'Critério de Fjortoft - Jato D')
	ax2.set_xlabel(u"Critério de Fjortoft - " + r"$(\bar{u} - u_o) \frac{dq}{dy}$ [$m s^{-1}$]")
	ax2.axvline(0., color='gray')
	ax2.legend()

	ax1.text(8.37751e-13, -357798, '(A)')
	ax2.text(5.31048e-12, -353211, '(B)')

	if figname!='': # então rola salvar
		plt.savefig(figname, dpi=100)
	else:
		plt.show()

def find_nearest(n, s, r, index=False):
	"""
		Find 'n' values in 's' with the lowest absolute
		difference from the number 'r'.

		Based on: https://stackoverflow.com/questions/24112259/finding-k-closest-numbers-to-a-given-number

		Use the heapq.nsmallest package

		input:
			n = how many numbers
			s = list/array with data
			r = reference value
			index = control if the returned value is the
			data itself or only the indexes

		output:
			values of s or indexes, based on index argument
	"""

	from heapq import nsmallest

	if index:
		# return the index of data
		ids = []
		# list with all data found
		val = nsmallest(n,s, key=lambda x: abs(x-r))
		# search for index
		for v in val:
			ids.append(np.where(s == v)[0])

		return ids
	else:
		# return the data itself
		return nsmallest(n,s, key=lambda x: abs(x-r))

def calcular_u0(dqdy, U, Lo, y):
	""" baseado no vetor de gradiente de VP basico,
	procurar onde dq/dy = 0 e utilizar este y para calcular u0

	"""
	# localizar o indice dos 4 pontos mais proximos a zero
	idy = find_nearest(4, dqdy, 0., index=True)

	# tratar um pouco os indices pq vem bagunçado e repetido
	tmp = []
	for k in [0,2]:
		a = idy[k]

		tmp.append(a[0])
		tmp.append(a[1])

	idy = np.asarray(tmp)

	# extrair somente os y de interesse
	ys = y[idy]

	if ys.size == 1:
		u = U * (1/cosh(ys/Lo))**2
	else:
		u = np.asarray([ U*1/cosh(i/Lo)**2 for i in ys ])

	return u, ys, idy

# CONSTANTES FORNECIDAS NO EXERCÍCIO
U  = 0.05 # m/s # velocidade maxima do jato
Lo = 5e4 # m # escala de decaimento do jato
L  = 2e5  # m -> metade do domínio meridional do canal idealizado no exercício

latitude = 30

# beta = 2*omega*sin(latitude central) : estamos considerando a aproximação do
# plano f
f_o   = 2 * (2*np.pi / 86400) * sin(latitude)

#
if PLOTandSAVE:
	figname = '../outputs/Fig2_1.png'
else:
	figname = ''

u, dqdy, y = plotar_ex2a(U, Lo, L, beta=0., figname=figname)

dqdy_min = dqdy.min()

# calcular e plotar ex 2b
if PLOTandSAVE:
	figname = '../outputs/Fig2_2.png'
else:
	figname = ''

# save the fjortoft curve to calculate the area in ex2D.III
fjor_jatoA = plotar_ex2b(y, dqdy, u, 0.03, Lo, figname=figname)

# calcular e plotar ex 3c
lat_radians = np.deg2rad(latitude) # converter a latitude para radianos
beta = 2 * sw.constants.OMEGA*cos(lat_radians) / sw.constants.earth_radius

if PLOTandSAVE:
	figname = '../outputs/Fig2_3.png'
else:
	figname = ''


u, dqdy, y = plotar_ex2c(U, Lo, L, beta, dqdy_min, figname=figname)

# calcular e plotar ex 3d - ../outputs/Fig2_4.png
U = 0.085

if PLOTandSAVE:
	figname = '../outputs/Fig2_4.png'
else:
	figname = ''


u, dqdy, y, u0 = plotar_ex2d(U, Lo, L, beta=beta, figname=figname)
print("As velocidades nos pontos de inflexão são: \n")
for x in u0:
	print(x)

#plotar o critério de fjortoft
# como são 4 raízes, mas somente duas velocidades, então envio somente duas
# velocidades para calcular e plotar o critério de Fjortoft
if PLOTandSAVE:
	figname = '../outputs/Fig2_4_2.png'
else:
	figname = ''

# save the fjortoft curve to calculate the area in ex2D.III
fjor_jatoD = plotar_ex2d_parteII(y, dqdy, u, u0[1:3], Lo, figname=figname)

#### complemento ex 2d
if PLOTandSAVE:
	figname = '../outputs/Fig2_4_3.png'
else:
	figname = ''

plotar_ex2d_parteIII(fjor_jatoA, fjor_jatoD, y, figname=figname)



import glob
import os

lstFile = glob.glob('../outputs/Fig*')
for fname in lstFile:
	if fname != '../outputs/Fig2_4_2.png':
		os.system('convert -trim %s %s' % (fname, fname))


###################################################
#### calcular d²u/dy² numericamente (challenge!) ##
###################################################
dudy = np.zeros(y.shape, np.float)				# criar matriz

dudy[0:-1] = np.diff(u)/np.diff(y)				# calcular du/dy
dudy[-1]   = (u[-1] - u[-2])/(y[-1] - y[-2])	# calcular du/dy nos extremos

dudySeg = np.zeros(y.shape, np.float)			# criar matriz
dudySeg[0:-1] = np.diff(dudy)/np.diff(y)		# calcular d²u/dy²
dudySeg[-1]   = (dudy[-1] - dudy[-2])/(y[-1] - y[-2]) # calcular d²u/dy² nos extremos

# plotar solução numérica e conferir que baterá com o
# critério de Rayleigh-Kuo plotado na Fig2.4
# plt.plot(beta - dudySeg, y)
