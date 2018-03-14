"""
    Codigo puxa os arquivo .mat salvos pela rotina gill74.m e realiza os plots
    considerando os 5 casos de jato para instabilidade baroclínica
"""
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.style.use('ggplot')
import glob
import os

### load and plot caso1
os.system('clear')
BASE_SAVE_DIR = '/home/tparente/danilo/mestrado/github/IOC5811/lista5/outputs/novoteste/'

def plotar_caso(data, caso, savefig=False):

    print('Plotting figure 1 - Buoyancy frequency profiles')
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15,12))

    ax[0][0].plot(data['N']*10e3, data['z'], 'r')
    ax[0][0].set_title('Buoyancy Frequency Profile', fontsize=15)
    ax[0][0].set_xlabel(r'N [$ 10^{-3} rad s^{-1}$]', fontsize=15)
    ax[0][0].set_ylabel('Depth [m]', fontsize=15)

    ax[0][0].grid()     # turn on grids

    ax[0][1].plot(data['Nz']*10e5, data['z'], 'b')
    ax[0][1].set_title('Buoyancy Frequency Vertical Gradient Profile', fontsize=15)
    ax[0][1].set_xlabel(r'dN/dz [$ 10^{-5} rad s^{-1}$]', fontsize=15)
    ax[0][1].set_ylabel('Depth [m]', fontsize=15)

    ax[0][1].grid()     # turn on grids

    ax[1][0].plot(data['U'], data['z'], 'r')
    # ax[0].plot([0,0], [0, -data['H']], 'k')
    ax[1][0].set_title('Background Velocity Vertical Profile', fontsize=15)
    ax[1][0].set_xlabel(r'U  [$ m s^{-1}$]', fontsize=15)
    ax[1][0].set_ylabel('Depth [m]', fontsize=15)

    ax[1][0].grid()     # turn on grids

    ax[1][1].plot(data['Uz']*10e3, data['z'], 'b')
    # ax[1].plot([0,0], [0, -data['H']], 'k')
    ax[1][1].set_title('Background Velocity Vertical Gradient Profile', fontsize=15)
    ax[1][1].set_xlabel(r'dU/dz  [$10^{-3} s^{-1}$]', fontsize=15)
    ax[1][1].set_ylabel('Depth [m]', fontsize=15)

    ax[1][1].grid()     # turn on grids

    if savefig == True:
        plt.savefig(BASE_SAVE_DIR+'caso'+caso+'_fig1.png')


    print('Plotting figure 2 - Potential Vorticity')
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8,8))

    ax.plot(data['Qy'], data['z'], 'g')
    # ax.plot([0,0], [0, -data['H']], 'k')
    ax.set_xlim([1.00000000e-11, 3e-11])
    # ax.set_ylim([-data['H'],0])
    ax.set_title('Potential Vorticity Meridional Gradient', fontsize=15)
    ax.set_xlabel(r'Qy in $(m s)^{-1}$', fontsize=15)
    ax.set_ylabel('Depth [m]', fontsize=15)
    ax.grid()

    if savefig == True:
        plt.savefig(BASE_SAVE_DIR+'caso'+caso+'_fig2.png')
    #
    # print('\n')
    # print('Topographic Effect over Instability')
    # print('\n')
    # print('Plotting figure 3 - phase speed and growth rate')
    #
    # fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15,12))
    #
    # ax[0][0].plot(data['k'].T*1000, -data['cr'], 'b')
    # ax[0][0].set_xlim([0, 0.05])
    # ax[0][0].set_ylim([0, 6])
    # ax[0][0].set_title('Phase Speed', fontsize=15)
    # ax[0][0].set_ylabel(r'Phase speed in $cm^{-1}$', fontsize=15)
    #
    # ax[0][0].grid()
    #
    # ax[0][1].plot(data['k'].T*1000, data['sig'], 'r')
    # ax[0][1].set_xlim([0, 0.05])
    # ax[0][1].set_ylim([0, 0.016])
    # ax[0][1].set_ylabel(r'Growth Rate in  $days^{-1}$', fontsize=15)
    # ax[0][1].set_xlabel(r'Wavenumber in $km^{-1}$', fontsize=15)
    # ax[0][1].grid()
    #
    # ax[1][0].plot(data['Pmax'], data['z'], 'r')
    # ax[1][0].set_title('Most Unstable Mode Amplitude', fontsize=15)
    # ax[1][0].set_xlabel('P mode amplitude', fontsize=15)
    # ax[1][0].set_ylabel('Depth [m]', fontsize=15)
    # ax[1][0].grid()
    #
    # nz = data['nz']-1
    # Pphmax = data['Pphmax']-data['Pphmax'][nz]
    # ax[1][1].plot(Pphmax[0,:,:], data['z'], 'b')
    # ax[1][1].set_title('Phase Relative to the Bottom Value', fontsize=15)
    # ax[1][1].set_xlabel('Phase in degress', fontsize=15)
    # ax[1][1].set_ylabel('Depth [m]', fontsize=15)
    # ax[1][1].grid()
    #
    # if savefig == True:
    #     plt.savefig(BASE_SAVE_DIR+'caso'+caso+'_fig3.png')

def plotar_topografia(caso, savefig=False):
    """
        plotar no mesmo grafico a diferença considerando s maior ou menos que 0

        plotar grafico com s = 0
    """
    # import .mat file and datas to plot
    snulo  = sio.loadmat('/home/tparente/danilo/mestrado/github/IOC5811/lista5/data/plotado/caso'+str(caso)+'_s0.mat')
    smaior = sio.loadmat('/home/tparente/danilo/mestrado/github/IOC5811/lista5/data/plotado/caso'+str(caso)+'_sPos.mat')
    smenor = sio.loadmat('/home/tparente/danilo/mestrado/github/IOC5811/lista5/data/plotado/caso'+str(caso)+'_sNeg.mat')

    print('\n')
    print('Topographic Effect over Instability')
    print('\n')
    print('Plotting figure 3 - phase speed and growth rate')

    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(15,12))

    ax[0][0].plot(snulo['k'].T*1000, -snulo['cr'], 'b', label='s = 0')
    ax[0][0].plot(smaior['k'].T*1000, -smaior['cr'], 'g', label='s > 0')
    ax[0][0].plot(smenor['k'].T*1000, -smenor['cr'], 'r', label='s < 0')
    ax[0][0].set_xlim([0, 0.05])
    ax[0][0].set_ylim([0, 6])
    ax[0][0].set_title('Phase Speed', fontsize=15)
    ax[0][0].set_ylabel(r'Phase speed in $cm^{-1}$', fontsize=15)

    ax[0][0].grid()

    ax[0][1].plot(snulo['k'].T*1000, snulo['sig'], 'b', label='s = 0')
    ax[0][1].plot(smaior['k'].T*1000, smaior['sig'], 'g', label='s > 0')
    ax[0][1].plot(smenor['k'].T*1000, smenor['sig'], 'r', label='s < 0')
    ax[0][1].set_xlim([0, 0.05])
    ax[0][1].set_ylim([0, 0.016])
    ax[0][1].set_title('Growth Rate', fontsize=15)
    ax[0][1].set_ylabel(r'Growth Rate in  $days^{-1}$', fontsize=15)
    ax[0][1].set_xlabel(r'Wavenumber in $km^{-1}$', fontsize=15)
    ax[0][1].grid()

    ax[1][0].plot(snulo['Pmax'], snulo['z'], 'b', label='s = 0')
    ax[1][0].plot(smaior['Pmax'], smaior['z'], 'g', label='s > 0')
    ax[1][0].plot(smenor['Pmax'], smenor['z'], 'r', label='s < 0')
    ax[1][0].set_title('Most Unstable Mode Amplitude', fontsize=15)
    ax[1][0].set_xlabel('P mode amplitude', fontsize=15)
    ax[1][0].set_ylabel('Depth [m]', fontsize=15)
    ax[1][0].grid()

    nz = snulo['nz']-1
    Pphmax = snulo['Pphmax']-snulo['Pphmax'][nz]
    ax[1][1].plot(Pphmax[0,:,:], snulo['z'], 'b', label='s = 0')

    Pphmax = smaior['Pphmax']-smaior['Pphmax'][nz]
    ax[1][1].plot(Pphmax[0,:,:], smaior['z'], 'g', label='s > 0')

    Pphmax = smenor['Pphmax']-smenor['Pphmax'][nz]
    ax[1][1].plot(Pphmax[0,:,:], smenor['z'], 'r', label='s < 0')

    ax[1][1].set_title('Phase Relative to the Bottom Value', fontsize=15)
    ax[1][1].set_xlabel('Phase in degress', fontsize=15)
    ax[1][1].set_ylabel('Depth [m]', fontsize=15)
    ax[1][1].grid()

    ax[1][1].legend(loc='lower right')

    if savefig == True:
        plt.savefig(BASE_SAVE_DIR+'caso'+str(caso)+'_fig3.png')

##############################################################
#                     Main Code                              #
##############################################################
caso = 5

while caso != 0:

    caso = input('Digite o caso que deseja plotar (1 a 5) ou digite 0 para encerrar: ')
    if caso == 0:
        break

    savefig =input('Deseja salvar (1) ou visualizar (0)? ')

    casoName = '/home/tparente/danilo/mestrado/github/IOC5811/lista5/data/plotado/caso' + str(caso) + '_s0.mat'

    data = sio.loadmat(casoName)
    plotar_caso(data, str(caso), savefig=savefig)
    plotar_topografia(caso, savefig=savefig)


print("\n")
print('Usando convert -trim para remover os espaços branco em excesso...')
lista = glob.glob(BASE_SAVE_DIR+'*.png')
for fname in lista:
    os.system('convert -trim %s %s' % (fname, fname))
