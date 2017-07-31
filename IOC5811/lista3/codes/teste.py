# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw
from mpl_toolkits.mplot3d import Axes3D as plt3
import matplotlib.gridspec as gridspec
import numpy.linalg as alg
from matplotlib import rc

def ekman_uv(rho,fo,he,z,taux,tauy):
    tau_u = (taux*np.cos((z/he)-np.pi/4) - tauy*np.sin((z/he)-np.pi/4))
    tau_v = (taux*np.sin((z/he)-np.pi/4) + tauy*np.cos((z/he)-np.pi/4))
    ue=(fo/np.abs(fo))*(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_u
    ve=(np.sqrt(2)/(rho*fo*he))*(np.exp(z/he))*tau_v
    return ue,ve

def ekman_plot(U,V,z,ue,ve,ze,Av,He,item='xx',Avt='xxxx',stp1=20,stp2=5):

    lw = 1.5

    U,V,ue,ve = 100*U,100*V,100*ue,100*ve

    X3 = np.vstack([np.zeros(U.shape),U])
    Y3 = np.vstack([np.zeros(V.shape),V])
    Z3 = np.vstack([z,z])

    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2)
    fig.suptitle('2-%s: Av(z) %s'%(item,Avt),fontsize=16)
    fig.text(0.5, 0.025,'Profundidade efetiva da Camada de Ekman = %.2f m'%(He),fontsize=18,
                horizontalalignment='center',verticalalignment='center')

    ax1 = fig.add_subplot(gs[0,:])
    #ax1.quiver(np.zeros(U[::stp1].shape),
    #            np.zeros(U[::stp1].shape),
    #            U[::stp1],V[::stp1],units='xy',
    #            angles='xy', scale=1,
    #            scale_units='xy')
    ax1.plot(X3[:,::stp1],Y3[:,::stp1],'g')
    ax1.plot(U,V,'g',label='Numerica')
    ax1.plot(np.vstack([np.zeros(ue.shape),ue])[:,::stp1],
            np.vstack([np.zeros(ve.shape),ve])[:,::stp1],'b',linewidth=lw)
    ax1.plot(ue,ve,'b',label='Analitica')
    ax1.grid('on')
    ax1.legend(loc='lower right')

    ax2 = fig.add_subplot(gs[1,0])
    ax2.plot(U,z,'g',label='U')
    ax2.plot(V,z,'--g',label='V')
    ax2.plot(ue,ze,'b',label='Ue')
    ax2.plot(ve,ze,'--b',label='Ve')
    ax2.legend(loc='lower right')


    ax3 = fig.add_subplot(gs[1,1], projection='3d')
    for x3,y3,z3 in zip(U[::stp2],V[::stp2],z[::stp2]):
        ax3.plot3D([0,x3],[0,y3],[z3,z3],'g')
    ax3.plot(U,V,z,'g')
    for x3,y3,z3 in zip(ue[::stp2],ve[::stp2],z[::stp2]):
        ax3.plot3D([0,x3],[0,y3],[z3,z3],'b')
    ax3.plot(ue,ve,ze,'b')
    ax3.view_init(elev=23., azim=85)
    plt.savefig('/Users/fernanda/Desktop/graficos_julio/vels_%s.png'%(item),dpi=250)

    # plot perfil Av
    fig2 = plt.figure(figsize=(12,8))
    fig2.suptitle('2-%s: Av(z) %s'%(item,Avt),fontsize=16)
    plt.plot(Av,z,'r')
    plt.xlabel(r'Coeficiente de Atrito Vertical (m$^{2}/$s$.10^{-4}$)')
    plt.ylabel('Profundidade (m)')
    plt.savefig('/Users/fernanda/Desktop/graficos_julio/perfil_%s.png'%(item),dpi=250)



def near(dat,val,how_many=1):
    dif = np.abs(dat-val)
    idx = np.argsort(dif)
    return dat[idx][:how_many]


##########################################
########### MUITO IMPORTANTE #############
##########################################
#Definindo todas as contantes
zl=0
l=2
delz = 0.2
ang=-23
f0 = sw.f(ang)
s = f0/np.abs(f0)
uar=10
Tau=1.225*0.0015*uar**2
rho=1027
hE = np.sqrt((2*0.5e-2)/np.abs(f0))
pE= np.pi*hE

#Definindo os vetores verticais, z e JJ
z = np.arange(0,-5*hE,-delz)
JJ = np.arange(z.shape[0])+1

#Definindo os vários Av(z)

#Av Constante
prop = ['c)','Constante em z']
Av = np.zeros(z.shape[0])+0.5e-2

#Av Linear
#prop = ['d)','varia linearmente']
#Av = (40+0.39*z)*10e-4

# Av Exponencial
#prop = ['e)','varia exponencialmente']
#Av = (2.5+40*np.exp(z/10.5))*10e-4

# Av Gaussiano
#prop = ['f)','modulado por envelope gaussiano de largura l']
#zl=0
#l=hE
#Av = (2.5+40*np.exp(-(z-zl)**2/(l**2)))*10e-4

# Av Real
#prop = ['g)','real [valores segundo Chao et al. (2001)]']
#z1 = np.array([0,10,15,27,35,50,58,65,78,85,100,150,200])
#Av1 = np.array([40,25,21,27,40,41,45,32,13,5,2.5,1.8,1])
#z=np.arange(0,200.2,0.2)
#Av=np.interp(z,z1,Av1)
#Av*=10e-4
#z=z*-1.
#JJ = np.arange(z.shape[0])+1

# Definindo dAv para cálculo dos Gamas
dAv = np.diff(Av)/np.diff(JJ)

# Definindo os tres gamas (diagonais da matriz A)
G1 = -1+((dAv/Av[:-1])*(delz/2))
G2 = (np.tile(2,z.shape[0]))+1j*((delz**2)*s*np.abs(f0))/Av
G3 = (-1-((dAv/Av[:-1])*(delz/2)))

# Definido o vetor solucao
S=np.append(np.zeros(z.shape[0]-1),(2*Tau*delz)/(rho*(Av[-1]+Av[-2])))

# Definindo a matriz A com os tres gamas e adicionando as condicoes
# de contorno na primeira e ultima linhas
A = np.diag(G2, 0)+np.diag(G1, -1)+np.diag(G3, 1)
A[0,0],A[-1,-1],A[0,1],A[-1,-2] = 1,1,0,-1

# Calculando as velocidades invertendo a matriz A e multiplicando
# matricialmente pelo vetor das solucoes S
A2=np.linalg.inv(A)
Vzao = np.dot(A2,S)

# Extraindo U e V das solucoes sabendo do enunciado que
# V = u+v*i
U,V = np.real(Vzao)[::-1],np.imag(Vzao)[::-1]

# Caluclo velocidades do modelo de ekman
ue,ve = ekman_uv(rho,f0,np.sqrt((2*0.5e-2)/np.abs(f0)),z,Tau,0)

# Truncando todas as solucoes de Av para os primeiros 66 metros de coluna de agua
U,V,ue,ve,z,Av=U[:332],V[:332],ue[:332],ve[:332],z[:332],Av[:332]

# Calculando as velocidades da solucao numerica e analitica (teorica)
SPD,SPDe = np.sqrt(U*U + V*V),np.sqrt(ue*ue + ve*ve)

# Encontrando a razao entre a velocidade de superficie e
# a velocidade na profundidade efetiva (pE)
rat=SPDe[0]/SPDe[np.argwhere(z==near(z,-pE,1))]

# Calculando velocidade teorica encontrada com a relacao de velocidades
# de Ekman para cada um dos Av(z)
vel_teoric=SPD[0]/rat

# Indice da profundidade da camada efetiva para cada um dos experimentos
prof_camada = np.argwhere(SPD==near(SPD,vel_teoric[0],1))

# Profundidade efetiva calculada com razao entre velocidades
He = z[prof_camada]
print 'Profundidade efetiva --> %f m'%(He)

# Plotando a porra toda!!
ekman_plot(U,V,z,ue,ve,z,Av,He,item=prop[0],Avt=prop[1],stp1=5,stp2=5)
plt.close('all')
