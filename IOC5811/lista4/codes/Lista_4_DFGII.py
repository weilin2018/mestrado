import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
#from mpl_toolkits.basemap import Basemap
from scipy.sparse         import linalg
from scipy                import interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#%*******************************************************
#%                    PROGRAM MPSI_S
#%
#% calcula numericamente a solução de Stommel para 
#% campos de rotacional do vento e de topografia 
#% de fundo teoricos especificados
#%
#% resolução por esquema iterativo em forma nao-dimensional
#% bacia oceanica quadrada de dimensao L=1e6
#%
#%
#%             por Ilson da Silveira 
#%
#%******************************************************


tol=1e-3 # Tolerância do erro médio para a iteração
L=1e6; 
H=5e3;     
fo=7.3e-5;
beta=2e-11;
he=np.sqrt(2e-2/fo) # espessura da camada de Ekman teórica [Av= 1e-2 m2/s]
r= fo*he/H; #Parametro de fricção
to=1e-4 #tensão de cisalhamento/densidade [m2/s2]

# A nao-dimensionalizacao é feita assumindo que a escala da
# funcao de corrente geostrofica psi é dada pela rel. de Sverdrup:
#                psi= tau0/(rho*H*beta)

r=r/(beta*L)
ff=fo/(beta*L) #Normalizando r e f

#Construindo grade

nmax=101
dx=0.01




x=np.arange(0,nmax)*dx
y=x.T

[xg,yg]=np.meshgrid(x,y)

xg=np.flipud(xg)
yg=np.flipud(yg)

yL=max(y)



#Campo de rotacional do vento nao-dimensionalizacao

#curl=np.zeros(xg.shape) #sem vento
curl=-np.pi/yL*np.sin(np.pi*yg/yL) #para giro único
#curl=-np.pi/yL*np.sin(2*np.pi*yg/yL) #para giros duplos
#curl=-np.pi/yL*np.exp(-(xg-yL/2)*(xg-yL/2));  #campo gaussiano
#curl=-np.pi/yL*np.ones((xg.shape));              # campo uniforme

#campo de topografia de fundo n-dimensional  assumindo H=5km
#ou seja hb=b/H

#hb=np.zeros(xg.shape) #fundo plano
#hb=-0.275*yg   #fundo inclinado linearmente em y;
#hb=-3*0.275*xg # fundo inclinado linearmente em x;
hb=0.1*np.exp(-30*(xg-yL/2)*(xg-yL/2)) # cordilheira meso-atlantica
#hb= 0.1*np.exp(-30*(xg-yL/2)*(xg-yL/2))*np.exp(-30*(yg-yL/2)*(yg-yL/2)); # banco submarino

#Campo de Vorticidade Potencial Ambiente

fac=1. #usar fac=0. para plano f e fac=1. para os outros casos

y0=0.5; #latitude central

B=fac*(yg-y0)+ff*hb #função VP ambiente 

[By,Bx]=np.gradient(B,dx,dx)

By=-By
####################################################################################################


#Plotagem

#Campo de VP ambiente

b1=0.1*np.floor(B.min()*10);
b2=0.1*np.ceil(B.max()*10);

lb=np.arange(b1,b2+0.05,0.05)

plt.contourf(x,y,np.flipud(B),levels=lb); plt.colorbar()

plt.axis('square')
plt.xlabel(u'x [em 10^6 m]')
plt.ylabel(u'y [em 10^6 m]')
plt.title(u'A Vorticidade Potencial Ambiente')




#Plotando topografia 

hbflip=np.flipud(hb)

lt=np.arange(hb.min(),hb.max()+0.005,0.005)

fig=plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2= fig.add_subplot(1, 2, 2, projection='3d')


plt.subplot(121)
plt.contourf(x,y,hbflip,levels=lt);
plt.xlabel(u'x [em 10^6 m]')
plt.ylabel(u'y [em 10^6 m]')
plt.title(u'A topografia de fundo em 2D (z[4*10²m])')
plt.colorbar()
plt.axis('square')


[xx,yy]=np.meshgrid(x,y)

plt.subplot(122)
ax = fig.add_subplot(122, projection='3d')
surf=ax.plot_surface(xx, yy, 4e2*hbflip, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.axis('square')
plt.title(u'A topografia de fundo em 3D (z[m])')
plt.xlabel(u'x [em 10^6 m]')
plt.ylabel(u'y [em 10^6 m]')


psi=np.zeros(xg.shape)


#nsel1 ~ 0:nmax-2
#nsel2 ~ 2:nmax
#nsel  ~ 1:nmax-1

#processo de iteração

for i in range(2000):

    psi_old=psi
        
    av=psi[0:nmax-2,1:nmax-1] + psi[2:nmax,1:nmax-1] +psi[1:nmax-1,0:nmax-2] + psi[1:nmax-1,2:nmax]
        
    F=1/r*(curl[1:nmax-1,1:nmax-1]*dx*dx - psi[1:nmax-1,2:nmax]*(By[1:nmax-1,1:nmax-1])*dx +psi[2:nmax,1:nmax-1]*(Bx[1:nmax-1,1:nmax-1])*dx)
    
    psi[1:nmax-1,1:nmax-1]=(av-F)/(4+(dx/r)*(By[1:nmax-1,1:nmax-1]-Bx[1:nmax-1,1:nmax-1]))
        
    crit=(np.abs(psi-psi_old));crit=crit.max()
        
    if crit <= tol:
        continue
    else:
        break


#calculando a vorticidade relativa
  

zeta=(psi[0:nmax-2,1:nmax-1]+psi[2:nmax,1:nmax-1]+psi[1:nmax-1,0:nmax-2]+psi[1:nmax-1,2:nmax]-4*psi[1:nmax-1,1:nmax-1])/dx/dx

#plotando função de corrente

p1=0.1*np.floor(psi.min()*10)
p2=0.1*np.ceil(psi.max()*10)

psiflip=np.flipud(psi)

plt.figure()

lp=np.arange(p1,p2,0.1)

plt.contourf(x,y,psiflip,levels=lp)
plt.colorbar()
plt.axis('square')
plt.xlabel(u'x [em 10^6 m]')
plt.ylabel(u'y [em 10^6 m]')
plt.title(u'Solução de Stommel')



 


