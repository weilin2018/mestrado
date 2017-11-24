import numpy as np
import matplotlib.pyplot as plt


L = 50 #metros
l = 10
H = 30 #metros
g = 9.8 # aceleracao
Am = 0.5 # amplitude da onda
m = [0,1,2] # modos
c = np.sqrt(g*H) # velocidade de fase

x = np.linspace(0,L,num=100) # espaco
t = np.linspace(0,100,num=100) # tempo
freq = 0.1 # por seg

def calc(modo):   
    raz = (g*Am*modo*np.pi)/(L*c*np.sqrt(((modo*np.pi)/(L))**2 + (3*np.pi/l)**2))
    arg = (modo*np.pi*x)/L

    eta = Am*np.cos((modo*np.pi*x)/L)*np.cos(freq*t)
    u = raz*np.sin(arg)*np.sin(freq*t)

    return u, eta

def plotar(axes,u,eta,style):

    axes.plot(x/50,eta/30,style,label=modo-1)
    axes.set_ylabel(r'$\frac{\eta}{H}$',fontsize=20)

    axes.set_xlabel(r'$\frac{x}{L}$',fontsize=20)

    axes.set_title(u'Modos de Oscilação Natural')

    axes.grid()

    plt.legend(loc='best')

fig, axes = plt.subplots()

modo = 2*m[0]+1

u,eta=calc(0)

plotar(axes,u,eta,style='--')

modo = m[1]+1

u,eta=calc(modo)

plotar(axes,u,eta,style='-.')

modo = m[2]+1

u,eta=calc(modo)

plotar(axes,u,eta,style='-')


plt.show()
