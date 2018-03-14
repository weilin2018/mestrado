# code made by Guilherme Seki

import numpy as np
import matplotlib.pyplot as plt
import sys
reload(sys)
sys.setdefaultencoding('utf8')
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

def cart2pol(x, y):
    """Convert from Cartesian to polar coordinates.

    Example
    -------
    >>> theta, radius = pol2cart(x, y)
    """
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)

    return theta, radius

def equation(Av=1e-2, ug=0.1, fo=-1e-4):
	he=np.sqrt(Av/abs(fo))
	z=np.linspace(0,6*he,80)
	u=ug*(1-np.exp(-z/he)*np.cos(z/he))
	s=fo/abs(fo)
	v=s*ug*np.exp(-z/he)*np.sin(z/he)
	return u,v

def compass( arrowprops=None):
    """
    Compass draws a graph that displays the vectors with
    components `u` and `v` as arrows from the origin.

    Examples
    --------
    >>> import numpy as np
    >>> u = [+0, +0.5, -0.50, -0.90]
    >>> v = [+1, +0.5, -0.45, +0.85]
    >>> compass(u, v)
    """
    u,v = equation()

    angles, radii = cart2pol(u, v)

    fig, ax = plt.subplots(subplot_kw=dict(polar=True))

#    kw = dict(arrowstyle="->", color='k')
#    if arrowprops:
#        kw.update(arrowprops)
#    [ax.annotate("", xy=(angles, radii), xytext=(0, 0),
#                 arrowprops=kw) for
#     angle, radius in zip(angles, radii)]
    ax.plot(angles, radii, color='k', lw=2)
    ax.set_ylim(0, np.max(radii))
    ax.set_title('Hodógrafo da Espiral de Ekman de fundo (HS) (m/s)')
    plt.show()

    return fig, ax

def perfil(Av=4*1e-2,theta=-np.pi/6,l=50*1e3,y=0):
	fo=2*7.2921*(1e-5)*np.sin(theta)
	he=np.sqrt(Av/abs(fo))
	ug=np.exp(-y**2/l**2)
	z=np.linspace(0,6*he,60)

	u=ug*(1-np.exp(-z/he)*np.cos(z/he))
	s=fo/abs(fo)
	v=s*ug*np.exp(-z/he)*np.sin(z/he)
	fig=plt.figure()
	ax = fig.add_subplot(111)
	ax.grid(True)
	ax.set_ylabel('Profundidaze z/Uo')
	ax.set_xlabel('Velocidades (u,v)/Uo (m/s) (u-preto,v-vermelho)')
	ax.plot(u,z/he,color='k')
	ax.plot(v,z/he,color='r')
	ax.set_title('Perfil de velocidades normalizadas')
	plt.show()
	return u,v

def exec6(l=100*1e3,s=-1,y=50*1e3,Av=4*1e-2,theta=-np.pi/6):#perfil e bombeamento
	fo=2*7.2921*(1e-5)*np.sin(theta)
	y=np.arange(-100*1e3,100*1e3,1)
	he=np.sqrt(Av/abs(fo))
	we=he*s*y*np.exp(-y**2/l**2)
	ug=np.exp(-y**2/l**2)

	fig, ax = plt.subplots()
	newax = ax.twiny()

	# Make some room at the bottom
	fig.subplots_adjust(bottom=0.20)

	# I'm guessing you want them both on the bottom...
	newax.set_frame_on(True)
	newax.patch.set_visible(False)
	newax.xaxis.set_ticks_position('bottom')
	newax.xaxis.set_label_position('bottom')
	newax.spines['bottom'].set_position(('outward', 40))
	ax.grid(True)
	ax.plot(we/1e6,y/1e3,color='r')
	newax.plot(ug,y/1e3,color='b')

	ax.set_xlabel('We (Sv)',color='r')
	newax.set_xlabel('ug/Uo',color='b')
	ax.set_ylabel('y (km)')
	ax.set_title('Perfil de corrente geotrófica (azul) e Bombeamento de Ekman (vermelho)')

	plt.show()

def exc7(s=-1,theta=2*np.pi*25/360):
	fo=2*7.2921*(1e-5)*np.sin(theta)
	H = 2*np.pi*np.sqrt(Av/abs(fo))
	w, h = l, c;
	Matrix = [[0 for x in range(w)] for y in range(h)]
	a1 = -1 + (Av[i]-Av[i-1])/Av[i]*2
	a2 = 2 + j(Z**2)*s*abs(fo)/Av[i]
	a3 = -1 - (Av[i]-Av[i-1])/Av[i]*2



compass()
