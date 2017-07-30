""" https://ocefpaf.github.io/python4oceanographers/blog/2015/02/09/compass/"""






"""
Funções:


    def exec2() - plotar hodógrafo
    def exec3() - plotar perfil
    def exec6() - perfil de distribuição meridional
        do bombeamento de ekman bentico

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib
matplotlib.style.use('ggplot')

import sys
reload(sys)
sys.setdefaultencoding('utf8')

def cart2pol(x, y):
    """Convert from Cartesian to polar coordinates.

    Example
    -------
    >>> theta, radius = pol2cart(x, y)
    """
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)

    return theta, radius

def calcular_UeV(Av=1e-2, ug=0.1, fo=-1e-4):
	he=np.sqrt(Av/abs(fo))
	z=np.linspace(0,6*he,80)
	u=ug*(1-np.exp(-z/he)*np.cos(z/he))
	s=fo/abs(fo)
	v=s*ug*np.exp(-z/he)*np.sin(z/he)
	return u,v,z

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
    u,v,z = calcular_UeV()

    angles, radii = cart2pol(u, v)

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

    # kw = dict(arrowstyle="->", color='k')
    # if arrowprops:
    #    kw.update(arrowprops)
    # [ax.annotate("", xy=(angles, radii), xytext=(0, 0),
    #             arrowprops=kw) for
    # angle, radius in zip(angles, radii)]

    ax.plot(angles, radii, color='k', lw=2)
    ax.set_ylim(0, np.max(radii))
    ax.set_title('Hodógrafo da Espiral de Ekman de fundo (HS) (m/s)')
    plt.show()


def plotar_3d(ax, u, v, Z):
    """ """

    for x,y,z in zip(u, v, Z):
        ax.plot3D(x,y,z,'k')



def exercicio2():
    """
        plotar hodógrafo do vetor velocidade em função da profundidade

        Av: coeficiente de viscosidade vertical
        ug: velocidade do interior geostrófico
        fo: parâmetro de coriolis
    """

    compass()
