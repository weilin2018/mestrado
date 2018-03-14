# -*- coding: utf-8 -*-

%reset -f
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw
import sys  
reload(sys)  
sys.setdefaultencoding('utf8')


def calculateUandV(razao, ax):
    Av     = 0.014
    a      = np.sqrt(sw.f(23.)/(2*Av))
    deltaE = np.pi/a
    Ho = razao * deltaE
    z = np.arange(0,Ho,0.1)
    tau = (deltaE)/(Av*np.pi)
    zeta = Ho - z
    # calculo dos coeficientes A e B segundo equacoes 
    A = tau*(np.cosh(a*Ho) * np.cos(a*Ho) + np.sinh(a*Ho)*np.sin(a*Ho))/(np.cosh(2*a*Ho) + np.cos(2*a*Ho))
    B = tau*(np.cosh(a*Ho) * np.cos(a*Ho) - np.sinh(a*Ho)*np.sin(a*Ho))/(np.cosh(2*a*Ho) + np.cos(2*a*Ho))

    u = A * np.sinh(a*zeta) * np.cos(a*zeta) - B * np.cosh(a*zeta) * np.sin(a*zeta)  
    v = A * np.cosh(a*zeta) * np.sin(a*zeta) + B * np.sinh(a*zeta) * np.cos(a*zeta)

    # calculo do angulo formado enter vento e corrente superficial, segundo equacao ()
    alpha = np.arctan((np.sinh(2*Ho*a) - np.sin(2*Ho*a))/(np.sinh(2*Ho*a) + np.sin(2*Ho*a)))
    degree = 180*alpha/np.pi

    return u,v,degree

def cart2pol(x, y):
    """Convert from Cartesian to polar coordinates.

    Example
    -------
    >>> theta, radius = pol2cart(x, y)
    """
    radius = np.hypot(x, y)
    theta = np.arctan2(y, x)

    return theta, radius

def compass(u,v, degree, razao, ax, color='b'):
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

    angles, radii = cart2pol(u, v)
    
#    kw = dict(arrowstyle="->", color='k')
#    if arrowprops:
#        kw.update(arrowprops)
#    [ax.annotate("", xy=(angles, radii), xytext=(0, 0),
#                 arrowprops=kw) for
#     angle, radius in zip(angles, radii)]
    label=r'$\frac{H_o}{\delta_E} = %1.2f$ / $\alpha = %1.2f$'%(razao,degree)
    ax.plot(angles, radii, color=color, lw=2,linewidth=1,label=label)
    ax.plot(u[-1],v[-1])
    ax.set_ylim(0, np.max(radii))
    ax.set_title('Hodógrafo da Espiral de Ekman no HN')
    #plt.show()

    return ax


razoes = [0.1, 0.25, 1.25]
colors = ['b', 'r', 'k']

fig, ax = plt.subplots(subplot_kw=dict(polar=True),figsize=(15,7.5))

for razao,cor in zip(razoes, colors):
    u,v,angle = calculateUandV(razao,ax)
    ax = compass(u,v,angle,razao,ax,color=cor)


ax.legend(loc='lower left')
plt.show()
