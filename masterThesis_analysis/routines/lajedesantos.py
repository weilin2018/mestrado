'''
Próximo passo:

'''
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import pandas as pd
import os
import pickle
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('masterThesisPack/')

import masterThesisPack as oceano

def clc():
    import os

    os.system('clear')

def intdir2uv(intensity,direction,dmag,angRot):
    """Convert intensity and direction to zonal and meridional velocities component.

    this function also corrects magnetic declination and, if wanted, also
    rotate the vectors from zonal/meridional to along/cross shore velocities.

    Parameters
    ----------
    intensity : array_like
        Description of parameter `intensity`.
    direction : array_like
        Description of parameter `direction`.
    dmag : float
        Description of parameter `dmag`.
    angRot : float
        Description of parameter `angRot`.

    Returns
    -------
    type
        Description of returned object.

    """

    d = direction
    i = intensity

    d = d + dmag

    d = d * (np.pi/180.)

    u = i * np.sin(d)
    v = i * np.cos(d)

    return u,v

def dirmag2uv(direction,speed,decl_mag,ref_direcao):
    """ Transforma direcao (referenciada no LESTE - circulo trigonometrico ou Referenciada no NORTE) e magnitude/speed/velocidade em componentes U
    (leste-oeste) e V (norte-sul) - CHECK! Tudo correto aqui.
    ATENCAO!!! Se a sua direcao esta referenciada em NORTE, use: u = seno(dir)*speed e v=cos(dir)*speed
    Se a sua direcao esta referenciada em LESTE (circulo trignometrico) use: u = cosseno(dir)*speed e v=seno(dir)*speed.
    ref_direcao = 'trignometrico' ou 'norte'. 'trignometrico' se o dir=0 eh no LESTE. 'norte' se dir=0 eh no NORTE.
    """
    direction = direction + decl_mag
    # print("direction = " + str(direction))
    # direction = np.mod(direction, 360)
    # print("direction = " + str(direction))
    direction = direction * np.pi / 180
    if ref_direcao=='trigonometrico':
        u = np.cos(direction)*speed
        v = np.sin(direction)*speed
    if ref_direcao=='norte':
        u = np.sin(direction)*speed # esta assim no BoB.
        v = np.cos(direction)*speed
    return u,v

# stickplot
def stickplot(df,ax):
    """Create a stickplot.

    With the u and v components given as argument in pd.DataFrame df,
    this function plot a stickplot, using MPL quiver.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe containing wu and wv components and datetimeindex.
    ax : matplotlib.axis
        axis to stickplot

    Example
    -------
    >>>
    >>>

    Credits
    -------
    Stephane Raynaud
    http://permalink.gmane.org/gmane.comp.python.matplotlib.general/24155
    """

    # creating the date axis
    # dates = pd.to_datetime(df.index)
    # extracting components from dataframe
    u = df['wu'].values
    v = df['wv'].values

    # calculating speed
    spd = np.sqrt(u**2 + v**2)
    maxSpd = np.nanmax(spd)

    # plotting
    # fig, ax = plt.subplots()

    qiv = ax.quiver(df.index, [[0]*len(df)], u, v, headlength=0, headwidth=0, headaxislength=0 )
    key = ax.quiverkey(qiv, 0.25, 0.75, maxSpd, "%0.2f $m^{2}s^{-1}$"%(maxSpd), labelpos='N', coordinates='axes' )

    # plot a horizontal line in y=0.0
    ax.axhline(y=0.0,xmin=df.index[0],xmax=df.index[-1],linewidth=1.,color='black')

    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis_date()

    # ax.set_xticks(['2012-01-07', '2012-01-21', '2012-02-04', '2012-02-18', '2012-03-03', '2012-03-17', '2012-03-31'])

    return ax

# importar pickles
PICKLE_DIR='/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/pickles/'

cfsv = pickle.load(open(PICKLE_DIR+'cfsv2_laje.pickle', 'r')) # convencao oceanog.
laje = pickle.load(open(PICKLE_DIR+'laje_Meteor.pickle', 'r')) # convencao meteor.

direction = laje['direction'] #+ 180 # converting from meteorological to mathematical convention
intensity = laje['int']

clc()

''' me interessa visualizar os dados com vetores apontando para onde o vento vai'''

wu,wv = dirmag2uv(direction,intensity,-21.03,'norte')
# criar dataframe
data =  pd.DataFrame({'wu':wu,'wv':wv},index=pd.DatetimeIndex(laje['dates']))

# os dados não estão compreendo o mesmo periodo. Desta forma, vamos selecionar
data = data['2015-04-22':]
cfsv = cfsv['2015-04-22':]

# importante: os dados da Laje são horários, então é necessário um resample
data = data.resample('6H').mean()

#########################################################
##### plotando e calculando os valores estatísticos #####
#########################################################

def laje_statisticalAnaysis(laje,cfsv2,whichSerie):
    ''' '''
    # statistical analysis for Laje de Santos
    laje_cut = laje
    cfsv_cut = cfsv2

    skillWu = oceano.skill_willmott(laje_cut.wu.values, cfsv_cut.wu.values)
    skillWv = oceano.skill_willmott(laje_cut.wv.values, cfsv_cut.wv.values)

    corrWu = calculateCorr(laje_cut.wu.values, cfsv_cut.wu.values)[0]
    corrWv = calculateCorr(laje_cut.wv.values, cfsv_cut.wv.values)[0]

    mseWu  = calculateMSE(laje_cut.wu.values, cfsv_cut.wu.values)
    mseWv  = calculateMSE(laje_cut.wv.values, cfsv_cut.wv.values)

    # plot data and skill
    fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True)

    ax[0].plot(laje_cut.wu,label='Laje')
    ax[0].plot(cfsv_cut.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-10,10)

    wuText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text('2015-04-07', 7.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend(loc='lower left')

    ax[1].plot(laje_cut.wv,label='Laje')
    ax[1].plot(cfsv_cut.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-10,10)

    wvText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text('2015-04-07', 7.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend(loc='lower left')

    ax[2] = stickplot(laje_cut,ax[2])

    ax[3] = stickplot(cfsv_cut,ax[3])

    plt.suptitle('Laje de Santos (blue) v CFSv2 (red) \n 2015-04 to 2015-12 \n %s' % (whichSerie),fontsize=26)
    plt.show()

laje_statisticalAnaysis(data,cfsv,'')


###
# converter CFSv2 para direcao e intensidade (convencao meteorologica)

from math import atan2

u_ms,v_ms = cfsv['wu'].values,cfsv['wv'].values
wind_abs = np.sqrt(u_ms**2 + v_ms**2)

wind_dir_trig_to = []

for u,v,spd in zip(u_ms,v_ms,wind_abs):
    wind_dir_trig_to.append(atan2(v/spd, u/spd))

wind_dir_trig_to = np.asarray(wind_dir_trig_to)

wind_dir_trig_to_degrees   = wind_dir_trig_to * 180/np.pi
wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180

wind_dir_cardinal          = 90 - wind_dir_trig_from_degrees




df1 = pd.DataFrame({'direction':direction},index=pd.DatetimeIndex(laje['dates']))
df2 = pd.DataFrame({'direction':wind_dir_trig_from_degrees},index=cfsv.index)

df1 = df1['2015-04-22':]
df2 = df2['2015-04-22':]

fig, ax = plt.subplots()

ax.plot(df1.index, df1.direction.values, label='laje')
ax.plot(df2.index, df2.direction.values, label='cfsv')
