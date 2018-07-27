# add some description here

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
from sklearn.metrics import mean_squared_error
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib import dates
# import datetime
from scipy.spatial import cKDTree

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('../artigoTGpack/')

import artigoTGpack as oceano

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################
# insert functions here
def locate_coord(ilat,ilon,lon,lat):
    """Find the nearest coordinate point to the [ilon,ilat] in [lon,lat].

    Parameters
    ----------
    ilat : float
        Latitude point.
    ilon : float
        Longitude point.
    lon : np.ndarray
        Matrix with all longitude coordinates.
    lat : np.ndarray
        Matrix with all latitude coordinates.

    Returns
    -------
    iss,jss : list
        Indexes of the closests points selected.

    """

    lo = lon.ravel()
    la = lat.ravel()

    coords = []

    for i,j in zip(la,lo):
        coords.append([i,j])

    coords = np.asarray(coords)

    locations_posi = [[ilat,ilon]]

    locs = np.asarray(locations_posi)

    tree = cKDTree(coords)

    dists,indexes = tree.query(locs,k=1)

    pontos = []

    for index in indexes:
        pontos.append(coords[index])

    pontos = np.asarray(pontos)

    ind = []

    for p in pontos:
        ind.append(np.where(lon == p[1]))

    ind = np.asarray(ind)

    iss, jss = [],[]

    for i,j in ind:
        iss.append(int(i))
        jss.append(int(j))

    return iss,jss

def load_bndo(DATA_DIR,kind='mean'):
    """Load and process data from BNDO dataset. If you set kind as 'mean',
    so this function also will remove the mean signal from the time series.

    Parameters
    ----------
    DATA_DIR : string
        Full path to the files.
    kind : string
        Which process perform, default is remove the mean signal.

    Returns
    -------
    observ : pd.DataFrame
        Pandas.DataFrame with all data included.

    """
    lfiles = glob.glob(DATA_DIR+'*')
    lfiles.sort()

    # ler, inicialmente, os dois primeiros arquivos para ampliar uma série
    files = lfiles[:2]

    file1 = pd.read_csv(files[0], skiprows=11, delimiter=';', names=['nivel', 'x'])
    file1.drop(file1.columns[len(file1.columns)-1], axis=1, inplace=True)

    file2 = pd.read_csv(files[1], skiprows=11, delimiter=';', names=['nivel', 'x'])
    file2.drop(file2.columns[len(file2.columns)-1], axis=1, inplace=True)

    # criar os dataframes
    dtRange = pd.date_range(start=file1.index[0], end=file1.index[-1], freq='H')
    df1 = pd.DataFrame({'nivel': file1['nivel'].values/100.}, index=dtRange)

    dtRange = pd.date_range(start="1997-02-02 00:00", end="1997-03-05 23:00", freq='H')
    df2 = pd.DataFrame({'nivel': file2['nivel'].values/100.}, index=dtRange)

    dtRange = pd.date_range(start='1997-02-01 00:00', end='1997-02-01 23:00', freq='H')
    df3 = pd.DataFrame({'nivel': np.ones(dtRange.shape[0])*np.nan}, index=dtRange)

    # concatenar as séries
    observ = pd.concat([df1, df3, df2])

    # controle de qualidade
    cond = observ['nivel'] > 4.
    observ[cond] = np.nan

    if kind=='mean':
        # removendo a média da série temporal
        observ['nivel'] = observ['nivel'] - observ['nivel'].mean()

    if kind=='original':
        # removendo a média da série temporal
        observ['nivel'] = observ['nivel']

    return observ

def load_ecom(SIMS_DIR,locs=[-23.,-44.01916667]):
    """Load ECOM product (elevation,time,lon and lat).

    Parameters
    ----------
    SIMS_DIR : string
        Full path to the directory where files are stored.
    locs : list
        Latitude and Longitude (must be in this order), to extract time series.

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with the elevation timeseries extracted at the given location.

    """

    ncin  = xr.open_dataset(SIMS_DIR)
    lon   = ncin['lon'].data
    lat   = ncin['lat'].data
    lon[lon == 0.] = np.nan
    lat[lat == 0.] = np.nan

    # locating the closest index from observed data location
    ilat = locs[0]
    ilon = locs[1]

    iss,jss = locate_coord(ilat,ilon,lon,lat)

    # import elevation
    elev = np.squeeze(ncin['elev'].data[:,iss,jss])
    time = ncin['time'].values

    df = pd.DataFrame({'elev':elev},index=pd.DatetimeIndex(time))

    return df

def plot(data1,data2,title):
    """Simple plot of two datasets.

    Parameters
    ----------
    data1 : pd.DataFrame
        Dataframe containing data from in situ observations.
    data2 : pd.DataFrame
        Dataframe containing product from ECOM simulation.
    title : string
        Title.
    """

    plt.ion()

    fig,ax = plt.subplots()

    ax.plot(data1,'k',label='bndo')
    ax.plot(data2,'r',label='ecom')

    ax.set_title(title,fontsize=14)

    plt.legend()

def statisticalParameters(data1,data2):
    """calculate skill, correlation coeficient and root mean squared error.

    Parameters
    ----------
    data1 : np.ndarray
        In situ observation timeseries.
    data2 : nd.ndarray
        Model product timeseries.

    Returns
    -------
    skill,corre,rmse : float
        Statistical parameters calculated.

    """

    skill = oceano.skill_willmott(data1,data2)
    corre = np.corrcoef(data1,data2)[0][1]
    rmse  = np.sqrt(mean_squared_error(data1,data2))

    return skill,corre,rmse

def removeTidalSignal(x,window_len,window='hanning'):

    # creating the window with maximum value normalized to one
    w = eval('np.'+window+'(window_len)')
    #
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    # apply filter
    y = np.convolve(w/w.sum(),s,mode='valid')

    return y

##############################################################################
#                               MAIN CODE                                    #
##############################################################################
# beginnig of the main code

BASE_DIR = oceano.make_dir()

DATA_DIR = BASE_DIR.replace('github/','artigo_data/BNDO/')
SIMS_DIR = BASE_DIR.replace('github/','artigo_data/simulacoes/2ndPhase/')

# define which simulation
run = 'pcse_elevation'

os.system('clear')

### ----------- reading BNDO data --------------- ###
observ = load_bndo(DATA_DIR+'1997/',kind='mean')

observ.plot(title='Serie completa de dados - BNDO')
plt.show()

# serie contem dados NaN. Do inicio da série até o dia 31/12, há Apenas
# um ponto sem dado, q pode ser interpolado de forma simples, para obter
# uma série temporal longa e válida para comparação. Esse procedimento é
# realizado abaixo:
first_period_observe = observ[:'1997-01-31'].copy()
first_period_observe.fillna(first_period_observe.mean(),inplace=True)

# ainda podemos pegar uma segunda parte da série temporal, sem nan values:
# second_period_observe = observ['1997-02-02':].copy()

### ---------- reading ECOM data ---------------- ###
fname  = '/'+run+'.cdf'
model  = load_ecom(SIMS_DIR+fname)

### --------- ajuste de maré: ecom x bndo ------- ###
time_data = model.index
time_product = time_data - pd.Timedelta(hours=3)

model = pd.DataFrame({'elev':model.elev.values},index=time_product)

# recortando os dados do modelo para bater com os dados de series temproais
# extraidos do BNDO
product_selected = model['1996-12-18':'1997-01-31'].copy()
# second_period_product = model['1997-02-02':].copy()

### ----------- plotar
plot(first_period_observe,product_selected,title='BNDO with hourly frequency')

### ---------- subsetting bndo ---------------- ###
# BNDO possui uma frequência horária, mas o ECOM possui uma saída de modelo
# a cada 6 horas. Desta forma, os vetores de dados terão tamanhos diferentes
# impossibilitando qualquer tipo de comparação (correlação, rsme, skill).
# Portanto, é necessário realizar uma subamostragem dos dados do BNDO para
# obtermos uma frequência a cada 6 horas.

# Analisando os shapes, vemos que o BNDO possui um vetor ~6x maior que o do ECOM:
first_period_observe.shape,product_selected.shape

# subamostrando os dados, pegando somente os dados a cada 6 horas
bndo_filtered = first_period_observe.iloc[::6,:]

# tomando os dados do BNDO a cada 6 horas [0,6,12,18], temos:
# se fizermos um teste tomando um resample de 6H do BNDO:
skill = oceano.skill_willmott(bndo_filtered.nivel.values,product_selected.elev.values)
corr  = np.corrcoef(bndo_filtered.nivel.values,product_selected.elev.values)[0][1]
skill,corr,rmse = statisticalParameters(bndo_filtered.nivel.values,product_selected.elev.values)

plot(bndo_filtered,product_selected,title='Both with 6-hourly frequency \n skill: %.2f corr: %.2f'%(skill,corr))
# figure saved in: /media/danilo/Danilo/mestrado/github/artigoTG/figures/2ndPhase/rerun00_v_bndo_.png

##### TESTING ZONE
"""
Realizei um teste aqui para verificar o que poderia estar reduzingo meu coeficientes
estatísticos: maré ou vento.

Modificacao:
Dottori sugeriu mudar a janela para 30 horas e usar um filtro Hanning ou Hamming.
A filtragem se mostrou muito melhor que apenas a média móvel, ao remover
totalmente o sinal da maré, obtendo somente o sinal de elevação do nível do mar
associado a maré meteorológica

Adicionei este sinal ao meu produto do modelo (reRUN00_ECOM) e calculei os
valores estatísticos novamente. Desta vez, pude aumentar os valores de skill para
0.88 e correlação para 0.89.

Com isso, chego a conclusão que, se eu conseguir melhorar como o modelo
representa o vento, eu talvez obtenha uma melhor validação.


"""
# removendo maré astronomica dos dados do bndo
mare_meteorologica = removeTidalSignal(bndo_filtered.nivel.values,window_len=5,window='hamming')
mare_meteorologica = pd.DataFrame({'MeteorTide_BNDO':mare_meteorologica[2:-2]},index=product_selected.index)


# maré meteorologica + ecom
met_ecom = mare_meteorologica.MeteorTide_BNDO.values + product_selected.elev.values
met_ecom = pd.DataFrame({'metBNDO_ECOM':met_ecom},index=product_selected.index)

# calcular novos parametros estatísticos
indNan = np.where(~np.isnan(met_ecom.metBNDO_ECOM.values))
skill2,corr2,rmse2 = statisticalParameters(met_ecom.metBNDO_ECOM.values[indNan],bndo_filtered.nivel.values[indNan])

"""
Um novo teste que pode ser feito é: remover o sinal da maré meteorológica da
série do BNDO (apenas por subtração) e comparar essa nova série com minha
saída de modelo.g

Isso é feito abaixo
"""

bndo_sem_mareMet = bndo_filtered.nivel.values - mare_meteorologica.MeteorTide_BNDO.values
bndo_semMet = pd.DataFrame({'BNDOsemMet':bndo_sem_mareMet},index=product_selected.index)

# usando o mesmo indNan já calculado acima
skill3,corr3,rmse3 = statisticalParameters(bndo_semMet.BNDOsemMet.values[indNan],product_selected.elev.values[indNan])

# calculando variancias
var_orig_bndo = np.nanvar(observ.nivel.values)
var_mete_bndo = np.nanvar(mare_meteorologica.MeteorTide_BNDO.values)
var_bndo_semM = np.nanvar(bndo_semMet.BNDOsemMet.values)

### ---------- plotting data ---------------- ###

fig,ax = plt.subplots(nrows=5,sharex=True)

# boxes parameters:
props = dict(boxstyle='round',facecolor='gray',alpha=0.5)
xbox = 0.03
ybox = 0.90

# set y-axis Limits
for i in range(ax.shape[0]):
    ax[i].set_ylim(-0.8,1.2)

bndo_filtered.plot(ax=ax[0],title='BNDO - Ilha Guaiba Terminal [6-hourly frequency]')
# inserting box with statistical informations
textstr = r'$\sigma$: %.4f'%(var_orig_bndo)
ax[0].text(xbox,ybox,textstr,transform=ax[0].transAxes,fontsize=14,verticalalignment='top',bbox=props)

mare_meteorologica.plot(ax=ax[1],title=u'Meteorological Tide extracted from BNDO [Hamming, window=30h]')
# inserting box with statistical informations
textstr = r'$\sigma$: %.4f'%(var_mete_bndo)
ax[1].text(xbox,ybox,textstr,transform=ax[1].transAxes,fontsize=14,verticalalignment='top',bbox=props)

product_selected.plot(ax=ax[2],title='ECOM - %s [x BNDO = Skill: %.2f, Corr: %.2f, RMSE: %.2fcm]'%(run,skill,corr,rmse/0.01))
met_ecom.plot(ax=ax[3],title='MetBNDO + %s [x BNDO = Skill: %.2f, Corr: %.2f, RMSE: %.2fcm]'%(run,skill2,corr2,rmse2/0.01))

bndo_semMet.plot(ax=ax[4],title='BNDO - Meteorological Tide. [x %s = Skill: %.2f, Corr: %.2f, RMSE: %.2fcm]'%(run,skill3,corr3,rmse3/0.01))
# inserting box with statistical informations
textstr = r'$\sigma$: %.4f'%(var_bndo_semM)
ax[4].text(xbox,ybox,textstr,transform=ax[4].transAxes,fontsize=14,verticalalignment='top',bbox=props)
