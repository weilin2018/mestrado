'''
Próximo passo:
    . gerar funcao para converter u,v em int,dir. Usar os dados da Laje
    como forma de validar essa funcao

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

##############################################################################
#                          [GEN] FUNCTIONS                                   #
##############################################################################

# read files from laje de santos
def readFiles_Laje(nfiles):
    '''
        read file by file, extracting each column and appending in
        lists.

        args
            nfiles  (list): list with all files

        returns
            a bunch of lists ...
    '''

    datas = []
    inten = []    # intensity
    direc = []    # directions

    for f in nfiles:
        # read file
        data = np.loadtxt(f)

        for i in np.arange(0,data.shape[0],1):
            dt = convert2datetime(data[i,0], data[i,1], data[i,2], data[i,3])

            datas.append(dt)
            inten.append(data[i,7])
            direc.append(data[i,8])

    return np.asarray(datas),np.asarray(inten),np.asarray(direc)

# read files from CFSv2
def readFiles_CFSv2(nfiles):
    """Read and extract data from CFSv2 files.

    Parameters
    ----------
    nfiles : list
        Sorted by name list.

    Returns
    -------
    wu, wv, time: arrays


    """
    # define some variable
    wu   = []
    wv   = []
    time = []

    # read files
    for f in nfiles:
        ncdata = xr.open_dataset(f)
        u     = ncdata['U_GRD_L103'].values
        v     = ncdata['V_GRD_L103'].values
        t     = ncdata['time'].values

        for k in [0,1,2,3]:
            wu.append(np.squeeze(u[k,:,:]))
            wv.append(np.squeeze(v[k,:,:]))
            time.append(np.squeeze(t[k]))

    wu   = np.asarray(np.squeeze(wu))
    wv   = np.asarray(np.squeeze(wv))
    time = np.asarray(time)

    return wu,wv,time


# read hourly product from CFSv2
def readHourly(nfiles):
    '''
    cada arquivo se refere a 1 mes
    cada arquivo, portanto, contem 7*dias_mes*4 dados
    '''

    return 0

# convert columns to datetime
def convert2datetime(year,month,day,hour):
    '''
        receive some year,months,days and hours
        convert to datetime
        return a valid datetime
    '''

    import datetime

    return datetime.datetime(*map(int, [year,month,day,hour]))

# read netcdf files for a given location
def readCDF(nfiles,ilat,ilon):
    '''
        retornar wu e wv num ponto especifico
    '''

    # define array to store data
    wu,wv,instant = [],[],[]

    # loop file by file
    for f in nfiles:
        ncdata = xr.open_dataset(f)
        # data for 4 different hours, so ...
        u    = ncdata['U_GRD_L103'].values[:,indLon,indLat]
        v    = ncdata['V_GRD_L103'].values[:,indLon,indLat]
        time = ncdata['time'].values

        for x,y,t in zip(u,v,time):
            wu.append(x)
            wv.append(y)
            instant.append(t)

    wu   = np.asarray(wu)
    wv   = np.asarray(wv)
    instant = np.asarray(instant)

    return wu,wv,instant

def load_CFSv2_grid():
    # puxar grade CFSv2
    grid = pickle.load(open(BASE_DIR+'masterThesis_analysis/routines/coordenadas_refinada.pickle','r'))
    glon = grid['lon']
    glat = grid['lat']

    return glon,glat

# read netcdf files and extracted data from a given location
def extractData_fromCFSv2(location=[-25.274, -44.93],summer='2014'):
    '''
        Function to extract data from a given location.

        args
            location   [array]: with latitude and longitude you want the data
            summer     [string]: string with the year (2014 or 2015)

        returns
            data       [DataFrame]: dataframe with wu and wv, sorted by
                        DateTimeIndex.

    '''

    BASE_DIR = oceano.make_dir()
    DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao%s/'%(summer))

    nfiles = glob.glob(DATA_DIR+"*.nc")
    nfiles.sort()

    # define lat and lon location
    ilat = location[0]
    ilon = location[1]

    # puxar grade CFSv2
    glon,glat = load_CFSv2_grid()

    # indices: ilat = 26 // ilon = 24
    #lon_ncep,lat_ncep = oceano.find_nearest(glon,glat,ilon,ilat)
    indLat,indLon = 25,25

    # read netcdf files
    nfiles = glob.glob(BASE_DIR.replace('github/', 'ventopcse/data/CFSv2/verao%s/*.nc'%(summer)))
    # extract wu and wv
    wu,wv,time = readCDF(nfiles,indLat,indLon)

    # convert all data to pd.DataFrame. To do that, we need a valid daterange to use as index
    dtRange = pd.date_range(start=time[0],end=time[-1],freq='6H')
    data = pd.DataFrame({'wu':wu,'wv':wv}, index=dtRange)

    return data

# sem uso
def magneticDeclination(lat,lon,date,h=32.8084):
    '''
        Function to calculate a magnetic declination, to use in observed data.
        args
            lat      (float): latitude from some coordinate
            lon      (float): longitude from some coordinate
            date      (list): date to calculate the declination, in this format:
                              (YYYY, MM, DD)
            h        (float): altitude in feet, default = 32.8084 (10m)
        returns
            magDeclination (float): magnetic declination
    '''

    import geomag

    return geomag.declination(dlat=lat,dlon=lon,h=h,time=date)

# corrects magnetic declination caused by the diffenrece between
# magnetic and true North
def correctionMagDec(direction,magDec):
    '''
    This function corrects this magnetic declination, by adding the calculated
    angle to direction data.

    The magnetic declination is the angle of difference between true North
    and magnetic North.

    For instance, if the declination at a certain point were 10° W, then a
    compass at that location pointing north (magnetic) would actually align
    10° W of true North. True North would be 10° E relative to the magnetic
    North direction given by the compass. Declination varies with location
    and slowly changes in time.

    Initially, this function should use the geomag package, but the package
    was not working properly. So, you can access this website:
        https://www.ngdc.noaa.gov/geomag-web/#declination
    and select the coordinate and date you want to know the magnetic declination
    and, finally, insert as an argument to this function (magDec).

    ps: in many cases the magnetic declination will be constant, considering
    a short timeserie. Check !

    args
        lat       (float):
        lon       (float):
        direciton (array):
        magDec    (float):

    returns
        correctedDirection (array):

    '''
    new_direction = [ang+abs(magDec) for ang in direction]

    return new_direction

# function only to calculate wu and wv using direction and intensity
def windComponent(intensity, direction, dates, decliMag, deg=False):
    '''
        apenas para facilitar quando eu quiser comparar a série original com
        a série com correção magnética.
    '''

    # keeping arguments as arrays
    direction = np.asarray(direction)
    intensity = np.asarray(intensity)

    # wu,wv = oceano.spdir2uv(intensity,direction,deg=False)
    # wu,wv = oceano.compass2uv(direction,intensity,kind='meteo2uv')
    # usando rotina de conversão da Paulinha
    wu,wv = oceano.dirmag2uv(direction,intensity,decliMag,'norte')

    # create pandas.dataframe
    i = pd.DatetimeIndex(dates)
    df = pd.DataFrame({'wu':wu, 'wv':wv},
                        index=i)

    # quality control
    df[df['wu'] > 20] = np.nan
    df[df['wu'] <-20] = np.nan

    df[df['wv'] > 20] = np.nan
    df[df['wv'] <-20] = np.nan

    '''
    Como os dados do NCEP são de 6 em 6 horas, a partir da 1h inicial, precisamos
    redimensionar os dados para essa frequencia
    '''
    df = df.resample('6H').mean()

    return df

# replace comma for dots
def replaceComma4dots(df):
    ''' '''

    values = []
    for x in df.values:
        x = str(x).replace(',','.')
        values.append(float(x))

    return np.asarray(values)

##############################################################################
#                          [GEN]PLOTS                                        #
##############################################################################

# plot two databases
def plot(data1,data2,data1Leg,data2Leg,suptitle='Comparando'):
    '''

    '''
    fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(16,8),sharex=True)

    ax[0].plot(data1.index, data1.wu.values,label='WU %s'%(data1Leg))
    ax[0].plot(data2.index, data2.wu.values,label='WU %s'%(data2Leg))

    ax[1].plot(data1.index, data1.wv.values,label='WV %s'%(data1Leg))
    ax[1].plot(data2.index, data2.wu.values,label='WV %s'%(data2Leg))

    ax[0].legend(loc='best')
    ax[1].legend(loc='best')

    plt.suptitle(suptitle, fontsize=24)

    plt.margins(0)          # code to remove blank spaces at the begin and end of graphic

    plt.show()

def plotAlldata(data1,data2,data3,data1Leg,data2Leg,suptitle='Comparando'):
    '''
        3 subplots

        1o: wu
        2o: wv
        3o: intensidade do vento in situ + diagrama direcional

        example:
        plotAlldata(santos['2013-11':'2014-03'],
                    cfsv2['2013-11':'2014-03'],
                    data['2013-11':'2014-03'],
                'PNBOIA','CFSv2','PNBOIA')
    '''

    fig, ax = plt.subplots(nrows=3,ncols=1,figsize=(16,8),sharex=True)

    ax[0].plot(data1.index, data1.wu.values,label='WU %s'%(data1Leg))
    ax[0].plot(data2.index, data2.wu.values,label='WU %s'%(data2Leg))
    ax[0].set_ylabel('WU intensity')

    ax[1].plot(data1.index, data1.wv.values,label='WV %s'%(data1Leg))
    ax[1].plot(data2.index, data2.wu.values,label='WV %s'%(data2Leg))
    ax[1].set_ylabel('WV intensity')

    ax[2].plot(data3.index, data3.Wdir.values,label='Direction PNBOIA/Santos')
    ax[2].set_ylabel('Wind\'s direction')
    ax[2].set_xlabel('Time')
    ax[2].axhline(y=90., color='k', linestyle='-',alpha=0.5)
    ax[2].axhline(y=270., color='k', linestyle='-',alpha=0.5)

    ax[2].axhspan(90, 270, facecolor='0.5', alpha=0.5)

    ax3 = ax[2].twinx()
    ax3.plot(data3.index, data3.Wspd.values,'g--',alpha=0.4,label='Wind\'s Speed PNBOIA/Santos')
    ax3.set_ylabel('Wind\'s Speed')
    ax3.tick_params('g',colors='g')

    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    # ax[2].legend(loc='best')

    plt.suptitle(suptitle, fontsize=24)

    plt.show()

# plot CFSv2 grid and some given location
def plotGrid(point):
    '''
        Functino to plot southeast brazilian coastline, as a background,
        and CFSv2 grid over the map. Also, this function plot some
        given location, to better visualize.

        args
            points      [tuple]: with array containing lat and lon

        return
            plt.show()

        example
            >>> points = ([-24.3194, -46.1803, 'Laje','s','b'], [-25.274, -44.93, 'PNBOIA','d','g'])
            >>> plotGrid(points)
    '''
    fig, ax = plt.subplots(figsize=(15,10))

    m = oceano.make_map(ax)

    # load grid
    glon,glat = load_CFSv2_grid()
    glon,glat = np.meshgrid(glon,glat)

    # plot grid and other features
    # m.plot(glon,glat,'k--',alpha=.4,latlon=True)
    # m.plot(glon.T,glat.T,'k--',alpha=.4,latlon=True)

    for p in point:
        m.scatter(p[1], p[0], c=p[4], marker=p[3], s=50, label=p[2], latlon=True)

    plt.legend(loc='best')
    plt.title('CFSv2 grid \n Observed data location and nearest grid points')
    plt.show()

def plotOriginalCorrected_against_Reanalysis(original,corrected,reanalysis,title):
    fig, ax = plt.subplots(nrows=2,ncols=1)

    ax[0].plot(original.index, original.wu.values, label='wu original')
    ax[0].plot(corrected.index, corrected.wu.values, label='wu corrected')
    ax[0].plot(reanalysis.index, reanalysis.wu.values, label='wu cfsv2')

    ax[0].margins(0)
    ax[0].legend(loc='lower right')

    ax[1].plot(original.index, original.wv.values, label='wv original')
    ax[1].plot(corrected.index, corrected.wv.values, label='wv corrected')
    ax[1].plot(reanalysis.index, reanalysis.wv.values, label='wv cfsv2')
    ax[1].margins(0)
    ax[1].legend(loc='lower right')

    plt.suptitle(title, fontsize=24)

    plt.show()

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

##############################################################################
#                          LAJE DE SANTOS                                    #
##############################################################################

def convertCFS_2_meteoConv(wu,wv):
    """Converting CFSv2 wind to meteorological convention.

    This function also convert the trigonometric angle to cardinal coordinate.

    Parameters
    ----------
    wu : array_like
        Array containing all values of u.
    wv : array_like
        Array containing all values of v.

    Returns
    -------
    wind_absolut : array_like
        Description of returned object.
    wind_direction_trig_from_degrees : array_like
        as
    wind_direction_cardinal : array_like
        as
    """

    from math import atan2

    # compute the wind's intensity
    wind_absolut = np.sqrt(wu**2 + wv**2)
    # normalizing vectors
    wun,wvn = [],[]
    for x,y,mean in zip(wu,wv,wind_absolut):
        wun.append(x/mean)
        wvn.append(y/mean)
    wun = np.asarray(wun)
    wnv = np.asarray(wvn)

    # compute the wind's direction
    wind_direction_trig_to = atan2(wvn,wun)
    # convert wind's direction tro degrees
    wind_direction_trig_to_degrees = (wind_direction_trig_to * 180)/np.pi

    # finally, converting wind's direction to meteorological convention
    # (wind is coming from)
    wind_direction_trig_from_degrees = wind_direction_trig_to_degrees + 180

    # Converting the trigonometrical angle to cardinal coordinantes
    wind_direction_cardinal = 90 - wind_direction_trig_from_degrees

    return wind_absolut, wind_direction_trig_from_degrees, wind_direction_cardinal

def lajeVcfsv2(DATA_DIR, NCEP_DIR):
    '''
        Function to read and plot data from meteorological station (laje de santos)
        against CFSv2 reanalysis.

        args
            DATA_DIR (string): full path to meteorological data
            NCEP_DIR (string): full path to CFSv2 data
    '''
    # select files
    nfiles = glob.glob(DATA_DIR + '*.txt')
    nfiles.sort()

    # extract information from files
    dates,intensity,direction = readFiles_Laje(nfiles)
    # correction of magnetic declination
    # ndirection = correctionMagDec(direction,-21.03)
    # here we need tranform to oceanographic convention, by adding 180 degrees
    # ndirection = [ang+180. for ang in ndirection]
    # ndirection = np.asarray(ndirection)
    # ndirection *= np.pi/180

    # convert intensity and direction to components (u,v)
    decliMag = -21.03
    laje_corrected = windComponent(intensity, direction, dates, decliMag)
    # laje_original  = windComponent(intensity, direction, dates, deg=True)

    # LOAD CFSv2

    # select files
    nfiles = glob.glob(NCEP_DIR + '*.nc')
    nfiles.sort()

    wu,wv,time = readFiles_CFSv2(nfiles)

    # convert time to pd.DateTimeIndex
    i = pd.DatetimeIndex(time)

    # create dataframe
    cfsv2 = pd.DataFrame({'wu':wu,'wv':wv}, index=i)

    # rotate both dataset to along and cross shore components
    angRot = (-18.*np.pi)/180
    wu_along, wv_across = oceano.rotaciona(laje_corrected.wu.values, laje_corrected.wv.values,angRot)
    laje_rotated = pd.DataFrame({'wu':wu_along,'wv':wv_across}, index=laje_corrected.index)
    wu_along, wv_across = oceano.rotaciona(cfsv2.wu.values, cfsv2.wv.values,angRot)
    cfsv2_rotated = pd.DataFrame({'wu':wu_along,'wv':wv_across}, index=cfsv2.index)

    return laje_corrected, cfsv2,laje_rotated, cfsv2_rotated

def laje_statisticalAnaysis(laje,cfsv2,whichSerie):
    ''' '''
    # statistical analysis for Laje de Santos
    laje_cut = laje.copy()
    cfsv_cut = cfsv2.copy()

    laje_cut.interpolate(inplace=True)
    cfsv_cut.interpolate(inplace=True)

    skillWu = oceano.skill_willmott(laje_cut.wu.values, cfsv_cut.wu.values)
    skillWv = oceano.skill_willmott(laje_cut.wv.values, cfsv_cut.wv.values)

    corrWu = calculateCorr(laje_cut.wu.values, cfsv_cut.wu.values)[0]
    corrWv = calculateCorr(laje_cut.wv.values, cfsv_cut.wv.values)[0]

    mseWu  = calculateMSE(laje_cut.wu.values, cfsv_cut.wu.values)
    mseWv  = calculateMSE(laje_cut.wv.values, cfsv_cut.wv.values)

    # plot data and skill
    fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True)

    ax[0].plot(laje.wu,label='Laje')
    ax[0].plot(cfsv2.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Along shore')

    wuText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text('2015-04-07', 7.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend(loc='lower left')

    ax[1].plot(laje.wv,label='Laje')
    ax[1].plot(cfsv2.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Cross shore')

    wvText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text('2015-04-07', 7.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend(loc='lower left')

    ax[2] = stickplot(laje*(-1),ax[2])
    ax[2].set_title('Laje')
    ax[3] = stickplot(cfsv2*(-1),ax[3])
    ax[3].set_title('CFSv2')

    plt.suptitle('Laje de Santos [%s] v CFSv2\n[2015-01-03 to 2015-01-25]' % (whichSerie),fontsize=26)
    plt.show()

##############################################################################
#                                PNBOIA                                      #
##############################################################################
def extractColumnsIndexByName(columns_names,filepath):
    '''
        Function to locate the number of columns by name, returning those indexes
        to read_csv extract only that we want.

        args
            columns_names    (array): array with the names
            filepath         (string): full path to the file we want to load

        returns
            indexes          (list): list with indexes related to columns_names
    '''

    # load all file
    data = pd.read_csv(filepath)
    # extract columns names to a array
    all_columns_names = data.columns.values

    # locate the position of each column passed as argument
    # [i for i,x in enumerate(testlist) if x == 1]
    indexes = []
    for column in columns_names:
        values = [i for i,x in enumerate(all_columns_names) if x == column]
        indexes.append(values[0])

    return indexes

def pnboiaVcfsv2(DATA_DIR, NCEP_DIR,region,magDec):
    """
    This function reads the data from moored buoy (PNBOIA), correcting them
    to magnetic declination.

    This function reads the .csv file in DATA_DIR (for Santos or Itajai buoys),
    convert all these data to a pandas.DataFrame. After, the function also
    corrects the magnetic declination, using some auxiliar functions.

    The functino also reads netCDFs files in NCEP_DIR, downloaded from
    CFSv2 dataset.

    Finally, is returned 3 pandas.DataFrame.

    Parameters
    ----------
    DATA_DIR : string
        Full path to pnboia files.
    NCEP_DIR : string
        Full ath to CFSv2 files.
    region : string
        Which region we want to import data (Santos or Itajai).
    magDec : float
        The magnetic declination, calculated using NASA's website.

    Returns
    -------
    pnboia_corrected : pandas.Dataframe
        Corrected data from moored buoy.
    pnboia_corrected : pandas.Dataframe
        Original data from moored buoy.
    pnboia_corrected : pandas.Dataframe
        Reanalysis data from CFSv2.
    """
    if region == 'Santos':
        nfiles = DATA_DIR+'santos.csv'

    if region == 'Itajai':
        nfiles = DATA_DIR+'itajai.csv'

    # columns positions
    cols_name = ['Wspd', 'Wdir','Wspdflag']
    columns_data = extractColumnsIndexByName(cols_name,nfiles)
    cols_date = ['Year', 'Month', 'Day', 'Hour', 'Minute']
    columns_date = extractColumnsIndexByName(cols_date,nfiles)

    data = pd.read_csv(nfiles,usecols=columns_date)

    dates = []

    for index,row in data.iterrows():
        # dt = convert2datetime(data[i,0], data[i,1], data[i,2], data[i,3])
        values = [row[0], row[1], row[2], row[3], row[4]]
        dt = datetime.datetime(*map(int, values))
        dates.append(dt)

    # convert dates to pd.DateTimeIndex
    i = pd.DatetimeIndex(dates)

    # read values, using dates as index
    data = pd.read_csv(nfiles,usecols=columns_data)
    data['datetime'] = i
    data.set_index('datetime',inplace=True)

    # basic quality control
    # using the flag to insert NaN values
    data[data['Wspdflag'] == 4] = np.nan

    # now we can remove the flag column
    data.drop(columns=['Wspdflag'], inplace=True)

    # replace comma by dots
    Wspd = replaceComma4dots(data['Wspd'])
    Wdir = replaceComma4dots(data['Wdir'])

    ndirection = correctionMagDec(Wdir,magDec)

    # converting from meteorological to oceanographic convention
    # ndirection = np.asarray(ndirection) + 180

    # create new dataframe with wu and wv values
    # pnboia_original  = windComponent(Wspd,Wdir,dates,deg=True)
    pnboia_corrected = windComponent(Wspd,ndirection,dates,magDec,deg=True)

    # rotating componentes to along and cross shore, only for corrected data
    # with an angle of 60degrees (pi/3)
    wu_along, wv_cross = oceano.rotaciona(pnboia_corrected.wu.values, pnboia_corrected.wv.values,np.pi/3)
    # create new dataframe, with data rotated
    pnboia_rotated = pd.DataFrame({'wu':wu_along, 'wv':wv_cross},index=pnboia_corrected.index)

    # LOAD CFSv2

    # select files
    # NCEP_DIR = NCEP_DIR.replace('LajeDeSantos/', 'pnboia/')
    nfiles = glob.glob(NCEP_DIR + '*.nc')
    nfiles.sort()

    wu,wv,time = readFiles_CFSv2(nfiles)

    # convert time to pd.DateTimeIndex
    i = pd.DatetimeIndex(time)

    # create dataframe
    wu_along, wv_cross = oceano.rotaciona(wu, wv,np.pi/3)
    cfsv2 = pd.DataFrame({'wu':wu_along,'wv':wv_cross}, index=i)

    return pnboia_corrected, pnboia_rotated, cfsv2

def statisticalAnalysis(observed,reanalysis,region,whichSerie,date='01-21'):
    """calculate statistical coefficients and plot a graphic with those values.

    This function analyze how reanalysis data represent the real information
    obtained in some point (observed or in situ data). Calculating the
    skill parameter, pearson's coefficient and mean squared error and, finally,
    plotting all data in 3 subplots:

        first:
            u component (observed and reanalysis) with statistical
            parameters in a box
        second:
            v component, in the same way as firts graphic
        third:
            a stickplot of observed data

    Parameters
    ----------
    observed : pd.DataFrame
        observed dataframe containing wu and wv from in situ observation.
    reanalysis : pd.DataFrame
        analysis dataframe containing wu and wv from reanalysis.
    region : string
        could be Santos, Itajai or Laje de Santos.
    whichSerie : string
        corrected or original observed data, where corrected is the data
        with magnetic declination correction.
    date  : string
        Which month-day we want to plot the box with statistical values. Default
        is 01-21.

    Returns
    -------
    Show a plot.

    """
    # short variale names
    obse_cut = observed
    cfsv_cut = reanalysis
    # calculate some statistical values
    skillWu = oceano.skill_willmott(obse_cut.wu.values, cfsv_cut.wu.values)
    skillWv = oceano.skill_willmott(obse_cut.wv.values, cfsv_cut.wv.values)

    corrWu = calculateCorr(obse_cut.wu.values, cfsv_cut.wu.values)[0]
    corrWv = calculateCorr(obse_cut.wv.values, cfsv_cut.wv.values)[0]

    mseWu  = calculateMSE(obse_cut.wu.values, cfsv_cut.wu.values)
    mseWv  = calculateMSE(obse_cut.wv.values, cfsv_cut.wv.values)

    # defining some values, based on data, to plot
    year      = observed.index[0].strftime('%Y')
    xPosition = '%s-%s' % (str(year),date)

    firstDate = observed.index[0].strftime('%Y-%m-%d')
    lastDate  = observed.index[-1].strftime('%Y-%m-%d')
    period    = '%s to %s' % (firstDate, lastDate)

    # plotting all information
    fig, ax = plt.subplots(nrows=4,ncols=1,sharex=True)

    ax[0].plot(obse_cut.wu,label='PNBOIA')
    ax[0].plot(cfsv_cut.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-15,15)
    ax[0].set_title('Along shore')

    wuText = r'Skill: %0.2f - Corr.: %0.2f - MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text(xPosition, -12.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend(loc='best')

    ax[1].plot(obse_cut.wv,label='PNBOIA')
    ax[1].plot(cfsv_cut.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-15,15)
    ax[1].set_title('Cross shore')

    wvText = r'Skill: %0.2f | Corr.: %0.2f | MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text(xPosition, -12.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend(loc='best')

    ax[2] = stickplot(obse_cut,ax[2])
    ax[2].set_title('PNBOIA')

    ax[3] = stickplot(cfsv_cut,ax[3])
    ax[3].set_title('CFSv2')

    plt.suptitle('PNBOIA/%s  v CFSv2 \n [%s]' % (whichSerie, period),fontsize=26)
    plt.show()

##############################################################################
#                                STATS                                       #
##############################################################################
def calculateMSE(x,y):
    """Calculate mean squared error.

     Calculate mean squared error using the formula:

            MSE = 1/n * sum((x-y)^2)

    Parameters
    ----------
    x : array
        real data as an array
    y : array
        modeled data as an array

    Returns
    -------
    mse : float
        the mean squared error calculated

    """

    return np.mean((x - y)**2)

def calculateCorr(x,y):
    """Calculate correlation between two series with same size.

    Using the Pearson's method, from scipy.stats.pearsonr, this function
    calculate the correlation between x and y, two array with same size and
    without missing values.

    Parameters
    ----------
    x : array
        real data.
    y : array
        modeled data.

    Returns
    -------
    float
        correlation coefficient.

    """
    import scipy.stats as stats

    return stats.pearsonr(x,y)

##############################################################################
#                          PLOTAR PONTOS                                     #
##############################################################################
points = (
    [-24.3194, -46.1803, 'Laje de Santos','s','b'],
    [-25.274, -44.93, 'PNBOIA/Santos',',','g'],
    [-28.489, -47.52, 'PNBOIA/Itajai',',','y'],
    # definir os pontos de grade mais próximo aos pontos de dados observados
    [-28.51790047, -47.45500183, '', 'o','k'],
    [ -25.24699974,-45, '','o','k'],
    [-24.2249, -46.226999999999975, '', 'o','k']
)

##############################################################################
#                               MAIN CODE                                    #
##############################################################################

os.system('clear')

# define some constants
BASE_DIR = oceano.make_dir()
NCEP_DIR = BASE_DIR.replace('github/', 'ventopcse/data/serie_cfsv2/LajeDeSantos/2015/')
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Est_lajeSantos/2015/atualizado/')
ARGO_DIR = BASE_DIR.replace('github/', 'ventopcse/data/')

os.system('clear')
print('Plotting map with locations')
plotGrid(points)

#-----------------------------------------------------------------------------#
# -----------------------------LAJE DE SANTOS---------------------------------#
#-----------------------------------------------------------------------------#
'''
    Para comparar os dados da Laje de Santos, selecionamos o período de
    Janeiro de 2015 (de 03 a 25), seguindo o seguinte critério:
        . período mais longo, dentro do verão de 2015, sem cortes por falta de
          dados
        . estar compreendido dentro de um período de episódio anomalo

    O primeiro gráfico mostra as componentes zonais e meridionais do vento
    antes e depois de corrigirmos a direção do vento com a declinação magnética.
    O importante desde gráfico é observar que, mesmo com baixa correção,
    a influência do ângulo na direção do vento acarreta em mudanças tanto na fase
    quanto na intensidade das componentes.

    O segundo gráfico apresenta os dados observados na estação, originais e
    corrigidos, contra os dados, para o mesmo período, das saídas da reanalises
    do CFSv2. Nele notamos que a componente zonal é bem representada pelo modelo,
    mas a componente meridional não.

    Já o terceiro gráfico, apresenta uma comparação entre o dado corrigido
    contra o dado de reanalise.

    Por fim, plotamos, no quarto gráfico, comparações estatísticas entre os
    dados observados originais e os dados de reanalise e, no quinto gráfico,
    o mesmo, mas para os dados corrigidos e os dados de reanalise.

    Importante: os dados brutos (que são carregados no programa) não foram
    corrigidos com a declinação magnética e, portanto, isso é feito no programa.
    E eles estão na convenção meteorológica.
'''
os.system('clear')
print('Plotting data from Laje de Santos with the first point from NCEP')
laje_corrected, cfsv2_laje, laje_rotated, cfsv2_rotated = lajeVcfsv2(DATA_DIR,NCEP_DIR)

# mesma coisa para os dados rotacionados
lajeCut = laje_rotated['2015-03-20':]
cfsvCut = cfsv2_rotated['2015-03-20':]
cfsvCut['wu'] *= -1
# interpolating
# lajeCut.interpolate(inplace=True)
# cfsvCut.interpolate(inplace=True)

laje_statisticalAnaysis(lajeCut,cfsvCut,'rotated/-18deg')

#-----------------------------------------------------------------------------#
# -----------------------------LAJE DE SANTOS - 2 ----------------------------#
#-----------------------------------------------------------------------------#
# os.system('clear')
print('Plotting data from Laje de Santos with the second point from NCEP')
laje_corrected, cfsv2_laje, laje_rotated, cfsv2_rotated = lajeVcfsv2(DATA_DIR,NCEP_DIR.replace('2015','laje2'))

lajeCut = laje_rotated['2015-03-20':]
cfsvCut = cfsv2_rotated['2015-03-20':]*(-1)
laje_statisticalAnaysis(lajeCut,cfsvCut,'rotated/60deg')

##################################################
# testando conversao de componentes para intensidade e direcao
##################################################

def con(wu,wv):
    import math

    ws = np.sqrt(wu**2 + wv**2)
    wd = []
    for i,j in zip(wu,wv):
        wd.append((180/np.pi)*math.atan2(-i,-j))

    wd = np.asarray(wd)

    return ws,wd

wu,wv = lajeCut.wu.values, lajeCut.wv.values
ws,wd = con(wu,wv)

polarPlot(ws,wd,'convertido')

polarPlot(intensity,direction,'original')
#-----------------------------------------------------------------------------#
# -----------------------------PNBOIA/SANTOS----------------------------------#
#-----------------------------------------------------------------------------#
'''
    Para analisar os dados da bóia do Programa Nacional de Boias (PNBOIA), procuramos
    visualizar, inicialmente, todo o período de dados obtidos. Desta forma,
    foi possível verificar qual o melhor período para realizar a análise e comparação
    dos dados, tomando como critério:
        . período de anomalia (verões de 2014 e 2015)
        . período sem anomalia

    Após análise inicial, selecionou-se os seguintes períodos para continuar:
        1o) Jan/2012 - Mar/2012
        2o) Jan/2014 - Mar/2014
        3o) Jan/2015 - Mar/2015

    Uma coisa importante observada é que, quando comparando os dados observados
    com os dados de reanálise, havia uma possível divergência de sistema de convenções.
    Desta forma, multiplicou-se toda a série de dados observados por -1, obtendo
    comportamentos mais condizentes. No entanto, para o ano de 2014, essa mudança
    de convenção não se mostrou satisfatória. Desta forma, esta série permaneceu
    sem alteração. Os valores estatísticos obtidos, após todo este tratamento de dados,
    são apresentados nas imagens correspondentes de cada ano.
'''
os.system('clear')
print('Plotting data from PNBOIA/Santos')
santos_corrected, santos_rotated, cfsv2_santos = pnboiaVcfsv2(DATA_DIR=ARGO_DIR, NCEP_DIR=NCEP_DIR.replace('LajeDeSantos/2015/', 'pnboia_santos/'),region='Santos',magDec=-21.38)

# converting values from oceanographic convention to meteorological convention
cfsv2_santos *= -1

#### statistical analysis

def stats_and_plot(santos, cfsv2):
    # slicing the periods we want
    san2012 = santos['2011-11-01 06:00:00':'2012-03']
    cfs2012 = cfsv2['2011-11-01 06:00:00':'2012-03']

    san2013 = santos['2012-11-01 06:00:00':'2013-03']
    cfs2013 = cfsv2['2012-11-01 06:00:00':'2013-03']

    san2014 = santos['2013-11':'2014-03']
    cfs2014 = cfsv2['2013-11':'2014-03']

    san2015 = santos['2014-11-17':'2015-03']
    cfs2015 = cfsv2['2014-11-17':'2015-03']

    san2016 = santos['2015-11':'2016-03']
    cfs2016 = cfsv2['2015-11':'2016-03']

    # performin an interpolation to remove nan values
    san2012.interpolate(inplace=True)
    san2013.interpolate(inplace=True)
    san2014.interpolate(inplace=True)
    san2015.interpolate(inplace=True)
    san2016.interpolate(inplace=True)

    # calculating statistical parameters and plotting
    statisticalAnalysis(san2012, cfs2012,'Rotated','Santos','12-01')
    statisticalAnalysis(san2013, cfs2013,'Rotated','Santos','12-01')
    statisticalAnalysis(san2014, cfs2014,'Rotated','Santos','12-01')
    statisticalAnalysis(san2015, cfs2015,'Rotated','Santos','12-01')
    statisticalAnalysis(san2016, cfs2016,'Rotated','Santos','12-01')

    # correlation moving window
    s2012 = pd.rolling_corr(san2012.wu.values,cfs2012.wu.values,window=10,pairwise=True)
    s2013 = []
    s2014 = []
    s2015 = []
    s2016 = []



# clc()
stats_and_plot(santos_rotated,cfsv2_santos)

#### plotting windroses

def plotsPolar(santos,cfsv2,year):
    if year == 2012:
        san_cut = santos['2011-12-01 06:00:00':'2012-03']
        cfs_cut = cfsv2['2011-12-01 06:00:00':'2012-03']

    if year == 2013:
        san_cut = santos['2012-11-01 06:00:00':'2013-03']
        cfs_cut = cfsv2['2012-11-01 06:00:00':'2013-03']

    if year == 2014:
        san_cut = santos['2013-12':'2014-03']
        cfs_cut = cfsv2['2013-12':'2014-03']*(-1)

    if year == 2015:
        san_cut = santos['2014-12':'2015-03']
        cfs_cut = cfsv2['2014-12':'2015-03']

    if year == 2016:
        san_cut = santos['2015-12':'2016-03']
        cfs_cut = cfsv2['2015-12':'2016-03']

    # performin an interpolation to remove nan values
    san_cut.interpolate(inplace=True)

    ws,wd = oceano.uv2intdir(san_cut.wu.values,san_cut.wv.values)
    polarPlot(ws,wd,'PNBOIA - %s'%(str(year)))
    ws,wd = oceano.uv2intdir(cfs_cut.wu.values,cfs_cut.wv.values)
    polarPlot(ws,wd,'CFSv2 - %s'%(str(year)))

    statisticalAnalysis(san_cut, cfs_cut,'Corrected','Santos','12-01')

for y in [2012,2013,2014,2015,2016]:
    plotsPolar(santos_corrected,cfsv2_santos,y)


#-----------------------------------------------------------------------------#
# -----------------------------LAJE x PNBOIA----------------------------------#
#-----------------------------------------------------------------------------#

def compareLajePNBOIA(laje,pnboia):
    fig, ax = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax[0].set_title('Along shore')
    ax[0].plot(laje.index,laje.wu.values,label='Laje')
    ax[0].plot(pnboia.index,pnboia.wu.values,label='PNBOIA')
    ax[0].legend()

    ax[1].set_title('Cross shore')
    ax[1].plot(laje.index,laje.wv.values,label='Laje')
    ax[1].plot(pnboia.index,pnboia.wv.values,label='PNBOIA')
    ax[1].legend()

    plt.show()

compareLajePNBOIA(lajeCut,santos_corrected['2015-03-20 00:00:00':'2015-12-31 18:00:00'])
