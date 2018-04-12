'''
Próximo passo:
    . arrumar a série de dados com nan values para poder ver o valor de
      correlação entre eles.

      line: 830

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
#                               FUNCTIONS                                    #
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
    '''
        Read file by file, extracting u,v and time from netcdf files.

        args
            nfiles    (list): list with all files

        returns
            wu,wv and time (arrays)
    '''
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
            wu.append(u[k,:,:])
            wv.append(v[k,:,:])
            time.append(t[k])

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
    m.plot(glon,glat,'k--',alpha=.4,latlon=True)
    m.plot(glon.T,glat.T,'k--',alpha=.4,latlon=True)

    for p in point:
        m.scatter(p[1], p[0], c=p[4], marker=p[3], s=50, label=p[2], latlon=True)

    plt.legend(loc='best')
    plt.title('Observed data location and nearest grid points')
    plt.show()

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
def windComponent(intensity, direction, dates, deg=False):
    '''
        apenas para facilitar quando eu quiser comparar a série original com
        a série com correção magnética.
    '''

    wu,wv = oceano.spdir2uv(intensity,direction,deg=True)

    # create pandas.dataframe
    i = pd.DatetimeIndex(dates)
    df = pd.DataFrame({'wu':wu, 'wv':wv},
                        index=i)

    # quality control
    df[df['wu'] > 10] = np.nan
    df[df['wu'] <-10] = np.nan

    df[df['wv'] > 10] = np.nan
    df[df['wv'] <-10] = np.nan

    '''
    Como os dados do NCEP são de 6 em 6 horas, a partir da 1h inicial, precisamos
    redimensionar os dados para essa frequencia
    '''
    df = df.resample('6H').mean()

    return df

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

# replace comma for dots
def replaceComma4dots(df):
    ''' '''

    values = []
    for x in df.values:
        x = str(x).replace(',','.')
        values.append(float(x))

    return np.asarray(values)

##############################################################################
#                          LAJE DE SANTOS                                    #
##############################################################################

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

    ndirection = correctionMagDec(direction,-21.02)

    # convert intensity and direction to components (u,v)
    laje_corrected = windComponent(intensity, ndirection, dates, deg=True)
    laje_original  = windComponent(intensity, direction, dates, deg=True)

    # LOAD CFSv2

    # select files
    nfiles = glob.glob(NCEP_DIR + '*.nc')
    nfiles.sort()

    wu,wv,time = readFiles_CFSv2(nfiles)

    # convert time to pd.DateTimeIndex
    i = pd.DatetimeIndex(time)

    # create dataframe
    cfsv2 = pd.DataFrame({'wu':wu,'wv':wv}, index=i)

    #  PLOT LAJEvCFSv2
    # plot(lajesantos,cfsv2,'Laje','CFSv2')

    return laje_corrected, laje_original, cfsv2

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
    fig, ax = plt.subplots(nrows=2,ncols=1)

    ax[0].plot(laje_cut.wu,label='Laje')
    ax[0].plot(cfsv_cut.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-10,10)

    wuText = r'Skill: %0.2f - Corr.: %0.2f - MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text('2015-01-07', 7.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend()

    ax[1].plot(laje_cut.wv,label='Laje')
    ax[1].plot(cfsv_cut.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-10,10)

    wvText = r'Skill: %0.2f - Corr.: %0.2f - MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text('2015-01-07', 7.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend()

    plt.suptitle('Laje de Santos [%s] v CFSv2 [2015-01-03 to 2015-01-25]' % (whichSerie),fontsize=26)
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
    '''
        Function to read and plot data from moored buoy (Santos)
        against CFSv2 reanalysis.
    args
        DATA_DIR  [string] = full path to pnboia file
        NCEP_DIR  [string] = full path to cfsv2 data
        region    [string] = which region to plot (santos or itajai)
        magDec     [float] = magnetic declination

    returns
        nothing!
    '''

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

    # create new dataframe with wu and wv values
    pnboia_original  = windComponent(Wspd,Wdir,dates,deg=True)
    pnboia_corrected = windComponent(Wspd,ndirection,dates,deg=True)

    # LOAD CFSv2

    # select files
    # NCEP_DIR = NCEP_DIR.replace('LajeDeSantos/', 'pnboia/')
    nfiles = glob.glob(NCEP_DIR + '*.nc')
    nfiles.sort()

    wu,wv,time = readFiles_CFSv2(nfiles)

    # convert time to pd.DateTimeIndex
    i = pd.DatetimeIndex(time)

    # create dataframe
    cfsv2 = pd.DataFrame({'wu':wu,'wv':wv}, index=i)

    #  PLOT LAJEvCFSv2
    # if region == 'Santos':
    #
    #     #  PLOT SANTOSvCFSv2
    #     suptitlte = 'PNBOIA/%s x CFSv2 [11/2013 - 03/2014]'%(region)
    #     plot(pnboia_original['2013-11':'2014-03'],cfsv2['2013-11':'2014-03'],'PNBOIA/%s'%(region),'CFSv2',suptitle=suptitlte)
    #
    #     suptitlte = 'PNBOIA/%s x CFSv2 [11/2014 - 03/2015]'%(region)
    #     plot(pnboia_original['2014-11':'2015-03'],cfsv2['2014-11':'2015-03'],'PNBOIA/%s'%(region),'CFSv2',suptitle=suptitlte)
    #
    # if region == 'Itajai':
    #     suptitlte = 'PNBOIA/%s x CFSv2 [11/2014 - 03/2015]'%(region)
    #     # skill     = oceano.skill_willmott(pnboia['2015-01-01 02:00':'2015-03'],cfsv2['2015-01-01 02:00':'2015-03'])
    #     plot(pnboia_original['2015-01-01 02:00':'2015-03'],cfsv2['2015-01-01 02:00':'2015-03'],'PNBOIA/%s'%(region),'CFSv2',suptitle=suptitlte)
    #     # plotAlldata(santos['2014-11':'2015-03'],cfsv2['2014-11':'2015-03'], data['2014-11':'2015-03'], 'PNBOIA','CFSv2','PNBOIA/Santos x CFSv2 [11/2014 - 03/2015]')

    return pnboia_corrected, pnboia_original, cfsv2

def statisticalAnalysis(observed,reanalysis,region,whichSerie):
    '''
        observed     (pd.DataFrame): observed dataframe containing wu and wv
                        from in situ observation
        analysis     (pd.DataFrame): analysis dataframe containing wu and wv
                        from CFSv2
        whichSerie   (string): corrected or original
        region       (region): could be Santos, Itajai or Laje de Santos
    '''
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
    xPosition = '%s-01-21' % (str(year))

    firstDate = observed.index[0].strftime('%Y-%m-%d')
    lastDate  = observed.index[-1].strftime('%Y-%m-%d')
    period    = '%s to %s' % (firstDate, lastDate)

    # plotting all information
    fig, ax = plt.subplots(nrows=2,ncols=1)

    ax[0].plot(obse_cut.wu,label='PNBOIA')
    ax[0].plot(cfsv_cut.wu,label='CFSv2')
    ax[0].margins(0)
    ax[0].set_ylim(-10,10)

    wuText = r'Skill: %0.2f - Corr.: %0.2f - MSE: %0.2f' % (skillWu,corrWu,mseWu)
    ax[0].text(xPosition, 7.5, wuText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[0].legend(loc='upper center')

    ax[1].plot(obse_cut.wv,label='PNBOIA')
    ax[1].plot(cfsv_cut.wv,label='CFSv2')
    ax[1].margins(0)
    ax[1].set_ylim(-10,10)

    wvText = r'Skill: %0.2f - Corr.: %0.2f - MSE: %0.2f' % (skillWv,corrWv,mseWv)
    ax[1].text(xPosition, 7.5, wvText, ha='center',va='center',bbox=dict(boxstyle='round', ec=(1.,0.5,0.5), fc=(1.,0.8,0.8)))
    ax[1].legend(loc='upper center')

    plt.suptitle('PNBOIA/%s  v CFSv2 \n [%s]' % (whichSerie, period),fontsize=26)
    plt.show()

##############################################################################
#                                STATS                                       #
##############################################################################
def calculateMSE(x,y):
    '''
        Calculate mean square error using the formula:

        args
            x   [array]: real data
            y   [array]: modeled data

        return
            mse [float];
    '''

    return np.mean((x - y)**2)

def calculateCorr(x,y):
    '''
        Calculate correlation between two timeseries

        args
            x   [array]: real data
            y   [array]: modeled data

        return
            corr[float];
    '''

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
NCEP_DIR = BASE_DIR.replace('github/', 'ventopcse/data/serie_cfsv2/LajeDeSantos/')
DATA_DIR = BASE_DIR.replace('github/', 'ventopcse/data/Est_lajeSantos/2015/')
ARGO_DIR = BASE_DIR.replace('github/', 'ventopcse/data/')

os.system('clear')
print('Plotting map with locations')
plotGrid(points)

#--------------------------------------------------------------------------------
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
'''
os.system('clear')
print('Plotting data from Laje de Santos')
laje_corrected, laje_original ,cfsv2_laje = lajeVcfsv2(DATA_DIR,NCEP_DIR)

# if you wanna to check the differences between original and corrected, just ...
plotTimeseries(laje_original['2015-01-03':'2015-01-25'],
               laje_corrected['2015-01-03':'2015-01-25'],
               'Comparation between original and corrected timeserie \n after magnetic declination correction - Laje de Santos Weather Station')

plotOriginalCorrected_against_Reanalysis(laje_original['2015-01-03':'2015-01-25'],
                                         laje_corrected['2015-01-03':'2015-01-25'],
                                         cfsv2_laje['2015-01-03':'2015-01-25'],
                                         'Comparacao entre dados observados [originais e corrigidos] com dados de reanalisis\nJaneiro de 2015 [Laje de Santos]')

# visualize data against CFSv2
plot(laje_original['2015-01-03':'2015-01-25'],cfsv2_laje['2015-01-03':'2015-01-25'],'Laje','CFSv2')
plot(laje_corrected,cfsv2_laje,'Laje','CFSv2')

# statistical analysis
laje_statisticalAnaysis(laje_original['2015-01-03':'2015-01-25'], cfsv2_laje['2015-01-03':'2015-01-25'],'Original')
laje_statisticalAnaysis(laje_corrected['2015-01-03':'2015-01-25'], cfsv2_laje['2015-01-03':'2015-01-25'],'Corrected')

#--------------------------------------------------------------------------------
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
santos_corrected, santos_original, cfsv2_santos = pnboiaVcfsv2(DATA_DIR=ARGO_DIR, NCEP_DIR=NCEP_DIR.replace('LajeDeSantos/', 'pnboia_santos/'),region='Santos',magDec=-21.38)

# converting values to meteorological convention
santos_original *= -1
santos_corrected *= -1

# if you wanna to check the differences between original and corrected, just ...
plotTimeseries(santos_original['2015-01':'2015-01'],
               santos_corrected['2015-01':'2015-01'],
               'Comparation between original and corrected timeserie \n after magnetic declination correction - PNBOIA/Santos')

# statistical analysis
statisticalAnalysis(santos_corrected['2012-01':'2012-03'], cfsv2_santos['2012-01':'2012-03'],'Corrected','Santos')
statisticalAnalysis(santos_corrected['2014-01':'2014-03']*(-1), cfsv2_santos['2014-01':'2014-03'],'Corrected','Santos')
statisticalAnalysis(santos_corrected['2015-01':'2015-03'], cfsv2_santos['2015-01':'2015-03'],'Corrected','Santos')

#--------------------------------------------------------------------------------
print('Plotting data from PNBOIA/Itajai')
itajai_corrected, itajai_original, cfsv2_itajai = pnboiaVcfsv2(DATA_DIR=ARGO_DIR, NCEP_DIR=NCEP_DIR.replace('LajeDeSantos/', 'pnboia_itajai/'),region='Itajai',magDec=-19.22)


# if you wanna to check the differences between original and corrected, just ...
plotTimeseries(itajai_original['2015-01-03':'2015-01-25'],
               itajai_corrected['2015-01-03':'2015-01-25'],
               'Comparation between original and corrected timeserie \n after magnetic declination correction - Laje de Santos Weather Station')

# statistical analysis
pnboia_statisticalAnaysis(itajai_original['2015-01':'2015-03'], cfsv2_itajai['2015-01':'2015-03'],'Original','Itajai')
pnboia_statisticalAnaysis(itajai_corrected['2015-01':'2015-03'], cfsv2_itajai['2015-01':'2015-03'],'Corrected','Itajai')

suptitlte = 'PNBOIA/%s x CFSv2 [2015/January - 2015/March]'%('Itajai')
# skill     = oceano.skill_willmott(pnboia['2015-01-01 02:00':'2015-03'],cfsv2['2015-01-01 02:00':'2015-03'])
plot(itajai_original['2015-01':'2015-01'],cfsv2['2015-01':'2015-01'],'PNBOIA/%s'%('Itajai'),'CFSv2',suptitle=suptitlte)
plot(itajai_corrected['2015-01':'2015-01'],cfsv2['2015-01':'2015-01'],'PNBOIA/%s'%('Itajai'),'CFSv2',suptitle=suptitlte)
