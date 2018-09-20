import os
from calendar import monthrange
import numpy as np
from datetime import date

def saveLog(log,message):
    log.write(message)

def select_user(dataset,dataset_login='/media/danilo/Danilo/mestrado/dataset_login'):
    dataset = "#%s"%dataset

    with open(dataset_login) as f:
        for line_i,line in enumerate(f.readlines()):
            if line.replace('\n','') == dataset:
                i = line_i
    with open(dataset_login) as f:
        user = f.readlines()[i+1]
    with open(dataset_login) as f:
        pasS = f.readlines()[i+2]

    return user.replace('\n',''),pasS.replace('\n','')


FTP_HEADER = 'ftp://my.cmems-du.eu/Core/GLOBAL_REANALYSIS_PHY_001_025/dataset-global-reanalysis-phy-001-025-ran-fr-glorys2v4-daily/'
DATA_DIR   = os.getcwd()+'/'

# selecting username and password
user,password = select_user('mercator',dataset_login='dataset_login')

# create log file
log = open(DATA_DIR+'log.system',mode='w')

# region to sellonlatbox
domain = '-50 -40 -30 -20'.split(' ')

YEARS = np.arange(2010,2016)
MONTH = [int(m) for m in '12,01,02,03'.split(',')]

# creating vector with dates
DATES = []
for y in YEARS:
    for m in MONTH:
        how_many_days = monthrange(y,m)[1]
        # print(how_many_days)
        for d in np.arange(1,how_many_days+1):
            # using datetime.date to create a date forma like 2014-01-01
            DATES.append(date(y,m,d))

# downloading
for date in DATES:
    y = str(date.year)
    m = str(date.month).zfill(2)
    d = str(date.isoformat().replace('-',''))

    # log info
    saveLog(log,'Downloading %s \n'%(d))

    url = FTP_HEADER + str(y) + '/' + str(m) + '/mercatorglorys2v4_gl4_mean_' + str(d) + '*'
    os.system("wget --user=%s --password='%s' " % (user,password)+ url)

# cutting region
for fname in glob.glob(DATA_DIR+'mercatorglory*'):
    outName = fname.split('/')[-1].split('_')[-2] + '.nc'
    os.system('cdo sellonlatbox,%s,%s,%s,%s %s %s' % (domain[0],domain[1],domain[2],domain[3],fname,outName))
    os.system('rm %s' % (fname))

# close log file
log.close()
