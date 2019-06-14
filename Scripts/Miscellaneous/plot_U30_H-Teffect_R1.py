"""
Script plots figure for manuscript of stratosphere variables of U10 and Z30

Notes
-----
    Author : Zachary Labe
    Date   : 17 October 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput_AllMembers as MO
import calc_Utilities as UT
import cmocean
import itertools
import scipy.stats as sts
from netCDF4 import Dataset

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO comparisons - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
qbophase = ['pos','non','neg']
period = 'D'

def calcVarResp(varnames,period,qbophase):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
                                                'surface')
    lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames,'FIT',
                                                'surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)
    
    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+101)
    
    ### Concatonate runs
    runs = [tashit,tasfict,tasfit]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'ON': 
        tas_mo = np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1) 
    elif period == 'DJ':     
        tas_mo = np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
                                                lon,'surface',1) 
    elif period == 'FM':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
    elif period == 'DJF':
        tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'surface',1)   
    elif period == 'M':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,2,:,:]
    elif period == 'D':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:]
    elif period == 'N':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-2,:,:]
    elif period == 'ND':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:]
    tas_mofitpos = tas_mo[2][pos_fit,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:]
    tas_mofitneg = tas_mo[2][neg_fit,:,:]
    
    
    ### Slice U10 at 65N
    latq = np.where((lat >= 64.5) & (lat <= 65.5))[0]
    lat = lat[latq].squeeze()
    tas_mohitpos = tas_mohitpos[:,latq,:].squeeze()
    tas_mohitneg = tas_mohitneg[:,latq,:].squeeze()
    tas_mofictpos = tas_mofictpos[:,latq,:].squeeze()
    tas_mofictneg = tas_mofictneg[:,latq,:].squeeze()
    tas_mofitpos = tas_mofitpos[:,latq,:].squeeze()
    tas_mofitneg = tas_mofitneg[:,latq,:].squeeze()
    
    ### Take zonal mean
    hitpos = np.nanmean(tas_mohitpos,axis=1)
    hitneg = np.nanmean(tas_mohitneg,axis=1)
    fictpos = np.nanmean(tas_mofictpos,axis=1)
    fictneg = np.nanmean(tas_mofictneg,axis=1)
    fitpos = np.nanmean(tas_mofitpos,axis=1)
    fitneg = np.nanmean(tas_mofitneg,axis=1)
    
    runs = [hitpos,hitneg,fictpos,fictneg,fitpos,fitneg]
    
    return runs,lat,lon

def readControl(period):
    directory = '/surtsey/zlabe/simu/CTLQ/monthly/'
    filename = directory + 'U30_1801-2000.nc'
    
    datac = Dataset(filename)
    tashitq = datac.variables['U30'][:]
    lat = datac.variables['latitude'][:]
    lon = datac.variables['longitude'][:]
    datac.close()
    
    ### Reshape array (200,12,96,144)
    tashit = np.reshape(tashitq,(tashitq.shape[0]//12,12,lat.shape[0],lon.shape[0]))
    
    ### Concatonate runs
    runs = [tashit]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'DJF':
        tas_mo,tas_mo = UT.calcDecJanFeb(runs[0],runs[0],lat,
                                              lon,'surface',1)   
    elif period == 'D':
        tas_mo = runs[0][:,-1,:,:]

    ### Read in QBO winters
    filenamefctlqp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
    filenamefctlqn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
    pos_ctlq = np.genfromtxt(filenamefctlqp,unpack=True,usecols=[0],dtype='int')
    neg_ctlq = np.genfromtxt(filenamefctlqn,unpack=True,usecols=[0],dtype='int')
    
    ### Composite by QBO phase
    if period == 'DJF':
        upos = tas_mo[pos_ctlq[:-1],:,:]
    else:
        upos = tas_mo[pos_ctlq,:,:]
    uneg = tas_mo[neg_ctlq,:,:]
    
    ### Slice U10 at 65N
    latq = np.where((lat >= 64.5) & (lat <= 65.5))[0]
    lat = lat[latq].squeeze()
    upos = upos[:,latq,:].squeeze()
    uneg = uneg[:,latq,:].squeeze()
    
    ### Take zonal mean 
    uupos = np.nanmean(upos,axis=1)
    uuneg = np.nanmean(uneg,axis=1)
    
    return uupos,uuneg

### Call variables
modu30,lat,lon = calcVarResp('U30',period,qbophase)
cpos,cneg = readControl(period)

### Calculate HT effect [hitpos,hitneg,fictpos,fictneg,fitpos,fitneg]
diffhit = np.nanmean(modu30[1],axis=0)-np.nanmean(modu30[0],axis=0)
difffict = np.nanmean(modu30[3],axis=0)-np.nanmean(modu30[2],axis=0)
difffit = np.nanmean(modu30[5],axis=0)-np.nanmean(modu30[4],axis=0)
diffcon = np.nanmean(cneg,axis=0) - np.nanmean(cpos,axis=0)
print('\nHIT -  HTeffect = %s m/s' % np.round(diffhit,3))
print('FICT - HTeffect = %s m/s' % np.round(difffict,3))
print('FIT -  HTeffect = %s m/s' % np.round(difffit,3))
print('CTLQ - HTeffect = %s m/s' % np.round(diffcon,3))

t,phit= sts.ttest_ind(modu30[1],modu30[0],equal_var=True) 
t,pfict= sts.ttest_ind(modu30[3],modu30[2],equal_var=True) 
t,pfit= sts.ttest_ind(modu30[5],modu30[4],equal_var=True) 
t,pcon= sts.ttest_ind(cneg,cpos,equal_var=True) 

t,panomcon = sts.ttest_ind(cneg,cpos,
                        equal_var=True)
t,panomfict = sts.ttest_ind(modu30[1] - modu30[0][:64],
                        modu30[3] - modu30[2][:64],
                        equal_var=True)
t,panomfit = sts.ttest_ind(modu30[1] - modu30[0][:64],
                        modu30[5] - modu30[4][:64],
                        equal_var=True)

###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJF
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

print('\nCompleted: Script done!')

