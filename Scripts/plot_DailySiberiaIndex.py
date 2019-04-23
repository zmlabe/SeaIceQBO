"""
Plot manuscript figures for daily indices over Siberia

Notes
-----
    Author : Zachary Labe
    Date   : 20 March 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import scipy.stats as sts
from netCDF4 import Dataset

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directorydata3 = '/home/zlabe/Documents/Research/SeaIceQBO/Data/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating Siberian Daily Indices - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
runnames = [r'HIT',r'FICT']
qbophase = ['pos','non','neg']
def readU10(varnames):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)

    tas_mo = [tashit,tasfict]

    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:]
    
    ### Compute comparisons for months - select region
    if varnames == 'U10':
        latq = np.where((lat >=65) & (lat <=90))[0]
        fictpos = tas_mofictpos
        fictneg = tas_mofictneg
        fictpos = fictpos[:,:,latq]
        fictneg = fictneg[:,:,latq]
        lat2s = lat2[latq,:]
        fictpos = UT.calc_weightedAve(fictpos,lat2s)
        fictneg = UT.calc_weightedAve(fictneg,lat2s)
        
        hitpos = tas_mohitpos
        hitneg = tas_mohitneg
        hitpos = hitpos[:,:,latq]
        hitneg = hitneg[:,:,latq]
        hitpos = UT.calc_weightedAve(hitpos,lat2s)
        hitneg = UT.calc_weightedAve(hitneg,lat2s)
        
    diffruns = [fictpos.squeeze(),fictneg.squeeze(),hitpos.squeeze(),hitneg.squeeze()]
    
    return diffruns

def readVariables(varnames):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)

    tas_mo = [tashit,tasfict]

    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:]
    
    ### Compute comparisons for months - select region
    if varnames == 'SLP':
        lonq = np.where((lon >=80) & (lon <=120))[0]
        fictpos = tas_mofictpos[:,:,:,lonq]
        fictneg = tas_mofictneg[:,:,:,lonq]
        latq = np.where((lat >=40) & (lat <=65))[0]
        fictpos = fictpos[:,:,latq]
        fictneg = fictneg[:,:,latq]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        fictpos = UT.calc_weightedAve(fictpos,lat2s)
        fictneg = UT.calc_weightedAve(fictneg,lat2s)
        
        hitpos = tas_mohitpos[:,:,:,lonq]
        hitneg = tas_mohitneg[:,:,:,lonq]
        hitpos = hitpos[:,:,latq]
        hitneg = hitneg[:,:,latq]
        hitpos = UT.calc_weightedAve(hitpos,lat2s)
        hitneg = UT.calc_weightedAve(hitneg,lat2s)
    elif varnames == 'T1000':
        lonq = np.where((lon >=70) & (lon <=140))[0]
        fictpos = tas_mofictpos[:,:,:,lonq]
        fictneg = tas_mofictneg[:,:,:,lonq]
        latq = np.where((lat >=35) & (lat <=60))[0]
        fictpos = fictpos[:,:,latq]
        fictneg = fictneg[:,:,latq]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        fictpos = UT.calc_weightedAve(fictpos,lat2s) - 273.15
        fictneg = UT.calc_weightedAve(fictneg,lat2s) - 273.15
        
        hitpos = tas_mohitpos[:,:,:,lonq]
        hitneg = tas_mohitneg[:,:,:,lonq]
        hitpos = hitpos[:,:,latq]
        hitneg = hitneg[:,:,latq]
        hitpos = UT.calc_weightedAve(hitpos,lat2s) - 273.15
        hitneg = UT.calc_weightedAve(hitneg,lat2s) - 273.15
        
    diffruns = [fictpos.squeeze(),fictneg.squeeze(),hitpos.squeeze(),hitneg.squeeze()]
    
    return diffruns

def readMHF100():
    ### Call functions for MHF data at 100 hPa
    lat,lon,time,lev,varhit = DO.readMeanExperiAll('MHF100','HIT','profile')
    lat,lon,time,lev,varfit = DO.readMeanExperiAll('MHF100','FIT','profile')
    lat,lon,time,lev,varfict = DO.readMeanExperiAll('MHF100','FICT','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)    
    ### Concatonate runs
    var_mo = [varhit,varfit,varfict]
    
    ### Composite by QBO phase and take zonal mean
    hitpos = var_mo[0][pos_hit,:].squeeze()
    fictpos = var_mo[2][pos_fict,:].squeeze()

    hitneg = var_mo[0][neg_hit,:].squeeze()
    fictneg = var_mo[2][neg_fict,:].squeeze()
    
    ### Take zonal mean
    hitpos = np.nanmean(hitpos,axis=2)
    fictpos = np.nanmean(fictpos,axis=2)
    
    hitneg = np.nanmean(hitneg,axis=2)
    fictneg = np.nanmean(fictneg,axis=2)
    
    ### Store data in a list
    diffruns = [fictpos,fictneg,hitpos,hitneg]
    
    return diffruns

def calcLinearTrend(data,length):
    """
    Calculates moving average for n number of months
    
    Parameters
    ----------
    data : 1d array
        [time series data]
    length : integer
        [n months]
        
    Returns
    -------
    ave : 1d array
        time series of smoothed data from moving average
    
    Usage
    -----
    ave = calcLinearTrend(data,years,length)
    """
    print('\n>>> Using calcMovingAverage function!')
    
    ### Calculate moving average for n months (length)
    aven = np.convolve(data, np.ones((length,))/length, mode='valid') 
    print('Completed: *%s MONTHS* averages!' % length)
    
    ### Append nans for consistent time
    empty = np.array([np.nan]*(length-1))
    ave = np.append(empty,aven,axis=0)
    
    print('*Completed: Finished calcMovingAverage function!\n')    
    return ave

### Call functions
diffmhf = readMHF100()
diffslp = readVariables('SLP')
difftemp = readVariables('T1000')
diffu10 = readU10('U10')
smooth = True
smoothq = 10
    
## Calculate MHF anomalies for neg
mhfneg = (np.nanmean(diffmhf[1],axis=0) - np.nanmean(diffmhf[3],axis=0))/ \
                    np.nanstd(diffmhf[3],axis=0)
slpneg = (np.nanmean(diffslp[1],axis=0) - np.nanmean(diffslp[3],axis=0))/ \
                    np.nanstd(diffslp[3],axis=0) 
tempneg = (np.nanmean(difftemp[1],axis=0) - np.nanmean(difftemp[3],axis=0))/ \
                    np.nanstd(difftemp[3],axis=0)                     
u10neg = (np.nanmean(diffu10[1],axis=0) - np.nanmean(diffu10[3],axis=0))/ \
                    np.nanstd(diffu10[3],axis=0)     
    
#### Calculate MHF anomalies for pos
mhfpos = (np.nanmean(diffmhf[0],axis=0) - np.nanmean(diffmhf[2],axis=0))/ \
                    np.nanstd(diffmhf[2],axis=0)
slppos = (np.nanmean(diffslp[0],axis=0) - np.nanmean(diffslp[2],axis=0))/ \
                    np.nanstd(diffslp[2],axis=0) 
temppos = (np.nanmean(difftemp[0],axis=0) - np.nanmean(difftemp[2],axis=0))/ \
                    np.nanstd(difftemp[2],axis=0)                     
u10pos = (np.nanmean(diffu10[0],axis=0) - np.nanmean(diffu10[2],axis=0))/ \
                    np.nanstd(diffu10[2],axis=0)  
                    
### Calculate smoothing                    
if smooth == True:
    mhfneg = calcLinearTrend(mhfneg,smoothq)
    slpneg = calcLinearTrend(slpneg,smoothq)
    tempneg = calcLinearTrend(tempneg,smoothq)
    u10neg = calcLinearTrend(u10neg,smoothq)
    mhfpos = calcLinearTrend(mhfpos,smoothq)
    slppos = calcLinearTrend(slppos,smoothq)
    temppos = calcLinearTrend(temppos,smoothq)
    u10pos = calcLinearTrend(u10pos,smoothq)
 
##############################################################################
##############################################################################
##############################################################################
###### Plot Daily Indices
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.axhline(0,linestyle=':',linewidth=2,color='dimgrey',dashes=(1,0.3))

plt.plot(mhfneg,color='darkgreen',label=r'\textbf{MHF100}',linewidth=2)
plt.plot(slpneg,color='darkred',label=r'\textbf{SLP}',linewidth=2)
plt.plot(tempneg,color='darkblue',label=r'\textbf{T1000}',linewidth=2)
plt.plot(u10neg,color='darkorange',label=r'\textbf{U10}',linewidth=2)

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{Normalized Indices}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,0.5),list(map(str,np.arange(-5,6,0.5))),fontsize=9)
plt.ylim([-0.5,0.5])

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

if smooth == True:
    plt.savefig(directoryfigure + 'SiberianTimeSeries_neg_smooth.png',dpi=300)
elif smooth == False:
    plt.savefig(directoryfigure + 'SiberianTimeSeries_neg.png',dpi=300)

###############################################################################

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.axhline(0,linestyle=':',linewidth=2,color='dimgrey',dashes=(1,0.3))

plt.plot(mhfpos,color='darkgreen',label=r'\textbf{MHF100}',linewidth=2)
plt.plot(slppos,color='darkred',label=r'\textbf{SLP}',linewidth=2)
plt.plot(temppos,color='darkblue',label=r'\textbf{T1000}',linewidth=2)
plt.plot(u10pos,color='darkorange',label=r'\textbf{U10}',linewidth=2)

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{Normalized Indices}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,0.5),list(map(str,np.arange(-5,6,0.5))),fontsize=9)
plt.ylim([-0.5,0.5])

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

if smooth == True:
    plt.savefig(directoryfigure + 'SiberianTimeSeries_pos_smooth.png',dpi=300)
elif smooth == False:
    plt.savefig(directoryfigure + 'SiberianTimeSeries_pos.png',dpi=300)

### Calculate correlations
if smooth == False:
    negcorr,negp = sts.pearsonr(slpneg,tempneg)
    poscorr,posp = sts.pearsonr(slppos,temppos)
elif smooth == True:
    negcorr,negp = sts.pearsonr(slpneg[smoothq-1:],tempneg[smoothq-1:])
    poscorr,posp = sts.pearsonr(slppos[smoothq-1:],temppos[smoothq-1:])
print('\nQBO-E---> %s correlation at p-val=%s' % (negcorr,negp)) 
print('QBO-W---> %s correlation at p-val=%s' % (poscorr,posp)) 