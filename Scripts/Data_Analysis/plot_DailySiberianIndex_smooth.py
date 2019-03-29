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
import cmocean

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
    if varnames == 'U10' or varnames == 'Z30':
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

def calcSmooth(data,length):
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

### Calculate smoothing                    
if smooth == True:
    fictmhfpos = np.empty((diffmhf[0].shape))
    fictmhfneg = np.empty((diffmhf[1].shape))
    fictslppos = np.empty((diffmhf[0].shape))
    fictslpneg = np.empty((diffmhf[1].shape))
    ficttemppos = np.empty((diffmhf[0].shape))
    ficttempneg = np.empty((diffmhf[1].shape))
    fictu10pos = np.empty((diffmhf[0].shape))
    fictu10neg = np.empty((diffmhf[1].shape))
    hitmhfpos = np.empty((diffmhf[0].shape))
    hitmhfneg = np.empty((diffmhf[1].shape))
    hitslppos = np.empty((diffmhf[0].shape))
    hitslpneg = np.empty((diffmhf[1].shape))
    hittemppos = np.empty((diffmhf[0].shape))
    hittempneg = np.empty((diffmhf[1].shape))
    hitu10pos = np.empty((diffmhf[0].shape))
    hitu10neg = np.empty((diffmhf[1].shape))
    for i in range(diffmhf[0].shape[0]):
        fictmhfpos[i,:] = calcSmooth(diffmhf[0][i,:],smoothq)
        hitmhfpos[i,:] = calcSmooth(diffmhf[2][i,:],smoothq)
        
        fictslppos[i,:] = calcSmooth(diffslp[0][i,:],smoothq)
        hitslppos[i,:] = calcSmooth(diffslp[2][i,:],smoothq)
        
        ficttemppos[i,:] = calcSmooth(difftemp[0][i,:],smoothq)
        hittemppos[i,:] = calcSmooth(difftemp[2][i,:],smoothq)
        
        fictu10pos[i,:] = calcSmooth(diffu10[0][i,:],smoothq)
        hitu10pos[i,:] = calcSmooth(diffu10[2][i,:],smoothq)
    for i in range(diffmhf[1].shape[0]):
        fictmhfneg[i,:] = calcSmooth(diffmhf[1][i,:],smoothq)
        hitmhfneg[i,:] = calcSmooth(diffmhf[3][i,:],smoothq)
        
        fictslpneg[i,:] = calcSmooth(diffslp[1][i,:],smoothq)
        hitslpneg[i,:] = calcSmooth(diffslp[3][i,:],smoothq)
        
        ficttempneg[i,:] = calcSmooth(difftemp[1][i,:],smoothq)
        hittempneg[i,:] = calcSmooth(difftemp[3][i,:],smoothq)
        
        fictu10neg[i,:] = calcSmooth(diffu10[1][i,:],smoothq)
        hitu10neg[i,:] = calcSmooth(diffu10[3][i,:],smoothq)
    
## Calculate MHF anomalies for neg
mhfneg = (np.nanmean(fictmhfneg,axis=0) - np.nanmean(hitmhfneg,axis=0))/ \
                    np.nanstd(hitmhfneg,axis=0)
slpneg = (np.nanmean(fictslpneg,axis=0) - np.nanmean(hitslpneg,axis=0))/ \
                    np.nanstd(hitslpneg,axis=0) 
tempneg = (np.nanmean(ficttempneg,axis=0) - np.nanmean(hittempneg,axis=0))/ \
                    np.nanstd(hittempneg,axis=0)                     
u10neg = (np.nanmean(fictu10neg,axis=0) - np.nanmean(hitu10neg,axis=0))/ \
                    np.nanstd(hitu10neg,axis=0)     
    
#### Calculate MHF anomalies for pos
mhfpos = (np.nanmean(fictmhfpos,axis=0) - np.nanmean(hitmhfpos,axis=0))/ \
                    np.nanstd(hitmhfpos,axis=0)
slppos = (np.nanmean(fictslppos,axis=0) - np.nanmean(hitslppos,axis=0))/ \
                    np.nanstd(hitslppos,axis=0) 
temppos = (np.nanmean(ficttemppos,axis=0) - np.nanmean(hittemppos,axis=0))/ \
                    np.nanstd(hittemppos,axis=0)                     
u10pos = (np.nanmean(fictu10pos,axis=0) - np.nanmean(hitu10pos,axis=0))/ \
                    np.nanstd(hitu10pos,axis=0)   
                    
### Calculate statistical significance
stat,pmhfneg = UT.calc_indttest(fictmhfneg,hitmhfneg)
stat,pmhfpos = UT.calc_indttest(fictmhfpos,hitmhfpos)
pmhfneg[np.isnan(pmhfneg)] = 0.0
pmhfpos[np.isnan(pmhfpos)] = 0.0
pvalsmhfneg = mhfneg * pmhfneg
pvalsmhfneg[pvalsmhfneg == 0.0] = np.nan
pvalsmhfpos = mhfpos * pmhfpos
pvalsmhfpos[pvalsmhfpos == 0.0] = np.nan

stat,pslpneg = UT.calc_indttest(fictslpneg,hitslpneg)
stat,pslppos = UT.calc_indttest(fictslppos,hitslppos)
pslpneg[np.isnan(pslpneg)] = 0.0
pslppos[np.isnan(pslppos)] = 0.0
pvalsslpneg = slpneg * pslpneg
pvalsslpneg[pvalsslpneg == 0.0] = np.nan
pvalsslppos = slppos * pslppos
pvalsslppos[pvalsslppos == 0.0] = np.nan

stat,ptempneg = UT.calc_indttest(ficttempneg,hittempneg)
stat,ptemppos = UT.calc_indttest(ficttemppos,hittemppos)
ptempneg[np.isnan(ptempneg)] = 0.0
ptemppos[np.isnan(ptemppos)] = 0.0
pvalstempneg = tempneg * ptempneg
pvalstempneg[pvalstempneg == 0.0] = np.nan
pvalstemppos = temppos * ptemppos
pvalstemppos[pvalstemppos == 0.0] = np.nan

stat,pu10neg = UT.calc_indttest(fictu10neg,hitu10neg)
stat,pu10pos = UT.calc_indttest(fictu10pos,hitu10pos)
pu10neg[np.isnan(pu10neg)] = 0.0
pu10pos[np.isnan(pu10pos)] = 0.0
pvalsu10neg = u10neg * pu10neg
pvalsu10neg[pvalsu10neg == 0.0] = np.nan
pvalsu10pos = u10pos * pu10pos
pvalsu10pos[pvalsu10pos == 0.0] = np.nan
 
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

plt.plot(mhfneg,color=cmocean.cm.thermal(0.1),
         linewidth=1)
plt.plot(slpneg,color=cmocean.cm.thermal(0.3),
         linewidth=1)
plt.plot(tempneg,color=cmocean.cm.thermal(0.61),
         linewidth=1)
plt.plot(u10neg,color=cmocean.cm.thermal(0.77),
         linewidth=1)

plt.plot(pvalsmhfneg,color=cmocean.cm.thermal(0.1),
         linewidth=3,label=r'\textbf{MHF100}')
plt.plot(pvalsslpneg,color=cmocean.cm.thermal(0.3),
         linewidth=3,label=r'\textbf{SHI}')
plt.plot(pvalstempneg,color=cmocean.cm.thermal(0.61),
         linewidth=3,label=r'\textbf{T1000}')
plt.plot(pvalsu10neg,color=cmocean.cm.thermal(0.77),
         linewidth=3,label=r'\textbf{U10}')

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{Normalized Indices}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,0.5),list(map(str,np.arange(-5,6,0.5))),fontsize=9)
plt.ylim([-0.8,0.8])

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

plt.savefig(directoryfigure + 'SiberianTimeSeries_neg_smooth2.png',dpi=300)

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

plt.plot(mhfpos,color=cmocean.cm.thermal(0.1),linewidth=1)
plt.plot(slppos,color=cmocean.cm.thermal(0.3),linewidth=1)
plt.plot(temppos,color=cmocean.cm.thermal(0.61),linewidth=1)
plt.plot(u10pos,color=cmocean.cm.thermal(0.77),linewidth=1)

plt.plot(pvalsmhfpos,color=cmocean.cm.thermal(0.1),
         linewidth=3,label=r'\textbf{MHF100}')
plt.plot(pvalsslppos,color=cmocean.cm.thermal(0.3),
         linewidth=3,label=r'\textbf{SHI}')
plt.plot(pvalstemppos,color=cmocean.cm.thermal(0.61),
         linewidth=3,label=r'\textbf{T1000}')
plt.plot(pvalsu10pos,color=cmocean.cm.thermal(0.77),
         linewidth=3,label=r'\textbf{U10}')

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{Normalized Indices}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,0.5),list(map(str,np.arange(-5,6,0.5))),fontsize=9)
plt.ylim([-0.8,0.8])

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

plt.savefig(directoryfigure + 'SiberianTimeSeries_pos_smooth2.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Calculate correlations
if smooth == False:
    negcorr,negp = sts.pearsonr(slpneg,tempneg)
    poscorr,posp = sts.pearsonr(slppos,temppos)
elif smooth == True:
    negcorr,negp = sts.pearsonr(slpneg[smoothq-1:],tempneg[smoothq-1:])
    poscorr,posp = sts.pearsonr(slppos[smoothq-1:],temppos[smoothq-1:])
print('\nQBO-E---> %s correlation at p-val=%s' % (negcorr,negp)) 
print('QBO-W---> %s correlation at p-val=%s' % (poscorr,posp)) 

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(212)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.axhline(0,linestyle=':',linewidth=2,color='dimgrey',dashes=(1,0.3))

plt.plot(mhfneg,color=cmocean.cm.thermal(0.1),
         linewidth=1)
plt.plot(slpneg,color=cmocean.cm.thermal(0.3),
         linewidth=1)
plt.plot(tempneg,color=cmocean.cm.thermal(0.61),
         linewidth=1)
plt.plot(u10neg,color=cmocean.cm.thermal(0.77),
         linewidth=1)

plt.plot(pvalsmhfneg,color=cmocean.cm.thermal(0.1),
         linewidth=3,label=r'\textbf{MHF100}')
plt.plot(pvalsslpneg,color=cmocean.cm.thermal(0.3),
         linewidth=3,label=r'\textbf{SHI}')
plt.plot(pvalstempneg,color=cmocean.cm.thermal(0.61),
         linewidth=3,label=r'\textbf{T1000}')
plt.plot(pvalsu10neg,color=cmocean.cm.thermal(0.77),
         linewidth=3,label=r'\textbf{U10}')

l = plt.legend(shadow=False,fontsize=7,loc='lower center',
           fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,-0.054),
           labelspacing=0.4,columnspacing=1.2,handletextpad=0.6)

plt.yticks(np.arange(-0.8,0.9,0.4),list(map(str,
           np.round(np.arange(-0.8,0.9,0.4),2))),fontsize=6)
plt.ylim([-0.8,0.8])

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

###############################################################################
ax = plt.subplot(211)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.annotate(r'\textbf{QBO-W}',xy=(30,0),xytext=(0.2,0.85),
             textcoords='figure fraction',color='k',
             fontsize=17,rotation=0,ha='center',va='center')
plt.annotate(r'\textbf{QBO-E}',xy=(30,0),xytext=(0.2,0.44),
             textcoords='figure fraction',color='k',
             fontsize=17,rotation=0,ha='center',va='center')
plt.annotate(r'\textbf{Normalized Indices}',xy=(30,0),xytext=(0.05,0.5),
             textcoords='figure fraction',color='k',
             fontsize=10,rotation=90,ha='center',va='center')

plt.axhline(0,linestyle=':',linewidth=2,color='dimgrey',dashes=(1,0.3))

plt.plot(mhfpos,color=cmocean.cm.thermal(0.1),linewidth=1)
plt.plot(slppos,color=cmocean.cm.thermal(0.3),linewidth=1)
plt.plot(temppos,color=cmocean.cm.thermal(0.61),linewidth=1)
plt.plot(u10pos,color=cmocean.cm.thermal(0.77),linewidth=1)

plt.plot(pvalsmhfpos,color=cmocean.cm.thermal(0.1),
         linewidth=3,label=r'\textbf{MHF100}')
plt.plot(pvalsslppos,color=cmocean.cm.thermal(0.3),
         linewidth=3,label=r'\textbf{SHI}')
plt.plot(pvalstemppos,color=cmocean.cm.thermal(0.61),
         linewidth=3,label=r'\textbf{T1000}')
plt.plot(pvalsu10pos,color=cmocean.cm.thermal(0.77),
         linewidth=3,label=r'\textbf{U10}')

plt.yticks(np.arange(-0.8,0.9,0.4),list(map(str,
           np.round(np.arange(-0.8,0.9,0.4),2))),fontsize=6)
plt.ylim([-0.8,0.8])

xlabels = [] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
plt.xlim([30,210])

plt.savefig(directoryfigure + 'SiberianTimeSeries_subplot_p95.png',dpi=300)