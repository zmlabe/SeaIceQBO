"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
calculates correlations between the strength of the Siberian High and 
CDI.

Notes
-----
    Author : Zachary Labe
    Date   : 16 August 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import calc_Utilities as UT
import read_DailyOutput_AllMembers as DO
import scipy.stats as sts

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directorydata3 = '/home/zlabe/Documents/Research/SeaIceQBO/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating CDI correlations - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = True
varnames = ['T1000']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'

### Calculate Siberian High 
def readVariables(varnames,period,location):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitno2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'D':
        tas_mo= np.empty((2,tashit.shape[0],31,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,90:121,:,:]
    else:
        ValueError('Wrong period selected!')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:]
    
    tas_mohitnon = tas_mo[0][non_hit,:,:,:]
    tas_mofictnon = tas_mo[1][non_fict,:,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:]
    
    ### Compute comparisons for months - select region
    if location == 'Siberia':
        lonq = np.where((lon >=80) & (lon <=120))[0]
        ficthitpos = tas_mofictpos[:,:,:,lonq]
        ficthitnon = tas_mofictnon[:,:,:,lonq] 
        ficthitneg = tas_mofictneg[:,:,:,lonq]
        latq = np.where((lat >=40) & (lat <=65))[0]
        ficthitpos = ficthitpos[:,:,latq]
        ficthitnon = ficthitnon[:,:,latq] 
        ficthitneg = ficthitneg[:,:,latq]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
        ficthitnon = UT.calc_weightedAve(ficthitnon,lat2s)
        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
        
    elif location == 'Aleutian': # North Pacific Index (Trenberth and Hurrell, 1994)
        lonq = np.where((lon >=160) & (lon <=220))[0]
        ficthitpos = tas_mofictpos[:,:,:,lonq]
        ficthitnon = tas_mofictnon[:,:,:,lonq] 
        ficthitneg = tas_mofictneg[:,:,:,lonq]
        latq = np.where((lat >=30) & (lat <=65))[0]
        ficthitpos = ficthitpos[:,:,latq]
        ficthitnon = ficthitnon[:,:,latq] 
        ficthitneg = ficthitneg[:,:,latq]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
        ficthitnon = UT.calc_weightedAve(ficthitnon,lat2s)
        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
    diffruns = [ficthitpos.squeeze(),ficthitnon.squeeze(),ficthitneg.squeeze()]
    
    return diffruns,lat,lon,lev

slp,lat,lon,lev = readVariables('SLP',period,'Siberia')

### Calculate December averages for SLP
slp_ficthitpos = np.nanmean(slp[0],axis=1)
slp_ficthitnon = np.nanmean(slp[1],axis=1)
slp_ficthitneg = np.nanmean(slp[2],axis=1)

### Meshgrid
lon2,lat2 = np.meshgrid(lon,lat)

### Read in CDI
data1 = Dataset(directorydata3 + 'ColdDayIndex_December.nc')
fitpos = data1.variables['cdi_fitpos'][:]
fitneg = data1.variables['cdi_fitneg'][:]
fictpos = data1.variables['cdi_fictpos'][:]
fictneg = data1.variables['cdi_fictneg'][:]
data1.close()

### Slice CDI's over region
lonq = np.where((lon >=70) & (lon <=140))[0]
fitposq = fitpos[:,:,lonq]
fitnegq = fitneg[:,:,lonq] 
fictposq = fictpos[:,:,lonq]
fictnegq = fictneg[:,:,lonq]
latq = np.where((lat >=35) & (lat <=60))[0]
fitposqq  = fitposq[:,latq,:]
fitnegqq = fitnegq[:,latq,:] 
fictposqq = fictposq[:,latq,:]
fictnegqq = fictnegq[:,latq,:]
lat2sq = lat2[latq,:]
lat2s = lat2sq[:,lonq]
fitposc = UT.calc_weightedAve(fitposqq,lat2s)
fitnegc = UT.calc_weightedAve(fitnegqq,lat2s)
fictposc = UT.calc_weightedAve(fictposqq,lat2s)
fictnegc = UT.calc_weightedAve(fictnegqq,lat2s)

###############################################################################
###############################################################################
###############################################################################    
### Plot scatter plots
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})

### Adjust axes in time series plots 
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

fig = plt.figure(figsize=(7,7))

ax = fig.add_subplot(111)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.scatter(slp_ficthitpos,fictposc,color='crimson',marker='o',label=r'\textbf{QBO-W}')
timex = np.linspace(1025,1040,len(slp_ficthitpos))
slope, intercept, r_value1, p_value1, std_err = sts.linregress(slp_ficthitpos,fictposc)
line1 = slope*timex + intercept
plt.plot(timex,line1,color='crimson',linewidth=3)

plt.scatter(slp_ficthitneg,fictnegc,color='deepskyblue',marker='o',label=r'\textbf{QBO-E}')
timex = np.linspace(1025,1040,len(slp_ficthitneg))
slope, intercept, r_value2, p_value2, std_err = sts.linregress(slp_ficthitneg,fictnegc)
line1 = slope*timex + intercept
plt.plot(timex,line1,color='deepskyblue',linewidth=3)

plt.yticks(np.arange(-80,61,10),list(map(str,np.arange(-80,61,10))),fontsize=10) 
plt.ylim([-50,50])
plt.xticks(np.arange(900,1041,5),list(map(str,np.arange(900,1041,5))),fontsize=10) 
plt.xlim([1025,1040])
plt.xlabel(r'\textbf{hPa}',fontsize=12)
plt.ylabel(r'\textbf{CDI}',fontsize=12)
        
plt.legend(shadow=False,fontsize=12,loc='upper center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.9,1),
           labelspacing=0.2,columnspacing=0.4,handletextpad=0.4)

plt.savefig(directoryfigure + 'CDISiberianHighRelationship_FICT.png',dpi=300)