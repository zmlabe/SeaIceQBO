"""
Plot computes histograms of stratospheric variability in December over the 
polar cap

Notes
-----
    Author : Zachary Labe
    Date   : 5 June 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import read_DailyOutput_AllRegional as DOR
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting U10 Stratosphere Histograms - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = 'U'
runnames = [r'\textbf{FIT}',r'\textbf{HIT}',r'\textbf{FICT}']
qbophase = ['pos','non','neg']

### Call functions for variable profile data for polar cap
lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,'HIT',
                                           'profile')
lat,lon,time,lev,tasfit = DO.readMeanExperiAll('%s' % varnames,'FIT',
                                           'profile')
lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,'FICT',
                                            'profile')

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

filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)

filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)  

### Concatonate runs with selected level
levq = np.where(lev == 10)[0] # selected at 10 hPa
var_mo = [tashit[:,:,levq],tasfit[:,:,levq],tasfict[:,:,levq]]

### Composite by QBO phase    
var_mofitpos = var_mo[1][pos_fit,:]
var_mohitpos = var_mo[0][pos_hit,:]
var_mofictpos = var_mo[2][pos_fict,:]

var_mofitneg = var_mo[1][neg_fit,:]
var_mohitneg = var_mo[0][neg_hit,:]
var_mofictneg = var_mo[2][neg_fict,:]

### Calculate over DJF (90-180)
timeq = np.arange(90,120)
monthqq = 'D'
var_wfitpos = var_mofitpos[:,timeq]
var_whitpos = var_mohitpos[:,timeq]
var_wfictpos = var_mofictpos[:,timeq]

var_wfitneg = var_mofitneg[:,timeq]
var_whitneg = var_mohitneg[:,timeq]
var_wfictneg = var_mofictneg[:,timeq]

### Calculate difference over time average
difffit = np.nanmean(var_wfitneg - var_whitneg,axis=1).squeeze()
difffict = np.nanmean(var_wfictneg - var_whitneg,axis=1).squeeze()

difffit = np.nanmean(var_wfitpos - var_whitpos,axis=1).squeeze()
difffict = np.nanmean(var_wfictpos - var_whitpos,axis=1).squeeze()
 
###############################################################################
###############################################################################
###############################################################################    
### Plot histogram distributions
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
        
fig = plt.figure()
ax = plt.subplot(111) 

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('w')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',which='both',bottom=False)

#plt.axhline(0,color='dimgrey',linestyle='--',dashes=(0.9,1),linewidth=2)
#bx = plt.boxplot(dataq,0,'',patch_artist=True,showmeans=True,meanline=True,
#                 whis=[5,95])
#
#for i in bx['caps']:
#    i.set(color='k',linewidth=0)
#for whisker in bx['whiskers']:
#    whisker.set(color='dimgrey',linestyle='-',linewidth=2)
#for box in bx['boxes']: 
#    box.set(color='deepskyblue')
#for box in bx['means']:
#    box.set(color='r',linewidth=2,linestyle='-')
#for box in bx['medians']:
#    box.set(linewidth=0)
#    
#for i in range(len(dataq)):
#    y = dataq[i]
#    x = np.random.normal(1+i,0.04,size=len(y))
#    plt.plot(x,y,'r.',alpha=0.3,zorder=5)
#
#plt.ylabel(r'\textbf{U30 [m/s]}',color='k',fontsize=12)
#
#plt.yticks(np.arange(-20,21,5),list(map(str,np.arange(-20,21,5))))
#plt.xticks(np.arange(1,9,1),runnames) 
#plt.ylim([-18,18])

plt.savefig(directoryfigure + 'stratosphereDistribution_NET_%s.png' % monthqq,
            dpi=300)
