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

#### Call functions for variable profile data for polar cap
#lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,'HIT',
#                                           'profile')
#lat,lon,time,lev,tasfit = DO.readMeanExperiAll('%s' % varnames,'FIT',
#                                           'profile')
#lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,'FICT',
#                                            'profile')
#
#### Create 2d array of latitude and longitude
#lon2,lat2 = np.meshgrid(lon,lat)
#
##### Read in QBO phases 
#filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
#neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
#
#filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
#neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)
#
#filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
#neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)  

#### Concatonate runs with selected level
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
timeq = np.arange(90,121)
monthqq = 'D'
var_wfitpos = var_mofitpos[:,timeq].squeeze()
var_whitpos = var_mohitpos[:,timeq].squeeze()
var_wfictpos = var_mofictpos[:,timeq].squeeze()

var_wfitneg = var_mofitneg[:,timeq].squeeze()
var_whitneg = var_mohitneg[:,timeq].squeeze()
var_wfictneg = var_mofictneg[:,timeq].squeeze()

### PDF of stratospheric winds
ficteast = np.ravel(var_wfictneg)
fictwest = np.ravel(var_wfictpos)

fiteast = np.ravel(var_wfitneg)
fitwest = np.ravel(var_wfitpos)

hiteast = np.ravel(var_whitneg)
hitwest = np.ravel(var_whitpos)
 
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
        
num_bins = np.arange(-20,60,2)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(221)
#ax = plt.subplot(221) 

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#ax.tick_params(axis='x',which='both',bottom=False)

ne,binse,patchese = plt.hist(ficteast,num_bins,density=1,alpha=0.7,
                          facecolor='deepskyblue',edgecolor='deepskyblue',
                          range=(-20,50),label=r'\textbf{QBO-E}')
nw,binsw,patchesw = plt.hist(fictwest,num_bins,density=1,alpha=0.5,
                          facecolor='crimson',edgecolor='crimson',
                          range=(-20,50),label=r'\textbf{QBO-W}')

plt.yticks(np.arange(0,0.09,0.01),list(map(str,np.arange(0,0.09,0.01))),
           fontsize=10)
plt.xticks(np.arange(-20,60,10),list(map(str,np.arange(-20,60,10))),
           fontsize=10) 
plt.xlim([-20,50])
plt.ylim([0,0.07])

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

plt.title(r'\textbf{FICT',color='dimgrey',
                    fontsize=30)
plt.xlabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)

###############################################################################

ax = fig.add_subplot(222)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#ax.tick_params(axis='x',which='both',bottom=False)

ne,binse,patchese = plt.hist(fiteast,num_bins,density=1,alpha=0.7,
                          facecolor='deepskyblue',edgecolor='deepskyblue',
                          range=(-20,50),label=r'\textbf{QBO-E}')
nw,binsw,patchesw = plt.hist(fitwest,num_bins,density=1,alpha=0.5,
                          facecolor='crimson',edgecolor='crimson',
                          range=(-20,50),label=r'\textbf{QBO-W}')

plt.yticks(np.arange(0,0.09,0.01),list(map(str,np.arange(0,0.09,0.01))),
           fontsize=10)
plt.xticks(np.arange(-20,60,10),list(map(str,np.arange(-20,60,10))),
           fontsize=10) 
plt.xlim([-20,50])
plt.ylim([0,0.07])

plt.title(r'\textbf{FIT',color='dimgrey',
                    fontsize=30)
plt.xlabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

###############################################################################

ax = fig.add_subplot(223)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#ax.tick_params(axis='x',which='both',bottom=False)

ne,binse,patchese = plt.hist(ficteast,num_bins,density=1,alpha=0.7,
                          facecolor='mediumseagreen',edgecolor='mediumseagreen',
                          range=(-20,50),label=r'\textbf{FICT}')
nw,binsw,patchesw = plt.hist(hiteast,num_bins,density=1,alpha=0.4,
                          facecolor='indigo',edgecolor='indigo',
                          range=(-20,50),label=r'\textbf{HIT}')

plt.yticks(np.arange(0,0.09,0.01),list(map(str,np.arange(0,0.09,0.01))),
           fontsize=10)
plt.xticks(np.arange(-20,60,10),list(map(str,np.arange(-20,60,10))),
           fontsize=10) 
plt.xlim([-20,50])
plt.ylim([0,0.07])

plt.title(r'\textbf{QBO-E',color='dimgrey',
                    fontsize=30)
plt.xlabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

###############################################################################

ax = fig.add_subplot(224)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#ax.tick_params(axis='x',which='both',bottom=False)

ne,binse,patchese = plt.hist(fictwest,num_bins,density=1,alpha=0.7,
                          facecolor='mediumseagreen',edgecolor='mediumseagreen',
                          range=(-20,50),label=r'\textbf{FICT}')
nw,binsw,patchesw = plt.hist(hitwest,num_bins,density=1,alpha=0.4,
                          facecolor='indigo',edgecolor='indigo',
                          range=(-20,50),label=r'\textbf{HIT}')

plt.yticks(np.arange(0,0.09,0.01),list(map(str,np.arange(0,0.09,0.01))),
           fontsize=10)
plt.xticks(np.arange(-20,60,10),list(map(str,np.arange(-20,60,10))),
           fontsize=10) 
plt.xlim([-20,50])
plt.ylim([0,0.07])

plt.title(r'\textbf{QBO-W}',color='dimgrey',
                    fontsize=30)
plt.xlabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

plt.tight_layout()

plt.savefig(directoryfigure + 'stratosphereDistribution_U10_%s.png' % monthqq,
            dpi=300)
