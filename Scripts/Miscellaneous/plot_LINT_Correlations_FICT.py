"""
Plot comparisons between WACCM4 sea ice experiments. Script plots a 
rolling pattern correlation between the climatologicaland forced waves 
(1, 2, and all)

Notes
-----
    Author : Zachary Labe
    Date   : 20 July 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/home/zlabe/Documents/Research/SeaIceQBO/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Linear Interference Correlations - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
runnames = [r'HIT',r'FICT']
qbophase = ['pos','non','neg']
levels = ['col','tropo','strato']
wave = ['wave1','wave2','wave_all']

corrs = np.empty((3,3,3,212))
for i in range(len(wave)):
    for j in range(len(levels)):
        corrs[i,j,:,:] = np.genfromtxt(directorydata + 'LINT_GEOPx-60N-%s_%s.txt' % \
                              (wave[i],levels[j]),skip_header=2,unpack=True,
                              delimiter=',',usecols=[0,1,2])

### Calculate moving linear trends
def calcMovingAve(data,length):
    """
    Calculates moving average for n number of days
    
    Parameters
    ----------
    data : 1d array
        [time series data]
    length : integer
        [n days]
        
    Returns
    -------
    ave : 1d array
        time series of smoothed data from moving average
    
    Usage
    -----
    ave = calcMovingAve(data,years,length)
    """
    print('\n>>> Using calcMovingAverage function!')
    
    ### Calculate moving average for n months (length)
    aven = np.convolve(data, np.ones((length,))/length, mode='valid') 
    print('Completed: *%s DAYS* averages!' % length)
    
    ### Append nans for consistent time
    empty = np.array([np.nan]*(length-1))
    ave = np.append(empty,aven,axis=0)
    
    print('*Completed: Finished calcMovingAve function!\n')    
    return ave

corrsmooth = np.empty(corrs.shape)
for w in range(len(wave)):
    for l in range(len(levels)):
       for p in range(len(qbophase)):
           corrsmooth[w,l,p,:] = calcMovingAve(corrs[w,l,p,:],10)
           
corrpos = corrsmooth[:,:,0,:]
corrnon = corrsmooth[:,:,1,:]
corrneg = corrsmooth[:,:,2,:]

#############################################################################
#############################################################################
#############################################################################
### Plot daily correlations
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
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
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.axhline(0,color='dimgrey',linewidth=1,linestyle='-')

plt.plot(corrneg[0,0,:],color='m',label=r'\textbf{Wave 1}',
         linewidth=3)
plt.plot(corrneg[1,0,:],color='deepskyblue',label=r'\textbf{Wave 2}',
         linewidth=3)
plt.plot(corrneg[2,0,:],color='darkgreen',label=r'\textbf{Wave All}',
         linewidth=3)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.ylim([-1,1])
plt.xlim([30,150])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'QBOExperiments_linearCorrelations_Daily_2.png',
            dpi=300)

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.axhline(0,color='dimgrey',linewidth=1,linestyle='-')

plt.plot(corrneg[0,0,:],color='m',label=r'\textbf{Column}',
         linewidth=3)
plt.plot(corrneg[0,1,:],color='deepskyblue',label=r'\textbf{Troposphere}',
         linewidth=3)
plt.plot(corrneg[0,2,:],color='darkgreen',label=r'\textbf{Stratosphere}',
         linewidth=3)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.ylim([-1,1])
plt.xlim([30,150])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'QBOExperiments_linearCorrelations_Wave1_DailyAtmo_2.png',
            dpi=300)

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.axhline(0,color='dimgrey',linewidth=1,linestyle='-')

plt.plot(corrneg[1,0,:],color='m',label=r'\textbf{Column}',
         linewidth=3)
plt.plot(corrneg[1,1,:],color='deepskyblue',label=r'\textbf{Troposphere}',
         linewidth=3)
plt.plot(corrneg[1,2,:],color='darkgreen',label=r'\textbf{Stratosphere}',
         linewidth=3)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.ylim([-1,1])
plt.xlim([30,150])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'QBOExperiments_linearCorrelations_Wave2_DailyAtmo_2.png',
            dpi=300)

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.axhline(0,color='dimgrey',linewidth=1,linestyle='-')

plt.plot(corrneg[2,0,:],color='m',label=r'\textbf{Column}',
         linewidth=3)
plt.plot(corrneg[2,1,:],color='deepskyblue',label=r'\textbf{Troposphere}',
         linewidth=3)
plt.plot(corrneg[2,2,:],color='darkgreen',label=r'\textbf{Stratosphere}',
         linewidth=3)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.ylim([-1,1])
plt.xlim([30,150])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'QBOExperiments_linearCorrelations_WaveAll_DailyAtmo_2.png',
            dpi=300)

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.axhline(0,color='dimgrey',linewidth=1,linestyle='-')

plt.plot(corrpos[0,0,:],color='m',label=r'\textbf{QBO-W}',
         linewidth=3)
plt.plot(corrneg[0,0,:],color='deepskyblue',label=r'\textbf{QBO-E}',
         linewidth=3)

xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.ylim([-1,1])
plt.xlim([30,150])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'QBOExperiments_linearCorrelations_QBO_Daily_2.png',
            dpi=300)

print('Completed: Script done!')