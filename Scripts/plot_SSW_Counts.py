"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
assess stratospheric sudden warming event counts.

Notes
-----
    Author : Zachary Labe
    Date   : 11 July 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput_SSW as SO
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
print('\n' '----Plotting SSW Counts for Simulations - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
qbophase = ['pos','non','neg']

### Call functions for SSW counts
sswhit = SO.readDailySSW('count','HIT')
sswfit = SO.readDailySSW('count','FIT')
sswfict = SO.readDailySSW('count','FICT')
sswcon = SO.readDailySSW_CTLQ('count','CTLQ')

hitsum = int(np.nansum(sswhit))
fitsum = int(np.nansum(sswfit))
fictsum = int(np.nansum(sswfict))
consum = int(np.nansum(sswcon))

print('\n----------------------------------------------------')
print('----------------------------------------------------')
print('---------------------SSW Counts---------------------')
print('\n HIT  --> %s' % hitsum)
print(' FIT  --> %s' % fitsum)
print(' FICT --> %s' % fictsum)
print(' CTLQ --> %s' % consum)

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

filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+100)
neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)

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

filenamefctlqp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
filenamefctlqno = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[1]
filenamefctlqn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
pos_ctlq = np.genfromtxt(filenamefctlqp,unpack=True,usecols=[0],dtype='int')
non_ctlq = np.genfromtxt(filenamefctlqno,unpack=True,usecols=[0],dtype='int')
neg_ctlq = np.genfromtxt(filenamefctlqn,unpack=True,usecols=[0],dtype='int')

### Slice by QBO phase
hitp = sswhit[pos_hit]
hitno = sswhit[non_hit]
hitn = sswhit[neg_hit]

fitp = sswfit[pos_fit]
fitno = sswfit[non_fit]
fitn = sswfit[neg_fit]

fictp = sswfict[pos_fict]
fictno = sswfict[non_fict]
fictn = sswfict[neg_fict]

conp = sswcon[pos_ctlq]
conno = sswcon[non_ctlq]
conn = sswcon[neg_ctlq]

### List simulations
hitq = [np.nansum(hitp),np.nansum(hitno),np.nansum(hitn)]
fitq = [np.nansum(fitp),np.nansum(fitno),np.nansum(fitn)]
fictq = [np.nansum(fictp),np.nansum(fictno),np.nansum(fictn)]
conq = [np.nansum(conp),np.nansum(conno),np.nansum(conn)]
vals = np.array([hitq,fitq,fictq,conq])

### Print SSW data
print('\n----------------------------------------------------')
print('----------------------------------------------------')
print('---------------------SSW WQBO Counts----------------')
print('\n HIT  --> %s' % int(np.nansum(hitp)))
print(' FIT  --> %s' % int(np.nansum(fitp)))
print(' FICT --> %s' % int(np.nansum(fictp)))
print(' CTLQ --> %s' % int(np.nansum(conp)))

print('\n----------------------------------------------------')
print('----------------------------------------------------')
print('---------------------SSW NQBO Counts----------------')
print('\n HIT  --> %s' % int(np.nansum(hitno)))
print(' FIT  --> %s' % int(np.nansum(fitno)))
print(' FICT --> %s' % int(np.nansum(fictno)))
print(' CTLQ --> %s' % int(np.nansum(conno)))

print('\n----------------------------------------------------')
print('----------------------------------------------------')
print('---------------------SSW EQBO Counts----------------')
print('\n HIT  --> %s' % int(np.nansum(hitn)))
print(' FIT  --> %s' % int(np.nansum(fitn)))
print(' FICT --> %s' % int(np.nansum(fictn)))
print(' CTLQ --> %s' % int(np.nansum(conn)))

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
        
### Plot bar graph
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(0)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')  
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=True)      
        
cc = ['crimson','darkgrey','deepskyblue']
labels = [r'\textbf{QBO-W}',r'\textbf{QBO-N}',r'\textbf{QBO-E}']
experiments = [r'\textbf{HIT}',r'\textbf{FIT}',r'\textbf{FICT}',
               r'\textbf{CTLQ}']
bins = np.arange(4)

plt.bar(bins+0.00, vals[:,0], color='crimson',width=0.25,
        label=labels[0])
plt.bar(bins+0.25, vals[:,1], color='darkgrey',width=0.25,
        label=labels[1])
plt.bar(bins+0.50, vals[:,2], color='deepskyblue',width=0.25,
        label=labels[2])

plt.yticks(np.arange(0,26,2),list(map(str,np.arange(0,26,2))),
           fontsize=14)
plt.ylim([0,10])

plt.xticks(np.arange(.25,4.25,1),experiments,fontsize=20)

plt.ylabel(r'\textbf{Number of SSW}',color='k',fontsize=20)  

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

plt.savefig(directoryfigure + 'ssw_counts.png',dpi=300)

print('\nCompleted: Script done!')