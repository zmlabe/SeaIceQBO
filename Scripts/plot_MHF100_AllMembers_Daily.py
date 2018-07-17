"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for daily meridional heat flux at 100 hPa.

Notes
-----
    Author : Zachary Labe
    Date   : 16 July 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
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
print('\n' '----Plotting Daily MHF100 for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
#MASK = False
#varnames = ['MHF100']
#runnames = [r'HIT',r'FIT',r'FICT']
#qbophase = ['pos','non','neg']
#phases = ['QBO-W','QBO-N','QBO-W']
#experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
#               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
#
#### Call functions for MHF data at 100 hPa
#lat,lon,time,lev,varhit = DO.readMeanExperiAll('MHF100','HIT','profile')
#lat,lon,time,lev,varfit = DO.readMeanExperiAll('MHF100','FIT','profile')
#lat,lon,time,lev,varfict = DO.readMeanExperiAll('MHF100','FICT','profile')
#
#### Create 2d array of latitude and longitude
#lon2,lat2 = np.meshgrid(lon,lat)
#
#### Read in QBO phases 
#filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitno2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
#non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
#neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
#
#filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
#filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
#filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
#non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+100)
#neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)
#
#filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
#non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
#neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)    
#### Concatonate runs
#var_mo = [varhit,varfit,varfict]
#
#### Composite by QBO phase    
#var_mofitpos = var_mo[1][pos_fit,:]
#var_mohitpos = var_mo[0][pos_hit,:]
#var_mofictpos = var_mo[2][pos_fict,:]
#
#var_mofitnon = var_mo[1][non_fit,:]
#var_mohitnon = var_mo[0][non_hit,:]
#var_mofictnon = var_mo[2][non_fict,:]
#
#var_mofitneg = var_mo[1][neg_fit,:]
#var_mohitneg = var_mo[0][neg_hit,:]
#var_mofictneg = var_mo[2][neg_fict,:]
#
#### Compute comparisons for months - taken ensemble average
#fithitpos = np.nanmean(var_mofitpos - var_mohitpos,axis=0).squeeze()
#fithitnon = np.nanmean(var_mofitnon - var_mohitnon,axis=0).squeeze()
#fithitneg = np.nanmean(var_mofitneg - var_mohitneg,axis=0).squeeze()
#
#ficthitpos = np.nanmean(var_mofictpos - var_mohitpos,axis=0).squeeze()
#ficthitnon = np.nanmean(var_mofictnon - var_mohitnon,axis=0).squeeze()
#ficthitneg = np.nanmean(var_mofictneg - var_mohitneg,axis=0).squeeze()
#
#diffruns = np.array([fithitpos,fithitnon,fithitneg,ficthitpos,ficthitnon,ficthitneg])
#
#### Compute zonal average
#diffrunsz = np.nanmean(diffruns,axis=2)
#
##############################################################################
##############################################################################
##############################################################################
#### Plot daily MHF anomalies
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#
#def adjust_spines(ax, spines):
#    for loc, spine in ax.spines.items():
#        if loc in spines:
#            spine.set_position(('outward', 2))
#        else:
#            spine.set_color('none')  
#    if 'left' in spines:
#        ax.yaxis.set_ticks_position('left')
#    else:
#        ax.yaxis.set_ticks([])
#
#    if 'bottom' in spines:
#        ax.xaxis.set_ticks_position('bottom')
#    else:
#        ax.xaxis.set_ticks([]) 
#        
#fig = plt.figure()
#ax = plt.subplot(2,1,1) 
#       
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('dimgrey')
#ax.spines['bottom'].set_color('dimgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)
#
#color=iter(cmocean.cm.phase(np.linspace(0.1,0.9,3)))
#for i in range(3):
#    c=next(color)
#    plt.plot(diffrunsz[i],linewidth=2,color=c,alpha=1,
#             label = r'\textbf{%s}' % phases[i],linestyle='-')
#    
#xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
#plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
#plt.yticks(np.arange(-10,11,2),map(str,np.arange(-10,11,2)),fontsize=8)
#plt.ylim([-4,7])
#plt.xlim([30,210])
#
#plt.ylabel(r'\textbf{K m/s}',color='k',fontsize=11)
#
#ax = plt.subplot(2,1,2) 
#       
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('dimgrey')
#ax.spines['bottom'].set_color('dimgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)
#
#color=iter(cmocean.cm.phase(np.linspace(0.1,0.9,3)))
#for i in range(3):
#    c=next(color)
#    plt.plot(diffrunsz[i+3],linewidth=2,color=c,alpha=1,
#             label = r'\textbf{%s}' % phases[i],linestyle='-')
#    
#plt.legend(shadow=False,fontsize=9,loc='upper center',
#           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,-0.13))
#    
#xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
#plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
#plt.yticks(np.arange(-10,11,2),map(str,np.arange(-10,11,2)),fontsize=8)
#plt.ylim([-4,7])
#plt.xlim([30,210])
#
#plt.ylabel(r'\textbf{K m/s}',color='k',fontsize=11)
#    
#plt.savefig(directoryfigure + '/QBO_Daily_2/QBOExperiments_MHF100_All.png',
#            dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot -QBO for FICT
fict = var_mo[2].squeeze()

### Calculate means
fictclimo = np.nanmean(fict,axis=0)
climoz = np.nanmean(fictclimo,axis=1)

### Negative and Positive QBO
var_mofictpos = np.nanmean(var_mo[2][pos_fict,:,:].squeeze(),axis=0)
var_mofictneg = np.nanmean(var_mo[2][neg_fict,:,:].squeeze(),axis=0)

### Calculate zonal mean
posz = np.nanmean(var_mofictpos,axis=1)
negz = np.nanmean(var_mofictneg,axis=1)

### Calculate anomalies
posanom = posz - climoz
neganom = negz - climoz

anomsfict = [posanom,neganom]
phasesq = ['QBO-W','QBO-E']

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

color=iter(cmocean.cm.balance(np.linspace(0.2,0.8,2)))
for i in range(2):
    c=next(color)
    plt.plot(anomsfict[i],linewidth=3,color=c,alpha=1,
             label = r'\textbf{%s}' % phasesq[i],linestyle='-')
    
xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),fontsize=8)
plt.ylim([-3,4])
plt.xlim([30,210])

plt.legend(shadow=False,fontsize=9,loc='upper center',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.1))

plt.ylabel(r'\textbf{K m/s}',color='k',fontsize=11)
    
plt.savefig(directoryfigure + '/QBO_Daily_2/QBOExperiments_MHF100_FICT.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot -QBO for FIT
fit = var_mo[1].squeeze()

### Calculate means
fitclimo = np.nanmean(fit,axis=0)
climoz = np.nanmean(fitclimo,axis=1)

### Negative and Positive QBO
var_mofitpos = np.nanmean(var_mo[1][pos_fit,:,:].squeeze(),axis=0)
var_mofitneg = np.nanmean(var_mo[1][neg_fit,:,:].squeeze(),axis=0)

### Calculate zonal mean
posz = np.nanmean(var_mofitpos,axis=1)
negz = np.nanmean(var_mofitneg,axis=1)

### Calculate anomalies
posanom = posz - climoz
neganom = negz - climoz

anomsfit = [posanom,neganom]
phasesq = ['QBO-W','QBO-E']

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

color=iter(cmocean.cm.balance(np.linspace(0.2,0.8,2)))
for i in range(2):
    c=next(color)
    plt.plot(anomsfit[i],linewidth=3,color=c,alpha=1,
             label = r'\textbf{%s}' % phasesq[i],linestyle='-')
    
xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),fontsize=8)
plt.ylim([-3,4])
plt.xlim([30,210])

plt.legend(shadow=False,fontsize=9,loc='upper center',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.1))

plt.ylabel(r'\textbf{K m/s}',color='k',fontsize=11)
    
plt.savefig(directoryfigure + '/QBO_Daily_2/QBOExperiments_MHF100_FIT.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################
### Plot -QBO for hit
hit = var_mo[0].squeeze()

### Calculate means
hitclimo = np.nanmean(hit,axis=0)
climoz = np.nanmean(hitclimo,axis=1)

### Negative and Positive QBO
var_mohitpos = np.nanmean(var_mo[0][pos_hit,:,:].squeeze(),axis=0)
var_mohitneg = np.nanmean(var_mo[0][neg_hit,:,:].squeeze(),axis=0)

### Calculate zonal mean
posz = np.nanmean(var_mohitpos,axis=1)
negz = np.nanmean(var_mohitneg,axis=1)

### Calculate anomalies
posanom = posz - climoz
neganom = negz - climoz

anomshit = [posanom,neganom]
phasesq = ['QBO-W','QBO-E']

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

color=iter(cmocean.cm.balance(np.linspace(0.2,0.8,2)))
for i in range(2):
    c=next(color)
    plt.plot(anomshit[i],linewidth=3,color=c,alpha=1,
             label = r'\textbf{%s}' % phasesq[i],linestyle='-')
    
xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
plt.xticks(np.arange(0,212,30),xlabels,fontsize=8)
plt.yticks(np.arange(-10,11,1),map(str,np.arange(-10,11,1)),fontsize=8)
plt.ylim([-3,4])
plt.xlim([30,210])

plt.legend(shadow=False,fontsize=9,loc='upper center',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0.5,1.1))

plt.ylabel(r'\textbf{K m/s}',color='k',fontsize=11)
    
plt.savefig(directoryfigure + '/QBO_Daily_2/QBOExperiments_MHF100_HIT.png',
            dpi=300)




print('Completed: Script done!')