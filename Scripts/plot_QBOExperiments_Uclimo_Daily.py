"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
plots daily U30 data over the polar cap

Notes
-----
    Author : Zachary Labe
    Date   : 31 July 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
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
print('\n' '----Plotting Daily U30 for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
varnames = ['U']
runnames = [r'HIT',r'FIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']

#### Call functions for variable profile data for polar cap
#for v in range(len(varnames)):
#    lat,lon,time,lev,varhitq = DO.readMeanExperiAll('%s' % varnames[v],
#                                                'HIT','profile')
#    lat,lon,time,lev,varfitq = DO.readMeanExperiAll('%s' % varnames[v],
#                                                'FIT','profile')
#    lat,lon,time,lev,varfictq = DO.readMeanExperiAll('%s' % varnames[v],
#                                                'FICT','profile')
#    
#    ### Create 2d array of latitude and longitude
#    lon2,lat2 = np.meshgrid(lon,lat)
#    
#    ### Read in QBO phases 
#    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#    filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#    filenamehitno2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
#    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
#    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
#
#    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#    filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
#    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#    filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
#    filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
#    filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
#    pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
#    non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+100)
#    neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)
#    
#    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
#    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
#    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100) 
#    
#    ### Slice Level
#    levq = np.where(lev == 30)[0]
#    varhit = varhitq[:,:,levq]
#    varfit = varfitq[:,:,levq]
#    varfict = varfictq[:,:,levq]
#    
#    ### Concatonate runs
#    var_mo = [varhit,varfit,varfict]
#    
#    ### Composite by QBO phase    
#    var_mofitpos = var_mo[1][pos_fit,:]
#    var_mohitpos = var_mo[0][pos_hit,:]
#    var_mofictpos = var_mo[2][pos_fict,:]
#    
#    var_mofitnon = var_mo[1][non_fit,:]
#    var_mohitnon = var_mo[0][non_hit,:]
#    var_mofictnon = var_mo[2][non_fict,:]
#    
#    var_mofitneg = var_mo[1][neg_fit,:]
#    var_mohitneg = var_mo[0][neg_hit,:]
#    var_mofictneg = var_mo[2][neg_fict,:]
#    
#    ### Calculate ensemble means
#    var_mofitposm = np.nanmean(var_mofitpos,axis=0).squeeze()
#    var_mohitposm = np.nanmean(var_mohitpos,axis=0).squeeze()
#    var_mofictposm = np.nanmean(var_mofictpos,axis=0).squeeze()
#    
#    var_mofitnonm = np.nanmean(var_mofitnon,axis=0).squeeze()
#    var_mohitnonm = np.nanmean(var_mohitnon,axis=0).squeeze()
#    var_mofictnonm = np.nanmean(var_mofictnon,axis=0).squeeze()
#    
#    var_mofitnegm = np.nanmean(var_mofitneg,axis=0).squeeze()
#    var_mohitnegm = np.nanmean(var_mohitneg,axis=0).squeeze()
#    var_mofictnegm = np.nanmean(var_mofictneg,axis=0).squeeze()
#    
#    ### Calculate percentile
#    var_mofitpos95 = np.percentile(var_mofitpos,95,axis=0).squeeze()
#    var_mohitposm95 = np.percentile(var_mohitpos,95,axis=0).squeeze()
#    var_mofictpos95 = np.percentile(var_mofictpos,95,axis=0).squeeze()
#    
#    var_mofitnon95 = np.percentile(var_mofitnon,95,axis=0).squeeze()
#    var_mohitnon95 = np.percentile(var_mohitnon,95,axis=0).squeeze()
#    var_mofictnon95 = np.percentile(var_mofictnon,95,axis=0).squeeze()
#    
#    var_mofitneg95 = np.percentile(var_mofitneg,95,axis=0).squeeze()
#    var_mohitneg95 = np.percentile(var_mohitneg,95,axis=0).squeeze()
#    var_mofictneg95 = np.percentile(var_mofictneg,95,axis=0).squeeze()
#    
#    ### Calculate percentile
#    var_mofitpos5 = np.percentile(var_mofitpos,5,axis=0).squeeze()
#    var_mohitpos5 = np.percentile(var_mohitpos,5,axis=0).squeeze()
#    var_mofictpos5 = np.percentile(var_mofictpos,5,axis=0).squeeze()
#    
#    var_mofitnon5 = np.percentile(var_mofitnon,5,axis=0).squeeze()
#    var_mohitnon5 = np.percentile(var_mohitnon,5,axis=0).squeeze()
#    var_mofictnon5 = np.percentile(var_mofictnon,5,axis=0).squeeze()
#    
#    var_mofitneg5 = np.percentile(var_mofitneg,5,axis=0).squeeze()
#    var_mohitneg5 = np.percentile(var_mohitneg,5,axis=0).squeeze()
#    var_mofictneg5 = np.percentile(var_mofictneg,5,axis=0).squeeze()
#    
#    ### Compute comparisons for months - taken ensemble average
#    fithitpos = var_mofitposm - var_mohitposm
#    fithitnon = var_mofitnonm - var_mohitnonm
#    fithitneg = var_mofitnegm - var_mohitnegm
#    
#    ficthitpos = var_mofictposm - var_mohitposm
#    ficthitnon = var_mofictnonm - var_mohitnonm
#    ficthitneg = var_mofictnegm - var_mohitnegm
#    
#    diffruns = [fithitpos,fithitnon,fithitneg,ficthitpos,ficthitnon,ficthitneg]
    
###############################################################################
###############################################################################
###############################################################################    
### Plot zero wind line
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
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.plot(var_mofictnegm,linewidth=3,
         color=cmocean.cm.thermal(0.3),label=r'\textbf{FICT}')
plt.plot(var_mofictneg5,linewidth=0.7,
         color=cmocean.cm.thermal(0.3))
plt.plot(var_mofictneg95,linewidth=0.7,
         color=cmocean.cm.thermal(0.3))
plt.fill_between(np.arange(0,212,1),var_mofictneg95,var_mofictneg5,
                 color=cmocean.cm.thermal(0.3),alpha=0.3)

plt.plot(var_mohitnegm,linewidth=3,
         color=cmocean.cm.thermal(0.7),label=r'\textbf{HIT}')
plt.plot(var_mohitneg5,linewidth=0.7,
         color=cmocean.cm.thermal(0.7))
plt.plot(var_mohitneg95,linewidth=0.7,
         color=cmocean.cm.thermal(0.7))
plt.fill_between(np.arange(0,212,1),var_mohitneg95,var_mohitneg5,
                 color=cmocean.cm.thermal(0.7),alpha=0.3)

#plt.plot(var_mofictposm,linewidth=3,
#         color=cmocean.cm.thermal(0.7))
#plt.plot(var_mofictpos5,linewidth=0.7,
#         color=cmocean.cm.thermal(0.7))
#plt.plot(var_mofictpos95,linewidth=0.7,
#         color=cmocean.cm.thermal(0.7))
#plt.fill_between(np.arange(0,212,1),var_mofictpos95,var_mofictpos5,
#                 color=cmocean.cm.thermal(0.7),alpha=0.3)

l = plt.legend(shadow=False,fontsize=8,loc='upper left',
           fancybox=True,frameon=False,ncol=2,bbox_to_anchor=(0, 1.0),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for text in l.get_texts():
    text.set_color('k')  

plt.yticks(np.arange(-30,41,5),list(map(str,np.arange(-30,41,5))),
           fontsize=9)
xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr']
plt.xticks(np.arange(0,212,30),xlabels,fontsize=9)
plt.ylabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)  
plt.ylim([-10,35])
plt.xlim([30,210])

plt.tight_layout()

plt.savefig(directoryfigure + 'U30_vari.png',dpi=300)