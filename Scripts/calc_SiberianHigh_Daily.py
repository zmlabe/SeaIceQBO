"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Composites are organized by QBO phase (positive, neutral, negative). 
Script calculates Siberian High.

Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2018
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
import scipy.stats as sts

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
print('\n' '----Calculating Siberian High - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
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
#    elif period == 'N':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = runs[i][:,-2,:,:]
#    elif period == 'ND':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
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
    if varnames == 'SLP':
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
    diffruns = [ficthitpos.squeeze(),ficthitnon.squeeze(),ficthitneg.squeeze()]
    
    return diffruns,lat,lon,lev

slp,lat,lon,lev = readVariables('SLP',period,'Siberia')

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
        
num_bins = np.arange(1010,1060.1,1)

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

ne,binse,patchese = plt.hist(slp[2].ravel(),num_bins,alpha=0.7,
                          facecolor='deepskyblue',edgecolor='deepskyblue',
                          label=r'\textbf{QBO-E}',normed=True)
nw,binsw,patchesw = plt.hist(slp[0].ravel(),num_bins,alpha=0.5,
                          facecolor='crimson',edgecolor='crimson',
                          label=r'\textbf{QBO-W}',normed=True)

### Fit a distribution
mneg,seng = sts.norm.fit(slp[2].ravel())
mpos,spos = sts.norm.fit(slp[0].ravel())

pdfneg = sts.norm.pdf(num_bins,mneg,seng)
pdfpos = sts.norm.pdf(num_bins,mpos,spos)

plt.plot(num_bins,pdfneg,color='darkblue',linewidth=3)
plt.plot(num_bins,pdfpos,color='darkred',linewidth=3)

plt.yticks(np.arange(0,0.25,0.02),list(map(str,np.arange(0,0.25,0.02))),
           fontsize=10)
plt.xticks(np.arange(1010,1061,5),list(map(str,np.arange(1010,1061,5))),
           fontsize=10) 
plt.xlim([1010,1055])
plt.ylim([0,0.08])

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

plt.ylabel(r'\textbf{Density}',color='dimgrey',fontsize=12)  
plt.xlabel(r'\textbf{SLP [hPa]}',color='dimgrey',fontsize=12)

plt.legend(shadow=False,fontsize=12,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)

plt.tight_layout()

### Calculate statistical significance of distributions
t,pvalue = sts.ks_2samp(slp[2].ravel(),slp[0].ravel())
print('\n \n' + 'p_value=' + str(pvalue))

plt.savefig(directoryfigure + 'SiberianHigh_Daily.png',dpi=300)