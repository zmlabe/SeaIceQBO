"""
Manuscript figure for regional locations of the jet

Notes
-----
    Author : Zachary Labe
    Date   : 19 October 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
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
print('\n' '----Plotting Daily Jet for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
N=14 # days
varnames = ['U700']
runnames = [r'HIT',r'FIT',r'FICT']
qbophase = ['pos','non','neg']
qbophaseq = [r'QBO-E',r'QBO-W',r'Difference']
experiments = [r'\textbf{ATLANTIC}',r'\textbf{ATLANTIC}',
               r'\textbf{ATLANTIC}',r'\textbf{PACIFIC}',
               r'\textbf{PACIFIC}',r'\textbf{PACIFIC}']

def readJet(varnames,qbophase,location):
    ### Call functions for variable profile data for polar cap
    for v in range(len(varnames)):
        lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames[v],
                                                    'HIT','surface')
        lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames[v],
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
        ### Concatonate runs
        var_mo = [varhit,varfict]
        
        ### Composite by QBO phase    
        var_mohitpos = var_mo[0][pos_hit,:]
        var_mofictpos = var_mo[1][pos_fict,:]
        var_mohitneg = var_mo[0][neg_hit,:]
        var_mofictneg = var_mo[1][neg_fict,:]
        
        ### Average over longitude
        if location == 'Atlantic':
            lonq = np.append(np.where((lon >=310) & (lon <=360))[0],
                             np.where((lon >=0) & (lon <=20))[0],axis=0)
        if location == 'Pacific':
            lonq = np.where((lon >= 120) & (lon <= 220))[0]
            
        var_mohitposz = np.nanmean(var_mohitpos[:,:,:,lonq],axis=3)
        var_mofictposz = np.nanmean(var_mofictpos[:,:,:,lonq],axis=3)
        var_mohitnegz = np.nanmean(var_mohitneg[:,:,:,lonq],axis=3)
        var_mofictnegz = np.nanmean(var_mofictneg[:,:,:,lonq],axis=3)
    
        ### Save climatologies
        climo = [var_mohitnegz,var_mohitposz,np.zeros(var_mohitposz.shape)]
        
        ### Compute comparisons for months - taken ensemble average 
        ficthitpos = np.nanmean(var_mofictposz - var_mohitposz,axis=0)
        ficthitneg = np.nanmean(var_mofictnegz - var_mohitnegz,axis=0)
        diffruns = [ficthitneg,ficthitpos,ficthitneg-ficthitpos]
            
        ### Smoothing for statistics
        varstat = [var_mohitnegz,var_mohitposz,var_mofictnegz,var_mofictposz]
        stats = []
        for yr in range(len(varstat)):
            varqstat = np.empty((varstat[yr].shape[0],212-(N-1),96))
            for ph in range(var_mohitnegz.shape[0]):
                for s in range(var_mohitnegz.shape[2]):
                    varqstat[ph,:,s] = np.convolve(varstat[yr][ph,:,s], np.ones((N,))/N, mode='valid') 
            empty = np.empty((varstat[yr].shape[0],N-1,varstat[yr].shape[2]))
            empty[:] = np.nan
            varnstat = np.append(np.asarray(varqstat),empty,axis=1)   
            stats.append(varnstat)
        
        ### Calculate significance
        stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(stats[3],stats[1])
        stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(stats[2],stats[0])
    
        pruns = [pvalue_FICTHITneg,pvalue_FICTHITpos,np.zeros(pvalue_FICTHITneg.shape)]
        
        ### Delete files for memory
        del var_mo
        del varhit
        del varfict
        
    return lat,lon,diffruns,climo,pruns

### Use function to read data
lat,lon,diffrunsa,climoa,prunsa = readJet(varnames,qbophase,'Atlantic')
lat,lon,diffrunsp,climop,prunsp = readJet(varnames,qbophase,'Pacific')
    
###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
for i in range(6):
    ax1 = plt.subplot(2,3,i+1)
    
    if i < 3:
        var = diffrunsa[i]
        pvar = prunsa[i].transpose()
        climovar = np.nanmean(climoa[i],axis=0)
    elif i >= 3:
        var = diffrunsp[i-3]
        pvar = prunsp[i-3].transpose()
        climovar = np.nanmean(climop[i-3],axis=0)
    
    ### Begin smoothing of N days
    varq = []
    for s in range(var.shape[1]):
        varqq=np.convolve(var[:,s], np.ones((N,))/N, mode='valid') 
        varq.append(varqq)
        empty = np.empty((96,N-1))
        empty[:] = np.nan
    varn = np.append(np.asarray(varq),empty,axis=1)
    
    ### Begin smoothing of N days
    climovarq = []
    for s in range(climovar.shape[1]):
        climovarqq=np.convolve(climovar[:,s], np.ones((N,))/N, mode='valid') 
        climovarq.append(climovarqq)
        empty = np.empty((96,N-1))
        empty[:] = np.nan
    climovarn = np.append(empty,np.asarray(climovarq),axis=1)
    
    ax1.spines['top'].set_color('dimgrey')
    ax1.spines['right'].set_color('dimgrey')
    ax1.spines['bottom'].set_color('dimgrey')
    ax1.spines['left'].set_color('dimgrey')
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['right'].set_linewidth(2)
    ax1.spines['top'].set_linewidth(2)
    ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')
    ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                    width=2,color='dimgrey')    
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    
    ### Set limits
    limit = np.arange(-4,4.1,0.1)
    barlim = np.arange(-4,4.1,2)
    if varnames[0] == 'U700':
        limit = np.arange(-2,2.1,0.1)
        barlim = np.arange(-2,2.1,1)
    time = np.arange(212)
    if i < 3:
        climlimit = np.arange(-80,81,2.5)
    elif i >= 3:
        climlimit = np.arange(0,81,2.1)
    
    cs = plt.contourf(time,lat,varn,limit,extend='both')
    if i == 0 or i==1 or i==3 or i==4:
        cs1 = plt.contour(time,lat,climovarn,climlimit,linewidths=2,colors='k')
        cs2 = plt.contourf(time,lat,pvar,colors='None',hatches=['////'])     
    cs.set_cmap(cmocean.cm.balance)
    
    if i < 3:
        ax1.text(0.5,1.08,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=15,
                 transform=ax1.transAxes)
    if i == 0 or i == 3:
        ax1.text(-0.27,0.5,r'%s' % experiments[i],color='dimgray',
                     fontsize=15,rotation=90,ha='center',va='center',
                     transform=ax1.transAxes)
    
    xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
    if i==3:
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=5)
        plt.yticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
        plt.minorticks_off()
        plt.xlim([60,120])
        plt.ylim([15,90])
    elif i==4:
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=5)
        plt.yticks([])
        plt.minorticks_off()
        plt.xlim([60,120])
        plt.ylim([15,90])
    elif i==5:
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=5)
        plt.yticks([])
        plt.minorticks_off()
        plt.xlim([60,120])
        plt.ylim([15,90])
    elif i==0:
        plt.xticks([])
        plt.yticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
        plt.minorticks_off()
        plt.xlim([60,120])
        plt.ylim([15,90])
    else:
        plt.xticks([])
        plt.yticks([])
        plt.xlim([60,120])
        plt.ylim([15,90])
    
cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01,labelsize=6)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(wspace=0.1,hspace=0.08)
plt.subplots_adjust(bottom=0.18)

ax1.text(-2.37,1,r'\textbf{Latitude ($^\circ$N)}',color='k',
         fontsize=6,rotation=90,ha='center',va='center',
         transform=ax1.transAxes)

plt.savefig(directoryfigure + 'Daily_Jet_All_%s.png' % (varnames[0]),dpi=300)
print('Completed: Script done!')