"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data for the jet averaged over longitude.

Notes
-----
    Author : Zachary Labe
    Date   : 3 August 2018
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
print('\n' '----Plotting Daily Jet for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
N=21 # days
varnames = ['U200']
runnames = [r'HIT',r'FIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']

### Call functions for variable profile data for polar cap
for v in range(len(varnames)):
    lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'HIT','surface')
    lat,lon,time,lev,varfit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'FIT','surface')
    lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames[v],
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
    ### Concatonate runs
    var_mo = [varhit,varfit,varfict]
    
    ### Composite by QBO phase    
    var_mofitpos = var_mo[1][pos_fit,:]
    var_mohitpos = var_mo[0][pos_hit,:]
    var_mofictpos = var_mo[2][pos_fict,:]
    
    var_mofitnon = var_mo[1][non_fit,:]
    var_mohitnon = var_mo[0][non_hit,:]
    var_mofictnon = var_mo[2][non_fict,:]
    
    var_mofitneg = var_mo[1][neg_fit,:]
    var_mohitneg = var_mo[0][neg_hit,:]
    var_mofictneg = var_mo[2][neg_fict,:]
    
    ### Average over longitude
    region = 'Regional'
    if region == 'Regional':
        location = 'Pacific'
        if location == 'Atlantic':
            lonq = np.append(np.where((lon >=310) & (lon <=360))[0],
                             np.where((lon >=0) & (lon <=20))[0],axis=0)
        if location == 'Pacific':
            lonq = np.where((lon >= 120) & (lon <= 220))[0]
        var_mofitposz = np.nanmean(var_mofitpos[:,:,:,lonq],axis=3)
        var_mohitposz = np.nanmean(var_mohitpos[:,:,:,lonq],axis=3)
        var_mofictposz = np.nanmean(var_mofictpos[:,:,:,lonq],axis=3)
        
        var_mofitnonz = np.nanmean(var_mofitnon[:,:,:,lonq],axis=3)
        var_mohitnonz = np.nanmean(var_mohitnon[:,:,:,lonq],axis=3)
        var_mofictnonz = np.nanmean(var_mofictnon[:,:,:,lonq],axis=3)
        
        var_mofitnegz = np.nanmean(var_mofitneg[:,:,:,lonq],axis=3)
        var_mohitnegz = np.nanmean(var_mohitneg[:,:,:,lonq],axis=3)
        var_mofictnegz = np.nanmean(var_mofictneg[:,:,:,lonq],axis=3)
        
    elif region == 'Global':
        location = 'Global'
        var_mofitposz = np.nanmean(var_mofitpos,axis=3)
        var_mohitposz = np.nanmean(var_mohitpos,axis=3)
        var_mofictposz = np.nanmean(var_mofictpos,axis=3)
        
        var_mofitnonz = np.nanmean(var_mofitnon,axis=3)
        var_mohitnonz = np.nanmean(var_mohitnon,axis=3)
        var_mofictnonz = np.nanmean(var_mofictnon,axis=3)
        
        var_mofitnegz = np.nanmean(var_mofitneg,axis=3)
        var_mohitnegz = np.nanmean(var_mohitneg,axis=3)
        var_mofictnegz = np.nanmean(var_mofictneg,axis=3)
        
    climo = [var_mohitnegz,var_mohitposz,var_mofictnegz,
             var_mohitnegz,var_mohitposz,var_mofictnegz]
    
    ### Compute comparisons for months - taken ensemble average
    fithitpos = np.nanmean(var_mofitposz - var_mohitposz,axis=0)
    fithitnon = np.nanmean(var_mofitnonz - var_mohitnonz,axis=0)
    fithitneg = np.nanmean(var_mofitnegz - var_mohitnegz,axis=0)
    
    ficthitpos = np.nanmean(var_mofictposz - var_mohitposz,axis=0)
    ficthitnon = np.nanmean(var_mofictnonz - var_mohitnonz,axis=0)
    ficthitneg = np.nanmean(var_mofictnegz - var_mohitnegz,axis=0)
    
    diffruns = [fithitneg,fithitpos,fithitneg-fithitpos,ficthitneg,ficthitpos,ficthitneg-ficthitpos]
        
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
    
###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
for i in range(3,len(diffruns)):
    ax1 = plt.subplot(1,3,i-2)
    
    var = diffruns[i]
    var2 = np.nanmean(var_mofictnegz,axis=0)
    pvar = pruns[i-3].transpose()
    climovar = np.nanmean(climo[i],axis=0)
    
    ### Begin smoothing of N days
    varq = []
    for s in range(var.shape[1]):
        varqq=np.convolve(var[:,s], np.ones((N,))/N, mode='valid') 
        varq.append(varqq)
        empty = np.empty((96,N-1))
        empty[:] = np.nan
    varn = np.append(np.asarray(varq),empty,axis=1)
    
    ### Begin smoothing of N days
    varq2 = []
    for s in range(var2.shape[1]):
        varqq2=np.convolve(var2[:,s], np.ones((N,))/N, mode='valid') 
        varq2.append(varqq2)
        empty = np.empty((96,N-1))
        empty[:] = np.nan
    varn2 = np.append(np.asarray(varq2),empty,axis=1)
    
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
    time = np.arange(212)
    if location == 'Pacific':
        climlimit = np.arange(-80,81,10)
    elif location == 'Atlantic':
        climlimit = np.arange(-80,81,8)
    
    cs = plt.contourf(time,lat,varn,limit,extend='both')
    if i == 3 or i==4:
        cs1 = plt.contour(time,lat,climovarn,climlimit,linewidths=2,colors='k')
#        cs1 = plt.contour(time,lat,varn2,np.arange(-80,81,2),linewidths=1,colors='dimgrey',linestyles='--',
#                          dashes=(1,0.3))
        cs2 = plt.contourf(time,lat,pvar,colors='None',hatches=['////'])     
    cs.set_cmap(cmocean.cm.balance)
    
    if i >= 3:
        qbophaseq = [r'QBO-E',r'QBO-W',r'Difference']
        ax1.text(0.5,1.05,r'\textbf{%s}' % qbophaseq[i-3],
                 ha='center',va='center',color='dimgray',fontsize=15,
                 transform=ax1.transAxes)
    
    xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
    plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
    plt.yticks(np.arange(-90,91,15),map(str,np.arange(-90,91,15)),fontsize=6)
    plt.xlim([60,120])
    plt.ylim([15,90])

    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(bottom=0.21)
    
    plt.savefig(directoryfigure + 'Daily_Jet_%s.png' % location,dpi=300)
    print('Completed: Script done!')