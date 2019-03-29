"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data for all variables.

Notes
-----
    Author : Zachary Labe
    Date   : 6 September 2017
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
print('\n' '----Plotting Daily Variables for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
varnames = ['GEOP','TEMP','U','V']
runnames = [r'HIT',r'FIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']

### Call functions for variable profile data for polar cap
for v in range(len(varnames)):
    lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'HIT','profile')
    lat,lon,time,lev,varfit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'FIT','profile')
    lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames[v],
                                                'FICT','profile')
    
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
    
    ### Compute comparisons for months - taken ensemble average
    fithitpos = np.nanmean(var_mofitpos - var_mohitpos,axis=0)
    fithitnon = np.nanmean(var_mofitnon - var_mohitnon,axis=0)
    fithitneg = np.nanmean(var_mofitneg - var_mohitneg,axis=0)
    
    ficthitpos = np.nanmean(var_mofictpos - var_mohitpos,axis=0)
    ficthitnon = np.nanmean(var_mofictnon - var_mohitnon,axis=0)
    ficthitneg = np.nanmean(var_mofictneg - var_mohitneg,axis=0)
    
    diffruns = [fithitpos,fithitnon,fithitneg,ficthitpos,ficthitnon,ficthitneg]
    
    ### Calculate significance for FM
    stat_FITHITpos,pvalue_FITHITpos = UT.calc_indttest(var_mo[1][pos_fit,:],
                                                       var_mo[0][pos_hit,:])
    stat_FITHITnon,pvalue_FITHITnon = UT.calc_indttest(var_mo[1][non_fit,:],
                                                       var_mo[0][non_hit,:])
    stat_FITHITneg,pvalue_FITHITneg = UT.calc_indttest(var_mo[1][neg_fit,:],
                                                       var_mo[0][neg_hit,:])
    
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(var_mo[2][pos_fict,:],
                                                         var_mo[0][pos_hit,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(var_mo[2][non_fict,:],
                                                         var_mo[0][non_hit,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(var_mo[2][neg_fict,:],
                                                         var_mo[0][neg_hit,:])

    pruns = [pvalue_FITHITpos,pvalue_FITHITnon,pvalue_FITHITneg,
                pvalue_FICTHITpos,pvalue_FICTHITnon,pvalue_FICTHITneg]
                                                 
    ############################################################################
    ############################################################################
    ############################################################################
    ##### Plot daily profile with height for selected variable
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'GEOP':
        limit = np.arange(-150,150.1,15)
        barlim = np.arange(-150,151,75)
    elif varnames[v] == 'TEMP':
        limit = np.arange(-3,3.1,0.2)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'U':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'V':
        limit = np.arange(-0.3,0.305,0.01)
        barlim = np.arange(-0.3,0.31,0.15)
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    timeq = np.arange(0,212,1)
    timeqq,levq = np.meshgrid(timeq,lev)
    
    fig = plt.figure()
    for i in range(len(experiments)):
        ax1 = plt.subplot(2,3,i+1)
        
        var = diffruns[i]
        pvar = pruns[i]
        
        if MASK == True:
            pvar2 = pvar.copy()
            pvar2[np.isnan(pvar2)]=0.0
            var = var*pvar2
        
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
        
        cs = plt.contourf(timeq,lev,var.transpose(),limit,extend='both')
        
        if MASK == False:                  
            plt.contourf(timeqq,levq,pvar.transpose(),colors='None',
                         hatches=['////'])                  
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
        plt.xlim([30,210])
        plt.ylim([1000,10])
        
        if varnames[v] == 'U':
            cmap = ncm.cmap('NCV_blu_red')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'TEMP':
            cmap = ncm.cmap('NCV_blu_red')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'GEOP':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'V':
            cmap = ncm.cmap('temp_diff_18lev')            
            cs.set_cmap(cmap) 

        ### Add experiment text to subplot
        if i < 3:
            qbophaseq = [r'QBO-W',r'QBO-N',r'QBO-E']
            ax1.text(0.5,1.1,r'\textbf{%s}' % qbophaseq[i],
                     ha='center',va='center',color='dimgray',fontsize=13,
                     transform=ax1.transAxes)
        if i == 0 or i == 3:
            ax1.text(-0.4,0.5,r'%s' % experiments[i],color='k',
                         fontsize=20,rotation=90,ha='center',va='center',
                         transform=ax1.transAxes)
    
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    if varnames[v] == 'U':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
    elif varnames[v] == 'TEMP':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')
    elif varnames[v] == 'GEOP':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')
    elif varnames[v] == 'V':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.28)
    plt.subplots_adjust(bottom=0.21)
    
    if MASK == True:
        plt.savefig(directoryfigure + 'allExperiments_%s_MASK_daily.png' % varnames[v],
                    dpi=300)
    else:
        plt.savefig(directoryfigure + 'allExperiments_%s_daily.png' % varnames[v],
            dpi=300)
    print('Completed: Script done!')                                