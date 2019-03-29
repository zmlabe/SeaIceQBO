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
import read_DailyOutput_AllRegional as DOR
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_Daily_2/Smooth/'
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
N = 7 # days
MASK = True
varnames = ['GEOP','TEMP','U','V']
runnames = [r'CIT',r'FSUB',r'FPOL']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',r'\textbf{FSUB-CIT}',
               r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}']

### Call functions for variable profile data for polar cap
for v in range(len(varnames)):
    lat,lon,time,lev,varcit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'CIT','profile')
    lat,lon,time,lev,varfsub = DOR.readMeanExperiAllRegional('%s' % varnames[v],
                                                'FSUB','profile')
    lat,lon,time,lev,varfpol = DOR.readMeanExperiAllRegional('%s' % varnames[v],
                                                'FPOL','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamecitp = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    filenamecitp2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    pos_cit = np.append(np.genfromtxt(filenamecitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitp2,unpack=True,usecols=[0],dtype='int')+100)
    non_cit = np.append(np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_cit = np.append(np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitn2,unpack=True,usecols=[0],dtype='int')+100)    

    filenamefsubp = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    filenamefsubp2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    pos_fsub = np.append(np.genfromtxt(filenamefsubp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fsub = np.append(np.genfromtxt(filenamefsubno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fsub = np.append(np.genfromtxt(filenamefsubn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubn2,unpack=True,usecols=[0],dtype='int')+100)
    
    filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    filenamefpolp2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    pos_fpol = np.append(np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fpol = np.append(np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fpol = np.append(np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpoln2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    var_mo = [varcit,varfsub,varfpol]
    
    ### Composite by QBO phase    
    var_mofsubpos = var_mo[1][pos_fsub,:,:]
    var_mocitpos = var_mo[0][pos_cit,:,:]
    var_mofpolpos = var_mo[2][pos_fpol,:,:]
    
    var_mofsubnon = var_mo[1][non_fsub,:,:]
    var_mocitnon = var_mo[0][non_cit,:,:]
    var_mofpolnon = var_mo[2][non_fpol,:,:]
    
    var_mofsubneg = var_mo[1][neg_fsub,:,:]
    var_mocitneg = var_mo[0][neg_cit,:,:]
    var_mofpolneg = var_mo[2][neg_fpol,:,:]
    varphasepos = [var_mofsubpos,var_mocitpos,var_mofpolpos]
    varphasenon = [var_mofsubnon,var_mocitnon,var_mofpolnon]
    varphaseneg = [var_mofsubneg,var_mocitneg,var_mofpolneg]
    
    poss = []
    for s in range(3):
        for y in range(var_mocitpos.shape[0]):
            for l in range(var_mocitpos.shape[2]):
                possq=np.convolve(varphasepos[s][y,:,l],np.ones((N,))/N,mode='valid') 
                poss.append(possq)
    poss = np.reshape(np.asarray(poss),(3,68,17,varcit.shape[1]-(N-1)))    
    
    nons = []
    for s in range(3):
        for y in range(var_mocitnon.shape[0]):
            for l in range(var_mocitpos.shape[2]):
                nonsq=np.convolve(varphasenon[s][y,:,l],np.ones((N,))/N,mode='valid') 
                nons.append(nonsq)
    nons = np.reshape(np.asarray(nons),(3,68,17,varcit.shape[1]-(N-1)))   
        
    negs = []
    for s in range(3):
        for y in range(var_mocitneg.shape[0]):
            for l in range(var_mocitpos.shape[2]):
                negsq=np.convolve(varphaseneg[s][y,:,l],np.ones((N,))/N,mode='valid') 
                negs.append(negsq)
    negs = np.reshape(np.asarray(negs),(3,64,17,varcit.shape[1]-(N-1)))   
    
    ### Compute comparisons for months - taken ensemble average
    fsubcitpos = np.nanmean(var_mofsubpos - var_mocitpos,axis=0)
    fsubcitnon = np.nanmean(var_mofsubnon - var_mocitnon,axis=0)
    fsubcitneg = np.nanmean(var_mofsubneg - var_mocitneg,axis=0)
    
    fpolcitpos = np.nanmean(var_mofpolpos - var_mocitpos,axis=0)
    fpolcitnon = np.nanmean(var_mofpolnon - var_mocitnon,axis=0)
    fpolcitneg = np.nanmean(var_mofpolneg - var_mocitneg,axis=0)
    diffruns = [fsubcitpos,fsubcitnon,fsubcitneg,
                   fpolcitpos,fpolcitnon,fpolcitneg]
    
    ### Calculate significance for time
    stat_fsubcitpos,pvalue_fsubcitpos = UT.calc_indttest(poss[0],poss[1])
    stat_fsubcitnon,pvalue_fsubcitnon = UT.calc_indttest(nons[0],nons[1])
    stat_fsubcitneg,pvalue_fsubcitneg = UT.calc_indttest(negs[0],negs[1])
    
    stat_fpolcitpos,pvalue_fpolcitpos = UT.calc_indttest(poss[2],poss[1])
    stat_fpolcitnon,pvalue_fpolcitnon = UT.calc_indttest(nons[2],nons[1])
    stat_fpolcitneg,pvalue_fpolcitneg = UT.calc_indttest(negs[2],negs[1])

    pruns = [pvalue_fsubcitpos,pvalue_fsubcitnon,pvalue_fsubcitneg,
                pvalue_fpolcitpos,pvalue_fpolcitnon,pvalue_fpolcitneg]
                                                 
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
    timeq = np.arange(N,213,1)
    timeqq,levq = np.meshgrid(timeq,lev)
    
    fig = plt.figure()
    for i in range(len(experiments)):
        ax1 = plt.subplot(2,3,i+1)
        
        varss = diffruns[i]
        pvar = pruns[i]
        
        ### Begin smoothing of N days
        varq = []
        for s in range(varss.shape[1]):
            varqq=np.convolve(varss[:,s], np.ones((N,))/N, mode='valid') 
            varq.append(varqq)
        var = np.asarray(varq)
        
        ### Mask out non-significant values
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
        
        cs = plt.contourf(timeq,lev,var,limit,extend='both')
        
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
        plt.savefig(directoryfigure + 'allExperiments_%s_MASK_daily_regionalSmooth.png' % varnames[v],
                    dpi=300)
    else:
        plt.savefig(directoryfigure + 'allExperiments_%s_daily_regionalSmooth.png' % varnames[v],
            dpi=300)
    print('Completed: Script done!')                                