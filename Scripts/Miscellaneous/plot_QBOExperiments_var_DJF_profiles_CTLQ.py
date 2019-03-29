"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, FICT. Profiles are 
organized by QBO phase (positive, neutral, negative) using indices from CTLQ

Notes
-----
    Author : Zachary Labe
    Date   : 2 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_MonthlyOutput_AllMembers as MO
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_DJF_2_CTLQ/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO profile comparisons - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['U','TEMP','GEOP','V','EGR']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'DJF'
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                               'profile')
    lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[v],'FIT',
                                               'profile')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FICT',
                                                'profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
    filenamehitno = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[1]
    filenamehitn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
    pos_hit = np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int')
    non_hit = np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int')
    neg_hit = np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int')  
    
    pos_hit[pos_hit>99] = pos_hit[pos_hit>99] + 1
    non_hit[non_hit>99] = non_hit[non_hit>99] + 1
    neg_hit[neg_hit>99] = neg_hit[neg_hit>99] + 1

    filenamefitp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
    filenamefitno = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[1]
    filenamefitn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
    pos_fit = np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int')
    non_fit = np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int')
    neg_fit = np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int')
    
    pos_fit[pos_fit>99] = pos_fit[pos_fit>99] + 1
    non_fit[non_fit>99] = non_fit[non_fit>99] + 1
    neg_fit[neg_fit>99] = neg_fit[neg_fit>99] + 1
    
    filenamefictp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
    filenamefictno = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[1]
    filenamefictn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
    pos_fict = np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int')
    non_fict = np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int')
    neg_fict = np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int')
    
    pos_fict[pos_fict>99] = pos_fict[pos_fict>99] + 1
    non_fict[non_fict>99] = non_fict[non_fict>99] + 1
    neg_fict[neg_fict>99] = neg_fict[neg_fict>99] + 1
    
    ### Concatonate runs
    runs = [tashit,tasfit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'DJF':
        tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3],
                          tashit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'profile',17)   
    else:
        ValueError('Wrong period selected! (DJF)')
        
    ### Composite by QBO phase    
    tas_mofitpos = tas_mo[1][pos_fit,:,:,:]
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:,:]
    
    tas_mofitnon = tas_mo[1][non_fit,:,:,:]
    tas_mohitnon = tas_mo[0][non_hit,:,:,:]
    tas_mofictnon = tas_mo[2][non_fict,:,:,:]
    
    tas_mofitneg = tas_mo[1][neg_fit,:,:,:]
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:,:]

    ### Compute climatology    
    climofitpos = np.nanmean(tas_mofitpos,axis=0)
    climohitpos = np.nanmean(tas_mohitpos,axis=0)
    climofictpos = np.nanmean(tas_mofictpos,axis=0)
    climofitnon = np.nanmean(tas_mofitnon,axis=0)
    climohitnon = np.nanmean(tas_mohitnon,axis=0)
    climofictnon = np.nanmean(tas_mofictnon,axis=0)
    climofitneg = np.nanmean(tas_mofitneg,axis=0)
    climohitneg = np.nanmean(tas_mohitneg,axis=0)
    climofictneg = np.nanmean(tas_mofictneg,axis=0)
    
    ### Take zonal mean of climatologies
    zclimohitpos = np.nanmean(climohitpos,axis=2)
    zclimohitnon = np.nanmean(climohitnon,axis=2)
    zclimohitneg = np.nanmean(climohitneg,axis=2)
    zclimofitpos = np.nanmean(climofitpos,axis=2)
    zclimofitnon = np.nanmean(climofitnon,axis=2)
    zclimofitneg = np.nanmean(climofitneg,axis=2)
    zclimofictpos = np.nanmean(climofictpos,axis=2)
    zclimofictnon = np.nanmean(climofictnon,axis=2)
    zclimofictneg = np.nanmean(climofictneg,axis=2)
    climo = [zclimohitpos,zclimohitnon,zclimohitneg,
             zclimohitpos,zclimohitnon,zclimohitneg]
    
    ### Compute comparisons for months - taken ensemble average
    fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)
    fithitnon = np.nanmean(tas_mofitnon - tas_mohitnon,axis=0)
    fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)
    
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitnon = np.nanmean(tas_mofictnon - tas_mohitnon,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    
    ### Take zonal mean of experiments
    zfithitpos = np.nanmean(fithitpos,axis=2)
    zfithitnon = np.nanmean(fithitnon,axis=2)
    zfithitneg = np.nanmean(fithitneg,axis=2)
    zficthitpos = np.nanmean(ficthitpos,axis=2)
    zficthitnon = np.nanmean(ficthitnon,axis=2)
    zficthitneg = np.nanmean(ficthitneg,axis=2)
    diffruns_mo = [zfithitpos,zfithitnon,zfithitneg,
                   zficthitpos,zficthitnon,zficthitneg]
    
    ### Calculate significance for FM
    stat_FITHITpos,pvalue_FITHITpos = UT.calc_indttest(
            np.nanmean(tas_mo[1][pos_fit,:,:,:],axis=3),
            np.nanmean(tas_mo[0][pos_hit,:,:,:],axis=3))
    stat_FITHITnon,pvalue_FITHITnon = UT.calc_indttest(
            np.nanmean(tas_mo[1][non_fit,:,:,:],axis=3),
            np.nanmean(tas_mo[0][non_hit,:,:,:],axis=3))
    stat_FITHITneg,pvalue_FITHITneg = UT.calc_indttest(
            np.nanmean(tas_mo[1][neg_fit,:,:,:],axis=3),
            np.nanmean(tas_mo[0][neg_hit,:,:,:]))
    
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(
            np.nanmean(tas_mo[2][pos_fict,:,:,:],axis=3),
            np.nanmean(tas_mo[0][pos_hit,:,:,:],axis=3))
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(
            np.nanmean(tas_mo[2][non_fict,:,:,:],axis=3),
            np.nanmean(tas_mo[0][non_hit,:,:,:],axis=3))
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(
            np.nanmean(tas_mo[2][neg_fict,:,:,:],axis=3),
            np.nanmean(tas_mo[0][neg_hit,:,:,:],axis=3))

    pruns_mo = [pvalue_FITHITpos,pvalue_FITHITnon,pvalue_FITHITneg,
                pvalue_FICTHITpos,pvalue_FICTHITnon,pvalue_FICTHITneg]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #### Plot U
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'U':
        limit = np.arange(-3,3.1,0.1)
        barlim = np.arange(-3,4,1)
    elif varnames[v] == 'TEMP':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,1)
    elif varnames[v] == 'GEOP':
        limit = np.arange(-60,61,1)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'V':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.1)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.08,0.081,0.005)
        barlim = np.arange(-0.08,0.09,0.04)
        
    zscale = np.array([1000,700,500,300,200,
                        100,50,30,10])
    latq,levq = np.meshgrid(lat,lev)
    
    fig = plt.figure()
    for i in range(len(diffruns_mo)):
        ax1 = plt.subplot(2,3,i+1)
        
        var = diffruns_mo[i]
        pvar = pruns_mo[i]
        climovar = climo[i]
        
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
        
        
        cs = plt.contourf(lat,lev,var,limit,extend='both')
        
        if varnames[v] == 'U': 
            cs2 = plt.contour(lat,lev,climovar,np.arange(-20,101,5),
                              linewidths=0.6,colors='dimgrey')
        plt.contourf(latq,levq,pvar,colors='None',hatches=['////'],
                     linewidth=5)   
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xlim([0,90])
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
        plt.minorticks_off()
        
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
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.curl           
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
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
        
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.5)
    plt.subplots_adjust(bottom=0.21)
    
    plt.savefig(directoryfigure + 'DJF_vertical_%s_QBO.png' % varnames[v],dpi=300)
    print('Completed: Script done!')
