"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes fsub, cit, fpol. Profiles are 
organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 26 March 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import cmocean

### Define directories
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_DJF_2/'
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
runnames = [r'CIT',r'FSUB',r'FPOL']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',
               r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}']
qbophase = ['pos','non','neg']
period = 'DJF'
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tascit = MO.readExperi(directorydata2,
                                            '%s' % varnames[v],'CIT','profile')
    lat,lon,time,lev,tasfsub = MO.readExperi(directorydata2,
                                            '%s' % varnames[v],'FSUB','profile')
    lat,lon,time,lev,tasfpol = MO.readExperi(directorydata2,
                                             '%s' % varnames[v],'FPOL','profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamefsubp = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    pos_fsub = np.genfromtxt(filenamefsubp,unpack=True,usecols=[0],dtype='int')
    non_fsub = np.genfromtxt(filenamefsubno,unpack=True,usecols=[0],dtype='int')
    neg_fsub = np.genfromtxt(filenamefsubn,unpack=True,usecols=[0],dtype='int')
    
    filenamecitp = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    pos_cit = np.genfromtxt(filenamecitp,unpack=True,usecols=[0],dtype='int')
    non_cit = np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int')
    neg_cit = np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int')
    
    filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    pos_fpol = np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int')
    non_fpol = np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int')
    neg_fpol = np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int')
    
    ### Concatonate runs
    runs = [tascit,tasfsub,tasfpol]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'DJF':
        tas_mo= np.empty((3,tascit.shape[0]-1,tascit.shape[2],tascit.shape[3],
                          tascit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'profile',17)   
    else:
        ValueError('Wrong period selected! (DJF)')
        
    ### Composite by QBO phase    
    tas_mofsubpos = tas_mo[1][pos_fsub,:,:,:]
    tas_mocitpos = tas_mo[0][pos_cit,:,:,:]
    tas_mofpolpos = tas_mo[2][pos_fpol,:,:,:]
    
    tas_mofsubnon = tas_mo[1][non_fsub,:,:,:]
    tas_mocitnon = tas_mo[0][non_cit,:,:,:]
    tas_mofpolnon = tas_mo[2][non_fpol,:,:,:]
    
    tas_mofsubneg = tas_mo[1][neg_fsub,:,:,:]
    tas_mocitneg = tas_mo[0][neg_cit,:,:,:]
    tas_mofpolneg = tas_mo[2][neg_fpol,:,:,:]

    ### Compute climatology    
    climofsubpos = np.nanmean(tas_mofsubpos,axis=0)
    climocitpos = np.nanmean(tas_mocitpos,axis=0)
    climofpolpos = np.nanmean(tas_mofpolpos,axis=0)
    climofsubnon = np.nanmean(tas_mofsubnon,axis=0)
    climocitnon = np.nanmean(tas_mocitnon,axis=0)
    climofpolnon = np.nanmean(tas_mofpolnon,axis=0)
    climofsubneg = np.nanmean(tas_mofsubneg,axis=0)
    climocitneg = np.nanmean(tas_mocitneg,axis=0)
    climofpolneg = np.nanmean(tas_mofpolneg,axis=0)
    
    ### Take zonal mean of climatologies
    zclimocitpos = np.nanmean(climocitpos,axis=2)
    zclimocitnon = np.nanmean(climocitnon,axis=2)
    zclimocitneg = np.nanmean(climocitneg,axis=2)
    zclimofsubpos = np.nanmean(climofsubpos,axis=2)
    zclimofsubnon = np.nanmean(climofsubnon,axis=2)
    zclimofsubneg = np.nanmean(climofsubneg,axis=2)
    zclimofpolpos = np.nanmean(climofpolpos,axis=2)
    zclimofpolnon = np.nanmean(climofpolnon,axis=2)
    zclimofpolneg = np.nanmean(climofpolneg,axis=2)
    climo = [zclimocitpos,zclimocitnon,zclimocitneg,
             zclimocitpos,zclimocitnon,zclimocitneg]
    
    ### Compute comparisons for months - taken ensemble average
    fsubcitpos = np.nanmean(tas_mofsubpos - tas_mocitpos,axis=0)
    fsubcitnon = np.nanmean(tas_mofsubnon - tas_mocitnon,axis=0)
    fsubcitneg = np.nanmean(tas_mofsubneg - tas_mocitneg,axis=0)
    
    fpolcitpos = np.nanmean(tas_mofpolpos - tas_mocitpos,axis=0)
    fpolcitnon = np.nanmean(tas_mofpolnon - tas_mocitnon,axis=0)
    fpolcitneg = np.nanmean(tas_mofpolneg - tas_mocitneg,axis=0)
    
    ### Take zonal mean of experiments
    zfsubcitpos = np.nanmean(fsubcitpos,axis=2)
    zfsubcitnon = np.nanmean(fsubcitnon,axis=2)
    zfsubcitneg = np.nanmean(fsubcitneg,axis=2)
    zfpolcitpos = np.nanmean(fpolcitpos,axis=2)
    zfpolcitnon = np.nanmean(fpolcitnon,axis=2)
    zfpolcitneg = np.nanmean(fpolcitneg,axis=2)
    diffruns_mo = [zfsubcitpos,zfsubcitnon,zfsubcitneg,
                   zfpolcitpos,zfpolcitnon,zfpolcitneg]
    
    ### Calculate significance for FM
    stat_fsubcitpos,pvalue_fsubcitpos = UT.calc_indttest(
            np.nanmean(tas_mo[1][pos_fsub,:,:,:],axis=3),
            np.nanmean(tas_mo[0][pos_cit,:,:,:],axis=3))
    stat_fsubcitnon,pvalue_fsubcitnon = UT.calc_indttest(
            np.nanmean(tas_mo[1][non_fsub,:,:,:],axis=3),
            np.nanmean(tas_mo[0][non_cit,:,:,:],axis=3))
    stat_fsubcitneg,pvalue_fsubcitneg = UT.calc_indttest(
            np.nanmean(tas_mo[1][neg_fsub,:,:,:],axis=3),
            np.nanmean(tas_mo[0][neg_cit,:,:,:]))
    
    stat_fpolcitpos,pvalue_fpolcitpos = UT.calc_indttest(
            np.nanmean(tas_mo[2][pos_fpol,:,:,:],axis=3),
            np.nanmean(tas_mo[0][pos_cit,:,:,:],axis=3))
    stat_fpolcitnon,pvalue_fpolcitnon = UT.calc_indttest(
            np.nanmean(tas_mo[2][non_fpol,:,:,:],axis=3),
            np.nanmean(tas_mo[0][non_cit,:,:,:],axis=3))
    stat_fpolcitneg,pvalue_fpolcitneg = UT.calc_indttest(
            np.nanmean(tas_mo[2][neg_fpol,:,:,:],axis=3),
            np.nanmean(tas_mo[0][neg_cit,:,:,:],axis=3))

    pruns_mo = [pvalue_fsubcitpos,pvalue_fsubcitnon,pvalue_fsubcitneg,
                pvalue_fpolcitpos,pvalue_fpolcitnon,pvalue_fpolcitneg]
    
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
    
    plt.savefig(directoryfigure + 'Regional/DJF_vertical_%s_QBO.png' % varnames[v],dpi=300)
    print('Completed: Script done!')
