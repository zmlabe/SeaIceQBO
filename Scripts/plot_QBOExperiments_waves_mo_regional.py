"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes CIT, FSUB, FPOL. Composites are 
organized by QBO phase (positive, neutral, negative). Plot shows 
climatological and forced wave-numbers 1 and 2.

Notes
-----
    Author : Zachary Labe
    Date   : 24 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput_AllMembers as MO
import read_MonthlyOutput_AllRegional as MOR
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
print('\n' '----Plotting QBO comparisons for waves - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['Z300']
runnames = [r'CIT',r'FSUB',r'FPOL']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',r'\textbf{FSUB-CIT}',
               r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}']
qbophase = ['pos','non','neg']
period = 'N'
wv = 'Z300xwave2'
MASK = False
for v in range(len(varnames)):
    ### Call function for data from reach run
    lat,lon,time,lev,tascit = MO.readExperiAll('%s' % varnames[v],
                                                         'CIT','surface')
    lat,lon,time,lev,tasfsub = MOR.readExperiAllRegional('%s' % varnames[v],
                                                         'FSUB','surface')
    lat,lon,time,lev,tasfpol = MOR.readExperiAllRegional('%s' % varnames[v],
                                                         'FPOL','surface')
    
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
                        np.genfromtxt(filenamecitp2,unpack=True,usecols=[0],dtype='int')+101)
    non_cit = np.append(np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_cit = np.append(np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitn2,unpack=True,usecols=[0],dtype='int')+101)    

    filenamefsubp = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    filenamefsubp2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    pos_fsub = np.append(np.genfromtxt(filenamefsubp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fsub = np.append(np.genfromtxt(filenamefsubno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fsub = np.append(np.genfromtxt(filenamefsubn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubn2,unpack=True,usecols=[0],dtype='int')+101)
    
    filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    filenamefpolp2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    pos_fpol = np.append(np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fpol = np.append(np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fpol = np.append(np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpoln2,unpack=True,usecols=[0],dtype='int')+101)
    
    ### Concatonate runs
    runs = [tascit,tasfsub,tasfpol]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'ON': 
        tas_mo = np.empty((3,tascit.shape[0],tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1) 
    elif period == 'DJ':     
        tas_mo = np.empty((3,tascit.shape[0]-1,tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
                                                lon,'surface',1) 
    elif period == 'FM':
        tas_mo= np.empty((3,tascit.shape[0],tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
    elif period == 'DJF':
        tas_mo= np.empty((3,tascit.shape[0]-1,tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'surface',1)   
    elif period == 'M':
        tas_mo= np.empty((3,tascit.shape[0],tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,2,:,:]
    elif period == 'D':
        tas_mo= np.empty((3,tascit.shape[0],tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:]     
    elif period == 'N':
        tas_mo= np.empty((3,tascit.shape[0],tascit.shape[2],tascit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-2,:,:]    
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mofsubpos = tas_mo[1][pos_fsub,:,:]
    tas_mocitpos = tas_mo[0][pos_cit,:,:]
    tas_mofpolpos = tas_mo[2][pos_fpol,:,:]
    
    tas_mofsubnon = tas_mo[1][non_fsub,:,:]
    tas_mocitnon = tas_mo[0][non_cit,:,:]
    tas_mofpolnon = tas_mo[2][non_fpol,:,:]
    
    tas_mofsubneg = tas_mo[1][neg_fsub,:,:]
    tas_mocitneg = tas_mo[0][neg_cit,:,:]
    tas_mofpolneg = tas_mo[2][neg_fpol,:,:]
    
    ### Compute comparisons for months - taken ensemble average
    fsubcitpos = np.nanmean(tas_mofsubpos - tas_mocitpos,axis=0)
    fsubcitnon = np.nanmean(tas_mofsubnon - tas_mocitnon,axis=0)
    fsubcitneg = np.nanmean(tas_mofsubneg - tas_mocitneg,axis=0)
    
    fpolcitpos = np.nanmean(tas_mofpolpos - tas_mocitpos,axis=0)
    fpolcitnon = np.nanmean(tas_mofpolnon - tas_mocitnon,axis=0)
    fpolcitneg = np.nanmean(tas_mofpolneg - tas_mocitneg,axis=0)
    diffruns_mo = [fsubcitpos,fsubcitnon,fsubcitneg,
                   fpolcitpos,fpolcitnon,fpolcitneg]
    
    ### Calculate significance for FM
    stat_fsubcitpos,pvalue_fsubcitpos = UT.calc_indttest(tas_mo[1][pos_fsub,:,:],
                                                       tas_mo[0][pos_cit,:,:])
    stat_fsubcitnon,pvalue_fsubcitnon = UT.calc_indttest(tas_mo[1][non_fsub,:,:],
                                                   tas_mo[0][non_cit,:,:])
    stat_fsubcitneg,pvalue_fsubcitneg = UT.calc_indttest(tas_mo[1][neg_fsub,:,:],
                                               tas_mo[0][neg_cit,:,:])
    
    stat_fpolcitpos,pvalue_fpolcitpos = UT.calc_indttest(tas_mo[2][pos_fpol,:,:],
                                                       tas_mo[0][pos_cit,:,:])
    stat_fpolcitnon,pvalue_fpolcitnon = UT.calc_indttest(tas_mo[2][non_fpol,:,:],
                                                   tas_mo[0][non_cit,:,:])
    stat_fpolcitneg,pvalue_fpolcitneg = UT.calc_indttest(tas_mo[2][neg_fpol,:,:],
                                               tas_mo[0][neg_cit,:,:])

    pruns_mo = [pvalue_fsubcitpos,pvalue_fsubcitnon,pvalue_fsubcitneg,
                pvalue_fpolcitpos,pvalue_fpolcitnon,pvalue_fpolcitneg]
    
    ### Read in waves
    latc,lonc,time,lev,waveh=  MO.readExperiAll('%s' % wv,'CIT','surface')
    
    if period == 'D':
        climowavehmo = waveh[:,-1,:,:]
    elif period == 'N':
        climowavehmo = waveh[:,-2,:,:]
        
    wavecitpos = np.nanmean(climowavehmo[pos_cit,:,:],axis=0)
    wavecitnon = np.nanmean(climowavehmo[non_cit,:,:],axis=0)
    wavecitneg = np.nanmean(climowavehmo[neg_cit,:,:],axis=0)
    climo = [wavecitpos,wavecitnon,wavecitneg,wavecitpos,wavecitnon,wavecitneg]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot variable data for QBO composites
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'Z300':
        limit = np.arange(-60,60.1,5)
        barlim = np.arange(-60,61,30) 
        
    if wv == 'Z300xwave1':
        limitc = np.arange(-200,201,50)
    elif wv == 'Z300xwave2':
        limitc = np.arange(-160,170,30)
    
    fig = plt.figure()
    for i in range(len(diffruns_mo)):
        var = diffruns_mo[i]
        pvar = pruns_mo[i]
        climoq = climo[i]
        
        if MASK == True:
            pvar2 = pvar.copy()
            pvar2[np.isnan(pvar2)]=0.0
            var = var*pvar2
        
        ax1 = plt.subplot(2,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        lon2c, lat2c = np.meshgrid(lonc,latc)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
               
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        if MASK == False:
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)
        cs2 = m.contour(lon2c,lat2c,climoq,limitc,colors='k',linewidths=1.5,
                        zorder=10,latlon=True)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'Z300':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap)  
                    
        ### Add experiment text to subplot
        if i < 3:
            qbophaseq = [r'QBO-W',r'QBO-N',r'QBO-E']
            ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')
        if i == 0 or i == 3:
            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(-0.1,0.5),
                         textcoords='axes fraction',color='k',
                         fontsize=20,rotation=90,ha='center',va='center')
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)

    if varnames[v] == 'Z300':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + '/QBO_%s_2/Regional2/QBOExperiments_%s_%s.png' % (period,period,wv),
                                                                  dpi=300)

print('Completed: Script done!')

