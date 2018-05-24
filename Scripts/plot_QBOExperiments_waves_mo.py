"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, FICT. Composites are 
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
varnames = ['Z300xwave2']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
wv = 'Z300xwave2'
MASK = False
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                               'surface')
    lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[v],'FIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FICT',
                                                'surface')
    
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
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    

    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
    filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+101)
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)
    
    ### Concatonate runs
    runs = [tashit,tasfit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'ON': 
        tas_mo = np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1) 
    elif period == 'DJ':     
        tas_mo = np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
                                                lon,'surface',1) 
    elif period == 'FM':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
    elif period == 'DJF':
        tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                  lon,'surface',1)   
    elif period == 'M':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,2,:,:]
    elif period == 'D':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:]
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofitpos = tas_mo[1][pos_fit,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:]
    
    tas_mohitnon = tas_mo[0][non_hit,:,:]
    tas_mofitnon = tas_mo[1][non_fit,:,:]
    tas_mofictnon = tas_mo[2][non_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofitneg = tas_mo[1][neg_fit,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:]
    
    ### Compute comparisons for months - taken ensemble average
    fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)
    fithitnon = np.nanmean(tas_mofitnon - tas_mohitnon,axis=0)
    fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)
    
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitnon = np.nanmean(tas_mofictnon - tas_mohitnon,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    diffruns_mo = [fithitpos,fithitnon,fithitneg,
                   ficthitpos,ficthitnon,ficthitneg]
    
    ### Calculate significance for FM
    stat_FITHITpos,pvalue_FITHITpos = UT.calc_indttest(tas_mo[1][pos_fit,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FITHITnon,pvalue_FITHITnon = UT.calc_indttest(tas_mo[1][non_fit,:,:],
                                                   tas_mo[0][non_hit,:,:])
    stat_FITHITneg,pvalue_FITHITneg = UT.calc_indttest(tas_mo[1][neg_fit,:,:],
                                               tas_mo[0][neg_hit,:,:])
    
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[2][pos_fict,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(tas_mo[2][non_fict,:,:],
                                                   tas_mo[0][non_hit,:,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[2][neg_fict,:,:],
                                               tas_mo[0][neg_hit,:,:])

    pruns_mo = [pvalue_FITHITpos,pvalue_FITHITnon,pvalue_FITHITneg,
                pvalue_FICTHITpos,pvalue_FICTHITnon,pvalue_FICTHITneg]
    
    ### Read in waves
    latc,lonc,time,lev,waveh=  MO.readExperiAll('%s' % wv,'HIT','surface')
    
    if period == 'D':
        climowavehmo = waveh[:,-1,:,:]
        
    wavehitpos = np.nanmean(climowavehmo[pos_hit,:,:],axis=0)
    wavehitnon = np.nanmean(climowavehmo[non_hit,:,:],axis=0)
    wavehitneg = np.nanmean(climowavehmo[neg_hit,:,:],axis=0)
    climo = [wavehitpos,wavehitnon,wavehitneg,wavehitpos,wavehitnon,wavehitneg]
    
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
    elif varnames[v] == 'Z300xwave1' or varnames[v] == 'Z300xwave2':
        limit = np.arange(-20,20.1,0.5)
        barlim = np.arange(-20,21,20) 
        
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
        
#        var, lons_cyclic = addcyclic(var, lon)
#        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
#        x, y = m(lon2d, lat2d)
        lon2c, lat2c = np.meshgrid(lonc,latc)
#        
#        pvar,lons_cyclic = addcyclic(pvar, lon)
#        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
               
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        cs = m.contourf(lon2c,lat2c,var,limit,extend='both',latlon=True)
        if MASK == False:
            cs1 = m.contourf(lon2c,lat2c,pvar,colors='None',hatches=['....'],
                             linewidths=0.4,latlon=True)
        cs2 = m.contour(lon2c,lat2c,climoq,limitc,colors='k',linewidths=1.5,
                        zorder=10,latlon=True)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'Z300':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z300xwave1' or varnames[v] == 'Z300xwave2':
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
    elif varnames[v] == 'Z300xwave1' or varnames[v] == 'Z300xwave2':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + '/QBO_%s_2/QBOExperiments_%s_%s.png' % (period,
                                                                        period,
                                                                  wv),
                                                                  dpi=300)

print('Completed: Script done!')

