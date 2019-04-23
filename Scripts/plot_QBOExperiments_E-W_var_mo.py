"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, FICT. Composites are 
organized by QBO-E - QBO-W

Notes
-----
    Author : Zachary Labe
    Date   : 31 January 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_MonthlyOutput as MO
import calc_Utilities as UT
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO comparisons - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['Z500','Z30','SLP','T2M','U10','U300','SWE','THICK','P','EGR']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'DJF'
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,tasfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,tasfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    pos_fit = np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int')
    non_fit = np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int')
    neg_fit = np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int')
    
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int')
    non_hit = np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int')
    neg_hit = np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int')
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int')
    non_fict = np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int')
    neg_fict = np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int')
    
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
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mofitpos = tas_mo[1][pos_fit,:,:]
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:]
    
    tas_mofitnon = tas_mo[1][non_fit,:,:]
    tas_mohitnon = tas_mo[0][non_hit,:,:]
    tas_mofictnon = tas_mo[2][non_fict,:,:]
    
    tas_mofitneg = tas_mo[1][neg_fit,:,:]
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:]

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
    climo = [climohitpos,climohitnon,climohitneg,
             climohitpos,climohitnon,climohitneg]
    
    ### Compute comparisons for months - taken ensemble average
    fithit = np.nanmean((tas_mofitneg-tas_mohitneg) - (tas_mofitpos[:32]-tas_mohitpos[:32]),axis=0)
    
    ficthit = np.nanmean((tas_mofictneg-tas_mohitneg) - (tas_mofictpos[:32]-tas_mohitpos[:32]),axis=0)
    diffruns_mo = [fithit,ficthit]
    
    ### Calculate significance for FM
    stat_FITHIT,pvalue_FITHIT = UT.calc_indttest(tas_mofitneg-tas_mohitneg,tas_mofitpos[:32]-tas_mohitpos[:32])
    
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(tas_mofictneg-tas_mohitneg,tas_mofictpos[:32]-tas_mohitpos[:32])

    pruns_mo = [pvalue_FITHIT,pvalue_FICTHIT]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot variable data for QBO composites
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'Z500':
        limit = np.arange(-60,60.1,1)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z30':
        limit = np.arange(-100,100.1,5)
        barlim = np.arange(-100,101,50)
    elif varnames[v] == 'SLP':
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,3)
    elif varnames[v] == 'U10' or varnames[v] == 'U300':
        limit = np.arange(-10,10.1,1)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    elif varnames[v] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,1) 
    elif varnames[v] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.2)
    
    fig = plt.figure()
    for i in range(len(diffruns_mo)):
        var = diffruns_mo[i]
        pvar = pruns_mo[i]
        
        ax1 = plt.subplot(1,2,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
        climoq,lons_cyclic = addcyclic(climo[i], lon)
        climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                         linewidths=0.4)
        if varnames[v] == 'Z30': # the interval is 250 m 
            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                            colors='k',linewidths=1.5,zorder=10)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'T2M':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
            cmap = ncm.cmap('nrl_sirkes')           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z30':
            cmap = ncm.cmap('nrl_sirkes')  
            cs.set_cmap(cmap)  
        elif varnames[v] == 'SLP':
            cmap = ncm.cmap('nrl_sirkes')           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U300':
            cmap = ncm.cmap('temp_diff_18lev')           
            cs.set_cmap(cmap)           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'SWE':
            cmap = cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
        elif varnames[v] == 'P':
            cmap = ncm.cmap('precip4_diff_19lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'THICK':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
                    
        ### Add experiment text to subplot
        ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(0.5,1.05),
                     textcoords='axes fraction',color='dimgrey',
                     fontsize=23,rotation=0,ha='center',va='center')
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'Z30':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + '/QBO_%s/QBOExperiments_E-W_%s_%s.png' % (period,
                                                                        period,
                                                                  varnames[v]),
                                                                  dpi=300)

print('Completed: Script done!')

