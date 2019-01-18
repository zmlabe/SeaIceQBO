"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, CIT, FIC, FICT. Composites are 
organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 28 November 2017
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
directoryfigure = '/home/zlabe/Desktop/NOQBO/'
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
varnames = ['U10']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
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
    elif period == 'N':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-2,:,:]
    elif period == 'ND':
        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofitpos = tas_mo[1][pos_fit,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofitneg = tas_mo[1][neg_fit,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:]

    diffruns_mo = [np.nanmean(tas_mofitneg,axis=0) - np.nanmean(tas_mofitpos,axis=0),
                   np.nanmean(tas_mofictneg,axis=0) - np.nanmean(tas_mofictpos,axis=0)]
    
    ### Calculate significance
    stat_FITHITdiff,pvalue_FITHITdiff = UT.calc_indttest(tas_mo[1][pos_fit,:,:],
                                                       tas_mo[1][pos_fit,:,:])
    stat_FICTHITdiff,pvalue_FICTHITdiff = UT.calc_indttest(tas_mo[2][pos_fict,:,:],
                                                   tas_mo[2][neg_fict,:,:])

    pruns_mo = [pvalue_FITHITdiff,pvalue_FICTHITdiff]
    
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
    elif varnames[v] == 'SLP':
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,3)
    elif varnames[v] == 'Z500':
        limit = np.arange(-50,50.1,0.1)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z30':
        limit = np.arange(-100,100.1,1)
        barlim = np.arange(-100,101,50) 
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        limit = np.arange(-6,6.1,0.1)
        barlim = np.arange(-6,7,2)
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
    elif varnames[v] == 'WAFZ850':
        limit = np.arange(-0.1,0.101,0.001)
        barlim = np.arange(-0.1,0.11,0.1)
    elif varnames[v] == 'WAFZ150':
        limit = np.arange(-0.01,0.0101,0.0001)
        barlim = np.arange(-0.01,0.011,0.01)
    elif varnames[v] == 'WAFY850':
        limit = np.arange(-5,5.1,0.1)
        barlim = np.arange(-5,6,5)
    elif varnames[v] == 'WAFY150':
        limit = np.arange(-2,2.01,0.05)
        barlim = np.arange(-2,3,2)
    
    fig = plt.figure()
    for i in range(2):
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
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                         linewidths=0.4)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'T2M':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        elif varnames[v] == 'SLP':
            cmap = cmocean.cm.balance          
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z30':
            cmap = cmocean.cm.balance  
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
            cmap = cmocean.cm.balance          
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
        elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
        elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
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
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
    elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')
    elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + 'QBOExperiments_difference_%s_%s.png' % (
                                                                        period,
                                                                  varnames[v]),
                                                                  dpi=300)

print('Completed: Script done!')

