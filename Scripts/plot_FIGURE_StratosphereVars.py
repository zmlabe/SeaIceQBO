"""
Script plots figure for manuscript of stratosphere variables of U10 and Z30

Notes
-----
    Author : Zachary Labe
    Date   : 17 October 2018
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
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['U10','Z30']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'

def readVar(varnames,qbophase,period):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
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
    runs = [tashit,tasfict]
    
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
    tas_mofictpos = tas_mo[1][pos_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:]

    ### Compute climatology    
    climohitpos = np.nanmean(tas_mohitpos,axis=0)
    climohitneg = np.nanmean(tas_mohitneg,axis=0)
    climo = [climohitpos,climohitneg]
    
    ### Compute comparisons for months - taken ensemble average
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    diffruns_mo = [ficthitpos,ficthitneg]
    
    ### Calculate significance for ND
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[1][pos_fict,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(tas_mo[1][non_fict,:,:],
                                                   tas_mo[0][non_hit,:,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[1][neg_fict,:,:],
                                               tas_mo[0][neg_hit,:,:])

    pruns_mo = [pvalue_FICTHITpos,pvalue_FICTHITneg]
    
    return lat,lon,climo,diffruns_mo,pruns_mo
    
### Read in data functions
lat,lon,climoU,diffU,pU = readVar('U10',qbophase,period)
lat,lon,climoZ,diffZ,pZ = readVar('Z30',qbophase,period)

### Set variables for plotting
varnames = [r'U10',r'U10',r'Z30',r'Z30']
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for QBO composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
for i in range(4):
    ### Set limits for contours and colorbars
    if varnames[i] == 'Z30':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,50) 
    elif varnames[i] == 'U10' :
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,2)
    
    if i <= 1:
        var = diffU[i]
        pvar = pU[i]
        climo = climoU[i]
    elif i > 1:
        var = diffZ[i-2]
        pvar = pZ[i-2]
        climo = climoZ[i-2]
    
    ax1 = plt.subplot(2,2,i+1)
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    climoq,lons_cyclic = addcyclic(climo, lon)
    climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    
    cs = m.contourf(x,y,var,limit,extend='both')
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                     linewidths=0.4)
    if varnames[i] == 'Z30': # the interval is 250 m 
        cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                        colors='k',linewidths=1.5,zorder=10)

    m.drawcoastlines(color='dimgray',linewidth=0.8)
    
    if varnames[i] == 'Z30':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    elif varnames[i] == 'U10':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap)  
                
    ### Add experiment text to subplot
    if i < 2:
        qbophaseq = [r'QBO-W',r'QBO-E']
        ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                     textcoords='axes fraction',color='dimgray',
                     fontsize=13,rotation=0,ha='center',va='center')
    if i == 0 or i == 2:
        ax1.annotate(r'\textbf{%s}' % varnames[i],xy=(0,0),xytext=(-0.1,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=20,rotation=90,ha='center',va='center')
        
###########################################################################
    if i == 1:
        cbar_ax = fig.add_axes([0.82,0.535,0.017,0.3])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnames[i] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnames[i] == 'Z30':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.2)   
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif i == 3:
        cbar_ax = fig.add_axes([0.82,0.158,0.017,0.3])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnames[i] == 'U10':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnames[i] == 'Z30':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=4.5)      
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)

plt.subplots_adjust(wspace=-0.4)
plt.subplots_adjust(hspace=0.001)

plt.savefig(directoryfigure + 'StratosphereVars.png', dpi=900)

print('Completed: Script done!')

