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
import itertools

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
experimentsn = [r'\textbf{$\Delta$NET}']
qbophaseq = [r'QBO-W',r'QBO-E']
qbophase = ['pos','non','neg']
letters = ["a","b","c","d","e","f","g","h","i"]
period = 'D'

def calcVarResp(varnames,period,qbophase):
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
    
    diffruns = [ficthitpos,ficthitneg]
    
    ### Calculate significance for ND
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[1][pos_fict,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(tas_mo[1][non_fict,:,:],
                                                   tas_mo[0][non_hit,:,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[1][neg_fict,:,:],
                                               tas_mo[0][neg_hit,:,:])

    pruns = [pvalue_FICTHITpos,pvalue_FICTHITneg]
    
    return diffruns,pruns,climo,lat,lon

### Call variables
diffu300,pu300,climou300,lat,lon = calcVarResp('U300',period,qbophase)
diffz500,pz500,climoz500,lat,lon = calcVarResp('Z500',period,qbophase)
diffz300,pz300,climoz300,latc,lonc = calcVarResp('Z300xwave1',period,qbophase)
diffz30,pz30,climoz30,lat,lon = calcVarResp('Z30',period,qbophase)

varnamesn = [r'Z500',r'Z500',r'U300',r'U300',r'Wave 1',r'Wave 1',r'Z30',r'Z30']
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJF
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for v in range(8):
    ax = plt.subplot(4,2,v+1)
    
    ### Retrieve variables and pvalues
    if v < 2: 
        var = diffz500[v]
        pvar = pz500[v]
    elif v >= 2 and v < 4: 
        var = diffu300[v-2]
        pvar = pu300[v-2]
    elif v >= 4 and v < 6:  
        var = diffz300[v-4]
        pvar = pz300[v-4]
        climos = climoz300[v-4]
    elif v >= 6 and v < 8:  
        var = diffz30[v-6]
        pvar = pz30[v-6]
        climos = climoz30[v-6]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'U300':
        limit = np.arange(-5,5.1,0.5)
        barlim = np.arange(-5,6,5)
    elif varnamesn[v] == 'Z500':
        limit = np.arange(-50,50.1,5)
        barlim = np.arange(-50,51,25) 
    elif varnamesn[v] == 'Wave 1':
        limit = np.arange(-50,51,5)
        barlim = np.arange(-50,51,25) 
    elif varnamesn[v] == 'Z30':
        limit = np.arange(-100,100.1,5)
        barlim = np.arange(-100,101,50) 
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
 
    if varnamesn[v] == 'U300' or varnamesn[v] == 'Z500' or varnamesn[v] == 'Z30':
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])
        
        if varnamesn[v] == 'Z30':
            climoq,lons_cyclic = addcyclic(climos, lon)
            climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                colors='k',linewidths=1.2,zorder=10)

    elif varnamesn[v] == 'Wave 1': # the interval is 250 m       
        lon2c, lat2c = np.meshgrid(lonc,latc)            
        m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
        
        cs = m.contourf(lon2c,lat2c,var,limit,extend='both',latlon=True)
        cs1 = m.contourf(lon2c,lat2c,pvar,colors='None',hatches=['....'],latlon=True)
        cs2 = m.contour(lon2c,lat2c,climos,np.arange(-200,201,50),colors='k',
                        linewidths=1.2,zorder=10,latlon=True)

    m.drawcoastlines(color='dimgrey',linewidth=0.6)
    
    if varnamesn[v] == 'U300':
        cmap = cmocean.cm.balance          
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'Z500':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'Wave 1':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'Z30':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    if any([v == 0,v == 2,v == 4,v==6]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=14,rotation=90,ha='center',va='center')
    if any([v == 0,v == 1]):
        ax.annotate(r'\textbf{%s}' % qbophaseq[v],xy=(0,0),xytext=(0.5,1.12),
                     textcoords='axes fraction',color='dimgrey',
                     fontsize=13,rotation=0,ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=7)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 1:
        cbar_ax = fig.add_axes([0.70,0.72,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U300':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')  
        elif varnamesn[v] == 'Wave 1':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=-0.4)       
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 3:
        cbar_ax = fig.add_axes([0.70,0.51,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U300':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.4)  
        elif varnamesn[v] == 'Wave 1':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')      
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 5:
        cbar_ax = fig.add_axes([0.70,0.29,0.013,0.182])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U300':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')  
        elif varnamesn[v] == 'Wave 1':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=2.4)       
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 7:
        cbar_ax = fig.add_axes([0.70,0.07,0.013,0.18])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'U300':
            cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k') 
        elif varnamesn[v] == 'Z500' or varnamesn[v] == 'Z30':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.5)  
        elif varnamesn[v] == 'Wave 1':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k')       
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)

plt.tight_layout()    
fig.subplots_adjust(wspace=-0.75,hspace=0)
    
plt.savefig(directoryfigure + 'LargeScaleVars.png',dpi=900)

print('Completed: Script done!')

