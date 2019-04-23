"""
Plot sea ice forcing files for WACCM4 experiments for concentration
and thickness. Data derives from the LENS ensemble mean. Also plot net
surface heat flux in QBO experiments.

Notes
-----
    Reference : Kay et al. [2014] - CESM Large Ensemble Project (LENS)
    Author : Zachary Labe
    Date   : 17 October 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import read_MonthlyOutput_AllMembers as MO
import cmocean
import palettable.cubehelix as cm
import datetime
import calc_Utilities as UT

### Define directories
directorydata = '/seley/zlabe/LENS/ForcingPerturb/' 
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting forcing files - %s----' % titletime)

months = [r'OCT',r'NOV',r'DEC']

### Read in HIT forcing file
def readSeaIce(directorydata):
    data = Dataset(directorydata + 'SST-SIT_lens_1976-2005.nc')
    lons = data.variables['lon'][:]
    lats = data.variables['lat'][:]
    sith = data.variables['ice_thick'][:,:,:]
    sich = data.variables['ice_cov'][:,:,:]
    data.close()
    
    ### Read in FIT forcing file
    data = Dataset(directorydata + 'SST-SIT_lens_2051-2080.nc')
    sitf = data.variables['ice_thick'][:,:,:]
    data.close()
    
    ### Read in FIC forcing file
    data = Dataset(directorydata + 'SST-SIC-SIT_lens_2051-2080_FIC.nc')
    sicf = data.variables['ice_cov'][:,:,:]
    data.close()
    
    lons2,lats2 = np.meshgrid(lons,lats)
    
    print('Completed: Data read!')
    
    ### Create 2d array of latitude and longitude
    lons2,lats2 = np.meshgrid(lons,lats)
    
    ### Average over DJF 
    sith[np.where(sith == 0)] = np.nan            
    varh = sith
    sitf[np.where(sitf == 0)] = np.nan            
    varf = sitf
    
    sich[np.where(sich == 0)] = np.nan            
    varch = sich*100 # convert SIC to 1-100%
    sicf[np.where(sicf == 0)] = np.nan            
    varcf = sicf*100 # convert SIC to 1-100%
    
    ### Use land/Arctic mask
    varh[np.where(varh == 2.)] = np.nan
    varhtemp = varh.copy()
    varhtemp[np.where(varhtemp > 0)] = 1.
    varf = varf*varhtemp
    varf[np.where(varf == 2.)] = np.nan
    
    ### Calculate differences 
    diffsitq = varf - varh
    diffsicq = varcf - varch
    
    diffsit = np.append(diffsitq[9:],diffsitq[:3],axis=0)
    diffsic = np.append(diffsicq[9:],diffsicq[:3],axis=0)
    
    diffs = np.append(diffsic[:3,:,:],diffsit[:3,:,:],axis=0)
    return diffs,lats2,lons2

def readHeatFlux():
    ### Call arguments
    varnames = 'RNET'
    experiments = [r'\textbf{FIT--HIT}',r'\textbf{FICT--HIT}']
    
    ### Call function for rnet data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
                                                'surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Calculate Oct-Mar
    tashitmo = np.append(tashit[:,9:,:,:],tashit[:,:3,:,:],axis=1)
    tasfictmo = np.append(tasfict[:,9:,:,:],tasfict[:,:3,:,:],axis=1)
    
    ### Take difference
    ficthit = np.nanmean(tasfictmo - tashitmo,axis=0)
    
    ### Statistical Test
    stat_FICTHIT,pvalue_FICTHIT = UT.calc_indttest(tasfictmo,tashitmo)

    return ficthit[:3],lat,lon,pvalue_FICTHIT[:3]

### Read in data functions
diffs,lat,lon = readSeaIce(directorydata)
heat,lat2,lon2,pvalues = readHeatFlux()

### Create variable names 
varnamesn = ['SIC','SIC','SIC','SIT','SIT','SIT','RNET','RNET','RNET']
months = [r'\textbf{OCT}',r'\textbf{NOV}',r'\textbf{DEC}']
letters = ["a","b","c","d","e","f","g","h","i"]

###########################################################################
###########################################################################
###########################################################################
### Plot variable data for Oct-Dec
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for v in range(9):
    ax = plt.subplot(3,3,v+1)
    
    ### Retrieve variables and pvalues
    if v <= 5:
        var = diffs[v]
    elif v >= 6:
        lat = lat2
        lon = lon2
        var = heat[v-6]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'SIC':
        limit = np.arange(-100,1,5)
        barlim = np.arange(-100,1,25)
        extendq = 'max'
    elif varnamesn[v] == 'SIT':
        limit = np.arange(-3,0.1,0.25)
        barlim = np.arange(-3,1,1)
        extendq = 'min'
    elif varnamesn[v] == 'RNET':
        limit = np.arange(-75,75.1,5)
        barlim = np.arange(-75,76,75) 
    
    m = Basemap(projection='npstere',boundinglat=51,lon_0=270,resolution='l',
                round =True,area_thresh=10000)
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    
    if v <= 5:
        cs = m.contourf(lon,lat,var,limit,extend=extendq,latlon=True)
    elif v >= 6:
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
    
        cs = m.contourf(x,y,var*-1,limit,extend='both')
    
    if v >= 6:
        pvar = pvalues[v-6]
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'])

    m.drawcoastlines(color='darkgray',linewidth=0.3)
    m.fillcontinents(color='dimgrey')
    
    if varnamesn[v] == 'SIC':
        cmap = cm.jim_special_16.mpl_colormap     
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'SIT':
        cmap = cm.jim_special_16.mpl_colormap       
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'RNET':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
            
    ### Add experiment text to subplot
    if any([v == 0,v == 3,v == 6]):
        ax.annotate(r'\textbf{%s}' % varnamesn[v],xy=(0,0),xytext=(-0.18,0.5),
                     textcoords='axes fraction',color='k',
                     fontsize=16,rotation=90,ha='center',va='center')
    if any([v == 0,v == 1,v == 2]):
        ax.annotate(r'\textbf{%s}' % months[v],xy=(0,0),xytext=(0.5,1.12),
                     textcoords='axes fraction',color='k',
                     fontsize=13,rotation=0,ha='center',va='center')
        
    ax.annotate(r'\textbf{[%s]}' % letters[v],xy=(0,0),
            xytext=(0.92,0.9),xycoords='axes fraction',
            color='dimgrey',fontsize=7)
        
    ax.set_aspect('equal')
            
    ###########################################################################
    if v == 2:
        cbar_ax = fig.add_axes([0.84,0.65,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=9,color='k') 
        elif varnamesn[v] == 'SIT':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.2)  
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=9,color='k')     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 5:
        cbar_ax = fig.add_axes([0.84,0.40,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=9,color='k') 
        elif varnamesn[v] == 'SIT':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=4.5)  
        elif varnamesn[v] == 'QNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=9,color='k')       
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif v == 8:
        cbar_ax = fig.add_axes([0.84,0.15,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnamesn[v] == 'SIC':
            cbar.set_label(r'\textbf{\%}',fontsize=9,color='k') 
        elif varnamesn[v] == 'SIT':
            cbar.set_label(r'\textbf{m}',fontsize=9,color='k',labelpad=1.2)  
        elif varnamesn[v] == 'RNET':
            cbar.set_label(r'\textbf{W/m$^{2}$}',fontsize=9,color='k')      
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
fig.subplots_adjust(wspace=-0.4,hspace=0.001)
       
plt.savefig(directoryfigure + 'SeaIce-HeatFlux.png',dpi=900)
print('Completed: Script done!')