"""
Plot compares FSUB-CIT and FPOL-CIT for the net surface heat flux from
October through March over all 200 ensemble members

Notes
-----
    Author : Zachary Labe
    Date   : 1 June 2018
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
print('\n' '----Plotting Heat Flux Comparison Regional - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = 'RNET'
runnames = [r'CIT',r'FSUB',r'FPOL']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FPOL--CIT}']

### Call function for rnet data from reach run
lat,lon,time,lev,tascit = MO.readExperiAll('%s' % varnames,'CIT',
                                           'surface')
lat,lon,time,lev,tasfsub = MOR.readExperiAllRegional('%s' % varnames,'FSUB',
                                           'surface')
lat,lon,time,lev,tasfpol = MOR.readExperiAllRegional('%s' % varnames,'FPOL',
                                            'surface')

### Create 2d array of latitude and longitude
lon2,lat2 = np.meshgrid(lon,lat)

### Calculate Oct-Mar
tascitmo = np.append(tascit[:,9:,:,:],tascit[:,:3,:,:],axis=1)
tasfsubmo = np.append(tasfsub[:,9:,:,:],tasfsub[:,:3,:,:],axis=1)
tasfpolmo = np.append(tasfpol[:,9:,:,:],tasfpol[:,:3,:,:],axis=1)

### Take difference
fsubcit = np.nanmean(tasfsubmo - tascitmo,axis=0)
fpolcit = np.nanmean(tasfpolmo - tascitmo,axis=0)

monthall = np.append(fsubcit,fpolcit,axis=0)

###########################################################################
###########################################################################
###########################################################################
### Plot heat flux
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

limit = np.arange(-75,75.1,1)
barlim = np.arange(-75,76,75)
months = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',
           r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']

fig = plt.figure(figsize=(8,3.5))
for i in range(len(monthall)):
        
    var = monthall[i]*-1

    ax1 = plt.subplot(2,6,i+1)
    m = Basemap(projection='npstere',boundinglat=51,lon_0=270,resolution='l',
                round =True,area_thresh=10000)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    m.drawcoastlines(color='darkgrey',linewidth=0.1)
    m.fillcontinents(color='dimgray')
        
    cs = m.contourf(x,y,var,limit,extend='both')
    
    cmap = ncm.cmap('NCV_blu_red')           
    cs.set_cmap(cmap) 
    
    ### Add experiment text to subplot
    ax1.annotate(r'\textbf{%s}' % months[i],xy=(0,0),xytext=(0.865,0.90),
                 textcoords='axes fraction',color='k',fontsize=11,
                 rotation=320,ha='center',va='center')
    
    if i==0:
        ax1.annotate(r'\textbf{$\Delta$Pacific}',xy=(0,0),xytext=(-0.15,0.5),
             textcoords='axes fraction',color='k',fontsize=12,
             rotation=90,ha='center',va='center')
    elif i==6:
        ax1.annotate(r'\textbf{$\Delta$Atlantic}',xy=(0,0),xytext=(-0.15,0.5),
             textcoords='axes fraction',color='k',fontsize=12,
             rotation=90,ha='center',va='center')

    ### Add concentration colorbar
    cbar_ax = fig.add_axes([0.412,0.13,0.2,0.033])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        drawedges=False)
    cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=12,
                              color='dimgray',labelpad=3)
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.001,labelsize=9)
    cbar.outline.set_edgecolor('dimgray')
    
plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0)
        
plt.savefig(directoryfigure + 'monthly_netSurfaceHeatFlux_regional.png',dpi=900)