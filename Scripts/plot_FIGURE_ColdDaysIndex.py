"""
Plot manuscript figure for cold days index (December)

Notes
-----
    Author : Zachary Labe
    Date   : 17 October 2018 
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directorydata3 = '/home/zlabe/Documents/Research/SeaIceQBO/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Cold Days Index - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = ['T1000']
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','neg']

### Read in CDI
data1 = Dataset(directorydata3 + 'ColdDayIndex_December.nc')
fictpos = data1.variables['cdi_fictpos'][:]
fictneg = data1.variables['cdi_fictneg'][:]
lat = data1.variables['lat'][:]
lon = data1.variables['lon'][:]
data1.close()
    
diffruns = [np.nanmean(fictpos,axis=0),np.nanmean(fictneg,axis=0)]
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for QBO composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
MASK = False

### Set limits for contours and colorbars
if varnames[0] == 'T1000':
    limit = np.arange(-20,21,1)
    barlim = np.arange(-20,21,10)
    
fig = plt.figure()
for i in range(len(diffruns)):
    var = diffruns[i]
    
    ax1 = plt.subplot(1,2,i+1)
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.fillcontinents(color='dimgrey',zorder=3)
    
    cs = m.contourf(x,y,var,limit,extend='both',zorder=4)

    m.drawcoastlines(color='dimgray',linewidth=0.8)
            
    cs.set_cmap(cmocean.cm.balance)   
        
    m.drawlsmask(land_color=(0,0,0,0),ocean_color='dimgrey',lakes=True,
                 resolution='c',zorder=5)
                
    ### Add experiment text to subplot
    if i < 2:
        qbophaseq = [r'QBO-W',r'QBO-E']
        ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                     textcoords='axes fraction',color='k',
                     fontsize=17,rotation=0,ha='center',va='center')
            
###########################################################################
cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{Cold Days Intensity}',fontsize=11,color='k',labelpad=1.4)  

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='x', size=.01,labelsize=8)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0.01)
plt.subplots_adjust(bottom=0.05)

plt.savefig(directoryfigure + 'ColdDaysIndex.png',
            dpi=900)

print('Completed: Script done!')
