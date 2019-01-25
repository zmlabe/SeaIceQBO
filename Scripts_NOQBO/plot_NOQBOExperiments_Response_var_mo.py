"""
Plot comparison for the NOQBO experiment

Notes
-----
    Author : Zachary Labe
    Date   : 17 January 2019
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
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting NOQBO Monthly Data - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 1941
years = np.arange(year1,year2+1,1)

period = 'DJF'

### Call arguments
def readVariables(varnames,runnames,period,directory):
    for v in range(len(varnames)):
        ### Call function for surface temperature data from reach run
        lat,lon,time,lev,varhit = MO.readExperi(directory,'%s' % varnames[v],
                                                '%s' % runnames,'surface')
        
        ### Create 2d array of latitude and longitude
        lon2,lat2 = np.meshgrid(lon,lat)
        
        ### Concatonate runs
        runs = [varhit]
        
        ### Separate per periods
        if period == 'ON': 
            tas_mo = np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = np.nanmean(runs[i][:,9:11,:,:],axis=1) 
        elif period == 'DJ':     
            tas_mo = np.empty((3,varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
                                                    lon,'surface',1) 
        elif period == 'FM':
            tas_mo= np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
        elif period == 'DJF':
            tas_mo= np.empty((varhit.shape[0]-1,varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:],tas_mo[:] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                      lon,'surface',1)   
        elif period == 'M':
            tas_mo= np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = runs[i][:,2,:,:]
        elif period == 'D':
            tas_mo= np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = runs[i][:,-1,:,:]
        elif period == 'N':
            tas_mo= np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = runs[i][:,-2,:,:]
        elif period == 'ND':
            tas_mo= np.empty((varhit.shape[0],varhit.shape[2],varhit.shape[3]))
            for i in range(len(runs)):
                tas_mo[:] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
        else:
            ValueError('Wrong period selected! (ON,DJ,FM)')
    
    return lat,lon,tas_mo
            
### Call functions
var1q = 'SST'
var2q = 'SIC'
lat,lon,var1r = readVariables([var1q],'NOQBO',period,directorydata)
lat,lon,var2r = readVariables([var2q],'NOQBO',period,directorydata)

lat,lon,var1c = readVariables([var1q],'CTLN',period,directorydata)
lat,lon,var2c = readVariables([var2q],'CTLN',period,directorydata)

### Fill in missing data
var1r[np.where(var1r < -1e10)] = np.nan
var2r[np.where(var2r < -1e10)] = np.nan
var2r[np.where(var2r > 1e10)] = np.nan
var1c[np.where(var1c < -1e10)] = np.nan
var2c[np.where(var2c < -1e10)] = np.nan
var2c[np.where(var2c > 1e10)] = np.nan

if var2q == 'RNET':
    var2c = var2c * -1
    var2r = var2r * -1

### Calculate difference 
var1a = var1r[:71,:,:] - var1c[:71,:,:] 
var1b = var1r[71:,:,:] - var1c[71:,:,:]

var2a = var2r[:71,:,:] - var2c[:71,:,:] 
var2b = var2r[71:,:,:] - var2c[71:,:,:]

### Calculate ensemble means
var1ma = np.nanmean(var1a,axis=0)
var1mb = np.nanmean(var1b,axis=0)
var2ma = np.nanmean(var2a,axis=0)
var2mb = np.nanmean(var2b,axis=0)

### Calculate statistical significance
statvar1,pvaluevar1 = UT.calc_indttest(var1a,var1b)
statvar2,pvaluevar2 = UT.calc_indttest(var2a,var2b)

###########################################################################
###########################################################################
###########################################################################
### Plot variable data 
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

if var1q == 'Z500':
    limitvar1c = np.arange(-50,51,5)
    barlimvar1c = np.arange(-50,51,25)
elif var1q == 'U10':
    limitvar1c = np.arange(-6,6.1,0.5)
    barlimvar1c = np.arange(-6,6.1,2)
elif var1q == 'U300':
    limitvar1c = np.arange(-6,6.1,0.25)
    barlimvar1c = np.arange(-6,6.1,2)
elif var1q == 'T2M':
    limitvar1c = np.arange(-3,3.1,0.25)
    barlimvar1c = np.arange(-3,3.1,1)
elif var1q == 'WAFZ850':
    limitvar1c = np.arange(-0.1,0.101,0.001)
    barlimvar1c = np.arange(-0.1,0.11,0.1)
elif var1q == 'P':
    limitvar1c = np.arange(-2,2.1,0.05)
    barlimvar1c = np.arange(-2,3,1)
elif var1q == 'SST':
    limitvar1c = np.arange(-10,10.1,0.25)
    barlimvar1c = np.arange(-10,11,5)
    
if var2q == 'SLP':
    limitvar2c = np.arange(-4,5,0.5)
    barlimvar2c = np.arange(-4,5,2)
elif var2q == 'Z30':
    limitvar2c = np.arange(-100,100.1,5)
    barlimvar2c = np.arange(-100,101,50)
elif var2q == 'Z500':
    limitvar2c = np.arange(-50,51,5)
    barlimvar2c = np.arange(-50,51,25)
elif var2q == 'RNET':
    limitvar2c = np.arange(-100,101,5)
    barlimvar2c = np.arange(-100,101,25)
elif var2q == 'WAFZ150':
    limitvar2c = np.arange(-0.01,0.0101,0.0001)
    barlimvar2c = np.arange(-0.01,0.011,0.01)
elif var2q == 'SWE':
    limitvar2c = np.arange(-25,25.1,1)
    barlimvar2c = np.arange(-25,26,25)
elif var2q == 'U700':
    limitvar2c = np.arange(-2,2.1,0.1)
    barlimvar2c = np.arange(-2,2.1,2)
elif var2q == 'SIC':
    limitvar2c = np.arange(-75,75.1,1)
    barlimvar2c = np.arange(-75,76,25)


ax = plt.subplot(2,3,1)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var1ma, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar1, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar1c,extend='both')
#cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
#                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar1c)
cbar.set_ticklabels(list(map(str,barlimvar1c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var1q == 'SST':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

ax.annotate(r'\textbf{1-70}',xy=(0,0),xytext=(0.5,1.15),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')
ax.annotate(r'\textbf{%s}' % var1q,xy=(0,0),xytext=(-0.1,0.5),
                         textcoords='axes fraction',color='k',
                         fontsize=20,rotation=90,ha='center',va='center')

###############################################################################
ax = plt.subplot(2,3,2)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var1mb, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar1, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar1c,extend='both')
#cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
#                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar1c)
cbar.set_ticklabels(list(map(str,barlimvar1c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var1q == 'SST':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

ax.annotate(r'\textbf{71-140}',xy=(0,0),xytext=(0.5,1.15),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')

###############################################################################
ax = plt.subplot(2,3,3)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var1ma-var1mb, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar1, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar1c,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar1c)
cbar.set_ticklabels(list(map(str,barlimvar1c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var1q == 'SST':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

ax.annotate(r'\textbf{Difference}',xy=(0,0),xytext=(0.5,1.15),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')

###############################################################################

ax = plt.subplot(2,3,4)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var2ma, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar2, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar2c,extend='both')
#cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
#                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar2c)
cbar.set_ticklabels(list(map(str,barlimvar2c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var2q == 'SIC':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

ax.annotate(r'\textbf{%s}' % var2q,xy=(0,0),xytext=(-0.1,0.5),
                         textcoords='axes fraction',color='k',
                         fontsize=20,rotation=90,ha='center',va='center')

###############################################################################
ax = plt.subplot(2,3,5)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var2mb, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar2, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar2c,extend='both')
#cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
#                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar2c)
cbar.set_ticklabels(list(map(str,barlimvar2c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var2q == 'SIC':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

###############################################################################
ax = plt.subplot(2,3,6)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var2ma-var2mb, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvaluevar2, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limitvar2c,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap(cmocean.cm.balance)
cbar = plt.colorbar(cs,fraction=0.035)
cbar.set_ticks(barlimvar2c)
cbar.set_ticklabels(list(map(str,barlimvar2c)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)

if var2q == 'SIC':
    m.fillcontinents(color='dimgrey')
    m.drawcoastlines(color='darkgray',linewidth=0.3)

plt.subplots_adjust(wspace=0.15)
plt.subplots_adjust(hspace=0.01)

plt.savefig(directoryfigure + '/NOQBO/Response/NOQBO_Response_%s_%s-%s.png' % (period,var1q,var2q),dpi=300)
