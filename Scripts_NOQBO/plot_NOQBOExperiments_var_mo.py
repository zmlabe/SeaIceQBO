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
                                                'NOQBO','surface')
        
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
var1q = 'T2M'
var2q = 'SLP'
lat,lon,var1 = readVariables([var1q],'NOQBO',period,directorydata)
lat,lon,var2 = readVariables([var2q],'NOQBO',period,directorydata)

var1[np.where(var1 < -1e10)] = np.nan
var2[np.where(var2 < -1e10)] = np.nan

### Slice time period in half
var1a = var1[:71,:,:]
var1b = var1[71:,:,:]
var2a = var2[:71,:,:]
var2b = var2[71:,:,:]

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
    limitvar1 = np.arange(5000,6105,5)
    barlimvar1 = np.arange(5000,6105,200)
    limitvar1c = np.arange(-50,51,1)
    barlimvar1c = np.arange(-50,51,25)
elif var1q == 'U10':
    limitvar1 = np.arange(-50,51,1)
    barlimvar1 = np.arange(-50,51,10)
    limitvar1c = np.arange(-6,6.1,0.1)
    barlimvar1c = np.arange(-6,6.1,2)
elif var1q == 'U300':
    limitvar1 = np.arange(-50,51,1)
    barlimvar1 = np.arange(-50,51,10)
    limitvar1c = np.arange(-6,6.1,0.1)
    barlimvar1c = np.arange(-6,6.1,2)
elif var1q == 'T2M':
    limitvar1 = np.arange(-30,21,1)
    barlimvar1 = np.arange(-20,20.1,10)
    limitvar1c = np.arange(-3,3.1,0.1)
    barlimvar1c = np.arange(-3,3.1,1)
    
if var2q == 'SLP':
    limitvar2 = np.arange(990,1031,1)
    barlimvar2 = np.arange(990,1031,10)
    limitvar2c = np.arange(-4,5,0.1)
    barlimvar2c = np.arange(-4,5,2)
elif var2q == 'Z30':
    limitvar2 = np.arange(21900,23500,5)
    barlimvar2 = np.arange(21900,23500,250)
    limitvar2c = np.arange(-100,100.1,1)
    barlimvar2c = np.arange(-100,101,50)
elif var2q == 'Z500':
    limitvar2 = np.arange(5000,6105,5)
    barlimvar2 = np.arange(5000,6105,200)
    limitvar2c = np.arange(-50,51,1)
    barlimvar2c = np.arange(-50,51,25)

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

cs = m.contourf(x,y,var,limitvar1,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap('cubehelix')
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar1)
cbar.set_ticklabels(list(map(str,barlimvar1)))

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

cs = m.contourf(x,y,var,limitvar1,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap('cubehelix')
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar1)
cbar.set_ticklabels(list(map(str,barlimvar1)))

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
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar1c)
cbar.set_ticklabels(list(map(str,barlimvar1c)))

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

cs = m.contourf(x,y,var,limitvar2,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap('cubehelix')
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar2)
cbar.set_ticklabels(list(map(str,barlimvar2)))

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

cs = m.contourf(x,y,var,limitvar2,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.5)

cs.set_cmap('cubehelix')
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar2)
cbar.set_ticklabels(list(map(str,barlimvar2)))

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
cbar = plt.colorbar(cs)
cbar.set_ticks(barlimvar2c)
cbar.set_ticklabels(list(map(str,barlimvar2c)))

plt.savefig(directoryfigure + '/NOQBO/NOQBO_%s_%s-%s.png' % (period,var1q,var2q),dpi=300)
    
#    ###########################################################################
#    ###########################################################################
#    ###########################################################################
#    ### Plot variable data for QBO composites
#    plt.rc('text',usetex=True)
#    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#    
#    ### Set limits for contours and colorbars
#    if varnames[v] == 'T2M':
#        limit = np.arange(-10,10.1,0.5)
#        barlim = np.arange(-10,11,5)
#    elif varnames[v] == 'SLP':
#        limit = np.arange(-6,6.1,0.5)
#        barlim = np.arange(-6,7,3)
#    elif varnames[v] == 'Z500':
#        limit = np.arange(-60,60.1,5)
#        barlim = np.arange(-60,61,30) 
#    elif varnames[v] == 'Z30':
#        limit = np.arange(-100,100.1,10)
#        barlim = np.arange(-100,101,50) 
#    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
#        limit = np.arange(-5,5.1,0.5)
#        barlim = np.arange(-5,6,1)
#    elif varnames[v] == 'SWE':
#        limit = np.arange(-25,25.1,1)
#        barlim = np.arange(-25,26,25)
#    elif varnames[v] == 'P':
#        limit = np.arange(-2,2.1,0.05)
#        barlim = np.arange(-2,3,1) 
#    elif varnames[v] == 'THICK':
#        limit = np.arange(-60,60.1,3)
#        barlim = np.arange(-60,61,30)
#    elif varnames[v] == 'EGR':
#        limit = np.arange(-0.2,0.21,0.02)
#        barlim = np.arange(-0.2,0.3,0.2)
#    elif varnames[v] == 'WAFZ850':
#        limit = np.arange(-0.1,0.101,0.001)
#        barlim = np.arange(-0.1,0.11,0.1)
#    elif varnames[v] == 'WAFZ150':
#        limit = np.arange(-0.01,0.0101,0.0001)
#        barlim = np.arange(-0.01,0.011,0.01)
#    elif varnames[v] == 'WAFY850':
#        limit = np.arange(-5,5.1,0.1)
#        barlim = np.arange(-5,6,5)
#    elif varnames[v] == 'WAFY150':
#        limit = np.arange(-2,2.01,0.05)
#        barlim = np.arange(-2,3,2)
#    
#    fig = plt.figure()
#    for i in range(len(diffruns_mo)):
#        var = diffruns_mo[i]
#        pvar = pruns_mo[i]
#        
#        ax1 = plt.subplot(2,3,i+1)
#        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
#                    area_thresh=10000.)
#        
#        var, lons_cyclic = addcyclic(var, lon)
#        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
#        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
#        x, y = m(lon2d, lat2d)
#        
#        pvar,lons_cyclic = addcyclic(pvar, lon)
#        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
#        climoq,lons_cyclic = addcyclic(climo[i], lon)
#        climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
#                  
#        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
#        
#        cs = m.contourf(x,y,var,limit,extend='both')
#        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
#                         linewidths=0.4)
#        if varnames[v] == 'Z30': # the interval is 250 m 
#            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
#                            colors='k',linewidths=1.5,zorder=10)
#
#        m.drawcoastlines(color='dimgray',linewidth=0.8)
#        
#        if varnames[v] == 'T2M':
#            cmap = ncm.cmap('NCV_blu_red')           
#            cs.set_cmap(cmap)   
#        elif varnames[v] == 'SLP':
#            cmap = cmocean.cm.balance          
#            cs.set_cmap(cmap)   
#        elif varnames[v] == 'Z500':
#            cmap = cmocean.cm.balance           
#            cs.set_cmap(cmap)  
#        elif varnames[v] == 'Z30':
#            cmap = cmocean.cm.balance  
#            cs.set_cmap(cmap)  
#        elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
#            cmap = ncm.cmap('NCV_blu_red')            
#            cs.set_cmap(cmap)  
#        elif varnames[v] == 'SWE':
#            cmap = cmap = cmocean.cm.balance
#            cs.set_cmap(cmap)
#        elif varnames[v] == 'P':
#            cmap = ncm.cmap('precip4_diff_19lev')            
#            cs.set_cmap(cmap) 
#        elif varnames[v] == 'THICK':
#            cmap = ncm.cmap('NCV_blu_red')           
#            cs.set_cmap(cmap) 
#        elif varnames[v] == 'EGR':
#            cmap = cmocean.cm.curl
#            cs.set_cmap(cmap)
#        elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
#            cmap = cmocean.cm.curl
#            cs.set_cmap(cmap)
#        elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
#            cmap = cmocean.cm.curl
#            cs.set_cmap(cmap)
#                    
#        ### Add experiment text to subplot
#        if i < 3:
#            qbophaseq = [r'QBO-W',r'QBO-N',r'QBO-E']
#            ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
#                         textcoords='axes fraction',color='dimgray',
#                         fontsize=13,rotation=0,ha='center',va='center')
#        if i == 0 or i == 3:
#            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(-0.1,0.5),
#                         textcoords='axes fraction',color='k',
#                         fontsize=20,rotation=90,ha='center',va='center')
#                
#    ###########################################################################
#    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
#    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                        extend='max',extendfrac=0.07,drawedges=False)
#    
#    if varnames[v] == 'T2M':
#        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
#    elif varnames[v] == 'Z500':
#        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
#    elif varnames[v] == 'Z30':
#        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
#    elif varnames[v] == 'SLP':
#        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='dimgray')  
#    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
#        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
#    elif varnames[v] == 'SWE':
#        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
#    elif varnames[v] == 'P':
#        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
#    elif varnames[v] == 'THICK':
#        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
#    elif varnames[v] == 'EGR':
#        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
#    elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
#        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')
#    elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
#        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')
#
#    cbar.set_ticks(barlim)
#    cbar.set_ticklabels(list(map(str,barlim)))
#    cbar.ax.tick_params(axis='x', size=.01)
#    cbar.outline.set_edgecolor('dimgrey')
#    
#    plt.subplots_adjust(wspace=0.01)
#    plt.subplots_adjust(hspace=0.01)
#    plt.subplots_adjust(bottom=0.15)
#    
#    plt.savefig(directoryfigure + '/NOQBO/NOQBO_%s_%s.png' % (period,varnames[v]),dpi=300)
#
#print('Completed: Script done!')
#
