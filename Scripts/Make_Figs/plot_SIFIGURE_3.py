"""
Plot manuscript figure of WAFz map at 850 hPa and 150 hPa

Notes
-----
    Author : Zachary Labe
    Date   : 11 December 2018
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
varnames = ['WAFZ850','WAFZ150']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'ND'

def readWAFz(directorydata,varnames,period):
    for v in range(1):
        ### Call function for surface temperature data from reach run
        lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
                                                   'surface')
        lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
                                                    'surface')
        
        ### Create 2d array of latitude and longitude
        lon2,lat2 = np.meshgrid(lon,lat)
        
        ### Read in QBO phases 
        filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
        filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
        filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
        filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
        pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                            np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
        neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                            np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
        
        filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
        filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
        filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
        filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
        pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                            np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
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
    
        ### Compute comparisons for months - taken ensemble average   
        ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
        ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
        diffruns_mo = ficthitneg
        
        ### Calculate significance for FM 
        stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[1][pos_fict,:,:],
                                                           tas_mo[0][pos_hit,:,:])
        stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[1][neg_fict,:,:],
                                                   tas_mo[0][neg_hit,:,:])
    
        pruns_mo = pvalue_FICTHITneg
        
        return diffruns_mo,pruns_mo,lat,lon
    
### Read in data    
data85,p85,lat,lon = readWAFz(directorydata,'WAFZ850',period)
data15,p15,lat,lon = readWAFz(directorydata,'WAFZ150',period)
data = [data85,data15]
pval = [p85,p15]
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for QBO composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit85 = np.arange(-0.08,0.081,0.001)
barlim85 = np.arange(-0.08,0.081,0.08)
limit15 = np.arange(-0.008,0.0081,0.0001)
barlim15 = np.arange(-0.008,0.0081,0.008)

fig = plt.figure()   
ax1 = plt.subplot(1,2,1)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(data[0], lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pval[0], lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs85 = m.contourf(x,y,var,limit85,extend='both')
cs85q = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)

m.drawcoastlines(color='dimgray',linewidth=0.8)
cmap = cmocean.cm.balance
cs85.set_cmap(cmap)

### Add experiment text to subplot
ax1.annotate(r'\textbf{850 hPa}',xy=(0,0),xytext=(0.5,1.08),
             textcoords='axes fraction',color='k',
             fontsize=17,rotation=0,ha='center',va='center')

###############################################################################
###############################################################################
###############################################################################
ax1 = plt.subplot(1,2,2)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(data[1], lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pval[1], lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs15 = m.contourf(x,y,var,limit15,extend='both')
cs15q = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)

m.drawcoastlines(color='dimgray',linewidth=0.8)
cmap = cmocean.cm.balance
cs15.set_cmap(cmap)
                   
### Add experiment text to subplot
ax1.annotate(r'\textbf{150 hPa}',xy=(0,0),xytext=(0.5,1.08),
             textcoords='axes fraction',color='k',
             fontsize=17,rotation=0,ha='center',va='center')
            
###########################################################################
cbar_ax = fig.add_axes([0.156,0.1,0.2,0.03])                
cbar = fig.colorbar(cs85,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim85)
cbar.set_ticklabels(list(map(str,barlim85)))
cbar.ax.tick_params(axis='x', size=.01,labelsize=8)
cbar.outline.set_edgecolor('dimgrey')

###########################################################################
cbar_ax = fig.add_axes([0.65,0.1,0.2,0.03])                
cbar = fig.colorbar(cs15,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim15)
cbar.set_ticklabels(list(map(str,barlim15)))
cbar.ax.tick_params(axis='x', size=.01,labelsize=8)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0.01)
plt.subplots_adjust(bottom=0.05)
plt.tight_layout()

plt.savefig(directoryfigure + 'WAFz_Map_ND_QBOE.png',dpi=900)

print('Completed: Script done!')

