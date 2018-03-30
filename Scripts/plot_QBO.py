"""
Script plots the monthly zonal wind in the equatorial stratosphere from 
-5S to 5N from an experiment in the WACCM4 simulations

Notes
-----
    Author : Zachary Labe
    Date   : 30 March 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_MonthlyOutput_AllMembers as MO
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_total/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO profile - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['U']
runnames = [r'HIT']
qbophase = ['pos','non','neg']
v=0
### Call function for surface temperature data from reach run
lat,lon,time,lev,uhitq= MO.readExperiAll('%s' % varnames[v],'HIT',
                                           'profile')
uhitq = np.asarray(uhitq)
uhitq[np.where(uhitq < -999)] = np.nan

### Create 2d array of latitude and longitude
monthq = [0,1,2,8,9,10,11]

latq = np.where((lat>=-5) & (lat<=5))[0]
uhit = np.nanmean(uhitq[:,:,:,latq,:],axis=3)
uhitzonal = np.nanmean(uhit[:,monthq,:,:],axis=3)

### Calculate anomalies
meanu = np.nanmean(uhitzonal,axis=0)
anom = uhitzonal - meanu

ustack = np.reshape(anom,(uhitzonal.shape[0]*uhitzonal.shape[1],
                               uhitzonal.shape[2]))

lon2,lat2 = np.meshgrid(lon,lat[latq])

ustack[np.where(ustack<-40)] = np.nan

###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
if varnames[v] == 'U':
    limit = np.arange(-20,20.1,0.1)
    barlim = np.arange(-20,30,10)
    
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
timeq,levq = np.meshgrid(np.arange(ustack.shape[0]),lev)

fig = plt.figure()

ax1 = plt.subplot(111)

var = ustack

ax1.spines['top'].set_color('dimgrey')
ax1.spines['right'].set_color('dimgrey')
ax1.spines['bottom'].set_color('dimgrey')
ax1.spines['left'].set_color('dimgrey')
ax1.spines['left'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)
ax1.spines['top'].set_linewidth(2)
ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                width=2,color='dimgrey')
ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                width=2,color='dimgrey')    
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')


cs = plt.contourf(timeq,levq,ustack.transpose(),limit,extend='both')

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

#plt.xticks(np.arange(0,96,15),map(str,np.arange(0,91,15)),fontsize=8)
plt.yticks(zscale,map(str,zscale),ha='right',fontsize=8)
plt.minorticks_off()

plt.ylabel(r'\textbf{Pressure (hPa)}',fontsize=13,color='dimgrey')
if varnames[v] == 'U':
    cmap = ncm.cmap('NCV_blu_red')            
    cs.set_cmap(cmap) 
    
plt.xlim([0,140])    
plt.ylim([100,10])

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

if varnames[v] == 'U':
    cbar.set_label(r'\textbf{[U] m/s}',fontsize=11,color='dimgray')
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(bottom=0.2)

plt.savefig(directoryfigure + 'QBO_HIT.png',dpi=300)
print('Completed: Script done!')
