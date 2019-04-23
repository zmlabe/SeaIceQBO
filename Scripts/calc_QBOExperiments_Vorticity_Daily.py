"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script 
calculates the gradient of vorticity at 30 hPa

Notes
-----
    Author : Zachary Labe
    Date   : 25 July 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import cmocean
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from mpl_toolkits.basemap import Basemap

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
print('\n' '----Calculating Vorticity for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = ['U','V']
runnames = [r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FICT--HIT}']

### Call functions for variable profile data for polar cap
lat,lon,time,lev,varuu = DO.readMeanExperiAll('U','FICT','profile2')
lat,lon,time,lev,varvv = DO.readMeanExperiAll('V','FICT','profile2')

lat,lon,time,lev,hituu = DO.readMeanExperiAll('U','HIT','profile2')
lat,lon,time,lev,hitvv = DO.readMeanExperiAll('V','HIT','profile2')

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
                    np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)   
  
filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)

### Slice by QBO phases
varu_pos = varuu[pos_fict,:,:]
varu_neg = varuu[neg_fict,:,:]

varv_pos = varvv[pos_fict,:,:]
varv_neg = varvv[neg_fict,:,:]

hitu_pos = hituu[pos_fict,:,:]
hitu_neg = hituu[neg_fict,:,:]

hitv_pos = hitvv[pos_fict,:,:]
hitv_neg = hitvv[neg_fict,:,:]

### Delete variables to get memory
del varuu,varvv,hituu,hitvv

## Calculate ensemble mean over QBO phases
varum_pos = np.nanmean(varu_pos,axis=0)
varum_neg = np.nanmean(varu_neg,axis=0)

varvm_pos = np.nanmean(varv_pos,axis=0)
varvm_neg = np.nanmean(varv_neg,axis=0)

hitum_pos = np.nanmean(hitu_pos,axis=0)
hitum_neg = np.nanmean(hitu_neg,axis=0)

hitvm_pos = np.nanmean(hitv_pos,axis=0)
hitvm_neg = np.nanmean(hitv_neg,axis=0)

### Prepare data for calculating
uwnd_pos, uwndinfo_pos = prep_data(varum_pos, 'tyx')
vwnd_pos, vwndinfo_pos = prep_data(varvm_pos, 'tyx')

uwnd_neg, uwndinfo_neg = prep_data(varum_neg, 'tyx')
vwnd_neg, vwndinfo_neg = prep_data(varvm_neg, 'tyx')

latn, uwnd_pos, vwnd_pos = order_latdim(lat, uwnd_pos, vwnd_pos)
latn, uwnd_neg, vwnd_neg = order_latdim(lat, uwnd_neg, vwnd_neg)

hituwnd_pos, hituwndinfo_pos = prep_data(hitum_pos, 'tyx')
hitvwnd_pos, hitvwndinfo_pos = prep_data(hitvm_pos, 'tyx')

hituwnd_neg, hituwndinfo_neg = prep_data(hitum_neg, 'tyx')
hitvwnd_neg, hitvwndinfo_neg = prep_data(hitvm_neg, 'tyx')

latn, hituwnd_pos, hitvwnd_pos = order_latdim(lat, hituwnd_pos, hitvwnd_pos)
latn, hituwnd_neg, hitvwnd_neg = order_latdim(lat, hituwnd_neg, hitvwnd_neg)

### Change lat/lon to mesh
lon2n,lat2n = np.meshgrid(lon,latn)

### Calculate VelocityWind object
wpos = VectorWind(uwnd_pos, vwnd_pos)
wneg = VectorWind(uwnd_neg, vwnd_neg)

hitwpos = VectorWind(hituwnd_pos, hitvwnd_pos)
hitwneg = VectorWind(hituwnd_neg, hitvwnd_neg)

### Calculate Absolute Velocity
avrt_pos = wpos.absolutevorticity()
avrt_neg = wneg.absolutevorticity()

hitavrt_pos = hitwpos.absolutevorticity()
hitavrt_neg = hitwneg.absolutevorticity()

### Change dimensions
fictapos = recover_data(avrt_pos,uwndinfo_pos)
fictaneg = recover_data(avrt_neg,uwndinfo_neg)

hitapos = recover_data(hitavrt_pos,hituwndinfo_pos)
hitaneg = recover_data(hitavrt_neg,hituwndinfo_neg)

### Calculate dy 
fictdyapos = np.gradient(fictapos,axis=2)
fictdyaneg = np.gradient(fictaneg,axis=2)

hitdyapos = np.gradient(hitapos,axis=2)
hitdyaneg = np.gradient(hitaneg,axis=2)

### Average over time
timeq = np.arange(90,120)
fictdyapos_t = np.nanmean(fictdyapos[timeq,:,:],axis=0)
fictdyaneg_t = np.nanmean(fictdyaneg[timeq,:,:],axis=0)

hitdyapos_t = np.nanmean(hitdyapos[timeq,:,:],axis=0)
hitdyaneg_t = np.nanmean(hitdyaneg[timeq,:,:],axis=0)

### Calculate QBO difference
fictdiff = fictdyaneg_t - fictdyapos_t
hitdiff = hitdyaneg_t - hitdyapos_t

### Calculate forcing differences
ficthit_pos = fictdyapos_t - hitdyapos_t
ficthit_neg = fictdyaneg_t - hitdyaneg_t

diff = np.nanmean(fictaneg[timeq,:,:],axis=0) - np.nanmean(hitaneg[timeq,:,:],axis=0)

### Plot testing
m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l')
m.drawcoastlines(color='dimgrey')

cs = m.contourf(lon2n,lat2n,diff*1e5,77,latlon=True)
cs.set_cmap(cmocean.cm.balance)

plt.colorbar(cs)
plt.savefig(directoryfigure+'test.png',dpi=300)