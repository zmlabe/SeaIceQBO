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

### Create 2d array of latitude and longitude
lon2,lat2 = np.meshgrid(lon,lat)

### Read in QBO phases    
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

### Calculate ensemble mean
varum_pos = np.nanmean(varu_pos,axis=0)
varum_neg = np.nanmean(varu_neg,axis=0)

varvm_pos = np.nanmean(varv_pos,axis=0)
varvm_neg = np.nanmean(varv_neg,axis=0)

### Calculate velocity 
uwnd_pos, uwndinfo_pos = prep_data(varum_pos, 'tyx')
vwnd_pos, vwndinfo_pos = prep_data(varvm_pos, 'tyx')

uwnd_neg, uwndinfo_neg = prep_data(varum_neg, 'tyx')
vwnd_neg, vwndinfo_neg = prep_data(varvm_neg, 'tyx')

latn, uwnd_pos, vwnd_pos = order_latdim(lat, uwnd_pos, vwnd_pos)
latn, uwnd_neg, vwnd_neg = order_latdim(lat, uwnd_neg, vwnd_neg)

### Calculate VelocityWind object
wpos = VectorWind(uwnd_pos, vwnd_pos)
wneg = VectorWind(uwnd_neg, vwnd_neg)

### Calculate Absolute Velocity
avrt_pos = wpos.absolutevorticity()
avrt_neg = wneg.absolutevorticity()

### Change dimensions
apos = recover_data(avrt_pos,uwndinfo_pos)
aneg = recover_data(avrt_neg,uwndinfo_neg)

### Calculate dy 
dyapos = np.gradient(apos,axis=2)
dyaneg = np.gradient(aneg,axis=2)

### Average over time
dyapos_t = np.nanmean(dyapos[90:120,:,:],axis=0)
dyaneg_t = np.nanmean(dyaneg[90:120,:,:],axis=0)