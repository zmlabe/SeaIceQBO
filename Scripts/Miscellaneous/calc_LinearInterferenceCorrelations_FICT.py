"""
Plot comparisons between WACCM4 sea ice experiments. Plot computes a 
rolling pattern correlation between the climatologicaland forced waves 
(1, 2, and all)

Notes
-----
    Author : Zachary Labe
    Date   : 20 July 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
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
print('\n' '----Calculating Linear Interference Correlations - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
vari = ['GEOPx-60N-wave1']
varidvar = ['GEOPxwave1']
levelq = ['tropo']
runnames = [r'HIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FICT--HIT}']

### Function to read in control
def readCorrs(file,varid,experi):
    """
    Function reads daily data from WACCM4 control for rolling correlations
    of the climatological waves

    Parameters
    ----------
    file : string
        name of file variable
    varid : string
        variable name to read (LINT60N)
    experi : string
        experiment name (FICT)
        

    Returns
    -------
    lon : 1d numpy array
        longitudes
    lev : 1d numpy array
        vertical levels
    var : 2d numpy array
        [time,level,latitude,longitude]

    Usage
    -----
    lon,lev,var = readCorrs(file,varid,experi)
    """
    print('\n>>> Using readCorrs function! \n')
    
    ### Directorys for ensemble members 1-200 (remote server - Surtsey/Green)
    directorydata1 = '/surtsey/zlabe/simu/'
    directorydata2 = '/home/zlabe/green/simu/'
    
    ### Number of ensembles (years)
    ENS1 = 100
       
    ### Call files
    var1 = np.empty((ENS1,212,17,1,145))
    var2 = np.empty((ENS1,212,17,1,145))
 
    ### Call files for directory 1 (1-100 members)       
    for i in range(0,ENS1,1):       
        totaldirectory1 = directorydata1 + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
        filename1 = totaldirectory1 + file + '_%s.nc' % (i+1)
        
        ### Read in Data
        data1 = Dataset(filename1,'r')
        var1[i,:,:,:,:] = data1.variables['%s' % (varid)][:,:,:,:]
        lev = data1.variables['level'][:]
        data1.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+1,varid))
    print('Completed: Read members 1-100!') 
    
    ### Read longitude data
    datalon = Dataset(totaldirectory1 + '%s_%s.nc' % (varid,i+1),'r')
    lon = datalon.variables['longitude'][:]
    datalon.close()
    
    lon[-1] = 360.
    
    ### Call files for directory 1 (101-200 members)       
    for i in range(0,ENS1,1):        
        totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                        '%s/' % (i+101)
        filename2 = totaldirectory2 + file + '_%s.nc' % (i+101)
        
        ### Read in Data
        data2 = Dataset(filename2,'r')
        var2[i,:,:,:,:] = data2.variables['%s' % (varid)][:,:,:,:]
        data2.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+101,varid))
    print('Completed: Read members 101-200!')  
    
    ### Append all data    
    var = np.append(var1,var2,axis=0)
    
    print('\n*Completed: Finished readCorrs function!')
    return lon,lev,var

lon,lev,varfict = readCorrs(vari[0],varidvar[0],'FICT')
lon,lev,varhit = readCorrs(vari[0],varidvar[0],'HIT')
    
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
### Concatonate runs
var_mo = [varhit,varfict]

### Composite by QBO phase    
var_mohitpos = var_mo[0][pos_hit].squeeze()
var_mohitnon = var_mo[0][non_hit].squeeze()
var_mohitneg = var_mo[0][neg_hit].squeeze()
var_mofictpos = var_mo[1][pos_fict].squeeze() 
var_mofictnon = var_mo[1][non_fict].squeeze()
var_mofictneg = var_mo[1][neg_fict].squeeze()

### Compute differences
fictpos = np.nanmean(var_mofictpos - var_mohitpos,axis=0)
fictnon = np.nanmean(var_mofictnon - var_mohitnon,axis=0)
fictneg = np.nanmean(var_mofictneg - var_mohitneg,axis=0)

diffruns = [fictpos,fictnon,fictneg]

### Historical means
hitpos = np.nanmean(var_mohitpos,axis=0)
hitnon = np.nanmean(var_mohitnon,axis=0)
hitneg = np.nanmean(var_mohitneg,axis=0)

### Calculate Correlations
days = 212
corpos = np.empty((days))
cornon = np.empty((days))
corneg = np.empty((days))
for i in range(days):
    corpos[i] = UT.calc_spatialCorrHeightLev(hitpos[i],fictpos[i],lev,lon,
          'yes',levelq[0])
    cornon[i] = UT.calc_spatialCorrHeightLev(hitnon[i],fictnon[i],lev,lon,
          'yes',levelq[0])
    corneg[i] = UT.calc_spatialCorrHeightLev(hitneg[i],fictneg[i],lev,lon,
          'yes',levelq[0])

corrs = np.array([corpos,cornon,corneg])
  
### Save as text files  
np.savetxt(directorydata3 + 'LINT_%s_%s.txt' % (vari[0],levelq[0]),
           np.round(corrs.transpose(),1),delimiter=',',
           fmt='%3.1f',header='  '.join(qbophase)+'\n',
           footer='\n File contains pattern correlations of linear' \
           'interference for each QBO phase')