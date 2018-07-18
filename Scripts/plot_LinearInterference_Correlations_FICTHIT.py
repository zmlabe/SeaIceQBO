"""
Plot computes a rolling pattern correlation between the climatological
and forced waves (1, 2, and all)

Notes
-----
    Author : Zachary Labe
    Date   : 18 July 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
import read_DailyOutput_AllMembers as DO
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
print('\n' '----Plotting GEOPxwave Correlations - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
varnames = 'U'
runnames = [r'\textbf{HIT',r'\textbf{FICT}']
qbophase = ['pos','non','neg']

### Function to read in control
def readCorrs(varid,wave,experi,lev):
    """
    Function reads daily data from WACCM4 control for rolling correlations
    of the climatological waves

    Parameters
    ----------
    varid : string
        variable name to read (LINT60N)
    wave : string
        which wave (1, 2, or all)
    experi : string
        experiment name (FICT)
    lev : string
        (all, tropo, strato)
        

    Returns
    -------
    lon : 1d numpy array
        longitudes
    var : 2d numpy array
        [time,longitude]

    Usage
    -----
    lon,var = readCorrs(varid,wave,experi,lev)
    """
    print('\n>>> Using readCorrs function! \n')
    
    ### Directorys for ensemble members 1-200 (remote server - Surtsey/Green)
    directorydata1 = '/surtsey/zlabe/simu/'
    directorydata2 = '/home/zlabe/green/simu/'
    
    ### Number of ensembles (years)
    ENS1 = 100
       
    ### Call files
    var1 = np.empty((ENS1,212,144))
    var2 = np.empty((ENS1,212,144))
 
    ### Call files for directory 1 (1-100 members)       
    for i in range(0,ENS1,1):       
        totaldirectory1 = directorydata1 + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
        filename1 = totaldirectory1 + varid + '_%s.nc' % (i+1)
        
        ### Read in Data
        data1 = Dataset(filename1,'r')
        var1[i,:,:] = data1.variables['cor_%s_wave%s' % (lev,wave)][:]
        data1.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+1,varid))
    print('Completed: Read members 1-100!') 
    
    ### Read longitude data
    datalon = Dataset(totaldirectory1 + 'Z500_%s.nc' % (i+1),'r')
    lon = datalon.variables['longitude'][:]
    datalon.close()
    
    ### Call files for directory 1 (101-200 members)       
    for i in range(0,ENS1,1):        
        totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                        '%s/' % (i+101)
        filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        
        ### Read in Data
        data2 = Dataset(filename2,'r')
        var2[i,:,:] = data2.variables['cor_%s_wave%s' % (lev,wave)][:]
        data2.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+101,varid))
    print('Completed: Read members 101-200!')  
    
    ### Append all data    
    var = np.append(var1,var2,axis=0)
    
    print('\n*Completed: Finished readCorrs function!')
    return lon,var

### Read in data
lon,all1 = readCorrs('LINT60N','1','FICT','all')
lon,all2 = readCorrs('LINT60N','2','FICT','all')
lon,allall = readCorrs('LINT60N','_all','FICT','all')