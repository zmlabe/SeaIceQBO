"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script 
plots the refractive index!

Notes
-----
    Author : Zachary Labe
    Date   : 27 July 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
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
print('\n' '----Plotting Refractive Index - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FICT--HIT}']

### Call functions for variable profile data for refractive index
def readProb(file,varid,experi):
    """
    Function reads daily data from WACCM4 control for probability of 
    wave propogation (refractive index)

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
    print('\n>>> Using readProb function! \n')
    
    ### Directorys for ensemble members 1-200 (remote server - Surtsey/Green)
    directorydata1 = '/surtsey/zlabe/simu/'
    directorydata2 = '/home/zlabe/green/simu/'
    
    ### Number of ensembles (years)
    ENS1 = 100
       
    ### Call files
    var1 = np.empty((ENS1,7,17,96))
    var2 = np.empty((ENS1,7,17,96))
 
    ### Call files for directory 1 (1-100 members)       
    for i in range(0,ENS1,1):       
        totaldirectory1 = directorydata1 + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
        filename1 = totaldirectory1 + file + '_%s.nc' % (i+1)
        
        ### Read in Data
        data1 = Dataset(filename1,'r')
        var1[i,:,:,:] = data1.variables['%s' % (varid)][:,:,:]
        lev = data1.variables['level'][:]
        data1.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+1,varid))
    print('Completed: Read members 1-100!') 

    ### Call files for directory 1 (101-200 members)       
    for i in range(0,ENS1,1):        
        totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                        '%s/' % (i+101)
        filename2 = totaldirectory2 + file + '_%s.nc' % (i+101)
        
        ### Read in Data
        data2 = Dataset(filename2,'r')
        var2[i,:,:,:] = data2.variables['%s' % (varid)][:,:,:]
        lat = data2.variables['latitude'][:]
        data2.close()
  
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+101,varid))
    print('Completed: Read members 101-200!')  
    
    ### Append all data    
    var = np.append(var1,var2,axis=0)
    
    print('\n*Completed: Finished readProb function!')
    return lat,lev,var

lat,lev,varfict = readProb('wave2prob','REFRACTIVE_PROB_MONTHLY','FICT')
lat,lev,varhit = readProb('wave2prob','REFRACTIVE_PROB_MONTHLY','HIT')

### Create 2d array of latitude and longitude
lev2,lat2 = np.meshgrid(lev,lat)

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

vardiff = np.nanmean(var_mofictneg[:,3,:,:] - var_mohitpos[:-4,3,:,:],axis=0)
var = np.nanmean(var_mofictneg[:,3,:,:],axis=0)

############################################################################
############################################################################
############################################################################
##### Plot refractive index
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax1 = plt.subplot(111)

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

zscale = np.array([1000,700,500,300,200,100,50,30,10])
limit = np.arange(0,101,5)
barlim = np.arange(0,101,50)

cs = plt.contourf(lat2,lev2,var.transpose(),limit)
cs.set_cmap('cubehelix')

cbar = plt.colorbar(cs)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
plt.minorticks_off()
plt.xticks(np.arange(0,91,15),map(str,np.arange(0,91,15)),fontsize=6)
plt.xlim([0,90])
plt.ylim([1000,10])

plt.savefig(directoryfigure + 'refractiveIndexTest.png',dpi=300)

###############################################################################

fig = plt.figure()
ax1 = plt.subplot(111)

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

zscale = np.array([1000,700,500,300,200,100,50,30,10])
limit = np.arange(-20,21,1)
barlim = np.arange(-20,21,10)

cs = plt.contourf(lat2,lev2,vardiff.transpose(),limit,extend='both')
cs.set_cmap(cmocean.cm.balance)

cbar = plt.colorbar(cs)
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')

plt.gca().invert_yaxis()
plt.yscale('log',nonposy='clip')

plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
plt.minorticks_off()
plt.xticks(np.arange(0,91,15),map(str,np.arange(0,91,15)),fontsize=6)
plt.xlim([0,90])
plt.ylim([1000,10])

plt.savefig(directoryfigure + 'refractiveIndexTest_Wave2_diff_FICTHIT_NEGPOS.png',dpi=300)
print('Completed: Script done!')    