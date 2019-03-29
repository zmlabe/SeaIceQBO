"""
Calculate linear interference for FICT-HIT

Notes
-----
    Author : Zachary Labe
    Date   : 29 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import datetime
import read_MonthlyOutput_AllMembers as MO
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_D_2/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting linear interference - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

historicalforced = []
futureforced = []
lonss = []
varnames = ['GEOPxwave_all','GEOPxwave1','GEOPxwave2']
period = 'D'
for v in range(len(varnames)):
    ### Call function for wave data from reach run
    lat,lon1,time,lev,varhit = MO.readExperiAll('%s' % varnames[v],
                                             'HIT','profile')
    lat,lon1,time,lev,varfict = MO.readExperiAll('%s' % varnames[v],
                                             'FICT','profile')
    
    ### Missing data
    varhit[np.where(varhit<-100000)]=np.nan

    ### Modify lons 
    lon1[-1] = 360

    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Concatonate runs
    runnames = [r'HIT',r'FICT']
    experiments = [r'\textbf{FICT--HIT}']
    runs = [varhit,varfict]
    
    ### Separate per periods (D)
    var_d = [varhit[:,-1,:,:,:],varfict[:,-1,:,:,:]]
    
    ### Take at 60N 
    latq = np.where((lat>59) & (lat<61))[0]
    
    ### Compute comparisons for M - taken ensemble average at 60N (index 79)
    diff_FICTHIT = np.nanmean(var_d[1] - var_d[0],axis=0)
    diffruns_d = diff_FICTHIT[:,latq,:].squeeze()
    historicalforcedq = np.nanmean(var_d[0][:,:,latq,:].squeeze(),axis=0)
    
    historicalforced.append(historicalforcedq)
    futureforced.append(diffruns_d)
    lonss.append(lon1)
    
###########################################################################
###########################################################################
###########################################################################
#### Plot climatological waves
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])

fig = plt.figure()
for i in range(3):
    ax1 = plt.subplot(3,1,i+1)
    
    ### Calculate correlations
    corr = UT.calc_spatialCorrHeight(historicalforced[i],futureforced[i],
                                lev,lonss[i],'yes')
    
    lonq,levq = np.meshgrid(lonss[i],lev)
    
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
            
    cs = plt.contourf(lonq,levq,futureforced[i],np.arange(-35,35.1,0.1),
                      extend='both') 
    if i==2:
        cs1 = plt.contour(lonq,levq,historicalforced[i],5,
                          colors='k',linewidths=1) 
    else:
        cs1 = plt.contour(lonq,levq,historicalforced[i],20,
                  colors='k',linewidths=1) 
    
    waveq = ['Total Wave','Wave 1','Wave 2']
    plt.text(384,80,r'\textbf{%s}' % waveq[i],color='k',
             fontsize=8,rotation=0,ha='center',va='center')
    plt.text(384,140,r'\textbf{R=%s}' % str(corr)[:4],color='k',
         fontsize=8,rotation=0,ha='center',va='center')
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    xxlabels = ['0','60E','120E','180','120W','60W','0']
    
    if i==2:
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,361,60),xxlabels,fontsize=8)
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
        plt.xlabel(r'\textbf{Longitude ($^\circ$)}',color='k',fontsize=8)
    else:
        plt.ylim([1000,10])
        plt.xticks([])
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
    
    cmap = cmocean.cm.balance            
    cs.set_cmap(cmap) 
    
plt.savefig(directoryfigure + 'linearInterference_December_FICTHIT_All.png',dpi=300)
print('Completed: Script done!')


print('Completed: Script done!')

