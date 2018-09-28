"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Composites are organized by QBO phase (positive, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 24 September 2018
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
import scipy.stats as sts

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
varnames = ['U10']
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
def readVar(varnames,runnames,qbophase,period):
    for v in range(len(varnames)):
        ### Call function for surface temperature data from reach run
        lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                                   'surface')
        lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[v],'FIT',
                                                   'surface')
        lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FIT',
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
        runs = [tashit,tasfit,tasfict]
        
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
            print(ValueError('Wrong period selected!'))
            
        ### Composite by QBO phase    
        vhitpos = tas_mo[0][pos_hit,:,:]
        vfictpos = tas_mo[2][pos_fict,:,:]
        
        vhitneg = tas_mo[0][neg_hit,:,:]
        vfictneg = tas_mo[2][neg_fict,:,:]
        
        ### Compute comparisons for months - taken ensemble average
        ficthitpos = np.nanmean(vfictpos - vhitpos,axis=0)
        ficthitneg = np.nanmean(vfictneg - vhitneg,axis=0)
    return vhitpos,vfictpos,vhitneg,vfictneg,ficthitpos,ficthitneg,lat,lon

### Read in Data
vhitpos,vfictpos,vhitneg,vfictneg,ficthitpos,ficthitneg,lat,lon = readVar(varnames,runnames,qbophase,period)   
    
def boot_matrix(z, B):
    """Bootstrap sample
    
    Returns all bootstrap samples in a matrix"""
    
    n = len(z)  # sample size
    idz = np.random.randint(0, n, size=(B, n))  # indices to pick for all boostrap samples
    return z[idz]

def bootstrap_mean(x, B=1000, alpha=0.05):
    """Bootstrap standard error and (1-alpha)*100% c.i. for the population mean
    
    Returns bootstrapped standard error and different types of confidence intervals"""
   
    # Deterministic things
    n = len(x)  # sample size
    orig = x.mean()  # sample mean
    se_mean = x.std()/np.sqrt(n) # standard error of the mean
    qt = sts.t.ppf(q=1 - alpha/2, df=n - 1) # Student quantile
    
    # Generate boostrap distribution of sample mean
    xboot = boot_matrix(x, B=B)
    sampling_distribution = xboot.mean(axis=1)
   
   # Standard error and sample quantiles
    se_mean_boot = sampling_distribution.std()
    quantile_boot = np.percentile(sampling_distribution, q=(100*alpha/2, 100*(1-alpha/2)))
 
    # RESULTS
    print("Estimated mean:", orig)
    print("Classic standard error:", se_mean)
    print("Classic student c.i.:", orig + np.array([-qt, qt])*se_mean)
    print("\nBootstrap results:")
    print("Standard error:", se_mean_boot)
    print("t-type c.i.:", orig + np.array([-qt, qt])*se_mean_boot)
    print("Percentile c.i.:", quantile_boot)
    print("Basic c.i.:", 2*orig - quantile_boot[::-1])

def bootstrap_t_pvalue(x, y, equal_var,B):
    """Bootstrap p values for two-sample t test
    
    Returns boostrap p value, test statistics and parametric p value"""
    
    # Original t test statistic
    orig = sts.ttest_ind(x, y,nan_policy='omit')
    
    # Generate boostrap distribution of t statistic
    xboot = boot_matrix(x - x.mean(), B=B) # important centering step to get sampling distribution under the null
    yboot = boot_matrix(y - y.mean(), B=B)
    
    sampling_distribution = sts.ttest_ind(xboot, yboot,axis=1,nan_policy='omit')[0]
    newstats = sampling_distribution[1]

    # Calculate proportion of bootstrap samples with at least as strong evidence against null    
    p = np.mean(sampling_distribution >= orig[0])
    
    # RESULTS
    print("p value for null hypothesis of equal population means:")
    print("Parametric:", orig[1])
    print("Bootstrap:", 2*min(p, 1-p))
    origg = orig[1]
    pv = 2*min(p, 1-p)
    print(p)

    return origg,pv

### Calculate bootstrap test
pv = np.empty((lat.shape[0],lon.shape[0]))
origg = np.empty((lat.shape[0],lon.shape[0]))
for i in range(lat.shape[0]):
    for j in range(lon.shape[0]):
        origg[i,j],pv[i,j] = bootstrap_t_pvalue(vhitpos[:,i,j],vfictpos[:,i,j],equal_var=False,B=10000)
    
testp = pv.copy()    
testp[np.where(testp >= 0.05)] = np.nan
testp[np.where(testp < 0.05)] = 1.

testpo = origg.copy()    
testpo[np.where(testpo >= 0.05)] = np.nan
testpo[np.where(testpo < 0.05)] = 1.

###########################################################################
###########################################################################
###########################################################################
### Plot variable data for QBO composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
if varnames[0] == 'U10':
    limit = np.arange(-5,5.1,0.5)
    barlim = np.arange(-5,6,1)
elif varnames[0] == 'Z30':
    limit = np.arange(-100,100.1,10)
    barlim = np.arange(-100,101,50) 

fig = plt.figure()
var = ficthitpos
pvar = testp

ax1 = plt.subplot(111)
m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
            area_thresh=10000.)

var, lons_cyclic = addcyclic(var, lon)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
x, y = m(lon2d, lat2d)

pvar,lons_cyclic = addcyclic(pvar, lon)
pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
          
m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)

cs = m.contourf(x,y,var,limit,extend='both')
cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                 linewidths=0.4)
m.drawcoastlines(color='dimgray',linewidth=0.8)

cs.set_cmap(cmocean.cm.balance)
plt.savefig(directoryfigure + 'testboot_FITHITpos.png',dpi=300)




















#    ###########################################################################
#    ###########################################################################
#    ###########################################################################
#    ### Plot variable data for QBO composites
#    plt.rc('text',usetex=True)
#    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#    
#    ### Set limits for contours and colorbars
#    if varnames[v] == 'U10':
#        limit = np.arange(-5,5.1,0.5)
#        barlim = np.arange(-5,6,1)
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
#
#        m.drawcoastlines(color='dimgray',linewidth=0.8)
#          
#
#        if varnames[v] == 'U10':
#            cmap = ncm.cmap('NCV_blu_red')            
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
#    if varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
#        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
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
#    plt.savefig(directoryfigure + '/QBOExperiments_%s_%s.png' % (period,
#                                                                        period,
#                                                                  varnames[v]),
#                                                                  dpi=300)
#
#print('Completed: Script done!')
#
