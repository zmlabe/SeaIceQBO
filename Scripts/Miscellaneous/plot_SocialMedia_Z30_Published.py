"""
Script plots figure for manuscript of stratosphere variables of U10 and Z30

Notes
-----
    Author : Zachary Labe
    Date   : 17 October 2018
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
import itertools

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
print('\n' '----Plotting QBO comparisons - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
experimentsn = [r'\textbf{$\Delta$NET}']
qbophaseq = [r'QBO-W',r'QBO-E',r'E--W']
qbophase = ['pos','non','neg']
letters = ["a","b","c","d","e","f","g","h","i","j","k","l","m"]
period = 'D'

def calc_indttest_90(varx,vary):
    """
    Function calculates statistical difference for 2 independent
    sample t-test

    Parameters
    ----------
    varx : 3d array
    vary : 3d array
    
    Returns
    -------
    stat = calculated t-statistic
    pvalue = two-tailed p-value

    Usage
    -----
    stat,pvalue = calc_ttest(varx,vary)
    """
    print('\n>>> Using calc_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_ind(varx,vary,nan_policy='omit')
    
    ### Significant at 90% confidence level
    pvalue[np.where(pvalue >= 0.1)] = np.nan
    pvalue[np.where(pvalue < 0.1)] = 1.
    pvalue[np.isnan(pvalue)] = 0.
    
    print('*Completed: Finished calc_ttest function!')
    return stat,pvalue

def calcVarResp(varnames,period,qbophase):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
                                                'surface')
    
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
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+101)
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

    ### Compute climatology    
    climohitpos = np.nanmean(tas_mohitpos,axis=0)
    climohitneg = np.nanmean(tas_mohitneg,axis=0)
    climo = [climohitpos,climohitneg,climohitneg]
    
    ### Compute comparisons for months - taken ensemble average
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    diffpos = tas_mofictpos - tas_mohitpos
    diffneg = tas_mofictneg - tas_mohitneg
    diffall = np.nanmean(diffneg,axis=0) - np.nanmean(diffpos,axis=0)
    
    diffruns = [ficthitpos,ficthitneg,diffall]
    
    ### Calculate significance for ND
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[1][pos_fict,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(tas_mo[1][non_fict,:,:],
                                                   tas_mo[0][non_hit,:,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[1][neg_fict,:,:],
                                               tas_mo[0][neg_hit,:,:])
    stat_diff,pvalue_diff = calc_indttest_90(diffneg,diffpos)

    pruns = [pvalue_FICTHITpos,pvalue_FICTHITneg,pvalue_diff]
    
    return diffruns,pruns,climo,lat,lon

### Call variables
#variables,pvalues,climos,lat,lon = calcVarResp('Z30',period,qbophase)
#varnamesn = [r'Z30',r'Z30']
    
###########################################################################
###########################################################################
###########################################################################
### Plot variable data for DJ
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('savefig',facecolor='black')
plt.rc('axes',edgecolor='darkgrey')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',labelcolor='white')
plt.rc('axes',facecolor='black')

fig = plt.figure()

for v in range(len(variables)-1):
    ax = plt.subplot(1,2,v+1)
    
    ### Retrieve variables and pvalues
    var = variables[v]
    pvar = pvalues[v]
    
    ### Set limits for contours and colorbars
    if varnamesn[v] == 'SLP':
        limit = np.arange(-4,4.1,0.05)
        barlim = np.arange(-4,5,4)
    elif varnamesn[v] == 'Z500':
        limit = np.arange(-60,60.1,5)
        barlim = np.arange(-60,61,30) 
    elif varnamesn[v] == 'Z30':
        limit = np.arange(-100,100.1,5)
        barlim = np.arange(-100,101,50) 
    
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
    
    pvar,lons_cyclic = addcyclic(pvar, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    climoq,lons_cyclic = addcyclic(climos[v], lon)
    climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
              
    m.drawmapboundary(fill_color='white',color='dimgrey',linewidth=0.7)
    
    cs = m.contourf(x,y,var,limit,extend='both')
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['...'])
    if varnamesn[v] == 'Z30': # the interval is 250 m 
        cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                        colors='k',linewidths=2,zorder=10)
    if varnamesn[v] == 'RNET':
        m.drawcoastlines(color='darkgray',linewidth=0.3)
        m.fillcontinents(color='dimgrey')
    else:
        m.drawcoastlines(color='dimgrey',linewidth=1.1)
    
    if varnamesn[v] == 'SLP':
        cmap = cmocean.cm.balance          
        cs.set_cmap(cmap)   
    elif varnamesn[v] == 'Z500':
        cmap = cmocean.cm.balance           
        cs.set_cmap(cmap)  
    elif varnamesn[v] == 'Z30':
        cmap = cmocean.cm.balance  
        cs.set_cmap(cmap)  
    
    if v == 0:        
        ax.annotate(r'\textbf{LOSS OF SEA ICE}',xy=(0,0),xytext=(0.5,1.16),
                     textcoords='axes fraction',color='darkgray',
                     fontsize=22,rotation=0,ha='center',va='center')
        ax.annotate(r'\textbf{\underline{UNDER QBO-W}}',xy=(0,0),xytext=(0.5,1.045),
             textcoords='axes fraction',color='w',
             fontsize=22,rotation=0,ha='center',va='center')
    elif v == 1:        
        ax.annotate(r'\textbf{LOSS OF SEA ICE}',xy=(0,0),xytext=(0.5,1.16),
                     textcoords='axes fraction',color='darkgray',
                     fontsize=22,rotation=0,ha='center',va='center')
        ax.annotate(r'\textbf{\underline{UNDER QBO-E}}',xy=(0,0),xytext=(0.5,1.045),
                 textcoords='axes fraction',color='w',
                 fontsize=22,rotation=0,ha='center',va='center')
        
    ax.set_aspect('equal')
            
    ###########################################################################
    cbar_ax = fig.add_axes([0.304,0.11,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    if varnamesn[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'Z30':
        cbar.set_label(r'\textbf{30-hPa Geopotential Height [m]}',fontsize=10,color='darkgray')  
    elif varnamesn[v] == 'SLP':
        cbar.set_label(r'\textbf{Sea Level Pressure [hPa]}',fontsize=11,color='darkgray',
                                 labelpad=1.4)  
    elif varnamesn[v] == 'U10' or varnamesn[v] == 'U300':
        cbar.set_label(r'\textbf{m/s}',fontsize=11,color='dimgray')  
    elif varnamesn[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='dimgray')
    elif varnamesn[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='dimgray') 
    elif varnamesn[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='dimgray') 
    elif varnamesn[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='dimgray')
    elif varnamesn[v] == 'RNET':
        cbar.set_label(r'\textbf{W/m$^{\bf{2}}$}',fontsize=11,color='dimgray') 
    
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))    
    cbar.ax.tick_params(axis='x', size=.001,labelcolor='darkgray',labelsize=7)
    cbar.outline.set_edgecolor('darkgray')
    plt.tight_layout()
    
    plt.annotate(r'\textbf{December}',xy=(0,0),xytext=(0.5,0.165),
         textcoords='figure fraction',color='darkgray',
         fontsize=8,rotation=0,ha='center',va='center')
 
    plt.text(-0.915,-2.50,r'\textbf{STATISTICS:} Stippling -- significant at p$<$0.05',
     fontsize=5,rotation='horizontal',ha='left',color='darkgrey')
    plt.text(-0.915,-1.90,r'\textbf{REFERENCE:} adapted from Labe et al. [2019, \textit{GRL}]',
         fontsize=5,rotation='horizontal',ha='left',color='darkgrey')
    plt.text(-0.915,-3.10,r'\textbf{GRAPHIC:} Zachary Labe (@ZLabe)',
         fontsize=5,rotation='horizontal',ha='left',color='darkgrey') 
    
    
plt.savefig(directoryfigure + 'SocialMedia_Z30_Published.png',dpi=1000)

print('Completed: Script done!')

