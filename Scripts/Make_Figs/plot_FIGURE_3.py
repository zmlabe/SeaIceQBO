"""
Plot manuscript figures for Siberian High/Temperature distributions

Notes
-----
    Author : Zachary Labe
    Date   : 29 October 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import scipy.stats as sts
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directorydata3 = '/home/zlabe/Documents/Research/SeaIceQBO/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calculating Siberian High - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
def readVariablesSLP(varnames,period,location):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'D':
        tas_mo= np.empty((2,tashit.shape[0],31,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,90:121,:,:]
    else:
        ValueError('Wrong period selected!')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:]
    
    ### Compute comparisons for months - select region
    if varnames == 'SLP':
        lonq = np.where((lon >=80) & (lon <=120))[0]
        fictpos = tas_mofictpos[:,:,:,lonq]
        fictneg = tas_mofictneg[:,:,:,lonq]
        latq = np.where((lat >=40) & (lat <=65))[0]
        fictpos = fictpos[:,:,latq]
        fictneg = fictneg[:,:,latq]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        fictpos = UT.calc_weightedAve(fictpos,lat2s)
        fictneg = UT.calc_weightedAve(fictneg,lat2s)
        
        hitpos = tas_mohitpos[:,:,:,lonq]
        hitneg = tas_mohitneg[:,:,:,lonq]
        hitpos = hitpos[:,:,latq]
        hitneg = hitneg[:,:,latq]
        hitpos = UT.calc_weightedAve(hitpos,lat2s)
        hitneg = UT.calc_weightedAve(hitneg,lat2s)
        
    diffruns = [fictpos.squeeze(),fictneg.squeeze(),hitpos.squeeze(),hitneg.squeeze()]
    
    return diffruns,lat,lon,lev

def readVariablesTemps(varnames,period,location):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'D':
        tas_mo= np.empty((2,tashit.shape[0],31,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,90:121,:,:]
    else:
        ValueError('Wrong period selected!')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:] - 273.15
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:] - 273.15
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:] - 273.15
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:] - 273.15
    
    ### Compute comparisons for months - select region
    if varnames == 'T1000':
        lonq = np.where((lon >=70) & (lon <=140))[0]
        ficthitpos = tas_mofictpos[:,:,:,lonq] - tas_mohitpos[:,:,:,lonq]
        ficthitneg = tas_mofictneg[:,:,:,lonq] - tas_mohitneg[:,:,:,lonq]
        latq = np.where((lat >=35) & (lat <=60))[0]
        ficthitpos = ficthitpos[:,:,latq,:]
        ficthitneg = ficthitneg[:,:,latq,:]
        lat2sq = lat2[latq,:]
        lat2s = lat2sq[:,lonq]
        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
    diffruns = [ficthitpos.squeeze(),ficthitneg.squeeze()]
    
    return diffruns,lat,lon,lev

### Read Data
temps,lat,lon,lev = readVariablesTemps('T1000',period,'Siberia')
slp,lat,lon,lev = readVariablesSLP('SLP',period,'Siberia')

### Read in CDI
data1 = Dataset(directorydata3 + 'ColdDayIndex_December.nc')
fictpos = data1.variables['cdi_fictpos'][:]
fictneg = data1.variables['cdi_fictneg'][:]
lat = data1.variables['lat'][:]
lon = data1.variables['lon'][:]
data1.close()

data2 = Dataset(directorydata3 + 'ColdDayIndex_pvals_December.nc')
pvalpos = data2.variables['p_W'][:]
pvalneg = data2.variables['p_E'][:]
data2.close()
    
diffruns = [np.nanmean(fictpos,axis=0),np.nanmean(fictneg,axis=0)]
pvarn = [pvalpos,pvalneg]

###############################################################################
###############################################################################
###############################################################################    
### Plot histogram distributions
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})

### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])

fig = plt.figure()
grid = plt.GridSpec(2,2)

for i in range(len(diffruns)):
    limit = np.arange(-20,21,1)
    barlim = np.arange(-20,21,10)
    var = diffruns[i]
    p = pvarn[i]
    
    ax1 = plt.subplot(grid[0,i])
    m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                area_thresh=10000.)
    
    var, lons_cyclic = addcyclic(var, lon)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    pvar,lons_cyclic = addcyclic(p, lon)
    pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
    x, y = m(lon2d, lat2d)
              
    circle = m.drawmapboundary(fill_color='white',color='dimgray',
                               linewidth=0.7)
    circle.set_clip_on(False)
    m.fillcontinents(color='dimgrey',zorder=3)
    
    cs = m.contourf(x,y,var,limit,extend='both',zorder=4)
    cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                     linewidths=0.4,zorder=5)
              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.fillcontinents(color='dimgrey',zorder=3)

    m.drawcoastlines(color='dimgray',linewidth=0.8)
            
    cs.set_cmap(cmocean.cm.balance)   
        
    m.drawlsmask(land_color=(0,0,0,0),ocean_color='dimgrey',lakes=True,
                 resolution='c',zorder=5)
                
    ### Add experiment text to subplot
    if i < 2:
        qbophaseq = [r'QBO-W',r'QBO-E']
        ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                     textcoords='axes fraction',color='k',
                     fontsize=17,rotation=0,ha='center',va='center')
            
###########################################################################
cbar_ax = fig.add_axes([0.5,0.58,0.02,0.3])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{$\Delta$Cold Days Intensity}',fontsize=8,color='k',rotation=0,
               va='center',ha='center',y=-0.1,labelpad=-16)

cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))
cbar.ax.tick_params(axis='y', size=.01,labelsize=6)
cbar.outline.set_edgecolor('dimgrey')

###############################################################################
###############################################################################
###############################################################################
ax = plt.subplot(grid[1,:1])

num_bins = np.arange(-20,10.1,0.5)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

### Fit a distribution
mneg,seng = sts.norm.fit(temps[1].ravel())
mpos,spos = sts.norm.fit(temps[0].ravel())

pdfneg = sts.norm.pdf(num_bins,mneg,seng)
pdfpos = sts.norm.pdf(num_bins,mpos,spos)

plt.plot(num_bins,pdfneg,color='darkblue',linewidth=2.0,label=r'\textbf{QBO-E}')
plt.plot(num_bins,pdfpos,color='darkred',linewidth=2.0,label=r'\textbf{QBO-W}')

plt.yticks(np.arange(0,0.25,0.02),list(map(str,np.arange(0,0.25,0.02))),
           fontsize=6)
plt.xticks(np.arange(-50,1061,5),list(map(str,np.arange(-50,1061,5))),
           fontsize=6) 
plt.xlim([-15,10])
plt.ylim([0,0.12])

l = plt.legend(shadow=False,fontsize=7,loc='upper left',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(0.72,1),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for text in l.get_texts():
    text.set_color('k')

plt.ylabel(r'\textbf{Density}',color='k',fontsize=8)  
plt.xlabel(r'\textbf{$\Delta$T1000 [$^{\circ}$C]}',color='k',fontsize=8)

ax = plt.subplot(grid[1,1:])
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(0)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
ax.set_yticks([])

plt.axvline(np.median(slp[1]),color='dimgrey',dashes=(1,0.3),linewidth=1.5,
            zorder=1)

pos = [1,2,3.5,4.5]
vp = plt.violinplot(slp,pos,showmeans=False,showmedians=True,vert=False,
                    widths=0.6)

plt.xlabel(r'\textbf{Sea Level Pressure [hPa]}',color='k',fontsize=8)
plt.xticks(np.arange(1010,1061,15),list(map(str,np.arange(1010,1061,15))),
           fontsize=6) 
plt.xlim([1010,1055])

for i in vp['bodies']:
    i.set_edgecolor('w') 
    i.set_alpha(0.7)
vp['cbars'].set_color('k')
vp['cmedians'].set_color('k')
vp['cmaxes'].set_color('w')
vp['cmins'].set_color('w')
vp['cmaxes'].set_linewidth(0.5)        
vp['cmins'].set_linewidth(0.5) 
vp['cmedians'].set_linewidth(1.5)
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')              
vp['bodies'][0].set_facecolor('darkred')
vp['bodies'][1].set_facecolor('darkblue')  
vp['bodies'][2].set_facecolor('darkred')
vp['bodies'][3].set_facecolor('darkblue')  
vp['bodies'][0].set_alpha(1)  
vp['bodies'][1].set_alpha(1)  
vp['bodies'][2].set_alpha(1)  
vp['bodies'][3].set_alpha(1)  
vp['bodies'][0].set_hatch('\\\\')
vp['bodies'][1].set_hatch('\\\\')

ax.annotate(r'\textbf{[c]}',xy=(1030,3),xytext=(0.98,1.0),
             textcoords='axes fraction',color='dimgray',
             fontsize=12,rotation=0,ha='center',va='center')
ax.annotate(r'\textbf{[b]}',xy=(1030,3),xytext=(-1.14,1.0),
             textcoords='axes fraction',color='dimgray',
             fontsize=12,rotation=0,ha='center',va='center')
ax.annotate(r'\textbf{[a]}',xy=(1030,3),xytext=(-1.14,2.3),
             textcoords='axes fraction',color='dimgray',
             fontsize=12,rotation=0,ha='center',va='center')

ax.annotate(r'\textbf{Historical}',xy=(1030,3),xytext=(1.14,0.77),
             textcoords='axes fraction',color='k',
             fontsize=10,rotation=270,ha='center',va='center')

ax.annotate(r'\textbf{Future}',xy=(1030,3),xytext=(1.14,0.215),
             textcoords='axes fraction',color='k',
             fontsize=10,rotation=270,ha='center',va='center')

ax.annotate(r'$\bf{\bullet\bullet}$',xy=(1030,3),xytext=(1.03,0.33),
             textcoords='axes fraction',color='k',
             fontsize=10,rotation=270,ha='center',va='center')
ax.annotate(r'$\bf{\bullet\bullet}$',xy=(1030,3),xytext=(1.03,0.108),
             textcoords='axes fraction',color='k',
             fontsize=10,rotation=270,ha='center',va='center')

plt.savefig(directoryfigure + 'SiberianHigh_Distributions.png',dpi=900)

### Calculate statistical significance of distributions
t,pvaluedist = sts.ks_2samp(temps[1].ravel(),temps[0].ravel())
print('\n \n' + 'Distribution p_value=' + str(pvaluedist))

### Calculate statistical text for violin plots
t,pvalpos = sts.ttest_ind(slp[0].ravel(),slp[2].ravel(),equal_var=True) 
t,pvalneg = sts.ttest_ind(slp[1].ravel(),slp[3].ravel(),equal_var=True) 

print('\n' + 'QBO-W FICT-HIT p_value=' + str(pvalpos))
print('QBO-E FICT-HIT p_value=' + str(pvalneg))

t,pvalfict = sts.ttest_ind(slp[0].ravel(),slp[1].ravel(),equal_var=True) 
t,pvalhit = sts.ttest_ind(slp[2].ravel(),slp[3].ravel(),equal_var=True) 

print('\n' + 'QBO-E & QBO-W FICT p_value=' + str(pvalfict))
print('QBO-E & QBO-W HIT p_value=' + str(pvalhit))

t,pvalanoms = sts.ttest_ind(slp[1].ravel()-slp[3].ravel(),
                            slp[0].ravel()-slp[2].ravel(),equal_var=True) 
print('\n' + 'Anomaly differences p_value=' + str(pvalanoms))