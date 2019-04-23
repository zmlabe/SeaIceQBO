"""
Plot temperature profile difference between HIT and FIT experiments. 
These are sea ice thickness perturbation experiments using WACCM4.

Notes
-----
    Author : Zachary Labe
    Date   : 14 August 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_MonthlyLatOutput as MO
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_epflux/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting EP flux for QBO %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)
qbophase = 'neg'

### Call function for vertical temperature data
lat,lon,time,lev,epy_h = MO.readExperi(directorydata,'EPY','HIT','profile')
lat,lon,time,lev,epz_h = MO.readExperi(directorydata,'EPZ','HIT','profile')
lat,lon,time,lev,div_h = MO.readExperi(directorydata,'DEPF','HIT','profile')

lat,lon,time,lev,epy_f = MO.readExperi(directorydata,'EPY','FIT','profile')
lat,lon,time,lev,epz_f = MO.readExperi(directorydata,'EPZ','FIT','profile')
lat,lon,time,lev,div_f = MO.readExperi(directorydata,'DEPF','FIT','profile')

### Separate per month
epy_moh = np.append(epy_h[:,9:,:,:],epy_h[:,:3,:,:],axis=1)
epz_moh = np.append(epz_h[:,9:,:,:],epz_h[:,:3,:,:],axis=1)
div_moh = np.append(div_h[:,9:,:,:],div_h[:,:3,:,:],axis=1)

epy_mof = np.append(epy_f[:,9:,:,:],epy_f[:,:3,:,:],axis=1)
epz_mof = np.append(epz_f[:,9:,:,:],epz_f[:,:3,:,:],axis=1)
div_mof = np.append(div_f[:,9:,:,:],div_f[:,:3,:,:],axis=1)

### Read in QBO phases 
filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % 'pos'
filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % 'non'
filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % 'neg'
pos_fit = np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int')
non_fit = np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int')
neg_fit = np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int')

### Calculate differences
if qbophase == 'pos':
    diff_epy = np.nanmean(epy_mof[pos_fit] - epy_moh[pos_fit],axis=0)
    diff_epz = np.nanmean(epz_mof[pos_fit] - epz_moh[pos_fit],axis=0)
    diff_div = np.nanmean(div_mof[pos_fit] - div_moh[pos_fit],axis=0)/30.
elif qbophase == 'non':
    diff_epy = np.nanmean(epy_mof[non_fit] - epy_moh[non_fit],axis=0)
    diff_epz = np.nanmean(epz_mof[non_fit] - epz_moh[non_fit],axis=0)
    diff_div = np.nanmean(div_mof[non_fit] - div_moh[non_fit],axis=0)/30.
elif qbophase == 'neg':
    diff_epy = np.nanmean(epy_mof[neg_fit] - epy_moh[neg_fit],axis=0)
    diff_epz = np.nanmean(epz_mof[neg_fit] - epz_moh[neg_fit],axis=0)
    diff_div = np.nanmean(div_mof[neg_fit] - div_moh[neg_fit],axis=0)/30.

##### Calculate significance
#stat_on,pvalue_on = UT.calc_indttest(np.nanmean(th_on,axis=3),
#                                     np.nanmean(tf_on,axis=3))
#stat_dj,pvalue_dj = UT.calc_indttest(np.nanmean(th_dj,axis=3),
#                                     np.nanmean(tf_dj,axis=3))
#stat_fm,pvalue_fm = UT.calc_indttest(np.nanmean(th_fm,axis=3),
#                                     np.nanmean(tf_fm,axis=3))

###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-2,2.1,0.1)
barlim = np.arange(-2,3,1)
    
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(6):
    ax1 = plt.subplot(2,3,i+1)
    
    clmq = i

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
            
    cs = plt.contourf(lat,lev,diff_div[i],limit,extend='both')
    cs1 = plt.quiver(lat[::4],lev,diff_epy[i][:,::4],diff_epz[i][:,::4],
                     pivot='mid',color='k',units='width',
                     scale=0.8e7,width=0.007)
    if i == 5:
        plt.quiverkey(cs1,0.34,-0.3,0.8e6,r'\textbf{0.8$\times$10$^{6}$}',
                      coordinates='axes',labelpos='E')

#        plt.contourf(latq,levq,pruns[i],colors='None',hatches=['////'],
#                     linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=7)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=7)
    plt.minorticks_off()
    
    plt.xlim([0,90])
    plt.ylim([1000,10])
        
    cmap = ncm.cmap('temp_diff_18lev')            
    cs.set_cmap(cmap) 

    labelmonths = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
    ax1.annotate(r'\textbf{%s}' % labelmonths[i],
                xy=(0, 0),xytext=(0.5,1.08),xycoords='axes fraction',
                fontsize=17,color='dimgrey',rotation=0,
                ha='center',va='center')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{m/s/day}',fontsize=11,color='dimgray',labelpad=1)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
           
plt.annotate(r'\textbf{FIT--HIT}',
        xy=(0, 0),xytext=(0.045,0.535),xycoords='figure fraction',
        fontsize=17,color='k',rotation=90,
        ha='center',va='center')        

plt.subplots_adjust(hspace=0.33)
plt.subplots_adjust(bottom=0.18)
plt.subplots_adjust(wspace=0.3)

plt.savefig(directoryfigure + 'ep_flux_FIT-HIT_QBO_%s.png' % qbophase,dpi=300)
print('Completed: Script done!')

