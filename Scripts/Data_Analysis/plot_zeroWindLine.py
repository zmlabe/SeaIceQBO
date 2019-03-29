"""
Calculates zero wind line per QBO phase in each experiment - zonal-mean zonal
wind at 30 hPa for DJF (meridional profile)

Notes
-----
    Author : Zachary Labe
    Date   : 21 June 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput_AllMembers as MO
import calc_Utilities as UT
import cmocean

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
print('\n' '----Plotting Zero Wind Line by QBO- %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['U']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'DJF'

def procData(varnames,period,qbophase):
    """
    Function reads in zonal wind data for listed experiments. It then 
    post-processes the data for calculating the zero-wind line

    Parameters
    ----------
    varnames : list of strings
        list of variables to download
    period : list of strings
        month(s) to calculate
    qbophase : list of strings
        list of qbo phases

    Returns
    -------
    avez : list of arrays
        arrays for post-processed data for each experiment
    lat : 1d array
        latitude

    Usage
    -----
    avez,lat = procData(varnames,period,qbophase)
    """
    print('\n>>> Using procData function!')
    
    ### Read in data
    for v in range(len(varnames)):
        ### Call function for surface temperature data from reach run
        lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                                   'profile')
        lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[v],'FIT',
                                                   'profile')
        lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FICT',
                                                    'profile')
        
        ### Create 2d array of latitude and longitude
        lon2,lat2 = np.meshgrid(lon,lat)
        
        ### Read in QBO phases 
        filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
        filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
        filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
        filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
        pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,
                                          usecols=[0],dtype='int'),
                            np.genfromtxt(filenamehitp2,unpack=True,
                                          usecols=[0],dtype='int')+101)
        neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,
                                          usecols=[0],dtype='int'),
                            np.genfromtxt(filenamehitn2,unpack=True,
                                          usecols=[0],dtype='int')+101)    
    
        filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
        filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
        filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
        filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
        pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,
                                          usecols=[0],dtype='int'),
                            np.genfromtxt(filenamefitp2,unpack=True,
                                          usecols=[0],dtype='int')+101)
        neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,
                                          usecols=[0],dtype='int'),
                            np.genfromtxt(filenamefitn2,unpack=True,
                                          usecols=[0],dtype='int')+101)
        
        filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
        filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
        filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
        filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
        pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,
                                           usecols=[0],dtype='int'),
                            np.genfromtxt(filenamefictp2,unpack=True,
                                          usecols=[0],dtype='int')+101)
        neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,
                                           usecols=[0],dtype='int'),
                            np.genfromtxt(filenamefictn2,unpack=True,
                                          usecols=[0],dtype='int')+101)    
        ### Concatonate runs
        runs = [tashit,tasfit,tasfict]
        
        ### Separate per periods (ON,DJ,FM)
        if period == 'DJF':
            tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],
                              tashit.shape[3],tashit.shape[4]))
            for i in range(len(runs)):
                tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                                      lon,'profile',17)   
        elif period == 'ND':
            tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],
                  tashit.shape[3],tashit.shape[4]))
            for i in range(len(runs)):
                tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
        else:
            ValueError('Wrong period selected! (DJF)')
            
        ### Concatonate runs with selected level
        levq = np.where(lev == 30)[0] # selected at 10 hPa
        tas_mo = [tas_mo[0][:,levq,:,:],tas_mo[1][:,levq,:,:],
                  tas_mo[2][:,levq,:,:]]
            
        ### Composite by QBO phase    
        tas_mofitpos = tas_mo[1][pos_fit,:,:].squeeze()
        tas_mohitpos = tas_mo[0][pos_hit,:,:].squeeze()
        tas_mofictpos = tas_mo[2][pos_fict,:,:].squeeze()
        
        tas_mofitneg = tas_mo[1][neg_fit,:,:].squeeze()
        tas_mohitneg = tas_mo[0][neg_hit,:,:].squeeze()
        tas_mofictneg = tas_mo[2][neg_fict,:,:].squeeze()
        
        ### Compute meridional average
        tas_mofitposz = np.nanmean(tas_mofitpos,axis=2)
        tas_mohitposz = np.nanmean(tas_mohitpos,axis=2)
        tas_mofictposz = np.nanmean(tas_mofictpos,axis=2)
        
        tas_mofitnegz = np.nanmean(tas_mofitneg,axis=2)
        tas_mohitnegz = np.nanmean(tas_mohitneg,axis=2)
        tas_mofictnegz = np.nanmean(tas_mofictneg,axis=2)
        
        avez = [tas_mofitposz,tas_mohitposz,tas_mofictposz,
                tas_mofitnegz,tas_mohitnegz,tas_mofictnegz]
        
        print('\n*Completed: Finished procData function!')
        return avez,lat

### Read in data function and post-process 
avez,lat = procData(varnames,period,qbophase)

### Calculate ensemble average
meanhitpos = np.nanmean(avez[1],axis=0)
meanhitneg = np.nanmean(avez[4],axis=0)

meanfitpos = np.nanmean(avez[0],axis=0)
meanfitneg = np.nanmean(avez[3],axis=0)

meanfictpos = np.nanmean(avez[2],axis=0)
meanfictneg = np.nanmean(avez[5],axis=0)

stdhitpos = np.nanstd(avez[1],axis=0)
stdhitneg = np.nanstd(avez[4],axis=0)

stdfitpos = np.nanstd(avez[0],axis=0)
stdfitneg = np.nanstd(avez[3],axis=0)

stdfictpos = np.nanstd(avez[2],axis=0)
stdfictneg = np.nanstd(avez[5],axis=0)

###############################################################################
###############################################################################
###############################################################################    
### Plot zero wind line
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
        
num_bins = np.arange(-20,60,2)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.axvline(0,linewidth=2,linestyle='--',dashes=(1,0.3),color='dimgrey')
plt.axhline(0,linewidth=2,linestyle='--',dashes=(1,0.3),color='dimgrey')

plt.plot(lat,meanfictpos+stdfictpos,linewidth=0.5,
         color=cmocean.cm.thermal(0.3))
plt.plot(lat,meanfictpos-stdfictpos,linewidth=0.5,
         color=cmocean.cm.thermal(0.3))
plt.fill_between(lat,meanfictpos+stdfictpos,meanfictpos-stdfictpos,
                 color=cmocean.cm.thermal(0.3),alpha=0.3)
plt.plot(lat,meanfictpos,linewidth=3,label=r'\textbf{QBO-W}',
         color=cmocean.cm.thermal(0.3))

plt.plot(lat,meanfictneg+stdfictneg,linewidth=0.5,
         color=cmocean.cm.thermal(0.7))
plt.plot(lat,meanfictneg-stdfictneg,linewidth=0.5,
         color=cmocean.cm.thermal(0.7))
plt.fill_between(lat,meanfictneg+stdfictneg,meanfictneg-stdfictneg,
                 color=cmocean.cm.thermal(0.7),alpha=0.3)
plt.plot(lat,meanfictneg,linewidth=3,label=r'\textbf{QBO-E}',
         color=cmocean.cm.thermal(0.7))

#plt.plot(lat,meanhitpos+stdhitpos,linewidth=0.5,
#         color=cmocean.cm.thermal(0.5))
#plt.plot(lat,meanhitpos-stdhitpos,linewidth=0.5,
#         color=cmocean.cm.thermal(0.5))
#plt.fill_between(lat,meanhitpos+stdhitpos,meanhitpos-stdhitpos,
#                 color=cmocean.cm.thermal(0.5),alpha=0.3)
#plt.plot(lat,meanhitpos,linewidth=3,label=r'\textbf{HIT-W}',
#         color=cmocean.cm.thermal(0.5))
#
#plt.plot(lat,meanhitneg+stdhitneg,linewidth=0.5,
#         color=cmocean.cm.thermal(0.7))
#plt.plot(lat,meanhitneg-stdhitneg,linewidth=0.5,
#         color=cmocean.cm.thermal(0.7))
#plt.fill_between(lat,meanhitneg+stdhitneg,meanhitneg-stdhitneg,
#                 color=cmocean.cm.thermal(0.7),alpha=0.3)
#plt.plot(lat,meanhitneg,linewidth=3,label=r'\textbf{HIT-E}',
#         color=cmocean.cm.thermal(0.7))

l = plt.legend(shadow=False,fontsize=8,loc='upper right',
           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(1, 1.0),
           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
for text in l.get_texts():
    text.set_color('k')  

plt.yticks(np.arange(-30,31,5),list(map(str,np.arange(-30,31,5))),
           fontsize=10)
plt.xticks(np.arange(-90,91,10),list(map(str,np.arange(-90,91,10))),
           fontsize=10) 
plt.ylabel(r'\textbf{U30 [m/s]}',color='dimgrey',fontsize=12)  
plt.xlabel(r'\textbf{Latitude [$^\circ$N]}',color='dimgrey',fontsize=12)
plt.xlim([40,-20])
plt.ylim([-30,20])

plt.tight_layout()

plt.savefig(directoryfigure + 'zeroWindLine_DJF.png',dpi=300)