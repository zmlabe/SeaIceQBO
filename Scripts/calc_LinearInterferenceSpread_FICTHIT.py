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
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting linear interference - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

historicalforcedpos = []
futureforcedpos = []
historicalforcednon = []
futureforcednon = []
historicalforcedneg = []
futureforcedneg = []
lonssq = []
varnames = ['GEOPxwave_all','GEOPxwave1','GEOPxwave2']
qbophase = ['pos','non','neg']
period = 'N'
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon1,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                               'profile')
    lat,lon1,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FICT',
                                               'profile')
    
    ### Modify lons 
    lon1[-1] = 360

    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
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
    
    ### Separate per periods (December)
    if period == 'D':
        tas_mo= np.empty((2,tashit.shape[0],tashit.shape[2],tashit.shape[3],
                          tashit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:,:]
    elif period == 'N':
        tas_mo= np.empty((2,tashit.shape[0],tashit.shape[2],tashit.shape[3],
                          tashit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-2,:,:,:]
    elif period == 'ND':
        tas_mo= np.empty((2,tashit.shape[0],tashit.shape[2],tashit.shape[3],
                          tashit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:,:],axis=1)
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofictpos = tas_mo[1][pos_fict,:,:]
    
    tas_mohitnon = tas_mo[0][non_hit,:,:]
    tas_mofictnon = tas_mo[1][non_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:]
    
    ### Take at 60N 
    latq = np.where((lat>59) & (lat<61))[0]
    
    ### Compute comparisons for fictHIT - taken ensemble average at 60N (index 79)
    ### +QBO
    diff_FICTHITpos = tas_mofictpos - np.nanmean(tas_mohitpos,axis=0)
    diffrunspos = diff_FICTHITpos[:,:,latq,:].squeeze() 
    historicalforcedqpos = np.nanmean(tas_mohitpos[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcedpos.append(historicalforcedqpos)
    futureforcedpos.append(diffrunspos)
    
    ### Neutral QBO
    diff_FICTHITnon = tas_mofictnon - np.nanmean(tas_mohitnon,axis=0)
    diffrunsnon = diff_FICTHITnon[:,:,latq,:].squeeze() 
    historicalforcedqnon = np.nanmean(tas_mohitnon[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcednon.append(historicalforcedqnon)
    futureforcednon.append(diffrunsnon)
    
    ### -QBO
    diff_FICTHITneg = tas_mofictneg - np.nanmean(tas_mohitneg,axis=0)
    diffrunsneg = diff_FICTHITneg[:,:,latq,:].squeeze() 
    historicalforcedqneg = np.nanmean(tas_mohitneg[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcedneg.append(historicalforcedqneg)
    futureforcedneg.append(diffrunsneg)
    
    ### Latitudes append
    lonssq.append(lon1)
    
### Calculate correlations
def corrCalc(his,fut,wavenumber,lev,lons):
    """
    Calculate correlations for each ensemble member to understand the spread
    in linear interference for each wave during the set month
    """
    if wavenumber == 'all':
        hisw = his[0]
        futw = fut[0]
        lon = lons[0]
    elif wavenumber == '1':
        hisw = his[1]
        futw = fut[1]
        lon = lons[1]
    elif wavenumber == '2':
        hisw = his[2]
        futw = fut[2]
        lon = lons[2]
    else:
        print(ValueError('WRONG WAVE NUMBER!'))
    
    corrs = np.empty((futw.shape[0]))
    for i in range(futw.shape[0]):
        corrs[i] = UT.calc_spatialCorrHeight(hisw,futw[i],lev,lon[i],'yes')
    
    return corrs

#corrsa = corrCalc(historicalforcedneg,futureforcedneg,'all',lev,lonssq)
#corrs1 = corrCalc(historicalforcedneg,futureforcedneg,'1',lev,lonssq)
#corrs2 = corrCalc(historicalforcedneg,futureforcedneg,'2',lev,lonssq)

corrsa = corrCalc(historicalforcedpos,futureforcedpos,'all',lev,lonssq)
corrs1 = corrCalc(historicalforcedpos,futureforcedpos,'1',lev,lonssq)
corrs2 = corrCalc(historicalforcedpos,futureforcedpos,'2',lev,lonssq)
    
###############################################################################
###############################################################################
###############################################################################    
### Plot box plot distributions
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
        
dataq = [corrsa,corrs1,corrs2]
wavenames = [r'\textbf{All Waves}',r'\textbf{Wave 1}', r'\textbf{Wave 2}']
        
fig = plt.figure()
ax = plt.subplot(111) 

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('w')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
ax.tick_params(axis='x',which='both',bottom=False)

plt.axhline(0,color='dimgrey',linestyle='--',dashes=(0.9,1),linewidth=2)
bx = plt.boxplot(dataq,0,'',patch_artist=True,whis=[5,95])

for i in bx['caps']:
    i.set(color='k',linewidth=0)
for whisker in bx['whiskers']:
    whisker.set(color='dimgrey',linestyle='-',linewidth=2)
for box in bx['boxes']: 
    box.set(color='deepskyblue')
for box in bx['means']:
    box.set(linewidth=0)
for box in bx['medians']:
    box.set(color='r',linewidth=2,linestyle='-')
    
for i in range(len(dataq)):
    y = dataq[i]
    x = np.random.normal(1+i,0.04,size=len(y))
    plt.plot(x,y,'r.',alpha=0.3,zorder=5)

plt.ylabel(r'\textbf{Correlation (R)}',color='k',fontsize=12)

plt.yticks(np.arange(-1,2,1),list(map(str,np.arange(-1,2,1))))
plt.xticks(np.arange(1,4,1),wavenames) 
plt.ylim([-1,1])

plt.savefig(directoryfigure + 'LinearInterferenceSpread_FICTHITpos_%s.png' % period,dpi=300)