"""
Plot plumb flux comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, FICT. Profiles are 
organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 28 September 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput_AllMembers as MO
import cmocean
import calc_Utilities as UT
import scipy.stats as sts

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
print('\n' '----Plotting QBO profile comparisons - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'ND'
varnames = ['WAFZ']

### Call function for WAFz data for each run
#lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[0],'HIT',
#                                           'profile')
#lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[0],'FIT',
#                                           'profile')
#lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[0],'FICT',
#                                            'profile')
#
#### Create 2d array of latitude and longitude
#lon2,lat2 = np.meshgrid(lon,lat)

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

filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+101)
non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+101)
neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+101)

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
runs = [tashit,tasfit,tasfict]

### Separate per periods (ON,DJ,FM)
if period == 'ON': 
    tas_mo = np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = np.nanmean(runs[i][:,9:11,:,:,:],axis=1) 
elif period == 'DJ':     
    tas_mo = np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
                                            lon,'profile',1) 
elif period == 'FM':
    tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = np.nanmean(runs[i][:,1:3,:,:,:],axis=1)
elif period == 'DJF':
    tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
                                              lon,'profile',1)   
elif period == 'M':
    tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = runs[i][:,2,:,:,:]
elif period == 'D':
    tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = runs[i][:,-1,:,:,:]
elif period == 'N':
    tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = runs[i][:,-2,:,:,:]
elif period == 'ND':
    tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3],tashit.shape[4]))
    for i in range(len(runs)):
        tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:,:],axis=1)
else:
    ValueError('Wrong period selected! (ON,DJ,FM)')
    
### Composite by QBO phase    
tas_mohitposq = tas_mo[0][pos_hit,:,:]
tas_mofitposq = tas_mo[1][pos_fit,:,:]
tas_mofictposq = tas_mo[2][pos_fict,:,:]

tas_mohitnonq = tas_mo[0][non_hit,:,:]
tas_mofitnonq = tas_mo[1][non_fit,:,:]
tas_mofictnonq = tas_mo[2][non_fict,:,:]

tas_mohitnegq = tas_mo[0][neg_hit,:,:]
tas_mofitnegq = tas_mo[1][neg_fit,:,:]
tas_mofictnegq = tas_mo[2][neg_fict,:,:]

### Take zonal average
tas_mohitpos = np.nanmean(tas_mohitposq,axis=3)
tas_mofitpos = np.nanmean(tas_mofitposq,axis=3)
tas_mofictpos = np.nanmean(tas_mofictposq,axis=3)

tas_mohitnon = np.nanmean(tas_mohitnonq,axis=3)
tas_mofitnon = np.nanmean(tas_mofitnonq,axis=3)
tas_mofictnon = np.nanmean(tas_mofictnonq,axis=3)

tas_mohitneg = np.nanmean(tas_mohitnegq,axis=3)
tas_mofitneg = np.nanmean(tas_mofitnegq,axis=3)
tas_mofictneg = np.nanmean(tas_mofictnegq,axis=3)

### Compute climatology    
climofitpos = np.nanmean(tas_mofitpos,axis=0)
climohitpos = np.nanmean(tas_mohitpos,axis=0)
climofictpos = np.nanmean(tas_mofictpos,axis=0)
climofitnon = np.nanmean(tas_mofitnon,axis=0)
climohitnon = np.nanmean(tas_mohitnon,axis=0)
climofictnon = np.nanmean(tas_mofictnon,axis=0)
climofitneg = np.nanmean(tas_mofitneg,axis=0)
climohitneg = np.nanmean(tas_mohitneg,axis=0)
climofictneg = np.nanmean(tas_mofictneg,axis=0)
climo = [climohitpos,climohitnon,climohitneg,
         climohitpos,climohitnon,climohitneg]

### Compute comparisons for months - taken ensemble average
#fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)
#fithitnon = np.nanmean(tas_mofitnon - tas_mohitnon,axis=0)
#fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)
#
#ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
#ficthitnon = np.nanmean(tas_mofictnon - tas_mohitnon,axis=0)
#ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)/np.nanstd(tas_mofitpos,axis=0)
fithitnon = np.nanmean(tas_mofitnon - tas_mohitnon,axis=0)/np.nanstd(tas_mofitnon,axis=0)
fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)/np.nanstd(tas_mofitneg,axis=0)

ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)/np.nanstd(tas_mofictpos,axis=0)
ficthitnon = np.nanmean(tas_mofictnon - tas_mohitnon,axis=0)/np.nanstd(tas_mofictnon,axis=0)
ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)/np.nanstd(tas_mofictneg,axis=0)
#fithitpos = np.nanmean(sts.zscore(tas_mofitpos,axis=0) - sts.zscore(tas_mohitpos,axis=0),axis=0)
#fithitnon = np.nanmean(sts.zscore(tas_mofitnon,axis=0) - sts.zscore(tas_mohitnon,axis=0),axis=0)
#fithitneg = np.nanmean(sts.zscore(tas_mofitneg,axis=0) - sts.zscore(tas_mohitneg,axis=0),axis=0)
#
#ficthitpos = np.nanmean(sts.zscore(tas_mofictpos,axis=0) - sts.zscore(tas_mohitpos,axis=0),axis=0)
#ficthitnon = np.nanmean(sts.zscore(tas_mofictnon,axis=0) - sts.zscore(tas_mohitnon,axis=0),axis=0)
#ficthitneg = np.nanmean(sts.zscore(tas_mofictneg,axis=0) - sts.zscore(tas_mohitneg,axis=0),axis=0)
diffruns_mo = [fithitneg,fithitpos,fithitneg-fithitpos,
               ficthitneg,ficthitpos,ficthitneg-ficthitpos]

#diffruns_mo = [sts.zscore(fithitneg),sts.zscore(fithitpos),sts.zscore(fithitneg)-sts.zscore(fithitpos),
#               sts.zscore(ficthitneg),sts.zscore(ficthitpos),sts.zscore(ficthitneg)-sts.zscore(ficthitpos)]

### Calculate significance for FM
stat_FITHITpos,pvalue_FITHITpos = UT.calc_indttest(tas_mohitpos,
                                                   tas_mofitpos)
stat_FITHITnon,pvalue_FITHITnon = UT.calc_indttest(tas_mohitnon,
                                                   tas_mofitnon)
stat_FITHITneg,pvalue_FITHITneg = UT.calc_indttest(tas_mohitneg,
                                                   tas_mofitneg)

stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mohitpos,
                                                     tas_mofictpos)
stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(tas_mohitnon,
                                                     tas_mofictnon)
stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mohitneg,
                                                     tas_mofictneg)

pruns_mo = [pvalue_FITHITneg,pvalue_FITHITpos,pvalue_FITHITneg,
            pvalue_FICTHITneg,pvalue_FICTHITpos,pvalue_FICTHITneg]

###########################################################################
###########################################################################
###########################################################################
#### Plot EP Flux
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
#limit = np.arange(-0.01,0.011,0.001)
#barlim = np.arange(-0.01,0.01,0.01)
limit = np.arange(-0.5,0.6,0.05)
barlim = np.arange(-0.5,1,0.5)
    
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(len(diffruns_mo)):
    ax1 = plt.subplot(2,3,i+1)
    
    var = diffruns_mo[i]
    pruns = pruns_mo[i]

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
            
    cs = plt.contourf(lat,lev,var,limit,extend='both')

    if any([i==0,i==1,i==3,i==4]):
        plt.contourf(lat,lev,pruns,colors='None',hatches=['////'],
                     linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=6)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    plt.minorticks_off()
    
    plt.xlim([30,90])
    plt.ylim([1000,10])
        
    cmap = cmocean.cm.balance          
    cs.set_cmap(cmap) 

    ### Add experiment text to subplot
    if i < 3:
        qbophaseq = [r'QBO-E',r'QBO-W',r'Difference']
        ax1.text(0.5,1.1,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes)
    if i == 0 or i == 3:
        ax1.text(-0.35,0.5,r'%s' % experiments[i],color='k',
                     fontsize=20,rotation=90,ha='center',va='center',
                     transform=ax1.transAxes)
    if i == 4:
        plt.xlabel(r'\textbf{Latitude ($^\circ$N)}',color='k',fontsize=8,
                             labelpad=1)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{Normalized WAFz}',fontsize=11,color='dimgray',labelpad=1)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')     

plt.subplots_adjust(wspace=0.24)
plt.subplots_adjust(bottom=0.21)

plt.savefig(directoryfigure + 'WAFz_QBO_%s.png' % period,dpi=300)
print('Completed: Script done!')

