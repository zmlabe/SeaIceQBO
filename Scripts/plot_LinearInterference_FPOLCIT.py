"""
Calculate linear interference for FPOL-CIT

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
import read_MonthlyOutput_AllRegional as MOR
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_N_2/Regional2/'
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
period = 'D'
for v in range(len(varnames)):
    ### Call function for wave data from reach run
    lat,lon1,time,lev,tascit = MO.readExperiAll('%s' % varnames[v],
                                                'CIT','profile')
    lat,lon1,time,lev,tasfpol = MOR.readExperiAllRegional('%s' % varnames[v],
                                                'FPOL','profile')
    
    ### Modify lons 
    lon1[-1] = 360

    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon1,lat)
    
    ### Read in QBO phases 
    filenamecitp = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    filenamecitp2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    pos_cit = np.append(np.genfromtxt(filenamecitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitp2,unpack=True,usecols=[0],dtype='int')+101)
    non_cit = np.append(np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_cit = np.append(np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitn2,unpack=True,usecols=[0],dtype='int')+101)    
    
    filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    filenamefpolp2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    pos_fpol = np.append(np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolp2,unpack=True,usecols=[0],dtype='int')+101)
    non_fpol = np.append(np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolno2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fpol = np.append(np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpoln2,unpack=True,usecols=[0],dtype='int')+101) 
    
    ### Concatonate runs
    runs = [tascit,tasfpol]
    
    ### Separate per periods (December)
    if period == 'D':
        tas_mo= np.empty((2,tascit.shape[0],tascit.shape[2],tascit.shape[3],
                          tascit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:,:]
    elif period == 'N':
        tas_mo= np.empty((2,tascit.shape[0],tascit.shape[2],tascit.shape[3],
                          tascit.shape[4]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,-1,:,:,:]
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mocitpos = tas_mo[0][pos_cit,:,:]
    tas_mofpolpos = tas_mo[1][pos_fpol,:,:]
    
    tas_mocitnon = tas_mo[0][non_cit,:,:]
    tas_mofpolnon = tas_mo[1][non_fpol,:,:]
    
    tas_mocitneg = tas_mo[0][neg_cit,:,:]
    tas_mofpolneg = tas_mo[1][neg_fpol,:,:]
    
    ### Take at 60N 
    latq = np.where((lat>59) & (lat<61))[0]
    
    ### Compute comparisons for fpolcit - taken ensemble average at 60N (index 79)
    ### +QBO
    diff_FPOLCITpos = np.nanmean(tas_mofpolpos - tas_mocitpos,axis=0)
    diffrunspos = diff_FPOLCITpos[:,latq,:].squeeze() 
    historicalforcedqpos = np.nanmean(tas_mocitpos[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcedpos.append(historicalforcedqpos)
    futureforcedpos.append(diffrunspos)
    
    ### Neutral QBO
    diff_FPOLCITnon = np.nanmean(tas_mofpolnon - tas_mocitnon,axis=0)
    diffrunsnon = diff_FPOLCITnon[:,latq,:].squeeze() 
    historicalforcedqnon = np.nanmean(tas_mocitnon[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcednon.append(historicalforcedqnon)
    futureforcednon.append(diffrunsnon)
    
    ### -QBO
    diff_FPOLCITneg = np.nanmean(tas_mofpolneg - tas_mocitneg,axis=0)
    diffrunsneg = diff_FPOLCITneg[:,latq,:].squeeze() 
    historicalforcedqneg = np.nanmean(tas_mocitneg[:,:,latq,:].squeeze(),axis=0)
    
    historicalforcedneg.append(historicalforcedqneg)
    futureforcedneg.append(diffrunsneg)
    
    ### Latitudes append
    lonssq.append(lon1)
    
###########################################################################
###########################################################################
###########################################################################
#### Plot climatological waves
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Append all data
historicalforced = historicalforcedpos + historicalforcednon + historicalforcedneg
futureforced = futureforcedpos + futureforcednon + futureforcedneg
lonss = lonssq + lonssq + lonssq

### Set limits for contours and colorbars
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])

fig = plt.figure()
for i in range(9):
    ax1 = plt.subplot(3,3,i+1)
    
    ### Calculate correlations
    corrpos = UT.calc_spatialCorrHeight(historicalforced[i],futureforced[i],
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
    if i==2 or i==5 or i==8:
        cs1 = plt.contour(lonq,levq,historicalforced[i],5,
                          colors='k',linewidths=1) 
    else:
        cs1 = plt.contour(lonq,levq,historicalforced[i],20,
                  colors='k',linewidths=1) 
    
    if i==0 or i==3 or i==6:
        qbophaseq = [r'QBO-W',r'QBO-W',r'QBO-W',r'QBO-N',r'QBO-N',r'QBO-N',
                     r'QBO-E',r'QBO-E',r'QBO-E']
        ax1.text(-0.28,0.5,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes,rotation=90)
    if i < 3:
        wavephaseq = ['Total Wave','Wave 1','Wave 2']
        ax1.text(0.5,1.2,r'\textbf{%s}' % wavephaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes)
    ax1.text(0.9,1.06,r'\textbf{R=%s}' % str(corrpos)[:4],color='k',
         fontsize=8,rotation=0,ha='center',va='center',
         transform=ax1.transAxes)
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    xxlabels = ['0','60E','120E','180','120W','60W','0']
    
    if i==6:
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,361,60),xxlabels,fontsize=6)
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
    elif i==7:
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,361,60),xxlabels,fontsize=6)
        plt.xlim([0,360])
        plt.yticks([])
        plt.minorticks_off()
        plt.xlabel(r'\textbf{Longitude ($^\circ$)}',color='k',fontsize=8)
    elif i==8:
        plt.ylim([1000,10])
        plt.xticks(np.arange(0,361,60),xxlabels,fontsize=6)
        plt.xlim([0,360])
        plt.yticks([])
        plt.minorticks_off()
    elif i==0 or i==3 or i==6:
        plt.ylim([1000,10])
        plt.xticks([])
        plt.xlim([0,360])
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
    else:
        plt.ylim([1000,10])
        plt.xticks([])
        plt.xlim([0,360])
        plt.yticks([])
        plt.minorticks_off()
    
    cmap = cmocean.cm.balance            
    cs.set_cmap(cmap) 
    
    plt.subplots_adjust(wspace=0.25)
    
plt.savefig(directoryfigure + 'linearInterference_November_FPOLCIT.png',dpi=300)
print('Completed: Script done!')


print('Completed: Script done!')

