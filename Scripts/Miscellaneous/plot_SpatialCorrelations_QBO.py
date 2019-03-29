"""
Plots OCT-MAR spatial correlations -- test script so far!

Notes
-----
    Author : Zachary Labe
    Date   : 16 November 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput as MO
import cmocean
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/Documents/Research/SITperturb/Data/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SITperturb/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting spatial correlations - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

months = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR']
varnames = ['U10','Z30','Z500','SLP','T2M','THICK']
qbophase = ['pos','non','neg']
corrvarpos = []
corrvarnon = []
corrvarneg = []
for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'HIT','surface')
    lat,lon,time,lev,tasfit = MO.readExperi(directorydata,
                                            '%s' % varnames[v],'FIT','surface')
    lat,lon,time,lev,tasfict = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'FIC','surface')
    lat,lon,time,lev,tasfic = MO.readExperi(directorydata,
                                             '%s' % varnames[v],'CIT','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    pos_fit = np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int')
    non_fit = np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int')
    neg_fit = np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int')
    
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int')
    non_hit = np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int')
    neg_hit = np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int')
    
    filenameficp = directorydata + 'FIC/monthly/QBO_%s_FIC.txt' % qbophase[0]
    filenameficno = directorydata + 'FIC/monthly/QBO_%s_FIC.txt' % qbophase[1]
    filenameficn = directorydata + 'FIC/monthly/QBO_%s_FIC.txt' % qbophase[2]
    pos_fic = np.genfromtxt(filenameficp,unpack=True,usecols=[0],dtype='int')
    non_fic = np.genfromtxt(filenameficno,unpack=True,usecols=[0],dtype='int')
    neg_fic = np.genfromtxt(filenameficn,unpack=True,usecols=[0],dtype='int')
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int')
    non_fict = np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int')
    neg_fict = np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int')
    
    ### Separate per months
    varmo_fit = np.append(tasfit[:,9:,:,:],tasfit[:,0:3,:,:],
             axis=1)
    varmo_hit = np.append(tashit[:,9:,:,:],tashit[:,0:3,:,:],
             axis=1)
    varmo_fict = np.append(tasfict[:,9:,:,:],tasfict[:,0:3,:,:],
              axis=1)
    varmo_fic = np.append(tasfic[:,9:,:,:],tasfic[:,0:3,:,:],
          axis=1)
        
    ### Composite by QBO phase    
    tas_mofitpos = varmo_fit[pos_fit,:,:,:]
    tas_mohitpos = varmo_hit[pos_hit,:,:,:]
    tas_moficpos = varmo_fic[pos_fic,:,:,:]
    tas_mofictpos = varmo_fict[pos_fict,:,:,:]
    
    tas_mofitnon = varmo_fit[non_fit,:,:,:]
    tas_mohitnon = varmo_hit[non_hit,:,:,:]
    tas_moficnon = varmo_fic[non_fic,:,:,:]
    tas_mofictnon = varmo_fict[non_fict,:,:,:]
    
    tas_mofitneg = varmo_fit[neg_fit,:,:,:]
    tas_mohitneg = varmo_hit[neg_hit,:,:,:]
    tas_moficneg = varmo_fic[neg_fic,:,:,:]
    tas_mofictneg = varmo_fict[neg_fict,:,:,:]
    
    ### Calculate differences [FIT-HIT and FICT - FIT]
    fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)
    fithitnon = np.nanmean(tas_mofitnon - tas_mohitnon,axis=0)
    fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)
    
    fictficpos = np.nanmean(tas_mofictpos - tas_moficpos,axis=0)
    fictficnon = np.nanmean(tas_mofictnon - tas_moficnon,axis=0)
    fictficneg = np.nanmean(tas_mofictneg - tas_moficneg,axis=0)
    
    corrsposn = []
    for i in range(fithitpos.shape[0]):
        corrsqpos = UT.calc_spatialCorr(fithitpos[i],fictficpos[i],lat,lon,'yes')
        corrsposn.append(corrsqpos)
    corrvarpos.append(corrsposn)
    
    corrsnonn = []
    for i in range(fithitnon.shape[0]):
        corrsqnon = UT.calc_spatialCorr(fithitnon[i],fictficnon[i],lat,lon,'yes')
        corrsnonn.append(corrsqnon)
    corrvarnon.append(corrsnonn)
    
    corrsnegn = []
    for i in range(fithitneg.shape[0]):
        corrsqneg = UT.calc_spatialCorr(fithitneg[i],fictficneg[i],lat,lon,'yes')
        corrsnegn.append(corrsqneg)
    corrvarneg.append(corrsnegn)

corrvarpos = np.asarray(corrvarpos)
corrvarnon = np.asarray(corrvarnon)
corrvarneg = np.asarray(corrvarneg)

corrvar = [corrvarpos,corrvarnon,corrvarneg]

#### Save file
np.savetxt(directorydata2 + 'patterncorr_qbo_pos.txt',
           corrvarpos.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FICT-FIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]. \n' \
           'This is for composites of QBO-W',newline='\n\n')
np.savetxt(directorydata2 + 'patterncorr_qbo_neg.txt',
           corrvarneg.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FICT-FIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]. \n' \
           'This is for composites of QBO-E',newline='\n\n')
np.savetxt(directorydata2 + 'patterncorr_qbo_non.txt',
           corrvarnon.transpose(),delimiter=',',
           fmt='%3.2f',header='  '.join(varnames)+'\n',
           footer='\n File contains pearsonr correlation coefficients' \
           '\n between FIT-HIT and FICT-FIT to get the relative \n' \
           ' contributions of SIT and SIC [monthly, OCT-MAR]. \n' \
           'This is for composites of QBO-Neutral',newline='\n\n')

###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

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
for r in range(len(corrvar)):
    ax = plt.subplot(3,1,r+1)
    
    adjust_spines(ax, ['left', 'bottom'])
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('dimgrey')
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
    
    if r < 2:
        ax.spines['bottom'].set_color('w')
        ax.tick_params('x',length=0,width=0,which='major',color='w')
    
    plt.plot([0]*len(corrvar[r]),linewidth=2,color='dimgrey',linestyle='--')
    
    color=iter(ncm.cmap('MPL_gnuplot2')(np.linspace(0,0.8,len(corrvar[r]))))
    for i in range(len(corrvar[r])):
        c=next(color)
        plt.plot(corrvar[r][i],linewidth=1.5,color=c,alpha=1,
                 label = r'\textbf{%s}' % varnames[i],linestyle='-',
                 marker='o',markersize=3)
        
    if r == 1:
        plt.ylabel(r'\textbf{Pattern Correlation [R]}',color='dimgrey',
                             fontsize=13)
    if r == 2:
        plt.legend(shadow=False,fontsize=9,loc='lower center',
                   fancybox=True,frameon=False,ncol=3,
                   bbox_to_anchor=(0.5,-1.04))
        
    qbophaseq = [r'QBO-W',r'QBO-N',r'QBO-E']
    ax.annotate(r'\textbf{%s}' % qbophaseq[r],xy=(0,0),xytext=(1.03,0.5),
                 textcoords='axes fraction',color='dimgray',
                 fontsize=14,rotation=270,ha='center',va='center')
    
    plt.yticks(np.arange(-1,1.1,0.5),list(map(str,np.arange(-1,1.1,0.5))),
               fontsize=8)
    plt.ylim([-1,1])
    
    xlabels = [r'OCT',r'NOV',r'DEC',r'JAN',r'FEB',r'MAR',r'APR']
    plt.xticks(np.arange(0,6,1),xlabels,fontsize=8)
    plt.xlim([0,5])
    
    if r<2:
        ax.tick_params(labelbottom='off') 
    
    ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.3)
    
    plt.subplots_adjust(bottom=0.2,hspace=0.25)

plt.savefig(directoryfigure + 'patterncorrs_monthly_qbo.png',dpi=300)
