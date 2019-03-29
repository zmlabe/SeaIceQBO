"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Composites are organized by QBO phase (positive, neutral, negative). 
Script calculates the zonal index for selection regions.

Notes
-----
    Author : Zachary Labe
    Date   : 15 August 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
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
print('\n' '----Calculating Zonal Index - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
def readVariables(varnames,period,region):
    lat,lon,time,lev,tashit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','surface')
    lat,lon,time,lev,tasfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','surface')
    
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
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
    ### Separate per periods (ON,DJ,FM)
    if period == 'D':
        tas_mo= np.empty((2,tashit.shape[0],90,tashit.shape[2],tashit.shape[3]))
        for i in range(len(runs)):
            tas_mo[i] = runs[i][:,60:150,:,:]
#            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
    else:
        ValueError('Wrong period selected!')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:,:] 
    tas_mofictpos = tas_mo[1][pos_fict,:,:,:] 
    
    tas_mohitnon = tas_mo[0][non_hit,:,:,:] 
    tas_mofictnon = tas_mo[1][non_fict,:,:,:] 
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:,:] 
    tas_mofictneg = tas_mo[1][neg_fict,:,:,:] 
    
    ### Compute comparisons for months - select region
    if region == 'Atlantic':
        lonq = np.append(np.where((lon >=310) & (lon <=360))[0],
                             np.where((lon >=0) & (lon <=20))[0],axis=0)
    elif region == 'Pacific':
        lonq = np.where((lon >= 120) & (lon <= 170))[0]
        
    ficthitpos = tas_mofictpos[:,:,:,lonq]
    ficthitnon = tas_mofictnon[:,:,:,lonq] 
    ficthitneg = tas_mofictneg[:,:,:,lonq]
    
    ### Calculate upper lats
    latqu = np.where((lat >=60) & (lat <=90))[0]
    ficthitposu = ficthitpos[:,:,latqu]
    ficthitnonu = ficthitnon[:,:,latqu] 
    ficthitnegu = ficthitneg[:,:,latqu]
    lat2squ = lat2[latqu,:]
    lat2su = lat2squ[:,lonq]
    ficthitposuu = UT.calc_weightedAve(ficthitposu,lat2su)
    ficthitnonuu = UT.calc_weightedAve(ficthitnonu,lat2su)
    ficthitneguu = UT.calc_weightedAve(ficthitnegu,lat2su)
    
    ### Calculate lower lats
    latql = np.where((lat >=20) & (lat <=50))[0]
    ficthitposl = ficthitpos[:,:,latql]
    ficthitnonl = ficthitnon[:,:,latql] 
    ficthitnegl = ficthitneg[:,:,latql]
    lat2sql = lat2[latql,:]
    lat2sl = lat2sql[:,lonq]
    ficthitposll = UT.calc_weightedAve(ficthitposl,lat2sl)
    ficthitnonll = UT.calc_weightedAve(ficthitnonl,lat2sl)
    ficthitnegll = UT.calc_weightedAve(ficthitnegl,lat2sl)
    
    ### Calculate Zonal Index (Z500u - Z500l)
    zdiffpos = ficthitposll - ficthitposuu
    zdiffnon = ficthitnonll - ficthitnonuu 
    zdiffneg = ficthitnegll - ficthitneguu
    diffruns_fict = [zdiffpos,zdiffnon,zdiffneg]
    
    ###########################################################################
    ### Calculate for HIT
    hitpos = tas_mohitpos[:,:,:,lonq]
    hitnon = tas_mohitnon[:,:,:,lonq] 
    hitneg = tas_mohitneg[:,:,:,lonq]
    
    ### Calculate upper lats
    latqu = np.where((lat >=60) & (lat <=90))[0]
    hitposu = hitpos[:,:,latqu]
    hitnonu = hitnon[:,:,latqu] 
    hitnegu = hitneg[:,:,latqu]
    lat2squ = lat2[latqu,:]
    lat2su = lat2squ[:,lonq]
    hitposuu = UT.calc_weightedAve(hitposu,lat2su)
    hitnonuu = UT.calc_weightedAve(hitnonu,lat2su)
    hitneguu = UT.calc_weightedAve(hitnegu,lat2su)
    
    ### Calculate lower lats
    latql = np.where((lat >=20) & (lat <=50))[0]
    hitposl = hitpos[:,:,latql]
    hitnonl = hitnon[:,:,latql] 
    hitnegl = hitneg[:,:,latql]
    lat2sql = lat2[latql,:]
    lat2sl = lat2sql[:,lonq]
    hitposll = UT.calc_weightedAve(hitposl,lat2sl)
    hitnonll = UT.calc_weightedAve(hitnonl,lat2sl)
    hitnegll = UT.calc_weightedAve(hitnegl,lat2sl)
    
    ### Calculate Zonal Index (Z500u - Z500l)
    zdiffposh = hitposll - hitposuu
    zdiffnonh = hitnonll - hitnonuu 
    zdiffnegh = hitnegll - hitneguu
    diffruns_hit = [zdiffposh,zdiffnonh,zdiffnegh]
    
    return diffruns_fict,diffruns_hit,lat,lon

#zdiff_ficta,zdiff_hita,lat,lon = readVariables('Z500',period,'Atlantic')
#zdiff_fictp,zdiff_hitp,lat,lon = readVariables('Z500',period,'Pacific')
#
#### Calculate ensemble means
#fictmposa = np.nanmean(zdiff_ficta[0],axis=0)
#fictmnona = np.nanmean(zdiff_ficta[1],axis=0)
#fictmnega = np.nanmean(zdiff_ficta[2],axis=0)
#
#hitmposa = np.nanmean(zdiff_hita[0],axis=0)
#hitmnona = np.nanmean(zdiff_hita[1],axis=0)
#hitmnega = np.nanmean(zdiff_hita[2],axis=0)
#
#fictmposp = np.nanmean(zdiff_fictp[0],axis=0)
#fictmnonp = np.nanmean(zdiff_fictp[1],axis=0)
#fictmnegp = np.nanmean(zdiff_fictp[2],axis=0)
#
#hitmposp = np.nanmean(zdiff_hitp[0],axis=0)
#hitmnonp = np.nanmean(zdiff_hitp[1],axis=0)
#hitmnegp = np.nanmean(zdiff_hitp[2],axis=0)

#############################################################################
#############################################################################
#############################################################################
### Plot daily zonal index
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 2))
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
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmposa,color='crimson',label=r'\textbf{QBO-W}',
         linewidth=3)
plt.plot(fictmnega,color='deepskyblue',label=r'\textbf{QBO-E}',
         linewidth=3)

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
#plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.xlim([0,90])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'ZonalIndex_H-Teffect_Atlantic.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmposp,color='crimson',label=r'\textbf{QBO-W}',
         linewidth=3)
plt.plot(fictmnegp,color='deepskyblue',label=r'\textbf{QBO-E}',
         linewidth=3)

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
#plt.yticks(np.arange(-1,1.1,0.5),map(str,np.round(np.arange(-1,1.1,0.5),2)),fontsize=8)
plt.xlim([0,90])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,0))

plt.savefig(directoryfigure + 'ZonalIndex_H-Teffect_Pacific.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmposa - hitmposa,color='crimson',label=r'\textbf{QBO-W}',
         linewidth=3)
plt.plot(fictmnega - hitmnega,color='deepskyblue',label=r'\textbf{QBO-E}',
         linewidth=3)

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
plt.yticks(np.arange(-100,101,10),map(str,np.arange(-100,101,10)),fontsize=8)
plt.xlim([0,90])
plt.ylim([-70,10])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,1))

plt.savefig(directoryfigure + 'ZonalIndex_H-Teffect_Atlantic_forced.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmposp - hitmposp,color='crimson',label=r'\textbf{QBO-W}',
         linewidth=3)
plt.plot(fictmnegp - hitmnegp,color='deepskyblue',label=r'\textbf{QBO-E}',
         linewidth=3)

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
#plt.yticks(np.arange(0,1001,10),map(str,np.arange(0,1001,10)),fontsize=8)
#plt.xlim([0,90])
#plt.ylim([320,550])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0.5,1))

plt.savefig(directoryfigure + 'ZonalIndex_H-Teffect_Pacific_forced.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################
fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmnega,color='deepskyblue',label=r'\textbf{FICTe}',
         linewidth=3)
plt.plot(hitmnega,color='deepskyblue',label=r'\textbf{HITe}',
         linewidth=2,linestyle='--',dashes=(1,0.3))

plt.plot(fictmposa,color='crimson',label=r'\textbf{FICTw}',
         linewidth=3)
plt.plot(hitmposa,color='crimson',label=r'\textbf{HITw}',
         linewidth=2,linestyle='--',dashes=(1,0.3))

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
plt.yticks(np.arange(0,1001,10),map(str,np.arange(0,1001,10)),fontsize=8)
plt.xlim([0,90])
plt.ylim([420,570])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,1))

plt.savefig(directoryfigure + 'ZonalIndex_Atlantic_forced2.png',
            dpi=300)

###############################################################################
###############################################################################
###############################################################################

fig = plt.figure()
ax = plt.subplot(1,1,1) 
       
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',pad=1)

plt.plot(fictmnegp,color='deepskyblue',label=r'\textbf{FICTe}',
         linewidth=3)
plt.plot(hitmnegp,color='deepskyblue',label=r'\textbf{HITe}',
         linewidth=2,linestyle='--',dashes=(1,0.3))

plt.plot(fictmposp,color='crimson',label=r'\textbf{FICTw}',
         linewidth=3)
plt.plot(hitmposp,color='crimson',label=r'\textbf{HITw}',
         linewidth=2,linestyle='--',dashes=(1,0.3))

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,91,30),xlabels,fontsize=8)
plt.yticks(np.arange(0,1001,10),map(str,np.arange(0,1001,10)),fontsize=8)
plt.xlim([0,90])
plt.ylim([320,550])

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=4,bbox_to_anchor=(0.5,1))

plt.savefig(directoryfigure + 'ZonalIndex_Pacific_forced2.png',
            dpi=300)