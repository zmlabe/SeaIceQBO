"""
Plot comparisons between SIT and SIC modeling experiments using 
WACCM4. Composites are organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2018
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
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'
#def readVariables(varnames,period,location):
#    ### Call function for surface temperature data from reach run
#    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames,'HIT',
#                                               'surface')
#    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames,'FICT',
#                                                'surface')
#    
#    ### Create 2d array of latitude and longitude
#    lon2,lat2 = np.meshgrid(lon,lat)
#    
#    ### Read in QBO phases 
#    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#    filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#    filenamehitno2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
#    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+101)
#    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
#    
#    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
#    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+101)
#    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
#                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)
#    
#    ### Concatonate runs
#    runs = [tashit,tasfict]
#    
#    ### Separate per periods (ON,DJ,FM)
#    if period == 'ON': 
#        tas_mo = np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = np.nanmean(runs[i][:,9:11,:,:],axis=1) 
#    elif period == 'DJ':     
#        tas_mo = np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i],tas_mo[i] = UT.calcDecJan(runs[i],runs[i],lat,
#                                                lon,'surface',1) 
#    elif period == 'FM':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = np.nanmean(runs[i][:,1:3,:,:],axis=1)
#    elif period == 'DJF':
#        tas_mo= np.empty((3,tashit.shape[0]-1,tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i],tas_mo[i] = UT.calcDecJanFeb(runs[i],runs[i],lat,
#                                                  lon,'surface',1)   
#    elif period == 'M':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = runs[i][:,2,:,:]
#    elif period == 'D':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = runs[i][:,-1,:,:]
#    elif period == 'N':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = runs[i][:,-2,:,:]
#    elif period == 'ND':
#        tas_mo= np.empty((3,tashit.shape[0],tashit.shape[2],tashit.shape[3]))
#        for i in range(len(runs)):
#            tas_mo[i] = np.nanmean(runs[i][:,-2:,:,:],axis=1)
#    else:
#        ValueError('Wrong period selected! (ON,DJ,FM)')
#        
#    ### Composite by QBO phase    
#    tas_mohitpos = tas_mo[0][pos_hit,:,:]
#    tas_mofictpos = tas_mo[1][pos_fict,:,:]
#    
#    tas_mohitnon = tas_mo[0][non_hit,:,:]
#    tas_mofictnon = tas_mo[1][non_fict,:,:]
#    
#    tas_mohitneg = tas_mo[0][neg_hit,:,:]
#    tas_mofictneg = tas_mo[1][neg_fict,:,:]
#    
#    ### Compute comparisons for months - select region
#    if location == 'Atlantic':
#        lonq = np.append(np.where((lon >=310) & (lon <=360))[0],
#                         np.where((lon >=0) & (lon <=20))[0],axis=0)
#        ficthitpos = np.nanmean(tas_mofictpos[:,:,lonq],axis=2) - np.nanmean(tas_mohitpos[:,:,lonq],axis=2)
#        ficthitnon = np.nanmean(tas_mofictnon[:,:,lonq],axis=2) - np.nanmean(tas_mohitnon[:,:,lonq],axis=2)
#        ficthitneg = np.nanmean(tas_mofictneg[:,:,lonq],axis=2) - np.nanmean(tas_mohitneg[:,:,lonq],axis=2)
#    elif location == 'Pacific':
#        lonq = np.where((lon >= 120) & (lon <= 170))[0]
#        ficthitpos = np.nanmean(tas_mofictpos[:,:,lonq],axis=2) - np.nanmean(tas_mohitpos[:,:,lonq],axis=2)
#        ficthitnon = np.nanmean(tas_mofictnon[:,:,lonq],axis=2) - np.nanmean(tas_mohitnon[:,:,lonq],axis=2)
#        ficthitneg = np.nanmean(tas_mofictneg[:,:,lonq],axis=2) - np.nanmean(tas_mohitneg[:,:,lonq],axis=2)
#    elif location == 'Global':
#        ficthitpos = np.nanmean(tas_mofictpos[:,:,:],axis=2) - np.nanmean(tas_mohitpos[:,:,:],axis=2)
#        ficthitnon = np.nanmean(tas_mofictnon[:,:,:],axis=2) - np.nanmean(tas_mohitnon[:,:,:],axis=2)
#        ficthitneg = np.nanmean(tas_mofictneg[:,:,:],axis=2) - np.nanmean(tas_mohitneg[:,:,:],axis=2)
#    if varnames == 'U10':
#        latq = np.where((lat >= 58) & (lat <= 61))[0]
#        ficthitpos = np.nanmean(tas_mofictpos[:,latq,:],axis=1) - np.nanmean(tas_mohitpos[:,latq,:],axis=1)
#        ficthitnon = np.nanmean(tas_mofictnon[:,latq,:],axis=1) - np.nanmean(tas_mohitnon[:,latq,:],axis=1)
#        ficthitneg = np.nanmean(tas_mofictneg[:,latq,:],axis=1) - np.nanmean(tas_mohitneg[:,latq,:],axis=1)
#        ficthitpos = np.nanmean(ficthitpos[:,:],axis=1)
#        ficthitnon = np.nanmean(ficthitnon[:,:],axis=1)
#        ficthitneg = np.nanmean(ficthitneg[:,:],axis=1)
#    elif varnames == 'U200':
#        latq = np.where((lat >= 30) & (lat <= 75))[0]
#        ficthitpos = np.nanmean(ficthitpos[:,latq],axis=1)
#        ficthitnon = np.nanmean(ficthitnon[:,latq],axis=1)
#        ficthitneg = np.nanmean(ficthitneg[:,latq],axis=1)
#    elif varnames == 'Z30':
#        latq = np.where((lat >= 65) & (lat <= 90))[0]
#        lat2s = lat2[latq,:]
#        ficthitpos = tas_mofictpos[:,latq,:] - tas_mohitpos[:,latq,:]
#        ficthitnon = tas_mofictnon[:,latq,:] - tas_mohitnon[:,latq,:]
#        ficthitneg = tas_mofictneg[:,latq,:] - tas_mohitneg[:,latq,:]
#        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
#        ficthitnon = UT.calc_weightedAve(ficthitnon,lat2s)
#        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
#    diffruns = [ficthitpos.squeeze(),ficthitnon.squeeze(),ficthitneg.squeeze()]
#    
#    return diffruns,lat,lon,lev

def readVariables(varnames,period,location):
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
    
    tas_mohitnon = tas_mo[0][non_hit,:,:]
    tas_mofictnon = tas_mo[1][non_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofictneg = tas_mo[1][neg_fict,:,:]
    
    ### Compute comparisons for months - select region
    if location == 'Atlantic':
        lonq = np.append(np.where((lon >=310) & (lon <=360))[0],
                         np.where((lon >=0) & (lon <=20))[0],axis=0)
        ficthitpos = tas_mohitpos[:,:,lonq]
        ficthitnon = tas_mohitnon[:,:,lonq] 
        ficthitneg = tas_mohitneg[:,:,lonq]
    elif location == 'Pacific':
        lonq = np.where((lon >= 120) & (lon <= 170))[0]
        ficthitpos = tas_mofictpos[:,:,lonq]
        ficthitnon = tas_mofictnon[:,:,lonq]
        ficthitneg = tas_mofictneg[:,:,lonq]
    if varnames == 'U10':
        latq = np.where((lat >= 58) & (lat <= 61))[0]
        ficthitpos = np.nanmean(tas_mofictpos[:,latq,:],axis=1) - np.nanmean(tas_mohitpos[:,latq,:],axis=1)
        ficthitnon = np.nanmean(tas_mofictnon[:,latq,:],axis=1) - np.nanmean(tas_mohitnon[:,latq,:],axis=1)
        ficthitneg = np.nanmean(tas_mofictneg[:,latq,:],axis=1) - np.nanmean(tas_mohitneg[:,latq,:],axis=1)
        ficthitpos = np.nanmean(ficthitpos[:,:],axis=1)
        ficthitnon = np.nanmean(ficthitnon[:,:],axis=1)
        ficthitneg = np.nanmean(ficthitneg[:,:],axis=1)
    elif varnames == 'U200':
        latq = np.where((lat >= 30) & (lat <= 75))[0]
        lat2s = lat2[latq,:]
        lat2s = lat2s[:,lonq]
        ficthitpos = ficthitpos[:,latq]
        ficthitnon = ficthitnon[:,latq]
        ficthitneg = ficthitneg[:,latq]
        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
        ficthitnon = UT.calc_weightedAve(ficthitnon,lat2s)
        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
    elif varnames == 'Z30':
        latq = np.where((lat >= 65) & (lat <= 90))[0]
        lat2s = lat2[latq,:]
        ficthitpos = tas_mohitpos[:,latq,:] 
        ficthitnon = tas_mohitnon[:,latq,:] 
        ficthitneg = tas_mohitneg[:,latq,:] 
        ficthitpos = UT.calc_weightedAve(ficthitpos,lat2s)
        ficthitnon = UT.calc_weightedAve(ficthitnon,lat2s)
        ficthitneg = UT.calc_weightedAve(ficthitneg,lat2s)
    diffruns = [ficthitpos.squeeze(),ficthitnon.squeeze(),ficthitneg.squeeze()]
    
    return diffruns,lat,lon,lev

u200a,lat,lon,lev = readVariables('U200',period,'Atlantic')
u200p,lat,lon,lev = readVariables('U200',period,'Pacific')
#u10g,lat,lon,lev = readVariables('U10',period,'Global')
z30g,lat,lon,lev = readVariables('Z30',period,'Polar')

################################################################################
################################################################################
################################################################################    
#### Plot scatter plots
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})
#
#### Adjust axes in time series plots 
#def adjust_spines(ax, spines):
#    for loc, spine in ax.spines.items():
#        if loc in spines:
#            spine.set_position(('outward', 5))
#        else:
#            spine.set_color('none')  
#    if 'left' in spines:
#        ax.yaxis.set_ticks_position('left')
#    else:
#        ax.yaxis.set_ticks([])
#
#    if 'bottom' in spines:
#        ax.xaxis.set_ticks_position('bottom')
#    else:
#        ax.xaxis.set_ticks([])
#
#fig = plt.figure(figsize=(12,7))
#
#ax = fig.add_subplot(121)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('dimgrey')
#ax.spines['bottom'].set_color('dimgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#
#plt.scatter(u10g[0],u200a[0],color='deepskyblue',marker='o',label=r'QBO-W')
#timex = np.linspace(-25,36,len(u10g[0]))
#slope, intercept, r_value, p_value1, std_err = sts.linregress(u10g[0],u200a[0])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='deepskyblue',linewidth=3)
#
#plt.scatter(u10g[1],u200a[1],color='gold',marker='o',label=r'QBO-N')
#timex = np.linspace(-25,36,len(u10g[1]))
#slope, intercept, r_value, p_value2, std_err = sts.linregress(u10g[1],u200a[1])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='gold',linewidth=3)
#
#plt.scatter(u10g[2],u200a[2],color='crimson',marker='o',label=r'QBO-E')
#timex = np.linspace(-25,36,len(u10g[2]))
#slope, intercept, r_value, p_value3, std_err = sts.linregress(u10g[2],u200a[2])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='crimson',linewidth=3)
#
#plt.yticks(np.arange(-10,11,2),list(map(str,np.arange(-10,11,2))),fontsize=10) 
#plt.ylim([-4,4])
#plt.xticks(np.arange(-25,40,5),list(map(str,np.arange(-25,40,5))),fontsize=10) 
#plt.xlim([-25,35])
#plt.title(r'\textbf{ATLANTIC}',fontsize=30,color='dimgrey')
#
#ax = fig.add_subplot(122)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('dimgrey')
#ax.spines['bottom'].set_color('dimgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
#
#plt.scatter(u10g[0],u200p[0],color='deepskyblue',marker='o',label=r'\textbf{QBO-W}')
#timex = np.linspace(-25,36,len(u10g[0]))
#slope, intercept, r_value, p_value4, std_err = sts.linregress(u10g[0],u200p[0])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='deepskyblue',linewidth=3)
#
#plt.scatter(u10g[1],u200p[1],color='gold',marker='o',label=r'\textbf{QBO-N}')
#timex = np.linspace(-25,36,len(u10g[1]))
#slope, intercept, r_value, p_value5, std_err = sts.linregress(u10g[1],u200p[1])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='gold',linewidth=3)
#
#plt.scatter(u10g[2],u200p[2],color='crimson',marker='o',label=r'\textbf{QBO-E}')
#timex = np.linspace(-25,36,len(u10g[2]))
#slope, intercept, r_value, p_value6, std_err = sts.linregress(u10g[2],u200p[2])
#line1 = slope*timex + intercept
#plt.plot(timex,line1,color='crimson',linewidth=3)
#
#plt.yticks(np.arange(-10,11,2),list(map(str,np.arange(-10,11,2))),fontsize=10) 
#plt.ylim([-4,4])
#plt.xticks(np.arange(-25,40,5),list(map(str,np.arange(-25,40,5))),fontsize=10) 
#plt.xlim([-25,35])
#plt.title(r'\textbf{PACIFIC}',fontsize=30,color='dimgrey')
#        
#plt.legend(shadow=False,fontsize=12,loc='lower center',
#           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(0, -0.15),
#           labelspacing=0.2,columnspacing=0.4,handletextpad=0.4)
#
#plt.savefig(directoryfigure + 'PolarVortexJetRelationship.png',dpi=300)

###############################################################################
###############################################################################
###############################################################################    
### Plot scatter plots
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

fig = plt.figure(figsize=(12,7))

ax = fig.add_subplot(121)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.scatter(z30g[0],u200a[0],color='deepskyblue',marker='o',label=r'QBO-W')
timex = np.linspace(22000,23201,len(z30g[0]))
slope, intercept, r_value1, p_value1, std_err = sts.linregress(z30g[0],u200a[0])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='deepskyblue',linewidth=3)

plt.scatter(z30g[1],u200a[1],color='gold',marker='o',label=r'QBO-N')
timex = np.linspace(22000,23201,len(z30g[1]))
slope, intercept, r_value2, p_value2, std_err = sts.linregress(z30g[1],u200a[1])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='gold',linewidth=3)

plt.scatter(z30g[2],u200a[2],color='crimson',marker='o',label=r'QBO-E')
timex = np.linspace(22000,23201,len(z30g[2]))
slope, intercept, r_value3, p_value3, std_err = sts.linregress(z30g[2],u200a[2])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='crimson',linewidth=3)

plt.yticks(np.arange(0,51,2),list(map(str,np.arange(0,51,2))),fontsize=10) 
plt.ylim([16,24])
plt.xticks(np.arange(22000,23201,200),list(map(str,np.arange(22000,23201,200))),fontsize=10) 
plt.xlim([22000,23200])
plt.title(r'\textbf{ATLANTIC}',fontsize=30,color='dimgrey')
plt.xlabel(r'\textbf{m}',fontsize=12)
plt.ylabel(r'\textbf{m/s}',fontsize=12)

ax = fig.add_subplot(122)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')

plt.scatter(z30g[0],u200p[0],color='deepskyblue',marker='o',label=r'\textbf{QBO-W}')
timex = np.linspace(22000,23201,len(z30g[0]))
slope, intercept, r_value4, p_value4, std_err = sts.linregress(z30g[0],u200p[0])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='deepskyblue',linewidth=3)

plt.scatter(z30g[1],u200p[1],color='gold',marker='o',label=r'\textbf{QBO-N}')
timex = np.linspace(22000,23201,len(z30g[1]))
slope, intercept, r_value5, p_value5, std_err = sts.linregress(z30g[1],u200p[1])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='gold',linewidth=3)

plt.scatter(z30g[2],u200p[2],color='crimson',marker='o',label=r'\textbf{QBO-E}')
timex = np.linspace(22000,23201,len(z30g[2]))
slope, intercept, r_value6, p_value6, std_err = sts.linregress(z30g[2],u200p[2])
line1 = slope*timex + intercept
plt.plot(timex,line1,color='crimson',linewidth=3)

plt.yticks(np.arange(0,51,2),list(map(str,np.arange(0,51,2))),fontsize=10) 
plt.ylim([25,35])
plt.xticks(np.arange(22000,23201,200),list(map(str,np.arange(22000,23201,200))),fontsize=10) 
plt.xlim([22000,23200])
plt.title(r'\textbf{PACIFIC}',fontsize=30,color='dimgrey')
plt.xlabel(r'\textbf{m}',fontsize=12)
plt.ylabel(r'\textbf{m/s}',fontsize=12)
        
plt.legend(shadow=False,fontsize=12,loc='lower center',
           fancybox=True,frameon=False,ncol=3,bbox_to_anchor=(-0.12, -0.15),
           labelspacing=0.2,columnspacing=0.4,handletextpad=0.4)

plt.savefig(directoryfigure + 'PolarVortexJetRelationship_HIT.png',dpi=300)