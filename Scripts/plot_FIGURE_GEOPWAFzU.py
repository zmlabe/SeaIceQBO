"""
Plot manuscript figure for GEOP, WAFz, and change in zonal wind

Notes
-----
    Author : Zachary Labe
    Date   : 26 March 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput_AllMembers as MO
import read_DailyOutput_AllMembers as DO
import cmocean
import calc_Utilities as UT

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO manuscript Figure 2- %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

runnames = [r'HIT',r'FICT']
experiments = [r'\textbf{$\Delta$NET}']
qbophase = ['pos','non','neg']
period = 'ND'
letters = ["a","b","c","d","e","f","g","h","i"]

### Call function for WAFz data for each run
def readWAFz(varnames,qbophase,period):
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[0],'HIT',
                                               'profile')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[0],'FICT',
                                                'profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
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
    tas_mofictposq = tas_mo[1][pos_fict,:,:]
    
    tas_mohitnegq = tas_mo[0][neg_hit,:,:]
    tas_mofictnegq = tas_mo[1][neg_fict,:,:]
    
    ### Take zonal average
    tas_mohitpos = np.nanmean(tas_mohitposq,axis=3)
    tas_mofictpos = np.nanmean(tas_mofictposq,axis=3)
    
    tas_mohitneg = np.nanmean(tas_mohitnegq,axis=3)
    tas_mofictneg = np.nanmean(tas_mofictnegq,axis=3)
    
    
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)/np.nanstd(tas_mofictpos,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)/np.nanstd(tas_mofictneg,axis=0)
    diffruns_mo = [ficthitpos,ficthitneg]
    
    ### Calculate significance 
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mohitpos,
                                                         tas_mofictpos)
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mohitneg,
                                                         tas_mofictneg)
    
    pruns_mo = [pvalue_FICTHITpos,pvalue_FICTHITneg]
    return diffruns_mo,pruns_mo,lat,lon,lev

def readU(varnames,qbophase,period):
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[0],'HIT',
                                               'profile')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[0],'FICT',
                                                'profile')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
    filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
    pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+101)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)    
    ### Concatonate runs
    runs = [tashit,tasfict]
    
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
    tas_mofictposq = tas_mo[1][pos_fict,:,:]
    
    tas_mohitnegq = tas_mo[0][neg_hit,:,:]
    tas_mofictnegq = tas_mo[1][neg_fict,:,:]
    
    ### Take zonal average
    tas_mohitpos = np.nanmean(tas_mohitposq,axis=3)
    tas_mofictpos = np.nanmean(tas_mofictposq,axis=3)
    
    tas_mohitneg = np.nanmean(tas_mohitnegq,axis=3)
    tas_mofictneg = np.nanmean(tas_mofictnegq,axis=3)
    
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    diffruns_mo = [ficthitpos,ficthitneg]
    
    ### Calculate significance 
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mohitpos,
                                                         tas_mofictpos)
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mohitneg,
                                                         tas_mofictneg)
    
    pruns_mo = [pvalue_FICTHITpos,pvalue_FICTHITneg]
    return diffruns_mo,pruns_mo,lat,lon,lev

### Call functions for variable profile data for polar cap
def readPolarCap(varnames,qbophase):
    for v in range(len(varnames)):
        lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames[v],
                                                    'HIT','profile')
        lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames[v],
                                                    'FICT','profile')
        
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
        var_mo = [varhit,varfict]
        
        ### Composite by QBO phase    
        var_mohitpos = var_mo[0][pos_hit,:]
        var_mofictpos = var_mo[1][pos_fict,:]
        
        var_mohitneg = var_mo[0][neg_hit,:]
        var_mofictneg = var_mo[1][neg_fict,:]
        
        ### Compute comparisons for months - taken ensemble average
        ficthitpos = np.nanmean(var_mofictpos - var_mohitpos,axis=0)
        ficthitneg = np.nanmean(var_mofictneg - var_mohitneg,axis=0)
        
        diffruns = [ficthitpos,ficthitneg]
        
        ### Calculate significance
        stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(var_mo[1][pos_fict,:],
                                                             var_mo[0][pos_hit,:])
        stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(var_mo[1][neg_fict,:],
                                                             var_mo[0][neg_hit,:])
    
        pruns = [pvalue_FICTHITpos,pvalue_FICTHITneg]
        return diffruns,pruns
        
### Read data
#diffruns_mo,pruns_mo,lat,lon,lev = readWAFz(['WAFZ'],qbophase,period)
#u_mo,u_p,lat,lon,lev = readU(['U'],qbophase,period)
#diffrungeop,prunsgeop = readPolarCap(['GEOP'],qbophase)

###########################################################################
###########################################################################
###########################################################################
#### Plot EP Flux
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-0.5,0.51,0.025)
barlim = np.arange(-0.5,0.6,0.5)
limitg = np.arange(-150,150.1,5)
barlimg = np.arange(-150,151,75)
timeq = np.arange(0,212,1)
timeqq,levqq = np.meshgrid(timeq,lev) 
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(4):
    ax1 = plt.subplot(2,2,i+1)
    
    if i >= 2:
        var = diffruns_mo[i-2]
        pruns = pruns_mo[i-2]
    
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
    
        plt.contourf(lat,lev,pruns,colors='None',hatches=['////'],
                     linewidth=5)       
        plt.contour(lat,lev,u_mo[i-2],np.arange(-90,91,0.25),colors='dimgrey',linewidths=0.7)
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        plt.xticks(np.arange(0,96,10),map(str,np.arange(0,91,10)),fontsize=6)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
        
        plt.xlim([10,90])
        plt.ylim([1000,10])
            
        cmap = cmocean.cm.balance          
        cs.set_cmap(cmap) 
    elif i < 2:
        var = diffrungeop[i]
        pvar = prunsgeop[i]
        
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
        
        csg = plt.contourf(timeq,lev,var.transpose(),limitg,extend='both')
        plt.contourf(timeqq,levqq,pvar.transpose(),colors='None',
                     hatches=['////'])                  
        
        plt.gca().invert_yaxis()
        plt.yscale('log',nonposy='clip')
        
        xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
        plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
        plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
        plt.minorticks_off()
        plt.xlim([30,210])
        plt.ylim([1000,10])
        
        cmap = cmocean.cm.balance           
        csg.set_cmap(cmap) 

    ### Add experiment text to subplot
    if i < 2:
        qbophaseq = [r'QBO-W',r'QBO-E']
        ax1.text(0.5,1.07,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes)
    if i == 0 or i == 2:
        plt.ylabel(r'\textbf{Pressure (hPa)}',color='k',fontsize=7,
                     labelpad=1)
    if i == 2 or i==3:
        plt.xlabel(r'\textbf{Latitude ($^\circ$N)}',color='k',fontsize=7,
                             labelpad=1)
    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
            xytext=(0.01,0.92),xycoords='axes fraction',
            color='k',fontsize=9)
        
cbar_axg = fig.add_axes([0.90,0.6,0.013,0.25])                
cbarg = fig.colorbar(csg,cax=cbar_axg,orientation='vertical',
                    extend='both',extendfrac=0.07,drawedges=False) 
cbarg.set_label(r'\textbf{m}',fontsize=8,color='k',labelpad=1.2)
cbarg.set_ticks(barlimg)
cbarg.set_ticklabels(list(map(str,barlimg))) 
cbarg.ax.tick_params(axis='y', size=.01,labelsize=8)
cbarg.outline.set_edgecolor('dimgrey')

cbar_ax = fig.add_axes([0.90,0.155,0.013,0.25])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                    extend='both',extendfrac=0.07,drawedges=False) 
cbar.set_label(r'\textbf{Normalized WAFz}',fontsize=8,color='k',labelpad=2.5)  
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='y', size=.01,labelsize=8)
cbar.outline.set_edgecolor('dimgrey')     

plt.tight_layout()    
fig.subplots_adjust(hspace=0.18,wspace=0.23,right=0.86)

plt.savefig(directoryfigure + 'Fig2.png',dpi=900)
print('Completed: Script done!')

