"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data for all wave activity flux (WAFz and WAFy).

Notes
-----
    Author : Zachary Labe
    Date   : 26 June 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import cmocean

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
print('\n' '----Plotting Daily Variables for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = False
varnames = ['WAFZ','WAFZ','WAFY','WAFY']
runnames = [r'HIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{WAFz}',r'\textbf{WAFz}',r'\textbf{WAFy}',
               r'\textbf{WAFy}']

def readWAF(varnames,runnames,experiments,qbophase):
    """
    Function reads in WAF data for listed experiments

    Parameters
    ----------
    varnames : string
        variable to download
    runnames : list of strings
        model experiments to read in
    experiments : list of strings
        model simulations to compare
    qbophase : list of strings
        list of qbo phases

    Returns
    -------
    diffruns : list of arrays
        arrays for each experiment variable
    pruns : list of arrays
        arrays of p-values for each experiment variable
    lev : 1d array
        leves

    Usage
    -----
    diffruns,pruns,lev = readWAF(varnames,runnames,experiments,qbophase)
    """
    print('\n>>> Using readWAF function!')
    
    ### Call functions for variable profile data for polar cap
    lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames,
                                                'HIT','profilesiberia')
    lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames,
                                                'FICT','profilesiberia')
    
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
                        np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],
                                      dtype='int')+100)
    non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],
                                      dtype='int')+100)
    neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],
                                      dtype='int')+100)    
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],
                                      dtype='int')+100)
    non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],
                                      dtype='int')+100)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],
                                      dtype='int')+100)   
    
    ### Concatonate runs
    var_mo = [varhit,varfict]
    
    ### Composite by QBO phase    
    var_mohitpos = var_mo[0][pos_hit,:]
    var_mofictpos = var_mo[1][pos_fict,:]
    
    var_mohitneg = var_mo[0][neg_hit,:]
    var_mofictneg = var_mo[1][neg_fict,:]
    
    ### Compute comparisons for months - taken ensemble average
    ficthitpos = np.nanmean(var_mofictpos - var_mohitpos,axis=0)/np.nanstd(var_mohitpos,axis=0)
    ficthitneg = np.nanmean(var_mofictneg - var_mohitneg,axis=0)/np.nanstd(var_mohitneg,axis=0)
    
    diffruns = [ficthitpos,ficthitneg]
    
    ### Calculate significance for days
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(var_mo[1][pos_fict,:],
                                                         var_mo[0][pos_hit,:])
    stat_FICTHITnon,pvalue_FICTHITnon = UT.calc_indttest(var_mo[1][non_fict,:],
                                                         var_mo[0][non_hit,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(var_mo[1][neg_fict,:],
                                                         var_mo[0][neg_hit,:])

    pruns = [pvalue_FICTHITpos,pvalue_FICTHITneg]
    
    print('\n*Completed: Finished readWAF function!')
    return diffruns,pruns,lev

### Read in data
diffwafz,pwafz,lev = readWAF('WAFZs',runnames,experiments,qbophase)
diffwafy,pwafy,lev = readWAF('WAFYs',runnames,experiments,qbophase)

### Assign to lists for plotting
diffruns = np.append(diffwafz,diffwafy,axis=0)
pruns = np.append(pwafz,pwafy,axis=0)
                                                 
############################################################################
############################################################################
############################################################################
##### Plot daily profile with height for selected variable
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
zscale = np.array([700,500,300,200,
                    100,50,30,10])
timeq = np.arange(0,212,1)
timeqq,levq = np.meshgrid(timeq,lev)

fig = plt.figure()
for i in range(4):
    ax1 = plt.subplot(2,2,i+1)
    
    var = diffruns[i]
    pvar = pruns[i]
    
    ### Set limits for contours and colorbars
#    if i <= 2:
#        limit = np.arange(-0.05,0.05001,0.001)
#        barlim = np.arange(-0.05,0.051,0.05)
#    elif i > 2:
#        limit = np.arange(-2,2.01,0.05)
#        barlim = np.arange(-2,3,2)
    limit = np.arange(-0.5,0.51,0.025)
    barlim = np.arange(-0.5,0.6,0.5)
    
    if MASK == True:
        pvar2 = pvar.copy()
        pvar2[np.isnan(pvar2)]=0.0
        var = var*pvar2
    
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
    
    cs = plt.contourf(timeq,lev,var.transpose(),limit,extend='both')
    
    if MASK == False:                  
        plt.contourf(timeqq,levq,pvar.transpose(),colors='None',
                     hatches=['////'])                  
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    xlabels = [r'Sep',r'Oct',r'Nov',r'Dec',r'Jan',r'Feb',r'Mar',r'Apr'] 
    plt.xticks(np.arange(0,212,30),xlabels,fontsize=6)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    plt.minorticks_off()
    plt.xlim([30,120])
    plt.ylim([700,10])
    
    cmap = cmocean.cm.balance
    cs.set_cmap(cmap)

    ### Add experiment text to subplot
    if i < 2:
        qbophaseq = [r'QBO-W',r'QBO-E']
        ax1.text(0.5,1.1,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes)
    if i == 0 or i == 2:
        ax1.text(-0.2,0.5,r'%s' % experiments[i],color='k',
                     fontsize=20,rotation=90,ha='center',va='center',
                     transform=ax1.transAxes)

    if i == 1:
        cbar_ax = fig.add_axes([0.92,0.60,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnames[i] == 'WAFZ':
            cbar.set_label(r'\textbf{Normalized}',fontsize=11,color='dimgray',
                           labelpad=0.5)
        elif varnames[i] == 'WAFY':
            cbar.set_label(r'\textbf{Normalized}',fontsize=11,color='dimgray',
                           labelpad=0.5)    
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=7) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
        
    elif i == 3:
        cbar_ax = fig.add_axes([0.92,0.19,0.015,0.2])                
        cbar = fig.colorbar(cs,cax=cbar_ax,orientation='vertical',
                            extend='both',extendfrac=0.07,drawedges=False)    
        if varnames[i] == 'WAFZ':
            cbar.set_label(r'\textbf{Normalized}',fontsize=11,color='dimgray',
                           labelpad=0.5)
        elif varnames[i] == 'WAFY':
            cbar.set_label(r'\textbf{Normalized}',fontsize=11,color='dimgray',
                           labelpad=0.5)     
        cbar.set_ticks(barlim)
        cbar.set_ticklabels(list(map(str,barlim)))
        cbar.ax.tick_params(labelsize=6,pad=8) 
        ticklabs = cbar.ax.get_yticklabels()
        cbar.ax.set_yticklabels(ticklabs,ha='center')
        cbar.ax.tick_params(axis='y', size=.001)
        cbar.outline.set_edgecolor('dimgrey')
        cbar.outline.set_linewidth(0.5)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')

plt.subplots_adjust(wspace=0.22)

if MASK == True:
    plt.savefig(directoryfigure + 'allExperiments_WAF_MASK_daily_FICTHIT_Siberia.png',
                dpi=300)
else:
    plt.savefig(directoryfigure + 'allExperiments_WAF_daily_FICTHIT_Siberia.png',
        dpi=300)
print('Completed: Script done!')                                