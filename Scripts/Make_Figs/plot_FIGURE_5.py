"""
Plot manuscript figure for U10 comparing SIT and NET experiments

Notes
-----
    Author : Zachary Labe
    Date   : 21 March 2019
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
varnames = ['U10']
runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{$\Delta$SIT}',r'\textbf{$\Delta$SIT}',
               r'\textbf{$\Delta$NET}',r'\textbf{$\Delta$NET}']
qbophase = ['pos','non','neg']
period = 'D'
letters = ["a","b","c","d","e","f","g","h","i"]

def calc_indttest_90(varx,vary):
    """
    Function calculates statistical difference for 2 independent
    sample t-test

    Parameters
    ----------
    varx : 3d array
    vary : 3d array
    
    Returns
    -------
    stat = calculated t-statistic
    pvalue = two-tailed p-value

    Usage
    -----
    stat,pvalue = calc_ttest(varx,vary)
    """
    print('\n>>> Using calc_ttest function!')
    
    ### Import modules
    import numpy as np
    import scipy.stats as sts
    
    ### 2-independent sample t-test
    stat,pvalue = sts.ttest_ind(varx,vary,nan_policy='omit')
    
    ### Significant at 90% confidence level
    pvalue[np.where(pvalue >= 0.1)] = np.nan
    pvalue[np.where(pvalue < 0.1)] = 1.
    pvalue[np.isnan(pvalue)] = 0.
    
    print('*Completed: Finished calc_ttest function!')
    return stat,pvalue

for v in range(len(varnames)):
    ### Call function for surface temperature data from reach run
    lat,lon,time,lev,tashit = MO.readExperiAll('%s' % varnames[v],'HIT',
                                               'surface')
    lat,lon,time,lev,tasfit = MO.readExperiAll('%s' % varnames[v],'FIT',
                                               'surface')
    lat,lon,time,lev,tasfict = MO.readExperiAll('%s' % varnames[v],'FICT',
                                                'surface')
    
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

    filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
    filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
    pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+101)
    
    filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
    filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
    pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+101)
    neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+101)
    
    ### Concatonate runs
    runs = [tashit,tasfit,tasfict]
    
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
    tas_mofitpos = tas_mo[1][pos_fit,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofitneg = tas_mo[1][neg_fit,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:]

    ### Compute climatology    
    climofitpos = np.nanmean(tas_mofitpos,axis=0)
    climohitpos = np.nanmean(tas_mohitpos,axis=0)
    climofictpos = np.nanmean(tas_mofictpos,axis=0)
    climofitneg = np.nanmean(tas_mofitneg,axis=0)
    climohitneg = np.nanmean(tas_mohitneg,axis=0)
    climofictneg = np.nanmean(tas_mofictneg,axis=0)
    climo = [climohitpos,climohitneg,climohitneg,
             climohitpos,climohitneg,climohitneg]
    
    ### Compute comparisons for months - taken ensemble average
    fithitpos = np.nanmean(tas_mofitpos - tas_mohitpos,axis=0)
    fithitneg = np.nanmean(tas_mofitneg - tas_mohitneg,axis=0)
    diffposfit = tas_mofitpos - tas_mohitpos
    diffnegfit = tas_mofitneg - tas_mohitneg
    diffallfit = np.nanmean(diffnegfit,axis=0) - np.nanmean(diffposfit,axis=0)
    
    ficthitpos = np.nanmean(tas_mofictpos - tas_mohitpos,axis=0)
    ficthitneg = np.nanmean(tas_mofictneg - tas_mohitneg,axis=0)
    diffposfict = tas_mofictpos - tas_mohitpos
    diffnegfict = tas_mofictneg - tas_mohitneg
    diffallfict = np.nanmean(diffnegfict,axis=0) - np.nanmean(diffposfict,axis=0)
    diffruns_mo = [fithitpos,fithitneg,diffallfit,ficthitpos,ficthitneg,diffallfict]
    
    ### Calculate significance for FM
    stat_FITHITpos,pvalue_FITHITpos = UT.calc_indttest(tas_mo[1][pos_fit,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FITHITneg,pvalue_FITHITneg = UT.calc_indttest(tas_mo[1][neg_fit,:,:],
                                               tas_mo[0][neg_hit,:,:])
    stat_fitdiff,pvalue_fitdiff = calc_indttest_90(diffnegfit,diffposfit)
    
    stat_FICTHITpos,pvalue_FICTHITpos = UT.calc_indttest(tas_mo[2][pos_fict,:,:],
                                                       tas_mo[0][pos_hit,:,:])
    stat_FICTHITneg,pvalue_FICTHITneg = UT.calc_indttest(tas_mo[2][neg_fict,:,:],
                                               tas_mo[0][neg_hit,:,:])
    stat_fictdiff,pvalue_fictdiff = calc_indttest_90(diffnegfict,diffposfict)

    pruns_mo = [pvalue_FITHITpos,pvalue_FITHITneg,pvalue_fitdiff,
                pvalue_FICTHITpos,pvalue_FICTHITneg,pvalue_fictdiff]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot variable data for QBO composites
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,5)
    elif varnames[v] == 'SLP':
        limit = np.arange(-6,6.1,0.5)
        barlim = np.arange(-6,7,3)
    elif varnames[v] == 'Z500':
        limit = np.arange(-60,60.1,5)
        barlim = np.arange(-60,61,30) 
    elif varnames[v] == 'Z30':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,50) 
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        limit = np.arange(-5,5.1,0.25)
        barlim = np.arange(-5,6,5)
    elif varnames[v] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    elif varnames[v] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,1) 
    elif varnames[v] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,30)
    elif varnames[v] == 'EGR':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.2)
    elif varnames[v] == 'WAFZ850':
        limit = np.arange(-0.1,0.101,0.001)
        barlim = np.arange(-0.1,0.11,0.1)
    elif varnames[v] == 'WAFZ150':
        limit = np.arange(-0.01,0.0101,0.0001)
        barlim = np.arange(-0.01,0.011,0.01)
    elif varnames[v] == 'WAFY850':
        limit = np.arange(-5,5.1,0.1)
        barlim = np.arange(-5,6,5)
    elif varnames[v] == 'WAFY150':
        limit = np.arange(-2,2.01,0.05)
        barlim = np.arange(-2,3,2)
    
    fig = plt.figure()
    for i in range(len(diffruns_mo)):
        var = diffruns_mo[i]
        pvar = pruns_mo[i]
        
        ax1 = plt.subplot(2,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
        climoq,lons_cyclic = addcyclic(climo[i], lon)
        climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        if any([i==2,i==5]):
            cs = m.contourf(x,y,var*pvar,limit,extend='both')
        else:
            cs = m.contourf(x,y,var,limit,extend='both')
        if any([i==0,i==1,i==3,i==4]):
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)
        if varnames[v] == 'Z30': # the interval is 250 m 
            cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                            colors='k',linewidths=1.5,zorder=10)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'T2M':
            cmap = cmocean.cm.balance         
            cs.set_cmap(cmap)   
        elif varnames[v] == 'SLP':
            cmap = cmocean.cm.balance          
            cs.set_cmap(cmap)   
        elif varnames[v] == 'Z500':
            cmap = cmocean.cm.balance           
            cs.set_cmap(cmap)  
        elif varnames[v] == 'Z30':
            cmap = cmocean.cm.balance  
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
            cmap = cmocean.cm.balance          
            cs.set_cmap(cmap)  
        elif varnames[v] == 'SWE':
            cmap = cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
        elif varnames[v] == 'P':
            cmap = ncm.cmap('precip4_diff_19lev')            
            cs.set_cmap(cmap) 
        elif varnames[v] == 'THICK':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap) 
        elif varnames[v] == 'EGR':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
        elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
        elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
            cmap = cmocean.cm.curl
            cs.set_cmap(cmap)
                    
        ### Add experiment text to subplot
        if i < 3:
            qbophaseq = [r'QBO-W',r'QBO-E',r'E--W']
            ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')
        if i == 0 or i == 3:
            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(-0.1,0.5),
                         textcoords='axes fraction',color='k',
                         fontsize=18,rotation=90,ha='center',va='center')
        ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
                               xytext=(0.92,0.9),xycoords='axes fraction',
                               color='dimgrey',fontsize=7)
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=11,color='k')  
    elif varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='k')  
    elif varnames[v] == 'Z30':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='k')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=11,color='k')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U500':
        cbar.set_label(r'\textbf{m/s}',fontsize=9,color='k')  
    elif varnames[v] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=11,color='k')
    elif varnames[v] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=11,color='k') 
    elif varnames[v] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=11,color='k') 
    elif varnames[v] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=11,color='k')
    elif varnames[v] == 'WAFZ850' or varnames[v] == 'WAFZ150':
        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='k')
    elif varnames[v] == 'WAFY850' or varnames[v] == 'WAFY150':
        cbar.set_label(r'\textbf{m$^{2}$/s$^{2}$}',fontsize=11,color='k')

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01,labelsize=7)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=-0.15,hspace=0.001,bottom=0.15)
    
    plt.savefig(directoryfigure + 'QBOExperiments_%s_%s_Mask.png' % (period,
                                                                  varnames[v]),
                                                                  dpi=900)

print('Completed: Script done!')

