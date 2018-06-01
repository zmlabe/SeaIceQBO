"""
Plot difference in climatologies of the polar vortex for HIT and FICT by 
compositing according to QBO phase

Notes
-----
    Author : Zachary Labe
    Date   : 1 June 2018
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
print('\n' '----Plotting QBO climatology differences - %s----' % titletime)

### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Call arguments
varnames = ['Z30','U10','Z500','SLP','T2M','U300','U30']
runnames = [r'HIT',r'FICT']
qbophase = ['pos','non','neg']
experiments = ['HIT','FICT']
period = 'DJF'

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
    else:
        ValueError('Wrong period selected! (ON,DJ,FM)')
        
    ### Composite by QBO phase    
    tas_mohitpos = tas_mo[0][pos_hit,:,:]
    tas_mofitpos = tas_mo[1][pos_fit,:,:]
    tas_mofictpos = tas_mo[2][pos_fict,:,:]
    
    tas_mohitnon = tas_mo[0][non_hit,:,:]
    tas_mofitnon = tas_mo[1][non_fit,:,:]
    tas_mofictnon = tas_mo[2][non_fict,:,:]
    
    tas_mohitneg = tas_mo[0][neg_hit,:,:]
    tas_mofitneg = tas_mo[1][neg_fit,:,:]
    tas_mofictneg = tas_mo[2][neg_fict,:,:]

    ### Compute climatology    
    climohitdiff = np.nanmean(tas_mohitneg-tas_mohitpos[:-4],axis=0)
    climofictdiff = np.nanmean(tas_mofictneg-tas_mofictpos[:-4],axis=0)
    climo = [climohitdiff,climofictdiff]
    
    stat_hitdiff,pvalue_hitdiff = UT.calc_indttest(tas_mohitneg,
                                                   tas_mohitpos[:-4])
    stat_fictdiff,pvalue_fictdiff = UT.calc_indttest(tas_mofictneg,
                                                     tas_mofictpos[:-4])

    pruns = [pvalue_hitdiff,pvalue_fictdiff]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot variable data for QBO composites
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'Z30':
        limit = np.arange(-100,100.1,5)
        barlim = np.arange(-100,101,50) 
    elif varnames[v] == 'U10' or varnames[v] == 'U30':
        limit = np.arange(-8,8.1,0.1)
        barlim = np.arange(-8,9,4)
    elif varnames[v] == 'Z500':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25) 
    elif varnames[v] == 'SLP':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,4) 
    elif varnames[v] == 'T2M':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,4) 
    elif varnames[v] == 'U300':
        limit = np.arange(-4,4.1,0.1)
        barlim = np.arange(-4,5,4) 
        
    fig = plt.figure()
    for i in range(len(climo)):
        var = climo[i]
        pvar = pruns[i]
        
        ax1 = plt.subplot(1,2,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=90,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        if varnames[v] == 'U10' or varnames[v] == 'U30':
            cs = m.contourf(x,y,var,limit,extend='both')
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)

        elif varnames[v] == 'Z30' or varnames[v] == 'Z500':
            cs = m.contourf(x,y,var,limit,extend='both')
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)
        elif varnames[v] == 'SLP':
            cs = m.contourf(x,y,var,limit,extend='both')
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)
        elif varnames[v] == 'T2M':
            cs = m.contourf(x,y,var,limit,extend='both')
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)
        elif varnames[v] == 'U300':
            cs = m.contourf(x,y,var,limit,extend='both')
            cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                             linewidths=0.4)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
         
        if varnames[v] == 'Z30' or varnames[v] == 'Z500':
            cmap = cmocean.cm.balance
            cs.set_cmap(cmap)  
        elif varnames[v] == 'U10' or varnames[v] == 'U30':   
            cmap = ncm.cmap('NCV_blu_red') 
            cs.set_cmap(cmap)  
        elif varnames[v] == 'SLP':
            cmap = cmocean.cm.balance
            cs.set_cmap(cmap)
        elif varnames[v] == 'T2M':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        elif varnames[v] == 'U300':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)  
                    
        ### Add experiment text to subplot
        ax1.annotate(r'\textbf{%s}' % experiments[i],xy=(0,0),xytext=(0.5,1.08),
                     textcoords='axes fraction',color='dimgray',
                     fontsize=21,rotation=0,ha='center',va='center')

    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'Z30' or varnames[v] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=13,color='dimgray')  
    elif varnames[v] == 'U10' or varnames[v] == 'U300' or varnames[v] == 'U30':
        cbar.set_label(r'\textbf{m/s}',fontsize=13,color='dimgray')  
    elif varnames[v] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=13,color='darkgray')
    elif varnames[v] == 'T2M':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=13,color='dimgray')  

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    
    plt.savefig(directoryfigure + '/QBO_Climo_2/QBOExperiments_CLIMODIFF_%s_%s.png' % (
                                                                        period,
                                                                  varnames[v]),
                                                                  dpi=300)

print('Completed: Script done!')