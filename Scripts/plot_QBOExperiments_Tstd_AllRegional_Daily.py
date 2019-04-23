"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data for temperature standard deviation at 1000 hPa (10%)

Notes
-----
    Author : Zachary Labe
    Date   : 22 May 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import read_DailyOutput_AllRegional as DOR
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
print('\n' '----Plotting Daily T1000 Sigma for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
MASK = True
varnames = ['T1000']
runnames = [r'CIT',r'FSUB',r'FPOL']
qbophase = ['pos','non','neg']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',r'\textbf{FSUB-CIT}',
               r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}']

### Call functions for variable profile data for polar cap
for v in range(len(varnames)):
    lat,lon,time,lev,varcit = DO.readMeanExperiAll('%s' % varnames[v],
                                                'CIT','surface')
    lat,lon,time,lev,varfsub = DOR.readMeanExperiAllRegional('%s' % varnames[v],
                                                'FSUB','surface')
    lat,lon,time,lev,varfpol = DOR.readMeanExperiAllRegional('%s' % varnames[v],
                                                'FPOL','surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Read in QBO phases 
    filenamecitp = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    filenamecitp2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
    filenamecitno2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
    filenamecitn2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
    pos_cit = np.append(np.genfromtxt(filenamecitp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitp2,unpack=True,usecols=[0],dtype='int')+100)
    non_cit = np.append(np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_cit = np.append(np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamecitn2,unpack=True,usecols=[0],dtype='int')+100)    

    filenamefsubp = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    filenamefsubp2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
    filenamefsubno2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
    filenamefsubn2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
    pos_fsub = np.append(np.genfromtxt(filenamefsubp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fsub = np.append(np.genfromtxt(filenamefsubno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fsub = np.append(np.genfromtxt(filenamefsubn,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefsubn2,unpack=True,usecols=[0],dtype='int')+100)
    
    filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    filenamefpolp2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
    filenamefpolno2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
    filenamefpoln2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
    pos_fpol = np.append(np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolp2,unpack=True,usecols=[0],dtype='int')+100)
    non_fpol = np.append(np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpolno2,unpack=True,usecols=[0],dtype='int')+100)
    neg_fpol = np.append(np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int'),
                        np.genfromtxt(filenamefpoln2,unpack=True,usecols=[0],dtype='int')+100)
    
    ### Concatonate runs
    var_mo = [varcit,varfsub,varfpol]
    
    ### Composite by QBO phase    
    var_mofsubpos = var_mo[1][pos_fsub,:,:]
    var_mocitpos = var_mo[0][pos_cit,:,:]
    var_mofpolpos = var_mo[2][pos_fpol,:,:]
    
    var_mofsubnon = var_mo[1][non_fsub,:,:]
    var_mocitnon = var_mo[0][non_cit,:,:]
    var_mofpolnon = var_mo[2][non_fpol,:,:]
    
    var_mofsubneg = var_mo[1][neg_fsub,:,:]
    var_mocitneg = var_mo[0][neg_cit,:,:]
    var_mofpolneg = var_mo[2][neg_fpol,:,:]
    
    ### Calculate over DJF (90-180)
    timeq = np.arange(120,150)
    monthqq = 'Jan'
    var_wfsubpos = var_mofsubpos[:,timeq,:,:]
    var_wcitpos = var_mocitpos[:,timeq,:,:]
    var_wfpolpos = var_mofpolpos[:,timeq,:,:]
    
    var_wfsubnon = var_mofsubnon[:,timeq,:,:]
    var_wcitnon = var_mocitnon[:,timeq,:,:]
    var_wfpolnon = var_mofpolpos[:,timeq,:,:]
    
    var_wfsubneg = var_mofsubneg[:,timeq,:,:]
    var_wcitneg = var_mocitneg[:,timeq,:,:]
    var_wfpolneg = var_mofpolneg[:,timeq,:,:]
    
    ### Compute standard deviation temperature at each grid point
    ten_wfsubpos = np.empty((var_wfsubpos.shape[0],lat.shape[0],lon.shape[0]))
    ten_wcitpos = np.empty((var_wcitpos.shape[0],lat.shape[0],lon.shape[0]))
    ten_wfpolpos = np.empty((var_wfpolpos.shape[0],lat.shape[0],lon.shape[0]))
    for ens in range(var_wfsubpos.shape[0]):
        for i in range(lat.shape[0]):
            for j in range(lon.shape[0]):
                ten_wfsubpos[ens,i,j] = np.nanstd(var_wfsubpos[ens,:,i,j],axis=0)
                ten_wcitpos[ens,i,j] = np.nanstd(var_wcitpos[ens,:,i,j],axis=0)
                ten_wfpolpos[ens,i,j] = np.nanstd(var_wfpolpos[ens,:,i,j],axis=0)
        print('Completed: Sigma of QBO-W for %s member!' % ens)
                
    ten_wfsubnon = np.empty((var_wfsubnon.shape[0],lat.shape[0],lon.shape[0]))
    ten_wcitnon = np.empty((var_wcitnon.shape[0],lat.shape[0],lon.shape[0]))
    ten_wfpolnon = np.empty((var_wfpolnon.shape[0],lat.shape[0],lon.shape[0]))
    for ens in range(var_wfsubnon.shape[0]):
        for i in range(lat.shape[0]):
            for j in range(lon.shape[0]):
                ten_wfsubnon[ens,i,j] = np.nanstd(var_wfsubnon[ens,:,i,j],axis=0)
                ten_wcitnon[ens,i,j] = np.nanstd(var_wcitnon[ens,:,i,j],axis=0)
                ten_wfpolnon[ens,i,j] = np.nanstd(var_wfpolnon[ens,:,i,j],axis=0)
        print('Completed: Sigma of QBO-N for %s member!' % ens)
                
    ten_wfsubneg = np.empty((var_wfsubneg.shape[0],lat.shape[0],lon.shape[0]))
    ten_wcitneg = np.empty((var_wcitneg.shape[0],lat.shape[0],lon.shape[0]))
    ten_wfpolneg = np.empty((var_wfpolneg.shape[0],lat.shape[0],lon.shape[0]))
    for ens in range(var_wfsubneg.shape[0]):
        for i in range(lat.shape[0]):
            for j in range(lon.shape[0]):
                ten_wfsubneg[ens,i,j] = np.nanstd(var_wfsubneg[ens,:,i,j],axis=0)
                ten_wcitneg[ens,i,j] = np.nanstd(var_wcitneg[ens,:,i,j],axis=0)
                ten_wfpolneg[ens,i,j] = np.nanstd(var_wfpolneg[ens,:,i,j],axis=0)
        print('Completed: Sigma of QBO-E for %s member!' % ens)
    
    ### Compute comparisons for days - taken ensemble average
    fsubcitpos = np.nanmean(ten_wfsubpos - ten_wcitpos,axis=0)
    fsubcitnon = np.nanmean(ten_wfsubnon - ten_wcitnon,axis=0)
    fsubcitneg = np.nanmean(ten_wfsubneg - ten_wcitneg,axis=0)
    
    fpolcitpos = np.nanmean(ten_wfpolpos - ten_wcitpos,axis=0)
    fpolcitnon = np.nanmean(ten_wfpolnon - ten_wcitnon,axis=0)
    fpolcitneg = np.nanmean(ten_wfpolneg - ten_wcitneg,axis=0)
    diffruns = [fsubcitpos,fsubcitnon,fsubcitneg,
                   fpolcitpos,fpolcitnon,fpolcitneg]
    
    ### Calculate significance for DJF
    stat_fsubcitpos,pvalue_fsubcitpos = UT.calc_indttest(ten_wfsubpos,ten_wcitpos)
    stat_fsubcitnon,pvalue_fsubcitnon = UT.calc_indttest(ten_wfsubnon,ten_wcitnon)
    stat_fsubcitneg,pvalue_fsubcitneg = UT.calc_indttest(ten_wfsubneg,ten_wcitneg)
    
    stat_fpolcitpos,pvalue_fpolcitpos = UT.calc_indttest(ten_wfpolpos,ten_wcitpos)
    stat_fpolcitnon,pvalue_fpolcitnon = UT.calc_indttest(ten_wfpolnon,ten_wcitnon)
    stat_fpolcitneg,pvalue_fpolcitneg = UT.calc_indttest(ten_wfpolneg,ten_wcitneg)

    pruns = [pvalue_fsubcitpos,pvalue_fsubcitnon,pvalue_fsubcitneg,
                pvalue_fpolcitpos,pvalue_fpolcitnon,pvalue_fpolcitneg]
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Plot variable data for QBO composites
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[v] == 'T1000':
        limit = np.arange(-2,2.1,0.1)
        barlim = np.arange(-2,3,2)
        
    fig = plt.figure()
    for i in range(len(diffruns)):
        var = diffruns[i]
        pvar = pruns[i]
        
        ax1 = plt.subplot(2,3,i+1)
        m = Basemap(projection='ortho',lon_0=0,lat_0=89,resolution='l',
                    area_thresh=10000.)
        
        var, lons_cyclic = addcyclic(var, lon)
        var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
        lon2d, lat2d = np.meshgrid(lons_cyclic, lat)
        x, y = m(lon2d, lat2d)
        
        pvar,lons_cyclic = addcyclic(pvar, lon)
        pvar,lons_cyclic = shiftgrid(180.,pvar,lons_cyclic,start=False)
                  
        m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
        
        cs = m.contourf(x,y,var,limit,extend='both')
        cs1 = m.contourf(x,y,pvar,colors='None',hatches=['....'],
                         linewidths=0.4)

        m.drawcoastlines(color='dimgray',linewidth=0.8)
        
        if varnames[v] == 'T1000':
            cmap = ncm.cmap('NCV_blu_red')           
            cs.set_cmap(cmap)   
        
        m.drawlsmask(land_color=(0,0,0,0),ocean_color='gainsboro',lakes=False,
                     resolution='c',zorder=5)
                    
        ### Add experiment text to subplot
        if i < 3:
            qbophaseq = [r'QBO-W',r'QBO-N',r'QBO-E']
            ax1.annotate(r'\textbf{%s}' % qbophaseq[i],xy=(0,0),xytext=(0.5,1.08),
                         textcoords='axes fraction',color='dimgray',
                         fontsize=13,rotation=0,ha='center',va='center')
        if i == 0 or i == 3:
            ax1.annotate(r'%s' % experiments[i],xy=(0,0),xytext=(-0.1,0.5),
                         textcoords='axes fraction',color='k',
                         fontsize=20,rotation=90,ha='center',va='center')
                
    ###########################################################################
    cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=False)
    
    if varnames[v] == 'T1000':
        cbar.set_label(r'\textbf{std. dev. [$^\circ$C]}',fontsize=11,color='dimgray')  

    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim)))
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    
    plt.subplots_adjust(wspace=0.01)
    plt.subplots_adjust(hspace=0.01)
    plt.subplots_adjust(bottom=0.15)
    
    plt.savefig(directoryfigure + '/QBO_std_2/QBOExperiments_%s_%sSTD_Regional.png' % (
                                                                  monthqq,
                                                                  varnames[v]),
                                                                  dpi=300)

print('Completed: Script done!')

