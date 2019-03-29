"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
for DAILY data to calculate a NAO index

Notes
-----
    Author : Zachary Labe
    Date   : 8 August 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
import datetime
import read_DailyOutput_AllMembers as DO
import calc_Utilities as UT
import cmocean
from eofs.standard import Eof
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
print('\n' '----Plotting Daily NAO for SeaIceQBO - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

#### Add parameters
#varnames = ['Z500']
#runnames = [r'HIT',r'FICT']
#qbophase = ['pos','non','neg']
#experiments = [r'\textbf{FICT--HIT}']
#
#### Call functions for variable profile data for polar cap
#lat,lon,time,lev,varhit = DO.readMeanExperiAll('%s' % varnames[0],
#                                            'HIT','surface')
#lat,lon,time,lev,varfict = DO.readMeanExperiAll('%s' % varnames[0],
#                                            'FICT','surface')
#
#### Create 2d array of latitude and longitude
#lon2,lat2 = np.meshgrid(lon,lat)
#
### Read in QBO phases 
#filenamehitp = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitno = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#filenamehitn = directorydata + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#filenamehitp2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[0]
#filenamehitno2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[1]
#filenamehitn2 = directorydata2 + 'HIT/monthly/QBO_%s_HIT.txt' % qbophase[2]
#pos_hit = np.append(np.genfromtxt(filenamehitp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitp2,unpack=True,usecols=[0],dtype='int')+100)
#non_hit = np.append(np.genfromtxt(filenamehitno,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitno2,unpack=True,usecols=[0],dtype='int')+100)
#neg_hit = np.append(np.genfromtxt(filenamehitn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamehitn2,unpack=True,usecols=[0],dtype='int')+100)    
#
#filenamefictp = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictno = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#filenamefictn = directorydata + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#filenamefictp2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[0]
#filenamefictno2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[1]
#filenamefictn2 = directorydata2 + 'FICT/monthly/QBO_%s_FICT.txt' % qbophase[2]
#pos_fict = np.append(np.genfromtxt(filenamefictp,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictp2,unpack=True,usecols=[0],dtype='int')+100)
#non_fict = np.append(np.genfromtxt(filenamefictno,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictno2,unpack=True,usecols=[0],dtype='int')+100)
#neg_fict = np.append(np.genfromtxt(filenamefictn,unpack=True,usecols=[0],dtype='int'),
#                    np.genfromtxt(filenamefictn2,unpack=True,usecols=[0],dtype='int')+100)    
#### Concatonate runs
#var_mo = [varhit,varfict]
#
#### Composite by QBO phase    
#var_mohitpos = var_mo[0][pos_hit,:]
#var_mofictpos = var_mo[1][pos_fict,:]
#
#var_mohitnon = var_mo[0][non_hit,:]
#var_mofictnon = var_mo[1][non_fict,:]
#
#var_mohitneg = var_mo[0][neg_hit,:]
#var_mofictneg = var_mo[1][neg_fict,:]
#                                                                    
#### Calculate ensemble mean
#z5_diffposq = var_mofictpos - var_mohitpos
#z5_diffnonq = var_mofictnon - var_mohitnon
#z5_diffnegq = var_mofictneg - var_mohitneg
#
#### Calculate for NDJ
#z5_diffpos = z5_diffposq[:,60:150,:,:]     
#z5_diffnon = z5_diffnonq[:,60:150,:,:]   
#z5_diffneg = z5_diffnegq[:,60:150,:,:]   
#
##### Slice over (20-90N) and (90W-40E)
#latq = np.where((lat>=20) & (lat<=90))[0]
#latnao = lat[latq]
#
#lonnew = np.mod(lon, 360.0) - 180.0
#lonq = np.where((lonnew>=-90) & (lonnew<=40))[0]
#lonnao = lonnew[lonq]
#
#z5_diffposq = z5_diffpos[:,:,latq,:]
#z5_diffnaopos = z5_diffposq[:,:,:,lonq]
#
#z5_diffnonq = z5_diffnon[:,:,latq,:]
#z5_diffnaonon = z5_diffnonq[:,:,:,lonq]
#
#z5_diffnegq = z5_diffneg[:,:,latq,:]
#z5_diffnaoneg = z5_diffnegq[:,:,:,lonq]
#
#### Calculate climatology
#z5n_hpos = np.nanmean(var_mohitpos[:,60:150,latq,:],axis=0)
#z5nao_hpos = z5n_hpos[:,:,lonq]
#
#z5n_hnon = np.nanmean(var_mohitnon[:,60:150,latq,:],axis=0)
#z5nao_hnon = z5n_hnon[:,:,lonq]
#
#z5n_hneg = np.nanmean(var_mohitneg[:,60:150,latq,:],axis=0)
#z5nao_hneg = z5n_hneg[:,:,lonq]
#
#### Calculate NAO
## Create an EOF solver to do the EOF analysis. Square-root of cosine of
## latitude weights are applied before the computation of EOFs.
#coslat = np.cos(np.deg2rad(latnao)).clip(0., 1.)
#wgts = np.sqrt(coslat)[..., np.newaxis]
#solverpos = Eof(z5nao_hpos, weights=wgts)
#solvernon = Eof(z5nao_hnon, weights=wgts)
#solverneg = Eof(z5nao_hneg, weights=wgts)
#
## Retrieve the leading EOF, expressed as the covariance between the leading PC
## time series and the input SLP anomalies at each grid point.
#eof1pos = solverpos.eofsAsCovariance(neofs=1).squeeze()*-1
#pc1pos = solverpos.pcs(npcs=1, pcscaling=1).squeeze()
#
#eof1non = solvernon.eofsAsCovariance(neofs=1).squeeze()*-1
#pc1non = solvernon.pcs(npcs=1, pcscaling=1).squeeze()
#
#eof1neg= solverneg.eofsAsCovariance(neofs=1).squeeze()*-1
#pc1neg = solverneg.pcs(npcs=1, pcscaling=1).squeeze()
#
#### Calculate NAO index
#def NAOIndex(anomz5,eofpattern,members):
#    """
#    Calculate NAO index by regressing Z500 onto the EOF1 pattern
#    """
#    print('\n>>> Using NAO Index function!')       
#    
#    if members == True:
#        nao = np.empty((anomz5.shape[0],anomz5.shape[1]))
#        for i in range(anomz5.shape[0]):
#            print('Regressing ensemble ---> %s!' % (i+1))
#            for j in range(anomz5.shape[1]):
#                varx = np.ravel(anomz5[i,j,:,:])
#                vary = np.ravel(eofpattern[:,:])
#                mask = np.isfinite(varx) & np.isfinite(vary)     
#                
#                nao[i,j],intercept,r,p_value,std_err = sts.stats.linregress(
#                                                                      varx[mask],
#                                                                      vary[mask]) 
#    elif members == False:   
#        nao = np.empty((anomz5.shape[0]))
#        for i in range(anomz5.shape[0]):
#            varx = np.ravel(anomz5[i,:,:])
#            vary = np.ravel(eofpattern[:,:])
#            mask = np.isfinite(varx) & np.isfinite(vary)     
#            
#            nao[i],intercept,r,p_value,std_err = sts.stats.linregress(
#                                                                  varx[mask],
#                                                                  vary[mask]) 
#        print('Completed: Regressed ensemble mean!')
#    else:
#        ValueError('Please select [True] or [False] for averageing!')
#        
#    print('*Completed: finished with NAO function!')
#    return nao
#            
#### Calculate NAO index
naoindexpos = NAOIndex(np.nanmean(z5_diffnaopos,axis=0),eof1pos,False)
pc1pos = (naoindexpos-np.mean(naoindexpos))/np.std(naoindexpos)

naoindexnon = NAOIndex(np.nanmean(z5_diffnaonon,axis=0),eof1non,False)
pc1non = (naoindexnon-np.mean(naoindexnon))/np.std(naoindexnon)

naoindexneg = NAOIndex(np.nanmean(z5_diffnaoneg,axis=0),eof1neg,False)
pc1neg = (naoindexneg-np.mean(naoindexneg))/np.std(naoindexneg)
    
##### Plot figure
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#
#fig = plt.figure()
#ax = plt.subplot(111)
#    
#varf = eof1neg[:,:]    
#
#m = Basemap(projection='ortho',lon_0=-20,lat_0=60,resolution='l',
#            area_thresh=10000.)
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='dimgray',linewidth=0.8)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
#
## Make the plot continuous
#barlim = np.arange(-8,9,4)
#values = np.arange(-8,8.1,0.25)
#
#lon2,lat2 = np.meshgrid(lonnao,latnao)
#
#cs = m.contourf(lon2,lat2,varf,45,
#                extend='both',latlon=True)
#cs1 = m.contour(lon2,lat2,varf,
#                linewidths=0.1,colors='darkgrey',
#                linestyles='-',latlon=True)
#        
#cmap = cmocean.cm.balance         
#cs.set_cmap(cmap) 
#
#cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)                    
##cbar.set_ticks(barlim)
##cbar.set_ticklabels(list(map(str,barlim)))  
#cbar.ax.tick_params(labelsize=8)   
#
#plt.savefig(directoryfigure + 'testeof1_neg.png',dpi=300)
#
############################################################################
############################################################################
############################################################################
#### Plot NAO index
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
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('dimgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')


plt.axhline(0,linestyle=':',linewidth=2,color='dimgrey',dashes=(1,0.3))
plt.plot(pc1pos,linewidth=3.5,color='crimson',alpha=1,linestyle='-',
         label=r'\textbf{W-QBO}')
plt.plot(pc1neg,linewidth=3.5,color='deepskyblue',alpha=1,linestyle='-',
         label=r'\textbf{E-QBO}')

plt.legend(shadow=False,fontsize=9,loc='lower center',
           fancybox=True,frameon=False,ncol=5)
plt.ylabel(r'\textbf{NAO Index (Z500)}',color='dimgrey',fontsize=13)

plt.yticks(np.arange(-5,6,1),list(map(str,np.arange(-5,6,1))),fontsize=9)
plt.ylim([-3,3])

xlabels = [r'Nov',r'Dec',r'Jan',r'Feb'] 
plt.xticks(np.arange(0,121,30),xlabels,fontsize=9)
plt.xlim([0,60])

plt.savefig(directoryfigure + 'NAOIndex_FICT_QBO.png',dpi=300)

##############################################################################
##############################################################################
##############################################################################
###### Plot NAO histogram
#naoindexpos = NAOIndex(z5_diffnaopos,eof1pos,True)
#pc1pos = (naoindexpos-np.mean(naoindexpos))/np.std(naoindexpos)
#
#naoindexnon = NAOIndex(z5_diffnaonon,eof1non,True)
#pc1non = (naoindexnon-np.mean(naoindexnon))/np.std(naoindexnon)
#
#naoindexneg = NAOIndex(z5_diffnaoneg,eof1neg,True)
#pc1neg = (naoindexneg-np.mean(naoindexneg))/np.std(naoindexneg)
#
#num_bins = np.arange(-3,3.1,0.1)
#
#fig = plt.figure()
#ax = plt.subplot(111) 
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('dimgrey')
#ax.spines['bottom'].set_color('dimgrey')
#ax.spines['left'].set_linewidth(2)
#ax.spines['bottom'].set_linewidth(2)
#ax.tick_params('both',length=4,width=2,which='major',color='dimgrey')
##ax.tick_params(axis='x',which='both',bottom=False)
#
#ne,binse,patchese = plt.hist(pc1neg[:,30:60].ravel(),num_bins,density=1,alpha=0.7,
#                          facecolor='deepskyblue',edgecolor='deepskyblue',
#                          range=(-3,3),label=r'\textbf{QBO-E}')
#nw,binsw,patchesw = plt.hist(pc1pos[:,30:60].ravel(),num_bins,density=1,alpha=0.5,
#                          facecolor='crimson',edgecolor='crimson',
#                          range=(-3,3),label=r'\textbf{QBO-W}')
#
#plt.yticks(np.arange(0,1.1,0.1),list(map(str,np.arange(0,1.1,0.1))),
#           fontsize=10)
#plt.xticks(np.arange(-3,4,1),list(map(str,np.arange(-3,4,1))),
#           fontsize=10) 
#plt.xlim([-3,3])
#plt.ylim([0,0.5])
#
#plt.legend(shadow=False,fontsize=7,loc='upper left',
#           fancybox=True,frameon=False,ncol=1,bbox_to_anchor=(-0.02, 1.015),
#           labelspacing=0.2,columnspacing=1,handletextpad=0.4)
#
#plt.ylabel(r'\textbf{Density}',color='dimgrey',fontsize=12)  
#plt.xlabel(r'\textbf{NAO Index (Z500)}',color='dimgrey',fontsize=12)
#
#plt.savefig(directoryfigure + 'NAOindex_histogram_FICT.png',dpi=300)