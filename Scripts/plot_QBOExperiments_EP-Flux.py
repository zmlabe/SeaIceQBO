"""
Plot EP Flux comparisons between SIT and SIC modeling experiments using 
WACCM4. Subplot includes FIT, HIT, FICT. Profiles are 
organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 6 June 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput_AllMembers as MO
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_D_2/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO profile comparisons - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

runnames = [r'HIT',r'FIT',r'FICT']
experiments = [r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',r'\textbf{FIT--HIT}',
               r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}',r'\textbf{FICT--HIT}']
qbophase = ['pos','non','neg']
period = 'D'

### Call function for vertical temperature data
lat,lon,time,lev,epy_h = MO.readExperiAll('EPY','HIT','profile')
lat,lon,time,lev,epz_h = MO.readExperiAll('EPZ','HIT','profile')
lat,lon,time,lev,div_h = MO.readExperiAll('DEPF','HIT','profile')

lat,lon,time,lev,epy_f = MO.readExperiAll('EPY','FIT','profile')
lat,lon,time,lev,epz_f = MO.readExperiAll('EPZ','FIT','profile')
lat,lon,time,lev,div_f = MO.readExperiAll('DEPF','FIT','profile')

lat,lon,time,lev,epy_fict = MO.readExperiAll('EPY','FICT','profile')
lat,lon,time,lev,epz_fict = MO.readExperiAll('EPZ','FICT','profile')
lat,lon,time,lev,div_fict = MO.readExperiAll('DEPF','FICT','profile')

### Separate per month
epy_moh = epy_h[:,-1,:,:] 
epz_moh = epz_h[:,-1,:,:] 
div_moh = div_h[:,-1,:,:]

epy_mof = epy_f[:,-1,:,:] 
epz_mof = epz_f[:,-1,:,:] 
div_mof = div_f[:,-1,:,:]

epy_mofict = epy_fict[:,-1,:,:] 
epz_mofict = epz_fict[:,-1,:,:] 
div_mofict = div_fict[:,-1,:,:]

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

### Calculate differences (per day so divide by 30)
# Positive
diff_epyp = np.nanmean(epy_mof[pos_fit] - epy_moh[pos_hit],axis=0)/30.
diff_epzp = np.nanmean(epz_mof[pos_fit] - epz_moh[pos_hit],axis=0)/30.
diff_divp = np.nanmean(div_mof[pos_fit] - div_moh[pos_hit],axis=0)/30.

diff_epyp_fict = np.nanmean(epy_mofict[pos_fict] - epy_moh[pos_hit],axis=0)/30.
diff_epzp_fict = np.nanmean(epz_mofict[pos_fict] - epz_moh[pos_hit],axis=0)/30.
diff_divp_fict = np.nanmean(div_mofict[pos_fict] - div_moh[pos_hit],axis=0)/30.

# Difference (QBO-E minus QBO-W)
diff_epyno = np.nanmean(epy_mof[neg_fit] - epy_mof[pos_fit][:-4],axis=0)/30.
diff_epzno = np.nanmean(epz_mof[neg_fit] - epz_mof[pos_fit][:-4],axis=0)/30.
diff_divno = np.nanmean(div_mof[neg_fit] - div_mof[pos_fit][:-4],axis=0)/30.

diff_epyno_fict = np.nanmean(epy_mofict[neg_fict] - epy_mofict[pos_fict][:-4],axis=0)/30.
diff_epzno_fict = np.nanmean(epz_mofict[neg_fict] - epz_mofict[pos_fict][:-4],axis=0)/30.
diff_divno_fict = np.nanmean(div_mofict[neg_fict] - div_mofict[pos_fict][:-4],axis=0)/30.
    
# Negative
diff_epyn = np.nanmean(epy_mof[neg_fit] - epy_moh[neg_hit],axis=0)/30.
diff_epzn = np.nanmean(epz_mof[neg_fit] - epz_moh[neg_hit],axis=0)/30.
diff_divn = np.nanmean(div_mof[neg_fit] - div_moh[neg_hit],axis=0)/30.

diff_epyn_fict = np.nanmean(epy_mofict[neg_fit] - epy_moh[neg_hit],axis=0)/30.
diff_epzn_fict = np.nanmean(epz_mofict[neg_fit] - epz_moh[neg_hit],axis=0)/30.
diff_divn_fict = np.nanmean(div_mofict[neg_fit] - div_moh[neg_hit],axis=0)/30.

diff_div = [diff_divn,diff_divp,diff_divno,diff_divn_fict,diff_divp_fict,diff_divno_fict]
diff_epy = [diff_epyn,diff_epyp,diff_epyno,diff_epyn_fict,diff_epyp_fict,diff_epyno_fict]
diff_epz = [diff_epzn,diff_epzp,diff_epzno,diff_epzn_fict,diff_epzp_fict,diff_epzno_fict]

###########################################################################
###########################################################################
###########################################################################
#### Plot EP Flux
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-1,1.1,0.05)
barlim = np.arange(-1,2,1)
    
zscale = np.array([700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure()
for i in range(6):
    ax1 = plt.subplot(2,3,i+1)
    
    epy = diff_epy[i] * np.sqrt(1000/levq) 
    epz = diff_epz[i] * np.sqrt(1000/levq) * 100 
    div = diff_div[i]

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
            
    cs = plt.contourf(lat,lev,div,limit,extend='both')
#    cs1 = plt.quiver(lat[::4],lev,diff_epy[i][:,::4],diff_epz[i][:,::4],
#                     pivot='mid',color='k',units='width',
#                     scale=0.8e7,width=0.007)
    cs1 = plt.quiver(lat[::5],lev,epy[:,::5],epz[:,::5],
                     pivot='mid',color='k',width=0.011,
                     headwidth=5,headlength=6,headaxislength=4)
#    if i == 5:
#        plt.quiverkey(cs1,0.34,-0.3,0.8e7,r'\textbf{0.8$\times$10$^{6}$}',
#                      coordinates='axes',labelpos='E')

#        plt.contourf(latq,levq,pruns[i],colors='None',hatches=['////'],
#                     linewidth=5)   
    
    plt.gca().invert_yaxis()
    plt.yscale('log',nonposy='clip')
    
    plt.xticks(np.arange(0,96,30),map(str,np.arange(0,91,30)),fontsize=6)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    plt.minorticks_off()
    
    plt.xlim([0,90])
    plt.ylim([700,10])
        
    cmap = cmocean.cm.balance          
    cs.set_cmap(cmap) 

    ### Add experiment text to subplot
    if i < 3:
        qbophaseq = [r'QBO-E',r'QBO-W',r'Difference']
        ax1.text(0.5,1.1,r'\textbf{%s}' % qbophaseq[i],
                 ha='center',va='center',color='dimgray',fontsize=13,
                 transform=ax1.transAxes)
    if i == 0 or i == 3:
        ax1.text(-0.35,0.5,r'%s' % experiments[i],color='k',
                     fontsize=20,rotation=90,ha='center',va='center',
                     transform=ax1.transAxes)
    if i == 4:
        plt.xlabel(r'\textbf{Latitude ($^\circ$N)}',color='k',fontsize=8,
                             labelpad=1)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)

cbar.set_label(r'\textbf{m/s/day}',fontsize=11,color='dimgray',labelpad=1)
    
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')     

plt.subplots_adjust(wspace=0.24)
plt.subplots_adjust(bottom=0.21)

plt.savefig(directoryfigure + 'epflux_QBO_all.png',dpi=300)
print('Completed: Script done!')

