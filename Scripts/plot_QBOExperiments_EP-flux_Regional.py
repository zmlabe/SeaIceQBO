"""
Plot EP Flux comparisons between SIT and SIC modeling experiments using 
regional WACCM4. Subplot includes fsub, cit, fpol. Profiles are 
organized by QBO phase (positive, neutral, negative)

Notes
-----
    Author : Zachary Labe
    Date   : 7 June 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_MonthlyOutput_AllMembers as MO
import read_MonthlyOutput_AllRegional as MOR
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/simu/'
directorydata2 = '/home/zlabe/green/simu/'
directoryfigure = '/home/zlabe/Desktop/QBO_D_2/Regional2/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceQBO/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting QBO profile regional comparisons - %s----' % titletime)

### Alott time series
year1 = 1900
year2 = 2000
years = np.arange(year1,year2+1,1)

runnames = [r'CIT',r'FSUB',r'FPOL']
experiments = [r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',r'\textbf{FSUB--CIT}',
               r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}',r'\textbf{FPOL--CIT}']
qbophase = ['pos','non','neg']
period = 'D'

### Call function for vertical temperature data
lat,lon,time,lev,epy_h = MO.readExperiAll('EPY','CIT','profile')
lat,lon,time,lev,epz_h = MO.readExperiAll('EPZ','CIT','profile')
lat,lon,time,lev,div_h = MO.readExperiAll('DEPF','CIT','profile')

lat,lon,time,lev,epy_f = MOR.readExperiAllRegional('EPY','FSUB','profile')
lat,lon,time,lev,epz_f = MOR.readExperiAllRegional('EPZ','FSUB','profile')
lat,lon,time,lev,div_f = MOR.readExperiAllRegional('DEPF','FSUB','profile')

lat,lon,time,lev,epy_fpol = MOR.readExperiAllRegional('EPY','FPOL','profile')
lat,lon,time,lev,epz_fpol = MOR.readExperiAllRegional('EPZ','FPOL','profile')
lat,lon,time,lev,div_fpol = MOR.readExperiAllRegional('DEPF','FPOL','profile')

### Separate per month
epy_moh = epy_h[:,-1,:,:] 
epz_moh = epz_h[:,-1,:,:] 
div_moh = div_h[:,-1,:,:]

epy_mof = epy_f[:,-1,:,:] 
epz_mof = epz_f[:,-1,:,:] 
div_mof = div_f[:,-1,:,:]

epy_mofpol = epy_fpol[:,-1,:,:] 
epz_mofpol = epz_fpol[:,-1,:,:] 
div_mofpol = div_fpol[:,-1,:,:]

### Read in QBO phases 
filenamecitp = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
filenamecitno = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
filenamecitn = directorydata + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
filenamecitp2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[0]
filenamecitno2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[1]
filenamecitn2 = directorydata2 + 'CIT/monthly/QBO_%s_CIT.txt' % qbophase[2]
pos_cit = np.append(np.genfromtxt(filenamecitp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamecitp2,unpack=True,usecols=[0],dtype='int')+101)
non_cit = np.append(np.genfromtxt(filenamecitno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamecitno2,unpack=True,usecols=[0],dtype='int')+101)
neg_cit = np.append(np.genfromtxt(filenamecitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamecitn2,unpack=True,usecols=[0],dtype='int')+101)    

filenamefsubp = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
filenamefsubno = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
filenamefsubn = directorydata2 + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
filenamefsubp2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[0]
filenamefsubno2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[1]
filenamefsubn2 = directorydata + 'FSUB/monthly/QBO_%s_FSUB.txt' % qbophase[2]
pos_fsub = np.append(np.genfromtxt(filenamefsubp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefsubp2,unpack=True,usecols=[0],dtype='int')+101)
non_fsub = np.append(np.genfromtxt(filenamefsubno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefsubno2,unpack=True,usecols=[0],dtype='int')+101)
neg_fsub = np.append(np.genfromtxt(filenamefsubn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefsubn2,unpack=True,usecols=[0],dtype='int')+101)

filenamefpolp = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
filenamefpolno = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
filenamefpoln = directorydata2 + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
filenamefpolp2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[0]
filenamefpolno2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[1]
filenamefpoln2 = directorydata + 'FPOL/monthly/QBO_%s_FPOL.txt' % qbophase[2]
pos_fpol = np.append(np.genfromtxt(filenamefpolp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefpolp2,unpack=True,usecols=[0],dtype='int')+101)
non_fpol = np.append(np.genfromtxt(filenamefpolno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefpolno2,unpack=True,usecols=[0],dtype='int')+101)
neg_fpol = np.append(np.genfromtxt(filenamefpoln,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefpoln2,unpack=True,usecols=[0],dtype='int')+101)

### Calculate differences (per day so divide by 30)
# Positive
diff_epyp = np.nanmean(epy_mof[pos_fsub] - epy_moh[pos_cit],axis=0)/30.
diff_epzp = np.nanmean(epz_mof[pos_fsub] - epz_moh[pos_cit],axis=0)/30.
diff_divp = np.nanmean(div_mof[pos_fsub] - div_moh[pos_cit],axis=0)/30.

diff_epyp_fpol = np.nanmean(epy_mofpol[pos_fpol] - epy_moh[pos_cit],axis=0)/30.
diff_epzp_fpol = np.nanmean(epz_mofpol[pos_fpol] - epz_moh[pos_cit],axis=0)/30.
diff_divp_fpol = np.nanmean(div_mofpol[pos_fpol] - div_moh[pos_cit],axis=0)/30.

# Difference (QBO-E minus QBO-W)
diff_epyno = np.nanmean(epy_mof[neg_fsub] - epy_mof[pos_cit][:-4],axis=0)/30.
diff_epzno = np.nanmean(epz_mof[neg_fsub] - epz_mof[pos_cit][:-4],axis=0)/30.
diff_divno = np.nanmean(div_mof[neg_fsub] - div_mof[pos_cit][:-4],axis=0)/30.

diff_epyno_fpol = np.nanmean(epy_mofpol[neg_fpol] - epy_mofpol[pos_fpol][:-4],axis=0)/30.
diff_epzno_fpol = np.nanmean(epz_mofpol[neg_fpol] - epz_mofpol[pos_fpol][:-4],axis=0)/30.
diff_divno_fpol = np.nanmean(div_mofpol[neg_fpol] - div_mofpol[pos_fpol][:-4],axis=0)/30.
    
# Negative
diff_epyn = np.nanmean(epy_mof[neg_fsub] - epy_moh[neg_cit],axis=0)/30.
diff_epzn = np.nanmean(epz_mof[neg_fsub] - epz_moh[neg_cit],axis=0)/30.
diff_divn = np.nanmean(div_mof[neg_fsub] - div_moh[neg_cit],axis=0)/30.

diff_epyn_fpol = np.nanmean(epy_mofpol[neg_fsub] - epy_moh[pos_fsub],axis=0)/30.
diff_epzn_fpol = np.nanmean(epz_mofpol[neg_fsub] - epz_moh[pos_fsub],axis=0)/30.
diff_divn_fpol = np.nanmean(div_mofpol[neg_fsub] - div_moh[pos_fsub],axis=0)/30.

diff_div = [diff_divn,diff_divp,diff_divno,diff_divn_fpol,diff_divp_fpol,diff_divno_fpol]
diff_epy = [diff_epyn,diff_epyp,diff_epyno,diff_epyn_fpol,diff_epyp_fpol,diff_epyno_fpol]
diff_epz = [diff_epzn,diff_epzp,diff_epzno,diff_epzn_fpol,diff_epzp_fpol,diff_epzno_fpol]

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

plt.savefig(directoryfigure + 'epflux_QBO_regional.png',dpi=300)
print('Completed: Script done!')

