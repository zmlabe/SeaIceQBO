"""
Plot comparisons between WACCM4 sea ice experiments. These are 
sea ice thickness and concentration perturbation experiments. This script is
assess stratospheric sudden warming event counts.

Notes
-----
    Author : Zachary Labe
    Date   : 11 July 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import nclcmaps as ncm
import datetime
import read_DailyOutput_SSW as SO
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
print('\n' '----Plotting SSW Counts for Simulations - %s----' % titletime)

#### Alott time series
year1 = 1800
year2 = 2000
years = np.arange(year1,year2+1,1)

### Add parameters
qbophase = ['pos','non','neg']

### Call functions for SSW counts
sswhit = SO.readDailySSW('count','HIT')
sswfit = SO.readDailySSW('count','FIT')
sswfict = SO.readDailySSW('count','FICT')
sswcon = SO.readDailySSW_CTLQ('count','CTLQ')

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

filenamefitp = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn = directorydata + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
filenamefitp2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[0]
filenamefitno2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[1]
filenamefitn2 = directorydata2 + 'FIT/monthly/QBO_%s_FIT.txt' % qbophase[2]
pos_fit = np.append(np.genfromtxt(filenamefitp,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitp2,unpack=True,usecols=[0],dtype='int')+100)
non_fit = np.append(np.genfromtxt(filenamefitno,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitno2,unpack=True,usecols=[0],dtype='int')+100)
neg_fit = np.append(np.genfromtxt(filenamefitn,unpack=True,usecols=[0],dtype='int'),
                    np.genfromtxt(filenamefitn2,unpack=True,usecols=[0],dtype='int')+100)

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

filenamefctlqp = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[0]
filenamefctlqno = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[1]
filenamefctlqn = directorydata + 'CTLQ/monthly/QBO_%s_CTLQ.txt' % qbophase[2]
pos_ctlq = np.genfromtxt(filenamefctlqp,unpack=True,usecols=[0],dtype='int')
non_ctlq = np.genfromtxt(filenamefctlqno,unpack=True,usecols=[0],dtype='int')
neg_ctlq = np.genfromtxt(filenamefctlqn,unpack=True,usecols=[0],dtype='int')