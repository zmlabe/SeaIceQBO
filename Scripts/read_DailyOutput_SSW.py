"""
Script reads in DAILY data from WACCM4 experiments (HIT2,HIT,FIT,FICT2,FICT)
for all 200 ensemble members! This data consists of stratospheric sudden
warming (SSW) events, which are simple ascii files.
 
Notes
-----
    Author : Zachary Labe
    Date   : 10 July 2018
    
Usage
-----
    [1] readDailySSW(vari,experi)
    [2] readDailySSW_CTLQ(vari,experi)
"""

def readDailySSW(vari,experi):
    """
    Function reads daily data from WACCM4 simulations that are statistics
    on SSW events

    Parameters
    ----------
    varid : string
        statistic type
    experi : string
        experiment name (CIT or HIT or FIT or FIC or FICT or FSUB or FPOL)

    Returns
    -------
    ssw : 1d numpy array
        [count]

    Usage
    -----
    ssw = readDailySSW(vari,experi)
    """
    print('\n>>> Using readDailySSW function! \n')
    
    ### Import modules
    import numpy as np
    
    if any([experi=='CIT',experi=='HIT',experi=='FIT',
            experi=='FIC',experi=='FICT']):
        ### Directory 1 for members 1-100 (remote server - Surtsey)
        directorydata1 = '/surtsey/zlabe/simu/'
        
        ### Directory 2 for  members 101-200 (mounted server - GreenPlanet)
        directorydata2 = '/home/zlabe/green/simu/'
        
        ### Number of ensembles (years)
        ENS1 = 100
        ENS2 = 200
    elif any([experi=='FSUB',experi=='FPOL']):
        ### Directory 1 for members 1-100 (mounted server - GreenPlanet)
        directorydata1 = '/home/zlabe/green/simu/'
        
        ### Directory 2 for  members 101-200 (remote server - Surtsey)
        directorydata2 = '/surtsey/zlabe/simu/'
        
        ### Number of ensembles (years)
        ENS1 = 100
        ENS2 = 200
        
    if vari=='count':
        ###################################################################
        ###################################################################
        ###################################################################
        ### Call files for directory 1 (1-100 members)   
        ssw1 = np.empty((ENS1))
        for i in range(0,ENS1,1):
            
            ### Set directories
            totaldirectory1 = directorydata1 + experi + '/daily/' + \
                              experi + '%s/' % (i+1)
            filename1 = totaldirectory1 + 'sswcount_%s.txt' % (i+1)
            
            ### Read in file
            ssw1[i] = np.genfromtxt(filename1,dtype='float',unpack=True)
        
        print('Completed: Read %s ensemble members 1-100!' % experi)
        ###################################################################
        ###################################################################
        ###################################################################
        ### Call files for directory 1 (1-100 members)   
        ssw2 = np.empty((ENS1))
        for i in range(ENS1,ENS2,1):
            
            ### Set directories
            totaldirectory2 = directorydata2 + experi + '/daily/' + \
                              experi + '%s/' % (i+1)
            filename2 = totaldirectory2 + 'sswcount_%s.txt' % (i+1)
            
            ### Read in file
            ssw2[i-100] = np.genfromtxt(filename2,dtype='int',unpack=True)
        
        print('Completed: Read %s ensemble members 101-200!'  % experi)
        
        ### Append lists
        ssw = np.append(ssw1,ssw2,axis=0)
      
    print('\n*Completed: Finished readDailySSW function!')
    return ssw

def readDailySSW_CTLQ(vari,experi):
    """
    Function reads daily data from WACCM4 control runs that are statistics
    on SSW events

    Parameters
    ----------
    varid : string
        statistic type
    experi : string
        experiment name (CTLQ)

    Returns
    -------
    ssw : 1d numpy array
        [count]

    Usage
    -----
    sswc = readDailySSW_CTLQ(vari,experi)
    """
    print('\n>>> Using readDailySSW_CTLQ function! \n')
    
    ### Import modules
    import numpy as np
    
    if any([experi=='CTLQ']):
        ### Directory for members 1-200 (remote server - Surtsey)
        directorydata = '/surtsey/zlabe/simu/'
        
        ### Number of ensembles (years)
        ENS = 200
    else:
        print(ValueError('Wrong experiment!!!'))
        
    if vari=='count':
        ###################################################################
        ###################################################################
        ###################################################################
        ### Call files for directory 1 (1-200 members)   
        sswc = np.empty((ENS))
        for i in range(0,ENS,1):
            
            ### Set directories
            totaldirectory = directorydata + experi + '/daily/' + \
                              experi + '%s/' % (i+1)
            filename = totaldirectory + 'sswcount_%s.txt' % (i+1)
            
            ### Read in file
            sswc[i] = np.genfromtxt(filename,dtype='float',unpack=True)
        
        print('Completed: Read %s ensemble members 1-200!' % experi)
        
    print('\n*Completed: Finished readDailySSW_CTLQ function!')
    return sswc
        
### Test function -- no need to use    
#vari = 'count'
#experi = 'CTLQ'
##ssw = readDailySSW(vari,experi)
##sswc = readDailySSW_CTLQ(vari,experi)