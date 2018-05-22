"""
Script reads in DAILY data from WACCM4 regional experiments (FSUB,FPOL)
for all 200 ensemble members! This data has already been averaged over 
latitude and longitude.
 
Notes
-----
    Author : Zachary Labe
    Date   : 16 May 2018
    
Usage
-----
    readMeanExperiAllRegional(varid,experi,level)
"""

def readMeanExperiAllRegional(varid,experi,level):
    """
    Function reads daily data from WACCM4 regional simulations that are 
    already averaged over latitude/longitude (See script called
    calc_global_ave.proc)

    Parameters
    ----------
    varid : string
        variable name to read
    experi : string
        experiment name (CIT or HIT or FIT or FIC or FICT)
    level : string
        Height of variable (surface or profile)
        

    Returns
    -------
    lat : 1d numpy array
        latitudes
    lon : 1d numpy array
        longitudes
    time : 1d numpy array
        standard time (days since 1870-1-1, 00:00:00)
    var : 1d or 2d numpy array
        [days] or [days,lev]

    Usage
    -----
    lat,lon,time,lev,var = readMeanExperiAll(varid,experi,level)
    """
    print('\n>>> Using readExperiAll function! \n')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Directory 1 for ensemble members 1-100 (mounted server - GreenPlanet)
    directorydata1 = '/home/zlabe/green/simu/'
    
    ### Directory 2 for ensemble members 101-200 (remote server - Surtsey)
    directorydata2 = '/surtsey/zlabe/simu/'
    
    ### Number of ensembles (years)
    ENS1 = 100
    
    ### Call files
    if level == 'surface': # 1d variable
        var = np.empty((ENS1,212,96,144))
        var2 = np.empty((ENS1,212,96,144))
    elif level == 'profile': # 2d variable
        var = np.empty((ENS1,212,17))
        var2 = np.empty((ENS1,212,17))

    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory 1 (1-100 members)       
    for i in range(0,ENS1,1):
        
        totaldirectory = directorydata1 + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
        filename = totaldirectory + varid + '_%s_' % (i+1) + 'mean.nc'
        
        if varid == 'Z500':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'T1000':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        
        ### Read in Data
        if level == 'surface': # 3d variables
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = 'surface'
            lat = data.variables['latitude'][:]
            lon = data.variables['longitude'][:]
            var[i,:,:,:] = data.variables['%s' % varid][:].squeeze()
            data.close()
        elif level == 'profile': # 4d variables
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = data.variables['level'][:]
            lat = data.variables['latitude'][:]
            lon = data.variables['longitude'][:]
            var[i,:,:] = data.variables['%s' % varid][:]
            data.close()
        else:
            print(ValueError('Wrong height - (surface or profile!)!'))    
        print('Completed: Read data for *%s%s* : %s!' % (experi[:3],
                                                         i+1,varid))
    print('Completed: Read members 1-200!')    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory 2 (101-200 members)       
    for i in range(0,ENS1,1):
        
        totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                        '%s/' % (i+101)
        filename2 = totaldirectory2 + varid + '_%s_' % (i+101) + 'mean.nc'
        
        if varid == 'Z500':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'T1000':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        
        ### Read in Data
        if level == 'surface': # 3d variables
            data2 = Dataset(filename2,'r')
            time2 = data.variables['time'][:]
            lev2 = 'surface'
            lat2 = data.variables['latitude'][:]
            lon2 = data.variables['longitude'][:]
            var2[i,:,:,:] = data.variables['%s' % varid][:].squeeze()
            data2.close()
        elif level == 'profile': # 4d variables
            data2 = Dataset(filename2,'r')
            time2 = data.variables['time'][:]
            lev2 = data.variables['level'][:]
            lat2 = data.variables['latitude'][:]
            lon2 = data.variables['longitude'][:]
            var2[i,:,:] = data.variables['%s' % varid][:]
            data2.close()
        else:
            print(ValueError('Wrong height - (surface or profile!)!'))    
        print('Completed: Read data for *%s%s* : %s!' % (experi[:3],
                                                         i+101,varid))
        
    print('Completed: Read members 101-200!')    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Concatenate 1-200 members
    varall = np.append(var,var2,axis=0)
        
    print('\n*Completed: Finished readExperiAllRegional function!')
    return lat,lon,time,lev,varall
        
### Test function -- no need to use    
#varid = 'GEOP'
#experi = 'FSUB'
#level = 'profile'  
#lat,lon,time,lev,var = readMeanExperiAllRegional(varid,experi,level)