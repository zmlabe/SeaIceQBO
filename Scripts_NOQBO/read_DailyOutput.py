"""
Script reads in daily data from WACCM4 experiment (NOQBO)
 
Notes
-----
    Author : Zachary Labe
    Date   : 17 January 2019
    
Usage
-----
    readMeanExperi(directory,varid,experi,level)
"""

def readMeanExperi(directory,varid,experi,level):
    """
    Function reads daily data from WACCM4 simulation that are 
    already averaged over latitude/longitude (See script called
    calc_global_ave.proc)

    Parameters
    ----------
    directory : string
        working directory for stored WACCM4 experiments (remote server)
    varid : string
        variable name to read
    experi : string
        experiment name (CIT or HIT or FIT)
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
    lat,lon,time,lev,var = readMeanExperi(directory,varid,experi,level)
    """
    print('\n>>> Using readExperi function! \n')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Number of ensembles (years)
    ENS = 141
    
    ### Call files
    if level == 'surface': # 1d variable
        var = np.empty((ENS,212,96,144))
    elif level == 'profile': # 2d variable
        var = np.empty((ENS,212,17))
        
    for i in range(0,ENS,1):        
        totaldirectory = directory + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
        filename = totaldirectory + varid + '_%s_' % (i+1) + 'mean.nc'
        
        if varid == 'Z500':
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
        
    print('\n*Completed: Finished readExperi function!')
    return lat,lon,time,lev,var
        
#### Test function -- no need to use    
#directory = '/seley/zlabe/simu/'
#varid = 'V'
#experi = 'NOQBO'
#level = 'profile'  
#lat,lon,time,lev,var = readMeanExperi(directory,varid,experi,level)