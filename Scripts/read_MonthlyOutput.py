"""
Script reads in monthly data from WACCM4 experiments (CIT,HIT,FIT)
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 August 2017
    
Usage
-----
    readExperi(directory,varid,experi,level)
"""

def readExperi(directory,varid,experi,level):
    """
    Function reads monthly data from WACCM4 simulations

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
    var : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]

    Usage
    -----
    lat,lon,time,lev,var = readExperi(directory,varid,experi,level)
    """
    print('\n>>> Using readExperi function! \n')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Call files
    totaldirectory = directory + experi + '/monthly/'
    filename = totaldirectory + varid + '_1900-2000.nc'
    
    if any([experi == 'FPOL',experi == 'FSUB']):
        directory = '/home/zlabe/green/simu/'
        totaldirectory = directory + experi + '/monthly/'
        filename = totaldirectory + varid + '_1900-2000.nc'
    
    if varid == 'EGR' and level == 'surface':
        filename = totaldirectory + varid + '_500_850.nc'
    
    ### Read in Data
    if level == 'surface': # 3d variables
        data = Dataset(filename,'r')
        time = data.variables['time'][:]
        lev = 'surface'
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    elif level == 'profile': # 4d variables
        data = Dataset(filename,'r')
        time = data.variables['time'][:]
        lev = data.variables['level'][:]
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
        varq = data.variables['%s' % varid][:]
        data.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:4],varid))
    
    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var = np.reshape(varq,(int(varq.shape[0]/12),months,
                              int(lat.shape[0]),int(lon.shape[0])))
    elif level == 'profile': # 4d variables
        var = np.reshape(varq,(int(varq.shape[0]/12),months,int(lev.shape[0]),
                      int(lat.shape[0]),int(lon.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var = var - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif varid == 'SWE':
        var = var*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')

    print('\n*Completed: Finished readExperi function!')
    return lat,lon,time,lev,var

### Test function -- no need to use    
#directory = '/surtsey/zlabe/simu/'
#varid = 'T2M'
##varid = 'TEMP'
#experi = 'HIT'
#level = 'surface'
#    
#lat,lon,time,lev,var = readExperi(directory,varid,experi,level)
