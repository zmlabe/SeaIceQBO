"""
Script reads in monthly data from WACCM4 experiments (HIT2,HIT,FIT,FICT2,FICT)
for all 200 ensemble members! Note that HIT2=CIT and FICT2=FIC for the file
naming conventions in the filename and function.
 
Notes
-----
    Author : Zachary Labe
    Date   : 15 March 2018
    
Usage
-----
    readExperiAll(directory,varid,experi,level)
"""

def readExperiAll(varid,experi,level):
    """
    Function reads monthly data from WACCM4 simulations for 200 members

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
    var : 4d numpy array or 5d numpy array 
        [year,month,lat,lon] or [year,month,level,lat,lon]

    Usage
    -----
    lat,lon,time,lev,var = readExperiAll(varid,experi,level)
    """
    print('\n>>> Using readExperiAll function! \n')
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Directory 1 for ensemble members 1-100 (remote server - Surtsey)
    directorydata1 = '/surtsey/zlabe/simu/'
    
    ### Directory 2 for ensemble members 101-200 (mounted server - GreenPlanet)
    directorydata2 = '/home/zlabe/green/simu/'
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory 1 (1-100 members)
    totaldirectory = directorydata1 + experi + '/monthly/'
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
        
    print('Completed: Read members 1-100!\n')
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Read in green planet data 
    
    ### Call files for directory 2 (101-200 members)
    totaldirectory2 = directorydata2 + experi + '/monthly/'
    filename2 = totaldirectory2 + varid + '_1900-2000.nc'
    
    if varid == 'EGR' and level == 'surface':
        filename2 = totaldirectory2 + varid + '_500_850.nc'

    ### Read in Data
    if level == 'surface': # 3d variables
        data2 = Dataset(filename2,'r')
        time2 = data.variables['time'][:]
        lev2 = 'surface'
        lat2 = data.variables['latitude'][:]
        lon2 = data.variables['longitude'][:]
        varq2 = data.variables['%s' % varid][:]
        data2.close()
    elif level == 'profile': # 4d variables
        data2 = Dataset(filename2,'r')
        time2 = data.variables['time'][:]
        lev2 = data.variables['level'][:]
        lat2 = data.variables['latitude'][:]
        lon2 = data.variables['longitude'][:]
        varq2 = data.variables['%s' % varid][:]
        data2.close()
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!'))    
    print('Completed: Read data for *%s* : %s!' % (experi[:4],varid))
    
    ### Reshape to split years and months
    months = 12
    if level == 'surface': # 3d variables
        var2 = np.reshape(varq2,(int(varq2.shape[0]/12),months,
                              int(lat2.shape[0]),int(lon2.shape[0])))
    elif level == 'profile': # 4d variables
        var2 = np.reshape(varq2,(int(varq2.shape[0]/12),months,int(lev2.shape[0]),
                      int(lat2.shape[0]),int(lon2.shape[0])))
    else:
        print(ValueError('Selected wrong height - (surface or profile!)!')) 
    print('Completed: Reshaped %s array!' % (varid))
    
    ### Convert units
    if varid in ('TEMP','T2M'):
        var2 = var2 - 273.15 # Kelvin to degrees Celsius 
        print('Completed: Changed units (K to C)!')
    elif varid == 'SWE':
        var2 = var2*1000. # Meters to Millimeters 
        print('Completed: Changed units (m to mm)!')
    
    print('Completed: Read members 101-200!')    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Concatenate 1-200 members
    varall = np.append(var,var2,axis=0)

    print('\n>>>Completed: Finished readExperiAll function!')
    return lat,lon,time,lev,varall

#### Test function -- no need to use    
#varid = 'T2M'
##varid = 'TEMP'
#experi = 'HIT'
#level = 'surface'
#    
#lat,lon,time,lev,var1 = readExperiAll(varid,experi,level)