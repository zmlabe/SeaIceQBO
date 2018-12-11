"""
Script reads in DAILY data from WACCM4 experiments (HIT2,HIT,FIT,FICT2,FICT)
for all 200 ensemble members! This data has already been averaged over 
latitude and longitude. Note that HIT2=CIT and FICT2=FIC for the file
naming conventions in the filename and function.
 
Notes
-----
    Author : Zachary Labe
    Date   : 15 March 2018
    
Usage
-----
    readMeanExperiAll(varid,experi,level)
"""

def readMeanExperiAll(varid,experi,level):
    """
    Function reads daily data from WACCM4 simulations that are 
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
    
    ### Directory 1 for ensemble members 1-100 (remote server - Surtsey)
    directorydata1 = '/surtsey/zlabe/simu/'
    
    ### Directory 2 for ensemble members 101-200 (mounted server - GreenPlanet)
    directorydata2 = '/home/zlabe/green/simu/'
    
    ### Number of ensembles (years)
    ENS = 100
    
    ### Call files
    if varid == 'MHF100':
        var = np.empty((ENS,212,1,1,144))
        var2 = np.empty((ENS,212,1,1,144))
    elif level == 'surface': # 1d variable
        var = np.empty((ENS,212,96,144))
        var2 = np.empty((ENS,212,96,144))
    elif level == 'profile': # 2d variable
        var = np.empty((ENS,212,17))
        var2 = np.empty((ENS,212,17))
    elif level == 'profile2': 
        var = np.empty((ENS,212,96,144))
        var2 = np.empty((ENS,212,96,144))
    elif level == 'profile3': 
        var = np.empty((ENS,182,17,21,49))
        var2 = np.empty((ENS,182,17,21,49))
    elif level == 'profilesiberia': # 2d variable
        var = np.empty((ENS,212,17))
        var2 = np.empty((ENS,212,17))

    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory 1 (1-100 members)       
    for i in range(0,ENS,1):
        
        totaldirectory = directorydata1 + experi + '/daily/' + experi + \
                        '%s/' % (i+1)
                        
        filename = totaldirectory + varid + '_%s_' % (i+1) + 'mean.nc'
        
        if varid == 'Z500':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'T1000':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'MHF100':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'U200' or varid == 'U300' or varid == 'U700':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'SLP':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif level == 'profile2':
            filename = totaldirectory + varid + '_%s.nc' % (i+1)
        elif varid == 'WAFYv':
            varidn = 'WAFY'
            filename = totaldirectory + varidn + '_%s.nc' % (i+1)
        elif varid == 'WAFYs':
            varidn = 'WAFY'
            filename = totaldirectory + varidn + '_%s_' % (i+1) + 'mean_Siberia.nc'
        elif varid == 'WAFZv':
            varidn = 'WAFZ'
            filename = totaldirectory + varidn + '_%s.nc' % (i+1)
        elif varid == 'WAFZs':
            varidn = 'WAFZ'
            filename = totaldirectory + varidn + '_%s_' % (i+1) + 'mean_Siberia.nc'
        
        ### Read in Data
        if varid == 'MHF100':
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = data.variables['level'][:]
            lat = data.variables['lat'][:]
            lon = data.variables['lon'][:]
            var[i,:,:] = data.variables['%s' % varid][:]
            data.close()  
        elif level == 'surface': # 3d variables
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
        elif level == 'profile2': # 4d variables
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = data.variables['level'][:]
            lat = data.variables['latitude'][:]
            lon = data.variables['longitude'][:]
            var[i,:,:] = data.variables['%s' % varid][:,7,:,:]
            data.close()   
        elif level == 'profile3': # 4d variables
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = data.variables['level'][:]
            lat = data.variables['latitude'][:]
            lon = data.variables['longitude'][:]
            latq = np.where((lat>40) & (lat<80))[0]
            lonq = np.where((lon>=120) & (lon<=240))[0]
            var[i,:,:,:,:] = data.variables['%s' % varidn][30:,:,latq,lonq]
            data.close()  
            lat = lat[latq]
            lon = lon[lonq]
        elif level == 'profilesiberia': # 4d variables
            data = Dataset(filename,'r')
            time = data.variables['time'][:]
            lev = data.variables['level'][:]
            lat = data.variables['latitude'][:]
            lon = data.variables['longitude'][:]
            var[i,:,:] = data.variables['%s' % varidn][:]
            data.close()
            
        else:
            print(ValueError('Wrong height - (surface or profile!)!'))    
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+1,varid))
    print('Completed: Read members 1-200!')    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    
    ### Call files for directory 1 (101-200 members)       
    for i in range(0,ENS,1):
        
        totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                        '%s/' % (i+101)
        filename2 = totaldirectory2 + varid + '_%s_' % (i+101) + 'mean.nc'
        
        if varid == 'Z500':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'T1000':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'MHF100':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'U200' or varid == 'U300' or varid == 'U700':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'SLP':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif level == 'profile2':
            filename2 = totaldirectory2 + varid + '_%s.nc' % (i+101)
        elif varid == 'WAFYv':
            varidn = 'WAFY'
            filename2 = totaldirectory2 + varidn + '_%s.nc' % (i+101)
        elif varid == 'WAFYs':
            varidn = 'WAFY'
            filename2 = totaldirectory2 + varidn + '_%s_' % (i+101) + 'mean_Siberia.nc'
        elif varid == 'WAFZv':
            varidn = 'WAFZ'
            filename2 = totaldirectory2 + varidn + '_%s.nc' % (i+101)
        elif varid == 'WAFZs':
            varidn = 'WAFZ'
            filename2 = totaldirectory2 + varidn + '_%s_' % (i+101) + 'mean_Siberia.nc'
        
        ### Read in Data
        if varid == 'MHF100':
            data2 = Dataset(filename2,'r')
            time2 = data2.variables['time'][:]
            lev2 = data2.variables['level'][:]
            lat2 = data2.variables['lat'][:]
            lon2 = data2.variables['lon'][:]
            var2[i,:,:] = data2.variables['%s' % varid][:]
            data2.close()  
        elif level == 'surface': # 3d variables
            data2 = Dataset(filename2,'r')
            time2 = data2.variables['time'][:]
            lev2 = 'surface'
            lat2 = data2.variables['latitude'][:]
            lon2 = data2.variables['longitude'][:]
            var2[i,:,:,:] = data2.variables['%s' % varid][:].squeeze()
            data2.close()
        elif level == 'profile': # 4d variables
            data2 = Dataset(filename2,'r')
            time2 = data2.variables['time'][:]
            lev2 = data2.variables['level'][:]
            lat2 = data2.variables['latitude'][:]
            lon2 = data2.variables['longitude'][:]
            var2[i,:,:] = data2.variables['%s' % varid][:]
            data2.close()   
        elif level == 'profile2': # 4d variables
            data2 = Dataset(filename,'r')
            time2 = data2.variables['time'][:]
            lev2 = data2.variables['level'][:]
            lat2 = data2.variables['latitude'][:]
            lon2 = data2.variables['longitude'][:]
            var2[i,:,:] = data2.variables['%s' % varid][:,7,:,:]
            data2.close() 
        elif level == 'profile3': # 4d variables
            data2 = Dataset(filename,'r')
            time2 = data2.variables['time'][:]
            lev2 = data2.variables['level'][:]
            lat2 = data2.variables['latitude'][:]
            lon2 = data2.variables['longitude'][:]
            var2[i,:,:,:,:] = data2.variables['%s' % varidn][30:,:,latq,lonq]
            data2.close() 
        elif level == 'profilesiberia': # 4d variables
            data2 = Dataset(filename2,'r')
            time2 = data2.variables['time'][:]
            lev2 = data2.variables['level'][:]
            lat2 = data2.variables['latitude'][:]
            lon2 = data2.variables['longitude'][:]
            var2[i,:,:] = data2.variables['%s' % varidn][:]
            data2.close()   
        else:
            print(ValueError('Wrong height - (surface or profile!)!'))    
        print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                         i+101,varid))
        
    print('Completed: Read members 101-200!')    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    ### Concatenate 1-200 members
    varall = np.append(var,var2,axis=0)
        
    print('\n*Completed: Finished readExperiAll function!')
    return lat,lon,time,lev,varall
        
### Test function -- no need to use    
#varid = 'MHF100'
#experi = 'HIT'
#level = 'none'  
#lat,lon,time,lev,var = readMeanExperiAll(varid,experi,level)