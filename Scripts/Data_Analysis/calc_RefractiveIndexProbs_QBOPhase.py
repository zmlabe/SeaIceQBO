"""
Script retrieves arrays for refractive indices composited by the QBO phase

Notes
-----
    Author : Zachary Labe
    Date   : 28 September 2018
"""

def callWaveIndex():
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    directorydata = '/surtsey/zlabe/simu/'
    directorydata2 = '/home/zlabe/green/simu/'
    qbophase = ['pos','non','neg']
    
    ### Call functions for variable profile data for refractive index
    def readProb(file,varid,experi):
        """
        Function reads daily data from WACCM4 control for probability of 
        wave propogation (refractive index)
    
        Parameters
        ----------
        file : string
            name of file variable
        varid : string
            variable name to read (LINT60N)
        experi : string
            experiment name (FICT)
            
    
        Returns
        -------
        lon : 1d numpy array
            longitudes
        lev : 1d numpy array
            vertical levels
        var : 2d numpy array
            [time,level,latitude,longitude]
    
        Usage
        -----
        lon,lev,var = readCorrs(file,varid,experi)
        """
        print('\n>>> Using readProb function! \n')
        
        ### Directorys for ensemble members 1-200 (remote server - Surtsey/Green)
        directorydata1 = '/surtsey/zlabe/simu/'
        directorydata2 = '/home/zlabe/green/simu/'
        
        ### Number of ensembles (years)
        ENS1 = 100
           
        ### Call files
        var1 = np.empty((ENS1,7,17,96))
        var2 = np.empty((ENS1,7,17,96))
     
        ### Call files for directory 1 (1-100 members)       
        for i in range(0,ENS1,1):       
            totaldirectory1 = directorydata1 + experi + '/daily/' + experi + \
                            '%s/' % (i+1)
            filename1 = totaldirectory1 + file + '_%s.nc' % (i+1)
            
            ### Read in Data
            data1 = Dataset(filename1,'r')
            var1[i,:,:,:] = data1.variables['%s' % (varid)][:,:,:]
            lev = data1.variables['level'][:]
            data1.close()
      
            print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                             i+1,varid))
        print('Completed: Read members 1-100!') 
    
        ### Call files for directory 1 (101-200 members)       
        for i in range(0,ENS1,1):        
            totaldirectory2 = directorydata2 + experi + '/daily/' + experi + \
                            '%s/' % (i+101)
            filename2 = totaldirectory2 + file + '_%s.nc' % (i+101)
            
            ### Read in Data
            data2 = Dataset(filename2,'r')
            var2[i,:,:,:] = data2.variables['%s' % (varid)][:,:,:]
            lat = data2.variables['latitude'][:]
            data2.close()
      
            print('Completed: Read data for *%s%s* : %s!' % (experi[:4],
                                                             i+101,varid))
        print('Completed: Read members 101-200!')  
        
        ### Append all data    
        var = np.append(var1,var2,axis=0)
        
        print('\n*Completed: Finished readProb function!')
        return lat,lev,var
    
    lat,lev,varfict = readProb('wave1prob','REFRACTIVE_PROB_MONTHLY','FICT')
    lat,lev,varhit = readProb('wave1prob','REFRACTIVE_PROB_MONTHLY','HIT')
    
    ### Create 2d array of latitude and longitude
    lev2,lat2 = np.meshgrid(lev,lat)
    
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
    
    ### Concatonate runs
    var_mo = [varhit,varfict]
    
    ### Composite by QBO phase    
    var_mohitpos = var_mo[0][pos_hit].squeeze()
    var_mohitnon = var_mo[0][non_hit].squeeze()
    var_mohitneg = var_mo[0][neg_hit].squeeze()
    var_mofictpos = var_mo[1][pos_fict].squeeze() 
    var_mofictnon = var_mo[1][non_fict].squeeze()
    var_mofictneg = var_mo[1][neg_fict].squeeze()
    
    return var_mohitpos,var_mohitneg,var_mofictpos,var_mofictneg
