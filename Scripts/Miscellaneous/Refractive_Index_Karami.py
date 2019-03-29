import warnings
import datetime
import math
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib as mpl
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import nclcmaps as ncm




# of years of data you have minus one
# Example: 1976 to 1978 requires loop_max of two
loop_max = 1

# Wave Information #
zonal_wave_number = 1
meridional_wave_number = 1

# Dec. Indexing #
start_day_dec = 365-31
end_day_dec = 365

# Jan/Feb Indexing #
start_day_janfeb = 0
end_day_janfeb = 31+28

# Degrees/radians Conversion #
deg_2_rad = math.pi/180.0

# Data Dimensions #
nlat = 73
nlev = 17
nlon = 144
nday = 31+31+28
day_num = np.linspace(1,nday,num=nday)
ilats = np.linspace(1,nlat,num=nlat)

# Plotting later #

# Latitude range of x-axis #

lat_north = 90
lat_south = 20

# Pressure levels of y-axis #

desired_levs = [1000,500,200,100,50,10]

# Years for title #

year_storage = []

# RI and Probability Data Storage #
Refractive_Index_Final = np.zeros((loop_max,nlat,nlev))
Probability_Final = []

# Create functions to reshape arrays #

# Modify dimensions of Brunt Vaisala and N2_z for later

def array_reshape(arr):

	arr_storage1 = []

	for i, val in enumerate(day_num):

		tmp = arr[:,:,i,:]
		arr_storage1.append(tmp)

	arr_storage1 = np.array(arr_storage1)
	
	arr_storage2 = []

	for i, val in enumerate(ilats):

		tmp = arr_storage1[:,:,i,:]
		arr_storage2.append(tmp)

	arr_storage2 = np.array(arr_storage2)
	return arr_storage2

def three_D_array(arr):

	tmp1 = ([arr,]*nday)
	tmp2 = np.array(tmp1).transpose()
	tmp3 = tmp2[:,:,np.newaxis]
	arr2 = np.repeat(tmp3,nlev,axis=2)
	return arr2

def u_bar_reshape(arr):

	arr_storage1 = []

	for i, val in enumerate(day_num):
		tmp = arr[:,i,:]
		arr_storage1.append(tmp)

	arr_storage1 = np.array(arr_storage1)

	arr_storage2 = []

	for i, val in enumerate(ilats):
		tmp = arr_storage1[:,:,i]
		arr_storage2.append(tmp)

	arr_storage2 = np.array(arr_storage2)
	return arr_storage2

# Begin Loop #
total_years = np.linspace(1976,2018,num=2018-1976+1)
max_year = total_years[0:loop_max+1]

for year_index, year in enumerate(max_year):

	# Load in all of the data #

	year = int(year)
	year_storage.append(year)

	year_next = year + 1
	if year_next not in max_year:
		print year_next
		break

	U1 = Dataset('U_%s.nc' % year)
	U1_data = U1['U']
	print U1_data
	stop

	T1 = Dataset('TEMP_%s.nc' % year)
	T1_data = T1['TEMP']

	U2 = Dataset('U_%s.nc' % year_next)
	U2_data = U2['U']

	T2 = Dataset('TEMP_%s.nc' % year_next)
	T2_data = T2['TEMP']

	# Get data dimensions #

	lon = U1['longitude']
	lat = U1['latitude']
	lev = U1['level']

	# Concatenate DEC of current year with JF of next year

	# Zonal Wind #

	U1_data_DEC = U1_data[start_day_dec:end_day_dec,:,:,:]
	U2_data_JF = U2_data[start_day_janfeb:end_day_janfeb,:,:,:]

	U_wind_DJF = np.concatenate((U1_data_DEC,U2_data_JF),axis=0)

	# Temperature #

	T1_data_DEC = T1_data[start_day_dec:end_day_dec,:,:,:]
	T2_data_JF = T2_data[start_day_janfeb:end_day_janfeb,:,:,:]

	T_DJF = np.concatenate((T1_data_DEC,T2_data_JF),axis=0)

	# Calculate the altitude and density at each isobar #

	surf_pres = 101300 # pascal
	scale_height = 8000 # meters
	density_surf = 1.2 # kg/m^3

	# Height - Confirmed #

	height = []

	for isobar in lev:

		# Equation = -scale height X log((100*isobar)/surface pressure)
		factor = (-1*scale_height*math.log(np.true_divide(np.multiply(100,isobar),surf_pres)))
		height.append(factor)

	# Density - Confirmed #

	density = []

	for alt in height:
		factor = np.multiply(np.exp(np.true_divide(-1*alt,scale_height)),density_surf)
		density.append(factor)

	# Potential Temperature - Confirmed #
	# shape = 17, 90, 73, 144 #

	omega = 7.2921E-5 # s^-1
	earths_radius = 6371000 # meters
	RCp = 0.286 #  R/Cp= 0.286

	# Reshape Temp array so that dimensions match up later #

	test_temp = []

	for ind, l in enumerate(lat):

		tmp = T_DJF[:,:,ind,:]
		tmp = np.array(tmp)
		test_temp.append(tmp)

	T_DJF = np.array(test_temp)
	print T_DJF.shape

	# Calculate potential temperature between isobars #

	potential_temperature = []

	for ind,isobar in enumerate(lev):

		T_lev = T_DJF[:,:,ind,:]
		factor = np.true_divide(surf_pres,100*isobar)**RCp
		pot_temp = np.multiply(T_lev,factor)
		potential_temperature.append(pot_temp)

	potential_temperature = np.array(potential_temperature)

	# One sided finite difference and center finite difference method

	potential_temperature = np.array(potential_temperature)

	TETAz = []

	for ind, alt in enumerate(height):
		if ind == 0:
			factor = np.true_divide(potential_temperature[1,:,:,:]-potential_temperature[0,:,:,:],(height[1]-height[0]))
			TETAz.append(factor)
		if ind >= 1 and ind <= 15:
			factor = np.true_divide(potential_temperature[ind+1,:,:,:]-potential_temperature[ind-1,:,:,:],(height[ind+1]-height[ind-1]))
			TETAz.append(factor)
		if ind == 16:
			factor = np.true_divide(potential_temperature[16,:,:,:]-potential_temperature[15,:,:,:],(height[16]-height[15]))
			TETAz.append(factor)

	TETAz = np.array(TETAz)

	# Brunt Vaisala - Confirmed # 
	# Angular frequency at which a vertically displaced parcel will oscillate
	# within a statically stable environment
	# Does it only apply in the stratopshere then?

	Brunt_vaisala = np.multiply(np.true_divide(9.8,potential_temperature),TETAz)
	Brunt_vaisala_square = np.multiply(Brunt_vaisala,Brunt_vaisala)

	# Average the zonal wind zonally - Confirmed #

	u_bar = np.mean(U_wind_DJF,axis=3)

	# Calculate change in zonal wind with z - between isobars #

	# One sided finite difference and center finite difference method

	u_bar_z = []

	for ind, alt in enumerate(height):
		if ind == 0:
			factor = np.true_divide(u_bar[:,1,:]-u_bar[:,0,:],(height[1]-height[0]))
			u_bar_z.append(factor)
		if ind >= 1 and ind <= 15:
			factor = np.true_divide(u_bar[:,ind+1,:]-u_bar[:,ind-1,:],(height[ind+1]-height[ind-1]))
			u_bar_z.append(factor)
		if ind == 16:
			factor = np.true_divide(u_bar[:,16,:]-u_bar[:,15,:],(height[16]-height[15]))
			u_bar_z.append(factor)

	u_bar_z = np.array(u_bar_z)

	# Take the derivative of u_bar with respect to heingt #

	u_bar_z_z = []

	for ind, alt in enumerate(height):
		if ind == 0:
			factor = np.true_divide(u_bar_z[1,:,:]-u_bar_z[0,:,:],(height[1]-height[0]))
			u_bar_z_z.append(factor)
		if ind >= 1 and ind <= 15:
			factor = np.true_divide(u_bar_z[ind+1,:,:]-u_bar_z[ind-1,:,:],(height[ind+1]-height[ind-1]))
			u_bar_z_z.append(factor)
		if ind == 16:
			factor = np.true_divide(u_bar_z[16,:,:]-u_bar_z[15,:,:],(height[16]-height[15]))
			u_bar_z_z.append(factor)	

	u_bar_z_z = np.array(u_bar_z_z)

	# Take the derive of density with respect to z - Confirmed #

	rho_z = []

	for ind, alt in enumerate(height):
		if ind == 0:
			factor = np.true_divide(density[1]-density[0],(height[1]-height[0]))
			rho_z.append(factor)
		if ind >= 1 and ind <= 15:
			factor = np.true_divide(density[ind+1]-density[ind-1],(height[ind+1]-height[ind-1]))
			rho_z.append(factor)
		if ind == 16:
			factor = np.true_divide(density[16]-density[15],(height[16]-height[15]))
			rho_z.append(factor)	

	rho_z = np.array(rho_z)

	# Take the derivative of Brunt Vaisala with respect to z #

	# One sided finite difference and center finite difference method

	N2_z = []

	for ind, alt in enumerate(height):
		if ind == 0:
			factor = np.true_divide(Brunt_vaisala_square[1,:,:,:]-Brunt_vaisala_square[0,:,:,:],(height[1]-height[0]))
			N2_z.append(factor)
		if ind >= 1 and ind <= 15:
			factor = np.true_divide(Brunt_vaisala_square[ind+1,:,:,:]-Brunt_vaisala_square[ind-1,:,:,:],(height[ind+1]-height[ind-1]))
			N2_z.append(factor)
		if ind == 16:
			factor = np.true_divide(Brunt_vaisala_square[16,:,:,:]-Brunt_vaisala_square[15,:,:,:],(height[16]-height[15]))
			N2_z.append(factor)

	N2_z = np.array(N2_z)

	# Modify dimensions of Brunt Vaisala and N2_z for later

	Brunt_vaisala = array_reshape(Brunt_vaisala)
	Brunt_vaisala_square = array_reshape(Brunt_vaisala_square)
	N2_z = array_reshape(N2_z)

	# Take the derivative of zonal wind with respect to change in latitude - Confirmed #

	# One sided finite difference and center finite difference method

	u_bar_phi = []

	for ind, val in enumerate(lat):
		if ind == 0:
			factor = np.true_divide(u_bar[:,:,1]-u_bar[:,:,0],0.043)
			u_bar_phi.append(factor)
		if ind >= 1 and ind <= 71:
			factor = np.true_divide(u_bar[:,:,ind+1]-u_bar[:,:,ind-1],0.087)
			u_bar_phi.append(factor)
		if ind == 72:
			factor = np.true_divide(u_bar[:,:,72]-u_bar[:,:,71],0.043)
			u_bar_phi.append(factor)

	u_bar_phi = np.array(u_bar_phi)

	# Weighting u_bar by the cosine of latitude - Confirmed #

	# In tmp1, I'm multiplying by the latitude index, rather than the latitude
	# I believe Karami meant to multiply by the latitude at the corresponding index
	# not just the index

	# Corrected in this script as of May 21, 2018 #

	u_bar_cosphi_phi = []

	for index, value in enumerate(lat):

		#tmp1 = np.multiply(deg_2_rad,ilats[index])
		tmp1 = np.multiply(deg_2_rad,value)
		factor1 = np.cos(tmp1)

		tmp2 = np.multiply(deg_2_rad,value)
		factor2 = -1 * np.sin(tmp2)

		factor3 = np.multiply(u_bar_phi[index],factor1)
		factor4 = np.multiply(u_bar[:,:,index],factor2)

		u_bar_cosphi_phi_tmp = np.add(factor3,factor4)
		u_bar_cosphi_phi.append(u_bar_cosphi_phi_tmp)

	u_bar_cosphi_phi = np.array(u_bar_cosphi_phi)

	# Taking the derivative of u_bar_cosphi_phi with respect to phi - Confirmed #

	u_bar_cosphi_phi_phi = []

	for ind, val in enumerate(lat):
		if ind == 0:
			factor = np.true_divide(u_bar_cosphi_phi[1,:,:]-u_bar_cosphi_phi[0,:,:],0.043)
			u_bar_cosphi_phi_phi.append(factor)
		if ind >= 1 and ind <= 71:
			factor = np.true_divide(u_bar_cosphi_phi[ind+1,:,:]-u_bar_cosphi_phi[ind-1,:,:],0.087)
			u_bar_cosphi_phi_phi.append(factor)
		if ind == 72:
			factor = np.true_divide(u_bar_cosphi_phi[72,:,:]-u_bar_cosphi_phi[71,:,:],0.043)
			u_bar_cosphi_phi_phi.append(factor)

	u_bar_cosphi_phi_phi = np.array(u_bar_cosphi_phi_phi)

	# Now u_bar_cosphi_phi_overphi_phi

	u_bar_cosphi_phi_overphi_phi = []

	for index, val in enumerate(lat):

		tmp1 = np.multiply(deg_2_rad,(val))
		tmp2 = np.cos(tmp1)**2
		factor1 = np.true_divide(1,tmp2)

		factor2 = np.cos(tmp1)
		factor3 = np.sin(tmp1)

		factor4 = u_bar_cosphi_phi_phi[index,:,:]
		factor5 = u_bar_cosphi_phi[index,:,:]

		factor6 = np.multiply(factor1,factor4)
		factor7 = np.multiply(factor6,factor2)

		factor8 = -1 * np.multiply(factor3,factor5)

		u_bar_cosphi_phi_overphi_phi_tmp = np.subtract(factor7,factor8)
		u_bar_cosphi_phi_overphi_phi.append(u_bar_cosphi_phi_overphi_phi_tmp)

	u_bar_cosphi_phi_overphi_phi = np.array(u_bar_cosphi_phi_overphi_phi)

	# Calculate the Coriolis parameter squared # 

	coriolis_parameter_square = []

	for index, value in enumerate(lat):

		tmp1 = np.multiply(deg_2_rad,value)
		tmp2 = np.sin(tmp1)
		tmp3 = 2 * np.multiply(tmp2,omega)
		coriolis_parameter_square_tmp = tmp3**2
		coriolis_parameter_square.append(coriolis_parameter_square_tmp)

	coriolis_parameter_square = np.array(coriolis_parameter_square)


	# There are 3 terms in meridional gradient of the zonal mean potential vorticity #

	# Part 1 #

	meridional_gradient_pv_1 = []

	for index, value in enumerate(lat):

		tmp1 = np.multiply(deg_2_rad,value)
		tmp2 = np.cos(tmp1)
		tmp3 = 2 * np.multiply(omega,tmp2)
		meridional_gradient_pv_1_tmp = np.true_divide(tmp3,earths_radius)
		meridional_gradient_pv_1.append(meridional_gradient_pv_1_tmp)

	meridional_gradient_pv_1 = np.array(meridional_gradient_pv_1)
	meridional_gradient_pv_1 = three_D_array(meridional_gradient_pv_1)


	# Part 2 #

	tmp1 = np.true_divide(-1,(earths_radius**2))
	meridional_gradient_pv_2 = np.multiply(tmp1,u_bar_cosphi_phi_overphi_phi)

	# Part 3 #

	# meridional_gradient_pv_3 consist 3 functions (density, u_bar_z and static
	# stability of atmosphere which are function of z); according to Li et al .
	# 2007. Other terms can be derived from sun et al 2014. 

	# Part 3.1 #

	u_bar_z = u_bar_reshape(u_bar_z)
	u_bar_z_z = u_bar_reshape(u_bar_z_z)


	tmp1 = scale_height * np.mean(Brunt_vaisala,axis=3)
	coriolis_parameter_square = coriolis_parameter_square[:,np.newaxis,np.newaxis]
	tmp2 = np.true_divide(coriolis_parameter_square,tmp1)
	meridional_gradient_pv_3_1 = np.multiply(tmp2,u_bar_z)

	meridional_gradient_pv_3_1 = np.array(meridional_gradient_pv_3_1)


	# meridional_gradient_pv_3_1(ilat,ilev,iday)=(coriolis_parameter_square(ilat)/
	# (scale_height*mean(Brunt_vaisala(:,ilat,ilev,iday))))*u_bar_z(ilat,ilev,iday);


	# Part 3.2 #

	tmp1 = np.mean(N2_z,axis=3)
	tmp2 = np.mean(Brunt_vaisala,axis=3)

	tmp3 = np.multiply(coriolis_parameter_square,tmp1)

	tmp4 = np.multiply(tmp2,tmp2)
	tmp5 = np.true_divide(tmp3,tmp4)

	meridional_gradient_pv_3_2 = np.multiply(tmp5,u_bar_z)

	# Part 3.3 #

	tmp1 = np.mean(Brunt_vaisala,axis=3)
	tmp2 = np.true_divide(coriolis_parameter_square,tmp1)

	meridional_gradient_pv_3_3 = -1 * np.multiply(tmp2,u_bar_z_z)

	# Part 3 Final #

	tmp1 = np.add(meridional_gradient_pv_3_1,meridional_gradient_pv_3_2)
	meridional_gradient_pv_3 = np.add(tmp1,meridional_gradient_pv_3_3)

	# Sum the three meridional gradient of the zonal mean potential vorticity terms #

	tmp1 = np.add(meridional_gradient_pv_1,meridional_gradient_pv_2)
	meridional_gradient_pv = np.add(tmp1,meridional_gradient_pv_3)



	# Refactive Index for stationary planetary waves have 4 terms #

	# RF 1 #

	u_bar_storage = []

	for i, val in enumerate(ilats):

		tmp = u_bar[:,:,i]
		u_bar_storage.append(tmp)

	u_bar_storage = np.array(u_bar_storage)
	u_bar = u_bar_storage

	Refractive_Index_1 = np.true_divide(meridional_gradient_pv,u_bar)

	# RF 2 #

	Refractive_Index_2 = []

	for index, value in enumerate(lat):

		tmp1 = np.multiply(deg_2_rad,value)
		tmp2 = np.cos(tmp1)
		tmp3 = np.multiply(earths_radius,tmp2)
		tmp4 = np.true_divide(zonal_wave_number,tmp3)**2
		tmp5 = np.multiply(-1,tmp4)
		Refractive_Index_2.append(tmp5)

	Refractive_Index_2 = np.array(Refractive_Index_2)
	Refractive_Index_2 = three_D_array(Refractive_Index_2)

	# RF 3 #

	Refractive_Index_3 = []

	for index, value in enumerate(lat):

		tmp1 = np.multiply(deg_2_rad,value)
		tmp2 = np.cos(tmp1)
		tmp3 = 2 * np.multiply(earths_radius,tmp2)
		tmp4 = np.true_divide((meridional_wave_number*math.pi),tmp3)**2
		tmp5 = np.multiply(-1,tmp4)
		Refractive_Index_3.append(tmp5)

	Refractive_Index_3 = np.array(Refractive_Index_3)
	Refractive_Index_3 = three_D_array(Refractive_Index_3)
	
	# RF 4 #

	tmp1 = 4 * np.multiply(np.mean(Brunt_vaisala,axis=3),(scale_height**2))
	Refractive_Index_4 = -1 * np.true_divide(coriolis_parameter_square,tmp1)

	# Calculate the final Refractive Index #

	tmp1 = np.add(Refractive_Index_1,Refractive_Index_2)
	tmp2 = np.add(Refractive_Index_3,Refractive_Index_4)
	tmp3 = np.add(tmp1,tmp2)
	Refractive_Index_Final1 = np.multiply((earths_radius**2),tmp3)

	# Fuzzy logic membership value function #

	ROCF = []

	a = -500
	b = 500

	some_linpspace = np.linspace(1,1301,num=1301)

	for index, ijob in enumerate(some_linpspace):
		ijob1 = ijob - 501
		tester = []
		if ijob1 <= a:
			MVF = 0
			ROCF.append(MVF)
		if ijob1 >= b:
			MVF = 1
			ROCF.append(MVF)
		if ijob1 > a and ijob1 < ((a+b)/2):
			MVF = 1/(1+math.exp(-0.1*(ijob1)))
			ROCF.append(MVF)
		if ijob1 >= ((a+b)/2) and (ijob1 < b):
			MVF = 1.0-2*(((ijob1-b)/(b-a))**2)
			ROCF.append(MVF)
	

	# Probability #

	Probability = []

	Refractive_Index_Final1_flat = Refractive_Index_Final1.flatten()

	for index, value in enumerate(Refractive_Index_Final1_flat):
		Prob = []
		if value >= 800:
			Prob = 1
			Probability.append(Prob)
		if value <=-500:
			Prob = 0
			Probability.append(Prob)
		if value > -499.999 and value < 800:
			ind = int((math.floor(value)+500))
			Prob = ROCF[ind]
			Probability.append(Prob)

	Probability = np.array(Probability)
	Probability = np.reshape(Probability, (nlat,nday,nlev))
	Probability_Final.append(Probability)

Probability_Final = np.array(Probability_Final)
Probability_Final1 = np.mean(Probability_Final,axis=2)
Probability_Final2 = np.multiply(Probability_Final1,100)
Probability_Final3 = np.mean(Probability_Final2,axis=0)

Probability_Plot = Probability_Final3.transpose()

ys = [0,2,4,6,8,10,12,14,15]
y_levs = lev[ys]

south_index,south = min(enumerate(lat), key=lambda x: abs(x[1]-lat_south))
north_index,north = min(enumerate(lat), key=lambda x: abs(x[1]-lat_north))

lat_range = lat[south_index:north_index+1]
lat_range2 = lat_range[0::4]
lat_range2 = [int(x) for x in lat_range2]

lats_plot = []

for val in lat_range2:
	label = str(val)+'N'
	lats_plot.append(label)

# Lat indices #

lat_indices = []

for val in lat_range2:
	ind, actual_value = min(enumerate(lat), key=lambda x: abs(x[1]-val))
	lat_indices.append(ind)

Probability_Plot = Probability_Plot[:,south_index:north_index+1]

lev_num, lat_num = Probability_Plot.shape
lat_spacing = np.linspace(0,lat_num,num=len(lat_indices))

plot_levs = []

for val in desired_levs:
	ind, actual_value = min(enumerate(lev), key=lambda x: abs(x[1]-val))
	plot_levs.append(int(actual_value))

x, y = np.meshgrid(lat_range,lev)

# Plotting #

fig = plt.figure(figsize=(10,8))
mpl.rcParams['font.sans-serif'].insert(0, 'Arial')

#cmap_control = ncm.cmap('MPL_cubehelix_r')
cmap_control = 'CMRmap'
levs1 = np.linspace(0,100,num=21)
ticks = [0,10,20,30,40,50,60,70,80,90,100]


plt.suptitle('%s-%s Probability of Propagation (k,l) = (%s,%s)' % (year_storage[0],year_storage[-1],zonal_wave_number,meridional_wave_number), fontsize=26, fontweight='semibold')


ax1 = fig.add_subplot(111)
plt.contourf(x,y,Probability_Plot,origin='upper',cmap=cmap_control,levels=levs1,extend='both',vmin=0.0,vmax=100.0)
ax1.invert_yaxis()
ax1.set_yscale('log')
plt.ylabel('Pressure (hPa)',fontsize=22,fontweight='semibold',labelpad=0.00001)
plt.yticks(plot_levs, plot_levs, fontsize=20, fontweight='semibold',rotation=0)
plt.xlabel('Latitude',fontsize=22,fontweight='semibold',labelpad=8)
plt.xticks(lat_range2, lats_plot, fontsize=20, fontweight='semibold',rotation=0)
ax1.xaxis.set_tick_params(pad=10)

cbar = plt.colorbar(pad=0.04,shrink=1.0,ticks=levs1[::2],orientation='vertical')
cbar.ax.tick_params(labelsize=22, width=2)
cbar.ax.set_xticklabels(ticks,weight='bold')
cbar.set_label('Probability (%)',fontsize=22)


plt.show()

































