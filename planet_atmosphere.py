"""
Name
----
planet_atmosphere.py

Description
-----------
RETrO: Refraction in Exoplanet Transit Observations

This script contains the functions that determine the properties of the planetary atmos-
phere. Each function is called from shoulder.py at different times. 

Right now, the atmosphere is set up as desribed in Dalba (2017), but the user can alter
the atmosphere in any way they see fit, as long as the functions returned what the 
shoulder.py code is expecting. This should be obvious by the names of the functions. 

Input
-----
Various information about the planetary atmosphere.

Output
------
Various information about the planetary atmosphere.

Author
------
P. A. Dalba --- Boston University
pdalba -at- bu -dot- edu
"""
#-----------------------------------------------------------------------------------------
#Import various math, science, and plotting packages.
import numpy
from numpy import *   #Note, numpy functions will not be explicitly called out
import scipy
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy import optimize
from scipy.optimize import minimize_scalar, fsolve
import os
import datetime
import pickle
#-----------------------------------------------------------------------------------------


#Define fundamental constants - All physical units in this script will be MKS.
#-----------------------------------------------------------------------------------------
k_B     = 1.38065e-23   #mks
m_H     = 1.67353e-27   #kg
G       = 6.67408e-11   #mks
AU      = 1.49598e11    #m
R_earth = 6.371e6       #m
M_earth = 5.9723e24     #kg
R_sun   = 6.957e8       #m
#-----------------------------------------------------------------------------------------


#The first function returns the mean molecular *mass* and the reference refractivity 
# based on the string qualifier of the atmosphere. This is hardcoded as described in 
# Dalba (2017). Note that for the H2/He atmosphere, a solar helium mass fraction is used.
def get_mu_nuref(atm_type):
	#Atmosphere type. 4 possible options: H2, N2, H2O, CO2.
	if atm_type == 'H2':
		#Set a helium mass fraction for this atmosphere. 
		Y = 0.25   #solar
		#Find mu assuming only helium and H2
		mu = 1./(Y/(4.*m_H)+(1.-Y)/(2.*m_H))   #kg
	
		#The reference refractivity must come from the mole fraction, which can be found
		# from the mass fraction. 
		f_He = 1./(2./Y - 1.)
		#The refractivities come from the NPL (Kaye & Laby) at STP, which means 1.01325 bar,
		# visible light, 273.15 K. Have to assume little change in refractivity solely due to 
		# temperature. This value will also have to be corrected for the fact that the STP
		# number density is not necessarily the ref number density in this atmosphere.
		nu_ref = (1.-f_He)*1.32e-4 +f_He*3.5e-5
		return mu, nu_ref
	if atm_type == 'N2':
		#Find mu assuming only N2
		mu = 28.*m_H
		#The refractivities come from the NPL (Kaye & Laby) at STP, which means 1.01325 bar,
		# visible light, 273.15 K. Have to assume little change in refractivity solely due to 
		# temperature. This value will also have to be corrected for the fact that the STP
		# number density is not necessarily the ref number density in this atmosphere.
		nu_ref = 2.98e-4
		return mu, nu_ref
	if atm_type == 'CO2':
		#Find mu assuming only CO2
		mu = 44.*m_H
		#The refractivities come from the NPL (Kaye & Laby) at STP, which means 1.01325 bar,
		# visible light, 273.15 K. Have to assume little change in refractivity solely due to 
		# temperature. This value will also have to be corrected for the fact that the STP
		# number density is not necessarily the ref number density in this atmosphere.
		nu_ref = 4.49e-4
		return mu, nu_ref	
	if atm_type == 'H2O':
		#Find mu assuming only H2O
		mu = 18.*m_H
		#The refractivities come from the NPL (Kaye & Laby) at STP, which means 1.01325 bar,
		# visible light, 273.15 K. Have to assume little change in refractivity solely due to 
		# temperature. This value will also have to be corrected for the fact that the STP
		# number density is not necessarily the ref number density in this atmosphere.
		nu_ref = 2.56e-4
		return mu, nu_ref
	print 'atm_type not specified correctly'
	return nan	


#The next set of function each retrieve an individual property of the atmosphere.
#-----------------------------------------------------------------------------------------
def get_temperature(z, lat):
	#Isothermal atmosphere
	return T_atm

def get_number_density(z, lat, nd_ref, z_ref, H):
	#Hydrostatic equilibrium and an isothermal atmosphere is assumed to recover the 
	# familiar exponential density profile.
	return nd_ref*exp(-(z-z_ref)/H)
	
def get_pressure(z, lat):
	#From the temperature and number density, use the Ideal Gas Law to find the pressure.
	return get_number_density(z=z,lat=lat, nd_ref=nd_ref, z_ref=z_ref, H=H)*k_B*\
		get_temperature(z=z,lat=lat)	
	
def get_refractivity(z, lat, z_top, nu_ref, z_ref, H):
	#The refractivity follows the number density (or pressure) profiles as an exponential  
	# profile within the atmosphere and free space outside of the atmosphere
	if z > z_top:
		return 0.
	return nu_ref*exp(-(z-z_ref)/H)
	
def get_drefractivity_dz(z, lat, z_top, nu_ref, z_ref, H):
	#Take the z-derivative of the refractivity profile within the atmosphere, 0 elsewhere
	if z > z_top:
		return 0.
	return -(nu_ref/H)*exp(-(z-z_ref)/H)

def get_drefractivity_dlat(z, lat): 
	#Assume no latitudinal variation - the refractivity variation is entirely radial.
	return 0.

def get_ray_curvature(z, lat, beta, dndz, dndlat, z_top, nu_ref, z_ref, H):  
	#Calculate the ray curvature using Eq. 2b of van der Werf (2008) within the atmoshpere
	# and assume it is zero outside of the atmosphere (straight-line rays).
	if z > z_top:
		return 0.
	ior = 1. + get_refractivity(z=z, lat=lat, z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, \
		H=H)
	return (1./ior)*(cos(beta)*dndz - (sin(beta)/z)*dndlat)
#-----------------------------------------------------------------------------------------





