"""
Name
----
retro_rk4.py

Description
-----------
RETrO: Refraction in Exoplanet Transit Observations

This script contains the RK4 integration scheme for the ray tracing portion of RETrO. 

It utilizes many of the funcitons from the planet atmosphere module. For the optical depth,
it only uses a constant cross section for absorption (sigma). The user is free to alter
this for more sophisticated opacity scheme (i.e. Rayleigh scattering). 

Input
-----
Various information about the planetary atmosphere and present ray position.

Output
------
The next ray position.

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

#Import RETrO modules here
import planet_atmosphere
from planet_atmosphere import *
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


#Define the 4-step Runge-Kutta ray integration procedure
#-----------------------------------------------------------------------------------------
#Pass in arrays with the current and parameters in this order: 
#
# Current:   [z, lat, beta, xi, tau, ds, ray] 
# Constant:  [sigma,nd_ref,z_ref,z_top,nu_ref,H]
#
# Remember: z, lat, beta, xi, & tau are all evaluated through DEs. The other parameters
# are just being past through to be found at every point along the integration.
def rk4_ray_trace(current,constant):
	#Unpack the 'constant' array for more clarity
	sigma,nd_ref,z_ref,z_top,nu_ref,H = constant
	
	#First step of rk4 - at the initial point. Calculate any non-integrated values.
	nd1 = get_number_density(z=current[0],lat=current[1],nd_ref=nd_ref,z_ref=z_ref,H=H)
	n1 = get_refractivity(z=current[0],lat=current[1],z_top=z_top,nu_ref=nu_ref, \
		z_ref=z_ref, H=H)
	dndz1 = get_drefractivity_dz(z=current[0],lat=current[1],z_top=z_top,nu_ref=nu_ref, \
		z_ref=z_ref, H=H)
	dndlat1 = get_drefractivity_dlat(z=current[0],lat=current[1])
	c1 = get_ray_curvature(z=current[0],lat=current[1], beta=current[2], dndz=dndz1,\
		dndlat=dndlat1,z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, H=H)
	#Now the integrated parameters
	z1 = sin(current[2])
	lat1 = cos(current[2])/current[0]
	beta1 = lat1 + c1
	xi1 = c1
	tau1 = constant[0]*nd1

	#Second step of rk4 - at the midpoint. Calculate any non-integrated values.
	nd2 = get_number_density(z=(current[0]+0.5*z1*current[5]),\
		lat=(current[1]+0.5*lat1*current[5]),nd_ref=nd_ref, z_ref=z_ref, H=H)
	n2 = get_refractivity(z=(current[0]+0.5*z1*current[5]),\
		lat=(current[1]+0.5*lat1*current[5]),z_top=z_top,nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndz2 = get_drefractivity_dz(z=(current[0]+0.5*z1*current[5]),\
		lat=(current[1]+0.5*lat1*current[5]),z_top=z_top,nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndlat2 = get_drefractivity_dlat(z=(current[0]+0.5*z1*current[5]),\
		lat=(current[1]+0.5*lat1*current[5]))
	c2 = get_ray_curvature(z=(current[0]+0.5*z1*current[5]),\
		lat=(current[1]+0.5*lat1*current[5]), beta=(current[2]+0.5*beta1*current[5]),\
			dndz=dndz2,dndlat=dndlat2,z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, H=H)
	#Now the integrated parameters
	z2 = sin(current[2]+0.5*beta1*current[5])
	lat2 = cos(current[2]+0.5*beta1*current[5])/(current[0]+0.5*z1*current[5])
	beta2 = lat2 + c2
	xi2 = c2
	tau2 = constant[0]*nd2

	#Third step of rk4 - at the midpoint. Calculate any non-integrated values.
	nd3 = get_number_density(z=(current[0]+0.5*z2*current[5]),\
		lat=(current[1]+0.5*lat2*current[5]),nd_ref=nd_ref, z_ref=z_ref, H=H)
	n3 = get_refractivity(z=(current[0]+0.5*z2*current[5]),\
		lat=(current[1]+0.5*lat2*current[5]),z_top=z_top,nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndz3 = get_drefractivity_dz(z=(current[0]+0.5*z2*current[5]),\
		lat=(current[1]+0.5*lat2*current[5]),z_top=z_top,nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndlat3 = get_drefractivity_dlat(z=(current[0]+0.5*z2*current[5]),\
		lat=(current[1]+0.5*lat2*current[5]))
	c3 = get_ray_curvature(z=(current[0]+0.5*z2*current[5]),\
		lat=(current[1]+0.5*lat2*current[5]), beta=(current[2]+0.5*beta2*current[5]),\
			dndz=dndz3,dndlat=dndlat3,z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, H=H)
	#Now the integrated parameters
	z3 = sin(current[2]+0.5*beta2*current[5])
	lat3 = cos(current[2]+0.5*beta2*current[5])/(current[0]+0.5*z2*current[5])
	beta3 = lat3 + c3
	xi3 = c3
	tau3 = constant[0]*nd3
	
	#Fourth step of rk4 - at the end point. Calculate any non-integrated values.
	nd4 = get_number_density(z=(current[0]+z3*current[5]),lat=(current[1]+lat3*current[5]),\
		nd_ref=nd_ref, z_ref=z_ref, H=H)
	n4 = get_refractivity(z=(current[0]+z3*current[5]),lat=(current[1]+lat3*current[5]),\
		z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndz4 = get_drefractivity_dz(z=(current[0]+z3*current[5]),\
		lat=(current[1]+lat3*current[5]),z_top=z_top,nu_ref=nu_ref, z_ref=z_ref, H=H)
	dndlat4 = get_drefractivity_dlat(z=(current[0]+z3*current[5]),\
		lat=(current[1]+lat3*current[5]))
	c4 = get_ray_curvature(z=(current[0]+z3*current[5]),lat=(current[1]+lat3*current[5]),\
		beta=(current[2]+beta3*current[5]),dndz=dndz4,dndlat=dndlat4,z_top=z_top, \
			nu_ref=nu_ref, z_ref=z_ref, H=H)
	#Now the integrated parameters
	z4 = sin(current[2]+beta3*current[5])
	lat4 = cos(current[2]+beta3*current[5])/(current[0]+z3*current[5])
	beta4 = lat4 + c4
	xi4 = c4
	tau4 = constant[0]*nd4

	#Find the "next" values of the integrated parameters.
	z_next = current[0] + (current[5]/6.)*(z1 + 2.*(z2 + z3) + z4) 
	lat_next = current[1] + (current[5]/6.)*(lat1 + 2.*(lat2 + lat3) + lat4)
	beta_next = current[2] + (current[5]/6.)*(beta1 + 2.*(beta2 + beta3) + beta4)
	xi_next = current[3] + (current[5]/6.)*(c1 + 2.*(c2 + c3) + c4)
	tau_next = current[4] + (current[5]/6.)*(tau1 + 2.*(tau2 + tau3) + tau4)

	#We must add a block that checks for a negative z-value. This is not realistic, but can
	# occur in the case of a radially traveling ray. In this case, the z, lat, beta, and 
	# xi values must be manually adjusted due to the singularity that exists when z goes
	# to zero.
	if z_next < 0.:
		#The z-value should actually just be the step size minus the current z-value.
		z_next = current[5] - current[0]
		#As the radial ray passes through zero, its latitude must jump by pi to the other
		# hemisphere. The beta value also jumps by pi (from negative pi/2 to positive 
		# pi/2) as the ray passes through the center. The ray curvature and therefore the 
		# value of xi is also effected by the blown up beta values, so we force there to
		# be no added curvature and just set xi_next as its previous value. 
		lat_next = current[1]+pi
		beta_next = current[2]+pi
		xi_next = current[3]		

	return array([z_next,lat_next,beta_next,xi_next,tau_next])
#-----------------------------------------------------------------------------------------

