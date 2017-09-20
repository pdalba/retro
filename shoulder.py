"""
Name
----
shoulder.py

Description
-----------
RETrO: Refraction in Exoplanet Transit Observations

This script models the secondary stellar image in an exoplanetary atmosphere due to 
refraction as a function of the angle the planet sweeps out in orbit from mid-transit.
This therefore simulates the full structure of the shoulder feature. It only does this
for one side (ingress or egress).

Input
-----
Requires input file to specify the stellar-orbital-planetary parameters.

Output
------
Various data files and pickle files describing the refraction effect for this particular
system.

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
import matplotlib    #Plotting only occurs if plotting is turned 'on'
import matplotlib.pyplot
from matplotlib.pyplot import *  
from matplotlib import colors, cm

#Import other RETrO modules
import planet_mass_radius
import planet_atmosphere
from planet_atmosphere import *
import retro_rk4
#-----------------------------------------------------------------------------------------
	

#-----------------------------------------------------------------------------------------
#Set path information
path = './'
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


#Read the inputs that define this atmosphere
#-----------------------------------------------------------------------------------------
with open(path+'input.txt','r') as input_file:
	while 1:
		#Skip any header lines
		line = input_file.readline()
		if '#' in line: continue
		split_line = line.split('=')
		if 'ID' in line: ID = int(split_line[-1])
		if 'T_ATM' in line: T_atm = float(split_line[-1])
		if 'SEMI' in line: semimajor = float(split_line[-1])
		if 'M_P' in line: M_p = float(split_line[-1])
		if 'R_STAR' in line: R_star = float(split_line[-1])
		if 'ATM' in line: atm_type = split_line[-1]
		if not line: break
input_file.close()
#-----------------------------------------------------------------------------------------


#Calculate other necessary aspects of the atmosphere based on the input parameters.
#-----------------------------------------------------------------------------------------
#Convert units of input parameters
semimajor *= AU
M_p *= M_earth
R_star *= R_sun

#Get the planet radius in m
z_ref = planet_mass_radius.get_zref(M_p)	
	 	
#Calculate the gravitational acceleration
g = G*M_p/z_ref**2   #m/s^2

#Get mu and nu_ref from atmosphere module
mu, nu_ref = get_mu_nuref(atm_type)

#Atmospheric pressure scale height. 
H = k_B*T_atm/(mu*g)    #m
#Atmospheric 1 bar (1e5 Pa) reference number density	
nd_ref = 1e5/(k_B*T_atm)  #m^-3
#STP number density (from old STP when the refractivities were measured)
nd_STP = 101325./(k_B*273.15)
#Since both the number density and nu profiles share the same exponential factor, the 
# ratio of the densities will equal the ratio of the refractivities. 
nu_ref *= nd_ref/nd_STP
#-----------------------------------------------------------------------------------------


#Other parameters required to run the model
#-----------------------------------------------------------------------------------------
#Other orbital parameters, either held constant or arbitrary
ecc = 0.
theta_peri = 0.*(pi/180.)  #longitude of periastron, rad
n_theta0 = 4.   #Number of theta0 angles to cover
n_orbit_steps = 50   #Number of orbital steps to take

#Other atmosphere parameters, either held constant or arbitrary.
abs_cross_sec = 0.   #m^2, ignore absorption

#Ray tracing
z_top = z_ref+20.*H  #altitude where ray path integration begins
ds = 0.1*H   #step size for ray path integration
n_rotations = 100  #number of times the ray tracing plane is rotated (half crescent)

#These parameters controlled the initial ray separation and how it changes. Each ray is 
# some delta-y below the one above it, and if a ray experiences critical refraction, the
# delta-y is reduced by some fraction; this is one iteration. This will assure that low
# altitudes are more finely sampled. The delta-y listed here is the starting value. It 
# will be altered. This will continue until the rays have reached bending angles great
# enough for all the interpolation to occur.
delta_y_init = H/4.
delta_y_cut = 4.

#Turn plotting on or off here.
plotting = 'off'

#For this parameter space search, always assume spherical symmetry. Much of this code
# relies upon this assumption.
sphere_sym = 'yes'
#-----------------------------------------------------------------------------------------

#Make a useful save file that has all of the parameters for this run
#-----------------------------------------------------------------------------------------
save_file1 = open(path+'param_list_ID'+str(ID)+'.txt','w')
while 1:
	save_file1.write('#This file created at '+str(datetime.datetime.now())+\
		'. All units are SI.'+os.linesep)
	save_file1.write('ID\t\t='+str(ID)+os.linesep)
	save_file1.write('T_atm\t\t='+str(T_atm)+os.linesep)
	save_file1.write('semimajor\t='+str(semimajor)+os.linesep)
	save_file1.write('M_p\t\t='+str(M_p)+os.linesep)
	save_file1.write('R_star\t\t='+str(R_star)+os.linesep)
	save_file1.write('atm_type\t='+atm_type+os.linesep)
	save_file1.write('z_ref\t\t='+str(z_ref)+os.linesep)
	save_file1.write('g\t\t='+str(g)+os.linesep)
	save_file1.write('mu\t\t='+str(mu)+os.linesep)
	save_file1.write('nu_ref\t\t='+str(nu_ref)+os.linesep)
	save_file1.write('H\t\t='+str(H)+os.linesep)
	save_file1.write('nd_ref\t\t='+str(nd_ref)+os.linesep)
	save_file1.write('ecc\t\t='+str(ecc)+os.linesep)
	save_file1.write('theta_peri\t='+str(theta_peri)+os.linesep)
	save_file1.write('n_theta0\t='+str(n_theta0)+os.linesep)
	save_file1.write('n_orbit_steps\t='+str(n_orbit_steps)+os.linesep)
	save_file1.write('abs_cross_sec\t='+str(abs_cross_sec)+os.linesep)
	save_file1.write('z_top\t\t='+str(z_top)+os.linesep)
	save_file1.write('ds\t\t='+str(ds)+os.linesep)
	save_file1.write('n_rotations\t='+str(n_rotations)+os.linesep)
	save_file1.write('delta_y_init\t='+str(delta_y_init)+os.linesep)
	save_file1.write('delta_y_cut\t='+str(delta_y_cut)+os.linesep)
	save_file1.write('plotting\t='+plotting+os.linesep)
	save_file1.write('sphere_sym\t='+sphere_sym+os.linesep)
	break
save_file1.close()
#-----------------------------------------------------------------------------------------



#Define the function to be used when determining the rays the "bracket" the star (deter- 
# mine the rays that just graze the near and far edges of the star).
# Inputs are dx: the x-axis distance between the ray exit point and the
# origin point on the stellar surface; beta_x0: impact angle with y-plane at D.
#-----------------------------------------------------------------------------------------
def minimize_edges(dx, beta_x0, x_exit, y_exit):
	dy = dx/tan(beta_x0)
	y_final, x_final = y_exit - dy, x_exit - dx
	#Now find the distance from that point and the center of the star
	return sqrt((x_final-star_loc[0])**2 + (y_final-star_loc[1])**2)
#-----------------------------------------------------------------------------------------


#Orbit determination
#-----------------------------------------------------------------------------------------
#The orbit is measured relative to theta, which is zero at mid-transit. First determine
# the theta where the max effect occurs, so the projected separation of the star-planet
# centers (X) is z_top+R_star. This is an implicit equation that will require a solver.
theta_start = lambda thet: thet - arcsin((z_top+R_star)*(1.+ecc*cos(thet-theta_peri))/\
	(semimajor*(1.-ecc**2)))
theta0 = fsolve(theta_start,(z_top+R_star)/semimajor)[0]

#With this initial theta, define a theta array that is n-theta0 from mid-transit, and the
# number of steps chosen above
theta = linspace(theta0,n_theta0*theta0,n_orbit_steps)
#Calculate the instantaneous orbital distances at these theta values (physical units)
a = semimajor*(1.-ecc**2)/(1.+ecc*cos(theta-theta_peri))
#Calculate the projected center separations (in physical units)
X = a*sin(theta)
#Array to store the flux signals
signal = zeros_like(X)

#Create a few storage arrays that will store crescent-related data for all orbit steps
all_x_inner, all_x_outer, all_y_inner, all_y_outer, all_chi = zeros([n_orbit_steps,\
	n_rotations]),zeros([n_orbit_steps,n_rotations]),zeros([n_orbit_steps,n_rotations]),\
		zeros([n_orbit_steps,n_rotations]),zeros([n_orbit_steps,n_rotations])
#-----------------------------------------------------------------------------------------


#Start the loop over theta (orbital steps). 
#-----------------------------------------------------------------------------------------
# First define where the star is wrt the planet, which will always be at the origin
for l in range(n_orbit_steps):
	#Instantaneous x-y position of the star assuming origin at the center of the planet.
	star_loc = array([-sqrt(a[l]**2 - X[l]**2),-X[l]]) 
	
	#Set the the value of delta_y
	delta_y = delta_y_init
	
	#Set up the rotation of the ray tracing plane and begin a loop  over the angle chi. 
	# Also set up arrays to store x,y positions of the crescent. If the atmosphere is
	# spherically symmetric, rays only need to be traced once, since each rotation will 
	# shrink the effective radius of the host star. 
	#-------------------------------------------------------------------------------------
	chi = linspace(0,(1.-1e-10)*arcsin(R_star/X[l]),n_rotations)
	y_inner = zeros([n_rotations])    #Y position of edge of crescent closest to planet
	y_outer = zeros([n_rotations])    #Y position of edge of crescent farthest from planet
	x_inner = zeros([n_rotations])    #X position of edge of crescent closest to planet
	x_outer = zeros([n_rotations])    #X position of edge of crescent closest to planet
	r_star_t = zeros([n_rotations])   #Distance of tangent line dropped from star center
	R_star_eff = zeros([n_rotations]) #Eff. radius of star after rotation
	for k in range(n_rotations):
		#First, using chi, determine the stellar tangent radius
		r_star_t[k] = sin(chi[k])*X[l]
		#Now the effective stellar radius
		R_star_eff[k] = sqrt(R_star**2 - r_star_t[k]**2)
		#---------------------------------------------------------------------------------

		#---------------------------------------------------------------------------------
		#Now, initialize the ray arrays. These will be amended as more rays are traced.  
		# This initialization process only needs to happen once if spherical symmetry is 
		# assumed. The first ray will always start at z_top-delta_y
		if (sphere_sym =='yes') & (k==0) & (l==0):	
			#Ray IDs and array for final bending angle, minimum approach distance, exit 
			# angles. Also start the ray counter and set the critical refraction variable 
			# to 0 (False)
			bending_angle = array([0.])  #Start this at 0. to initiate the while loop
			z_min = array([])
			x_exit = array([])
			y_exit = array([])
			beta_x0 = array([])
			tau_final = array([])
			ray = array([0])
			crit_ref = 0
			
			#Also create a few storage lists that will save all of the low level ray
			# tracing details
			all_z, all_lat, all_beta, all_xi, all_tau = [],[],[],[],[]
			#-----------------------------------------------------------------------------

			#Begin the loop that will trace all rays.
			#-----------------------------------------------------------------------------
			#This is a while loop that forces ray tracing until the final accepted ray has
			# bent enough to allow all interpolations to occur. This critical bending is
			# approximated using the number of orbital steps, the size of the first step,
			# and the orbital parameter (e.g., simple bending angle is ~R_star/semi).
			bending_angle_thresh = (n_theta0+1.)*(2.*z_top+R_star)/semimajor
			n_iter = 0
			while abs(bending_angle[-1]) < bending_angle_thresh:
				#Redefine the bending angle array to remove the starting condition
				if n_iter == 0: bending_angle = array([])
				#Use a while loop to trace rays so long as critical refraction (or extreme 
				# bending beyond 90 deg) is not met.
				while crit_ref == 0:
					#Very first ray, start at the z_top minus delta_y
					if (ray[-1]==0) & (n_iter==0):
						y_initial = array([z_top - delta_y])       
						x_initial = sqrt(z_top**2 - y_initial**2)
						z_initial = ones_like(y_initial)*z_top
						lat_initial = arctan(y_initial/x_initial)
						beta_initial = lat_initial - pi/2.
						xi_initial = zeros_like(beta_initial)
						tau_initial = zeros_like(beta_initial)
					#All other rays, based the y_initial off the previous value, and then 
					# find the other initial values normally.
					else:
						y_initial = append(y_initial,y_initial[-1]-delta_y)	
						x_initial = append(x_initial,sqrt(z_top**2 - y_initial[-1]**2))
						z_initial = append(z_initial,z_top)
						lat_initial = append(lat_initial,arctan(y_initial[-1]/\
							x_initial[-1]))
						beta_initial = append(beta_initial, lat_initial[-1]-pi/2.)
						xi_initial = append(xi_initial,0.)
						tau_initial = append(tau_initial,0.)
							
					#Initialize the arrays for the integrated quantities, the constants,   
					# and the refractive invariant (which is just the impact param, or  
					# starting y-values).
					z, lat, beta, xi, tau = array([z_initial[-1]]), \
						array([lat_initial[-1]]),array([beta_initial[-1]]), \
							array([xi_initial[-1]]), array([tau_initial[-1]])
					constant = array([abs_cross_sec,nd_ref,z_ref,z_top,nu_ref,H]) 
	
					#Now begin the ray stepping and continue it until some termination 
					# condition is met.  
					while z[-1] <= z_top:
						#Create the current array
						current=array([z[-1],lat[-1],beta[-1],xi[-1],tau[-1],ds,ray[-1]])	
						#Call the RK4 function
						next = retro_rk4.rk4_ray_trace(current,constant)
						#Update the integrated quantity arrays
						z = append(z,next[0])
						lat = append(lat,next[1])
						beta = append(beta,next[2])
						xi = append(xi,next[3])
						tau = append(tau,next[4])
				
						#Monitor that the ray has not experienced critical refraction. If  
						# it has, break out of the ray stepping while loop. This condition 
						# assumes spherical symmetry. Also have the check to make sure the 
						# overall termination condition for the ray stepping while loop 
						# has not been met, otherwise the ray curvature will be zero since  
						# the ray is outside of the atmosphere.
						if z[-1] <= z_top:
							if abs(1./get_ray_curvature(z=z[-1],lat=lat[-1],beta=beta[-1],\
								dndz=get_drefractivity_dz(z=z[-1],lat=lat[-1],z_top=z_top,\
									nu_ref=nu_ref,z_ref=z_ref,H=H), dndlat = \
										get_drefractivity_dlat(z=z[-1],lat=lat[-1]),\
											z_top=z_top, nu_ref=nu_ref, z_ref=z_ref, \
												H=H)) < z[-1]: 
													crit_ref = 1
													break
						#Check to see if the bending has extended beyond 90 deg. In some
						# atmospheres, rays can not critically refract but bend >90 deg.
						# Such a ray will never reach the star and will not be considered.
						# I will assume this is essentially critical refraction, and 
						# terminate the ray in the same fashion as above.				
							if abs(xi[-1]) > pi/2.:
								crit_ref = 1
								break			

				
					#Only continue analyzing this ray if it did not experience critical 
					# refraction. 
					if crit_ref == 0:
						#Update the array keeping track of the bending angle, z_min, etc
						bending_angle = append(bending_angle,xi[-1])
						tau_final = append(tau_final,tau[-1])
						z_min = append(z_min,min(z))
						x_exit = append(x_exit,z[-1]*cos(lat[-1]))
						y_exit = append(y_exit,z[-1]*sin(lat[-1]))
						
						#Also update the storage arrays that are holding all the low-level
						# ray path integration info. For now, only keep a subset of all
						# the data, otherwise the data storage could become unreasonable.
						# The amount to keep should correlate to a step size of around a
						# scale height.
						all_z.append(z[::int(ceil(H/ds))])
						all_lat.append(lat[::int(ceil(H/ds))])
						all_beta.append(beta[::int(ceil(H/ds))])
						all_xi.append(xi[::int(ceil(H/ds))])
						all_tau.append(tau[::int(ceil(H/ds))])
						
						#Now for the projected geometry. It projects the ray back to the    
						# x=0 axis to find an alternate interior angle that is then the    
						# impact angle at the plane behind the Sun.
						beta_x0 = append(beta_x0,pi + beta[-1] - lat[-1])
	
						#The following block is used for plotting the path of the rays. 
						# The if statement will initialize all the plots for the first ray
						if (ray[-1]==0) & (n_iter==0) & (plotting=='on'): 
							#This plot for the path of the rays
							fig1=figure(1)
							circle = linspace(0,2.*pi,1000)  #radians
							ax1=fig1.add_subplot(111,aspect='equal')
							ax1.plot(z_top*sin(circle),z_top*cos(circle),c='k',lw=2)
							ax1.plot(z_ref*sin(circle),z_ref*cos(circle),c='k',ls='--')
							ax1.plot(R_star_eff[k]*sin(circle)+star_loc[0],R_star_eff[k]*\
								cos(circle)+star_loc[1],c='r',lw=2)
							ax1.axhline(0.,c='k',ls='--')
							ax1.axvline(0.,c='k',ls='--')
							xticks(fontsize='large')
							yticks(fontsize='large')
		
						#Update the plots for each ray
						if plotting=='on':
							ax1.plot(cos(lat)*z,sin(lat)*z)
						print 'Completed ray '+str(ray[-1])+', with xi= '+str(xi[-1])
					
						#Update the ray counter
						ray = append(ray,ray[-1]+1)	
						
			
				#If critical refraction was met, the while loop that keeps ray tracing
				# going will break and the code will come here.  
				# Alter the delta_y to be some fraction of what it was previously
				delta_y /= delta_y_cut
			
				#In some cases the delta_y becomes so small that it does not numerically
				# distinguish the next y_initial value from the previous one. If this 
				# happens, the ray tracing must end. Any crescent resulting from rays at
				# these altitudes will have widths of ~zero. 
				if (y_initial[-2]-delta_y) == y_initial[-2]:
					#Remove the most recent entries in the ray arrays. The most recent led
					# to critical refraction. The one before that represents the lowest
					# possible y_initial, and the maximally bent ray. Keep that one
					y_initial = delete(y_initial,-1)
					x_initial = delete(x_initial,-1)
					z_initial = delete(z_initial,-1)
					lat_initial = delete(lat_initial,-1)
					beta_initial = delete(beta_initial,-1)
					xi_initial = delete(xi_initial,-1)
					tau_initial = delete(tau_initial,-1)
					ray = delete(ray,-1)	
					#At this point, the ray tracing portion of the code needs to end. A
					# break here will terminate the ray tracing while loop. 
					print ''
					print 'Ending ray tracing early due to small delta_y limit issue'
					print ''
					break
			
				#Otherwise, remove the two most recent entries in the ray arrays. The most 
				# recent led to critical refraction, and the one before that usually bends 
				# far from the others. This prevents large gaps in sampling at the star 
				# plane.
				for cut in range(2):
					y_initial = delete(y_initial,-1)
					x_initial = delete(x_initial,-1)
					z_initial = delete(z_initial,-1)
					lat_initial = delete(lat_initial,-1)
					beta_initial = delete(beta_initial,-1)
					xi_initial = delete(xi_initial,-1)
					tau_initial = delete(tau_initial,-1)
					#If this is not the final iteration, then also remove the last 
					# successful ray to make the sampling more even
					if (cut==0) & (abs(bending_angle[-2])<bending_angle_thresh):
						bending_angle = delete(bending_angle,-1)
						z_min = delete(z_min,-1)
						x_exit = delete(x_exit,-1)
						y_exit = delete(y_exit,-1)
						beta_x0 = delete(beta_x0,-1)
						ray = delete(ray,-1)
						del all_z[-1]
						del all_lat[-1]
						del all_beta[-1]
						del all_xi[-1]
						del all_tau[-1]
					#If this is the final iteration, we can keep the last ray, and break
					# out of this for loop. The overall while loop tracing the rays will
					# end because its condition will be met. In this case, the small
					# delta y issue did not occur
					if abs(bending_angle[-2])>bending_angle_thresh: 
						ray = delete(ray,-1)
						break
				
				#Reset the critical refraction variable
				crit_ref = 0
				
				#print a message to mark the end of this iteration
				print ''
				print 'Delta_y iteration '+str(n_iter)+' has completed.'
				print ''	
	
				#Step up the iteration count
				n_iter += 1
			#At this point, all rays are traced given the desired number of iterations.
			n_rays = size(y_initial)
		#---------------------------------------------------------------------------------

		#House-keeping commands for the plot	
		#---------------------------------------------------------------------------------
		if plotting=='on':
			ax1.minorticks_on()
			ax1.set_xlabel('Distance from Origin [m]',fontsize='x-large')
			ax1.set_ylabel('Distance from Origin [m]',fontsize='x-large')
		#---------------------------------------------------------------------------------


		#Determine the "edge" of the crescent. 
		#---------------------------------------------------------------------------------
		#Determine the rays that are tangential to the star. These determine the 
		# boundaries of the secondary image of the star out-of-transit. This will only 
		# run if it is called. 
		stellar_min_z = zeros([n_rays]) #Distance b/t closest approach and star center
		x_final = zeros([n_rays])
		y_final = zeros([n_rays])
		#The x-range to be tested here are +/- a stellar radius. The point is to find  
		# the closest approach distance of the ray to the center of the star. For the   
		# two rays where this distance is R_star_eff, they bound the system.
		for i in range(n_rays):
			#Yplane is the plane of the star
			yplane_lims = array([star_loc[0]-R_star_eff[k],star_loc[0]+R_star_eff[k]])
			#dx_lims should be positive and in correct order for limits
			dx_lims = abs((yplane_lims - x_exit[i])[::-1]) 
			#Now run the optimize function
			min_res_edge = minimize_scalar(minimize_edges, bounds=(dx_lims[0],\
				dx_lims[1]),method='bounded', args=(beta_x0[i],x_exit[i],y_exit[i]))
			#Make sure the optimization worked
			if min_res_edge.message != 'Solution found.':
				print ''
				stop = input('Minimization unsuccessful - please intervene')
			else:
				x_final[i] = x_exit[i] - min_res_edge.x
				y_final[i] = y_exit[i] - min_res_edge.x/tan(beta_x0[i])	
				stellar_min_z[i] = sqrt((star_loc[0]-x_final[i])**2 + (star_loc[1]-\
					y_final[i])**2)

		#Include the minimum points on the plot of the rays. This is a busy plot, but it
		# show where all the rays end up at the plane of the host star. Keep this off 
		# unless diagnosing ray tracing issues.		
		#if plotting=='on':
		#	ax1.scatter(x_final,y_final,facecolor='none',edgecolor='r',lw=2,s=30)		
		
		#With the y_final values of all the rays, there needs to be a check that the star
		# central yplane is fully sampled. This  may very well be the case if the ray 
		# tracing ended due to the small delta_y problem. But that may only affect distant
		# orbit steps and not the close in ones. If this is the very first orbital step,
		# assume that the interpolation will work
		if l==0: fully_sampled = 'yes'
		 
		# Since the bending_angle_thresh was just an approximate threshold value of the
		# padding, this is a double check that the interpolation below won't break.
		# If the problem occurs on this step, then the crescent is assumed (justifiably
		# so) to have a width of zero. Do not use effective stellar radius, if the central
		# width is zero, all the rotations' widths should be zero. 
		if (y_final[-1] > (star_loc[1]-R_star)) | (fully_sampled == 'no'):
			#Skip the interpolation and manually set the y_inner and outer values
			y_initial_outer = y_initial[-1]
			y_initial_inner = y_initial[-1]
			if k==0: print 'Forcing the next crescent width to 0 due to the small '+\
				'delta_y limit'
			#If this is not already set, set it now to make sure the other edges of the 
			# crescent are not calculated. After this is set, no other rotations or 
			# orbital steps will switch it back, which is as expected.
			fully_sampled = 'no'
			
		#Stellar min z is a distance, so it is positive. However, now give it a sign 
		# whether or not the final ray point is above or below the y=star_loc[1] 
		# plane. In other words, this stellar min z is positive if the ray finished  
		# 'above' the equator of the star, and negative if it finished below. This way 
		# we can find the two bounding points of the star as +/- R_star_eff. Above the
		# equator is the outer ray, below is the inner ray. Only run this interp funciton
		# if the problem described just above did not occur.
		else:
			y_minz_interp = interp1d(stellar_min_z*sign(y_final-star_loc[1]),y_initial)
			y_initial_outer = y_minz_interp(R_star_eff[k])
			y_initial_inner = y_minz_interp(-R_star_eff[k])

		#Final x,y calculation
		#Based on the y inner and outer values found in the rotated plane, we can 
		# transform back to the unrotated plane and plot the points
		y_outer[k] = y_initial_outer*sin(chi[k])
		y_inner[k] = y_initial_inner*sin(chi[k])
		x_outer[k] = y_initial_outer*cos(chi[k])
		x_inner[k] = y_initial_inner*cos(chi[k])

		#Plotting
		#Create a new plot that gives the observers perspective
		if (k == 0) & (l==0) & (plotting=='on'):
			fig2 = figure(2)
			circle = linspace(0,2.*pi,1000)  #radians
			ax2=fig2.add_subplot(111,aspect='equal')
			ax2.plot(z_top*sin(circle),z_top*cos(circle),c='k',ls='--')
			ax2.plot(z_ref*sin(circle),z_ref*cos(circle),c='k')
			ax2.plot(R_star_eff[k]*sin(circle)+star_loc[1],R_star_eff[k]*cos(circle),\
				c='r',lw=2)
			xticks(fontsize='large')
			yticks(fontsize='large')	
			ax2.minorticks_on()
			ax2.set_xlabel('Distance from Origin [m]',fontsize='x-large')
			ax2.set_ylabel('Distance from Origin [m]',fontsize='x-large')
		
		#This marks the end of the kth ray-tracing plane rotation.
		#---------------------------------------------------------------------------------

	#Put the crescent data into its storage array
	#-------------------------------------------------------------------------------------
	all_x_inner[l,:], all_x_outer[l,:] = x_inner, x_outer
	all_y_inner[l,:], all_y_outer[l,:] = y_inner, y_outer
	all_chi[l,:] = chi
	#-------------------------------------------------------------------------------------

	#Plot the crescent by simply combining the edge points in a closed curve order
	#-------------------------------------------------------------------------------------
	x_crescent = append(append(x_inner,x_outer[::-1]),append(x_outer,x_inner[::-1]))
	y_crescent = append(append(y_inner,y_outer[::-1]),append(-y_outer,-y_inner[::-1]))

	#Add them to the plot, with the appropriate project separation
	if plotting=='on':
		ax2.plot(x_crescent+(X[l]-X[0]),y_crescent)
	#-------------------------------------------------------------------------------------


	#Determine the area of the crescent by integrating the equation for circle sector area
	# Do the integral for inner and outer, then subtract them to account for the switch in
	# integral limits. This integral will require the radius of the x, y points as a func-
	# tion of chi (rotation angle). Use interpolation to find this. Obviously, more rays  
	# will yield a more precise area.
	#-------------------------------------------------------------------------------------
	#First, find the radii of the inner and outer points
	radii_inner = sqrt(x_inner**2 + y_inner**2)
	radii_outer = sqrt(x_outer**2 + y_outer**2)

	#Interpolate these as functions of chi
	radii_inner_interp = interp1d(chi,radii_inner)
	radii_outer_interp = interp1d(chi,radii_outer)

	#Now compute the integral of a circle section using quad from scipy using throw-away 
	# lambda functions
	circle_sec_inner = lambda alpha: 0.5*radii_inner_interp(alpha)**2
	circle_sec_outer = lambda alpha: 0.5*radii_outer_interp(alpha)**2 

	#Final run the integration, subtract, and double (symmetry about x-axis) to get area
	outer_area = quad(circle_sec_outer,0.,max(chi))
	inner_area = quad(circle_sec_inner,0.,max(chi))
	crescent_area = 2.*(outer_area[0] - inner_area[0])
	signal[l] = crescent_area/(pi*R_star**2)
	print 'Ray tracing signal = '+str(signal[l]*1e6)+' ppm'
	print 'End of orbital step '+str(l+1)
	print ''
	#-------------------------------------------------------------------------------------

#Make a final plot of the shoulder
#-----------------------------------------------------------------------------------------
if plotting=='on':
	fig3 = figure(3)
	ax3 = fig3.add_subplot(111)
	ax3.plot(X/X[0],signal*1e6,c='k',lw=2)
	xticks(fontsize='large')
	yticks(fontsize='large')	
	ax3.minorticks_on()
	ax3.set_xlabel('Normalized X Separation',fontsize='x-large')
	ax3.set_ylabel('Relative Flux Increase [ppm]',fontsize='x-large')
	ax3.set_xlim(0.9,)
#-----------------------------------------------------------------------------------------


#Create save files in the form of pickles to save the information for this run
#-----------------------------------------------------------------------------------------
#Save quantities related to the ray tracing. 
with open(path+'ray_ID'+str(ID)+'.pickle','wb') as ray_pickle:
	pickle.dump((ray, x_initial, y_initial, z_initial, beta_initial, lat_initial,\
		bending_angle,z_min,x_exit,y_exit,beta_x0,tau_final),ray_pickle)
ray_pickle.close()	

#Save quantities related to the crescents
with open(path+'crescent_ID'+str(ID)+'.pickle','wb') as crescent_pickle:
	pickle.dump((theta,a,X,signal,all_x_inner,all_x_outer,all_y_inner,all_y_outer,\
		all_chi),crescent_pickle)
crescent_pickle.close()	

#Save quantities related to ray path integration
with open(path+'raypath_ID'+str(ID)+'.pickle','wb') as raypath_pickle:
	pickle.dump((all_z,all_lat,all_beta,all_xi,all_tau),raypath_pickle)
raypath_pickle.close()	
#-----------------------------------------------------------------------------------------



#End of script indicator
#-----------------------------------------------------------------------------------------	
print ''
print ''
print 'Ray tracing script for run ID '+str(ID)+' has completed.'
print ''
print ''
#-----------------------------------------------------------------------------------------
