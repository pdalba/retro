"""
Name
----
shoulder_visualization.py

Description
-----------
RETrO: Refraction in Exoplanet Transit Observations

This script creates several plots to visualize the output of the shoulder.py. It can 
work on a single run or multiple runs (in a loop, where the figures close after each).

Input
-----
Requires the path to the particular parameter space run and the ID(s) of the runs to be
visualized.

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
import matplotlib.pyplot as plt
from matplotlib.pyplot import *  
from matplotlib import colors, cm
plt.ion()
close('all')
#-----------------------------------------------------------------------------------------


#Path info
#-----------------------------------------------------------------------------------------
#Provide the path to the directory with the pickle files
param_search_path = './'

#Also provide a list of IDs to run.
ID_list = [1]
#-----------------------------------------------------------------------------------------


#Begin visualization loop
#-----------------------------------------------------------------------------------------
print ''
print 'Beginning loop of '+str(size(ID_list))+' cases.'

for i in range(size(ID_list)):
	close('all')
	#Load in all pickled data from this run
	with open(param_search_path+'/ray_ID'+str(ID_list[i])+'.pickle','rb') as ray_pickle:
		ray, x_initial, y_initial, z_initial, beta_initial, lat_initial,\
			bending_angle,z_min,x_exit,y_exit,beta_x0,tau_final = \
				pickle.load(ray_pickle)
	ray_pickle.close()
	with open(param_search_path+'/crescent_ID'+str(ID_list[i])+'.pickle','rb') as \
		crescent_pickle:
			theta,a,X,signal,all_x_inner,all_x_outer,all_y_inner,all_y_outer,all_chi = \
				pickle.load(crescent_pickle)
	crescent_pickle.close()
	with open(param_search_path+'/raypath_ID'+str(ID_list[i])+'.pickle','rb') as \
		raypath_pickle:
			all_z,all_lat,all_beta,all_xi,all_tau = pickle.load(raypath_pickle)
	raypath_pickle.close()									

	#Also read the param list file
	with open(param_search_path+'/param_list_ID'+str(ID_list[i])+'.txt','r') as param_file:
		while 1:
			#Skip any header lines
			line = param_file.readline()
			if '#' in line: continue
			line_split = line.split('=')
			if 'z_top' in line: z_top = float(line_split[-1])
			if 'z_ref' in line: z_ref = float(line_split[-1])
			if 'R_star' in line: R_star = float(line_split[-1])
			if not line: break
	param_file.close()
	
	#Calculate the star loc at the moment of max effect (theta0)
	star_loc = array([-sqrt(a[0]**2 - X[0]**2),-X[0]])
	
	
	#Make a plot showing the ray paths in the planetary atmosphere
	fig1=figure(1)
	circle = linspace(0,2.*pi,1000)  #radians
	ax1=fig1.add_subplot(111,aspect='equal')
	ax1.plot(z_top*sin(circle),z_top*cos(circle),c='k',lw=2)
	ax1.plot(z_ref*sin(circle),z_ref*cos(circle),c='k',ls='--')
	ax1.axhline(0.,c='k',ls='--')
	ax1.axvline(0.,c='k',ls='--')
	xticks(fontsize='large')
	yticks(fontsize='large')
	for j in range(size(all_z)):
		ax1.plot(cos(all_lat[j])*all_z[j],sin(all_lat[j])*all_z[j])
	ax1.minorticks_on()
	ax1.set_xlabel('Distance from Origin [m]',fontsize='x-large')
	ax1.set_ylabel('Distance from Origin [m]',fontsize='x-large')	
	ax1.set_title('ID '+str(ID_list[i]),fontsize='x-large')
	#fig1.savefig('raypaths.png')
	
	#Plot the crescents from the observers perspective
	fig2 = figure(2)
	ax2=fig2.add_subplot(111,aspect='equal')
	ax2.plot(z_top*sin(circle),z_top*cos(circle),c='k',ls='--')
	ax2.plot(z_ref*sin(circle),z_ref*cos(circle),c='k')
	ax2.plot(R_star*sin(circle)+star_loc[1],R_star*cos(circle),c='r',lw=2)
	xticks(fontsize='large')
	yticks(fontsize='large')	
	ax2.minorticks_on()
	ax2.set_xlabel('Distance from Origin [m]',fontsize='x-large')
	ax2.set_ylabel('Distance from Origin [m]',fontsize='x-large')
	for k in range(shape(all_x_inner)[0]):
		x_crescent = append(append(all_x_inner[k,:],all_x_outer[k,:][::-1]),\
			append(all_x_outer[k,:],all_x_inner[k,:][::-1]))
		y_crescent = append(append(all_y_inner[k,:],all_y_outer[k,:][::-1]),\
			append(-all_y_outer[k,:],-all_y_inner[k,:][::-1]))
		ax2.plot(x_crescent+(X[k]-X[0]),y_crescent)	
	ax2.set_title('ID '+str(ID_list[i]),fontsize='x-large')
	#fig2.savefig('crescents.png')
	
	#Make a plot of the shoulder vs. X separation
	fig3 = figure(3)
	ax3 = fig3.add_subplot(111)
	ax3.plot(X/X[0],signal*1e6,c='k',lw=2)
	xticks(fontsize='large')
	yticks(fontsize='large')	
	ax3.minorticks_on()
	ax3.set_xlabel('Normalized X Separation',fontsize='x-large')
	ax3.set_ylabel('Relative Flux Increase [ppm]',fontsize='x-large')
	ax3.set_xlim(0.9,)
	ax3.set_title('ID '+str(ID_list[i]),fontsize='x-large')
	#fig3.savefig('signal_vs_X.png')
	
	#Make a plot of the shoulder vs. phase
	fig4 = figure(4)
	ax4 = fig4.add_subplot(111)
	ax4.plot(theta,signal*1e6,c='k',lw=2)
	xticks(fontsize='large')
	yticks(fontsize='large')	
	ax4.minorticks_on()
	ax4.set_xlabel('Orbital Phase',fontsize='x-large')
	ax4.set_ylabel('Relative Flux Increase [ppm]',fontsize='x-large')
	#ax4.set_xlim(0.9,)
	ax4.set_title('ID '+str(ID_list[i]),fontsize='x-large')
	#fig4.savefig('signal_vs_phase.png')
	
	print ''
	stop = raw_input('Showing ID '+str(ID_list[i])+'. Press enter to continue.')
#-----------------------------------------------------------------------------------------







