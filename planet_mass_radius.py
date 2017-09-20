"""
Name
----
planet_mass_radius.py

Description
-----------
RETrO: Refraction in Exoplanet Transit Observations

This script returns the radius (z_ref) of a planet give its mass. The default is to rely
only the Chen & Kipping (2017) mass-radius relation, but the user can input any relation
they see fit.

Chen & Kipping (2017): http://adsabs.harvard.edu/abs/2017ApJ...834...17C

Input
-----
Planet mass in kg

Output
------
Planet radius in m

Author
------
P. A. Dalba --- Boston University
pdalba -at- bu -dot- edu
"""
#-----------------------------------------------------------------------------------------
#Import various math, science, and plotting packages.
import numpy
from numpy import *   #Note, numpy functions will not be explicitly called out
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


#Broken power law of R vs. M relation from Chen & Kipping (2017)
C = array([1.008,0.,0.])
S = array([0.279, 0.589, -0.044])
mass_trans = array([2.04,0.414*317.83,0.08*333e3])   #Earth masses

#Must determine the transition radii for the terran-neptunian planets
r_trans = zeros_like(mass_trans)
r_trans[0] = C[0]*mass_trans[0]**S[0]

#Solve for C for the neptunian power law relation using this transition point and then get
# the next transition point
C[1] = r_trans[0]/(mass_trans[0]**S[1])
r_trans[1] = C[1]*mass_trans[1]**S[1]
	
#Solve for C for the jovian power law relation using this transition point and then get
# the next transition point
C[2] = r_trans[1]/(mass_trans[1]**S[2])
r_trans[2] = C[2]*mass_trans[2]**S[2]	

#Determine the radius using the correct power law. z_ref in units of meters
def get_zref(M_p):
	M = M_p/M_earth
	if M < mass_trans[0]: return R_earth*C[0]*M**S[0]	
	if mass_trans[0] < M < mass_trans[1]: return R_earth*C[1]*M**S[1]
	if mass_trans[1] < M < mass_trans[2]: return R_earth*C[2]*M**S[2]	
	print 'Mass outside of deterministic relation. No z_ref specified.'
	return nan
	