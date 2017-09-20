# retro

Refraction in Exoplanet Transit Observations

RETrO is a tool to simulate atmospheric refraction phenomena in the atmospheres of transiting and non-transiting exoplanets. RETrO is written in Python v2.7 and described in an article that is accepted for publication by The Astrophysical Journal. 

ApJ link: (to be linked after publication by ApJ)
arXiv link: (to be linked on 9/21/2017)

If you make use of this tool in your research, please cite the above article.


## Getting Started

I recommend that you read the descriptions of the codes below. In the demo directory, I have an example input file (input.txt), which shoulder.py will look by name when it is run. The outputs from shoulder.py and shoulder_visualization.py are also stored in the demo directory. Just download all of the .py files, place input.txt in the same directory, and run shoulder.py. Then run shoulder_visualization.py. The plots and .pickle files that are created should match those that are stored in the demo directory.

If you would prefer to skip the demo check, follow these steps to get up and running quickly:

1. Create an input.txt file with the necessary parameters including unique run ID, atmospheric temperature in Kelvin, semi-major axis in AU, planet mass in Earth masses, stellar radius in solar radii, and atmosphere type (i.e., H2, H2O, N2, CO2). See input.txt in the demo directory for an example file that is properly formatted. Store input.txt in the same directory as all of the .py files, or else specify within shoulder.py where to look for the input.txt.

2. Run shoulder.py. It should create a full parameter list "param_list_ID##.py" and three .pickle files. The "##" is the ID set in the input.txt file.

3. Edit the loop in shoulder_visualization.py to include the ID or IDs you would like to visualize. 

4. Run shoulder_visualization.py to see plots of the data output by shoulder.py.

5. The most desired output exists in the "crescent_ID##.pickle" file. "Theta" is the orbital phase, and "signal" is the relative flux increase for egress (only) caused by atmospheric refraction assuming a clear atmosphere. 

6. Re-run again with other input.txt files, but be careful to specify new IDs. Right now, there is no safeguard against overwriting output files. There are also a number of other (less fundamental, less physical) parameters listed in shoulder.py. If you want very specialized runs, I recommend reading that code and amending those parameters as you see fit.


## Brief description of codes

### shoulder.py

shoulder is the primary tool that is used by retro. It takes in a few basic parameters that describe a particular exoplanet system and---using a few other functions that define the orbit, atmosphere, etc. (see below)---it produces the relative flux increase due to atmospheric refraction that an observer would see before or after a transit. It does this by looping through the planet's position in its orbit at a resolution set by the user. At each orbital position, the ray tracing plane is rotated in a loop in order to map out the entire image of the stellar mirage (also at a resolution set by the user). 

As it is written in its present form, it keeps track of optical depth, but does not incorporate it in the actual refracted light signal calculated at the end. So the result produce by shoulder.py is for a *clear* atmosphere. 

Backward ray tracing is used, meaning that rays are initialized on the observer side of the exoplanet atmosphere. It is difficult to know which rays will effectively bound the stellar mirage at the beginning of the ray tracing. Therefore, I use an automatic scheme that traces one ray at a time and varies their initial impact parameters as their total bending angles increase. As the bending increases, the rays are spaced more closely when initiated. Ray tracing continues until either the rays fully sample the stellar disk or the spacing of the rays decreases below the bit precision. This is the "small delta_y issue" that may occur for certain exoplanets. This just means that the width of the crescent has shrunk to essentially zero, so no flux increase will be observed for this exoplanet at this point in its orbit. For certain systems, this is the case at all points in the orbit. The sensitivity of the automatic ray tracing scheme can be adjusted by the user based on several of the parameters defined in shoulder.py including the step size.
 
shoulder.py outputs a parameter list for the exoplanet system of interest and a few .pickle files containing information about the refracted light crescent, light curve information, and ray path information. The ray path pickle file contains a portion of the ray path information of each ray that was traced. Here, "portion" means every ceil(H/ds) step in the ray path. Search for "ceil()" in shoulder.py to see how this is specified. You can choose to save less or more information, but be warned that certain atmospheres (i.e., large radii, small scale heights) can lead to very large pickle files.

### planet_mass_radius.py

This code contains the details to convert planet mass, which is directly included in the input.txt file, to planet radius. Currently, it is set up to use the deterministic portion of the Chen & Kipping (2017) mass-radius relations. You are free to include your own M-R relation. The get_mass() function is called in shoulder.py, where Mp is sent in kilograms and the return (z_ref) is expected in meters.

### planet_atmosphere.py

This code contains many functions that define the atmospheric properties of the exoplanet that is being investigated. The first function returns the mean molecular *mass* and the relative refractivity based on the four atmosphere types specified in Dalba (2017). The other functions define the structure of the atmosphere and are regularly called by retro_rk4.py. This can be changed by the user as desired. Note that, as of this version, all atmospheres are assumed to be isothermal, so get_temperature() is never actually called.

### retro_rk4.py

This code contains the fourth-order Runge-Kutta integration scheme that is at the core at RETrO. I urge caution in altering this code unless you are familiar with RK4 schemes.

The one alteration that may be required involves the optical depth. retro_rk4 expects a "sigma" to be an input, where sigma is the absorption cross section that will be multiplied by the number densities to get dtau/ds along the path. This is a highly simplified calculation of optical depth. The motivated user could expand upon this. However, as mentioned above, the optical depth is not factored into the relative flux increase for this version of shoulder.py.

### shoulder_visualization.py

shoulder_visualization reads in the output from shoulder.py and creates several useful figures displaying the data. It works based on the ID of the specific shoulder.py run(s) as specified in the input.txt file(s). You must go into shoulder_visualization.py and specify which IDs (in a list) you would like to visualize. 

###

Good luck, and enjoy RETrO!

Paul A. Dalba, PhD Candidate
Boston University, Department of Astronomy
https://blogs.bu.edu/pdalba
pdalba@bu.edu
