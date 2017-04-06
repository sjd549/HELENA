# HELENA
Hpem ELectronic ENgine Analysis: data plotting and analysis software for TECPLOT output files.

Files recognised: 

TECPLOT2D.PDT, kin.pdt, movie1.pdt, movie_icp.pdt


Requires the following software/modules to function:

python-pip, python-numpy, python-matplotlib, ffmpeg, findtools.




HELENA is designed for use in it's own seperate folder. When executed it will search for directories and sub-directories which contain HPEM output files and catagorize these into seprate 'simulations' based on which directory they came from. 
Each simulation will then be processed in turn with the output from requested diagnostics being saved in seperate directories within the simulation directory. Diagnostics which perform comparisons between simulations will be saved in the upper HELENA directory.

E.g. to compare a set of 5 simulations:

Copy the 5 directories containing the output files into your HELENA directory.
Edit the HELENA switchboard, requesting your diagnostics, parameters and image settings.
Run HELENA and allow the diagnostics to be completed.
Obtain your processed images and graphs from within the simulation directories.





Users can manipulate diagnostics from the switchboard. Multiple diagnostics may be run at the same time although users should be aware that plotting phase resolved or iterative data will take a long time to process depending on the number of variables or simulations requested.

IterVariables = ['AR+','E','PPOT','TE','TG-AVE']	#Requested Movie_icp (iteration) Variables.
PhaseVariables = ['AR+','E','PPOT','TE','TG-AVE']	#Requested Movie1 (phase) Variables.
electrodeloc = [0,0] 						#Centre Cell of powered electrode [R,Z]. (T,B,L,R)
phasecycles = 1								#Number of phase cycles to be plotted.

Variables = ['E','TE','P-POT','TG-AVE','PRESSURE'] #Requested TECPLOT2D variables.
MultiVar = ['AR+']						#Additional variables plotted ontop of [Variables] 
radialineouts = [] 						#Radial 1D-Profiles to be plotted (fixed (Z-mesh) height)
heightlineouts = [0]					#Axial 1D-Profiles to be plotted (fixed (R-mesh) radius)


savefig routines involve the saving of 2D or 1D images.
They take requested variables from the listings above and plot 2D images, and/or 1D images with the profile locations being dictated by radial and hight lineouts.

savefig_itermovie = False				#Produces 2D images, 1D graph and mp3 of simulation convergence.
savefig_plot2D = False					#Produces 2D images of requested variables.

savefig_radialines = False				#Produces 1D graphs of radial profiles for requested variables.
savefig_heightlines = False				#Produces 1D graphs of axial profiles for requested variables.
savefig_multiprofiles = False			#Produces 1D graphs of all variables with "multivar" plotted ontop
savefig_comparelineouts = True			#Produces 1D graph of each variable for all TECPLOT files found.

savefig_phaseresolvelines = False		#1D Phase Resolved Images with viewpoint along 'lineouts'.
savefig_phaseresolve2D = False			#2D Phase Resolved Images.
savefig_sheathdynamics = False			#Produces PROES style images with viewpoint along 'lineouts'.




Trend comparison reduces each simulation down to a single datapoint for requested variables and compares trends across a batch of simulations. print_"string" switches allow for more detailed information to be printed to terminal for each type of trend analysis.
Trends will be taken from the radial and height lineouts requested and produced plots will overlay trends from multiple lineouts by default.

savefig_trendcomparison = True			#Runs all trend analysis, reduced verbosity.
print_meshconvergence = False			#Plots all trend data against number of mesh cell points.
print_generaltrends = False				#Plots maximum value of each variable for each simulation.
print_KnudsenNumber = False				#Plots 2D knudsen number images and 1D knudsen trend. 
print_totalpower = False				#Calculates total power deposited for ions and electrons.
print_DCbias = False					#Calculates DC bias at given electrodeloc, allows for DBD,CCP.
print_thrust = False					#Calculates thrust and Isp at given axial location.





image_"string" options allow for modification of plotted images. They apply the changes to all images produced from all of the requested diagnostics. 

image_aspectratio = [10,10]					#[x,y] in inches
image_plotsymmetry = True					#Plots full image if mesh symmetry was used
image_contourplot = True					#Plots contour lines over 2D images.
image_normalize = False						#Normalizes values to each image seperately.
image_plotgrid = False						#Plots outline of grid vertices
image_logplot = False						#Logs data in images -- NB. legends not currently updated.
image_rotate = True							#Rotates 2D images 90 degrees anti-clockwise.








