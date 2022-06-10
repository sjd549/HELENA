#!/usr/bin/env python3

#################################
#		Point of Contact		#
#								#
#	    Dr. Scott J. Doyle      #
#	  University of Michigan    #
#	  Electrical Engineering    # 
#	 & Computer Science Dept.   #
#    1301 Beal Ave, Ann Arbor,  #
#        MI 48109-2122 USA      #
#	  Scott.Doyle@Physics.org   #
#                               #
#################################
#           'HELENA'            #
#Hpem ELectronic ENgine Analysis#
#################################


#====================================================================#
				 #PROGRAM FLAGS AND MODULE IMPORTS#
#====================================================================#

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--first", action="store_true", dest="install", default=False, help="Install prompt for required python modules")
(options, args) = parser.parse_args()

if 'True' in str(options): 
	import os, sys
	import os.path

	print( '')
	print( 'First time use requires installation of additional python modules')
	print( 'Please type your password when prompted to allow installation:')
	print( '')
	try:
		os.system('sudo apt-get install python-pip')
		os.system('sudo apt-get install python-matplotlib')
		os.system('sudo apt-get install python-numpy')
		os.system('sudo apt-get install python-scipy')
		os.system('sudo apt-get install ffmpeg')
		os.system('pip install tqdm')
	except:
		print( '')
		print( 'Error installing required packages')
		print( 'Please attempt manual installation')
		print( '')
	#endtry
	print( '')
	print( '')
#endif

#==============#

#Import core modules
import matplotlib.cm as cm
import numpy as np
import scipy as sp
import math as m
import subprocess
import os, sys
import os.path
import re

#Enforce matplotlib to avoid instancing undisplayed windows
#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib			#matplotlib.use('Agg')

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import colors as col
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *
#from Read_data_functions.py import overlay_GEC_geometry
 

#====================================================================#
				  		 #LOW LEVEL INPUTS#
#====================================================================#

#Various debug and streamlining options.
Magmesh = 1							#initmesh.exe magnification factor. (Obsolete - legacy)
QuickConverge = False				#Supresses 2D Convergence images in savefig_convergence
ffmpegMovies = False				#If False: Suppresses ffmpeg routines, saves RAM.
IDEBUG = False						#Produces debug outputs for most diagnostics.

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#Numerical Calculation Methods:
GlobSheathMethod = 'AbsDensity'			#Set Global Sheath Calculation Method.
#Choices: ('AbsDensity','IntDensity')
GlobThrustMethod = 'AxialMomentum'		#Set Global Thrust Calculation Method. 
#Choices:('ThermalVelocity','AxialMomentum')
GlobMeanCalculation = 'MeanFraction'	#Definition of 'mean' EDF value
#Choices: ('MeanEnergy','MeanFraction')

#Overrides or 'fudge factors' for diagnostics
DCbiasaxis = 'Auto'							#Force Direction Over Which DCBias is Calculated
#Choices:('Axial','Radial','Auto')
SheathIonSpecies = ['AR+']					#Force Sheath Ion Species (blank for auto)
#['AR+'] #['O+']

#Data Filtering and Smoothing Methods:
KineticFiltering = True						#Pre-fit kinetic data employing a SavGol filter
PlotKineticFiltering = False				#Plot Filtered Profiles, or employ only in trends.
Glob_SavWindow, Glob_SavPolyOrder = 25, 3	#Window > FeatureSize, Polyorder ~= Smoothness

#Apply azimuthal direction (phase) to relevant variables if true, else plot magnitude only
PlotAzimuthalDirection = True

#Minimum plotted EDF energy fraction, cuts x-axis at index where :: f(e) = f(e)*EDF_Threshold
#Note: IEDF/EEDF trends are only taken within range :: EDF_threshold < f(e) < 1.0 
EDF_Threshold = 0.01						# i.e. = 0.0 to plot all

#Define units for particular variables
# THIS NEEDS REWORKING:
# 1) CONVERT UNITS TO SI FOR MATHS ON INPUT		(Current variable function but renamed?)
# 2) APPLY PLOTTED UNITS SEPERATELY				(New variable function (plot variables?)
# 3) APPLY SI and CGS PRESETS
#		- PROBABLY REMOVE INDIVIDUAL UNIT SETTINGS (too messy)
PressureUnit = 'Torr'						#'Torr','mTorr','Pa'			#RM OUTDATED, TO REMOVE
Units = 'SI'								#'SI','CGS'

####################

#Commonly used variable sets.
Phys = ['P-POT','TE','EF-TOT','ERADIAL','ETHETA','EAXIAL','PHASER','PHASE','PHASEZ','EAMB-Z','EAMB-R','RHO','BR','BRS','BZ','BZS','BT','BTS','BRF','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','EFLUX-R','EFLUX-Z','JZ-NET','JR-NET','J-THETA','J-THETA-5','J-TH(MAG)','J-TH(PHA)','TG-AVE','PRESSURE','POW-RF','POW-RF-E','POW-ICP','EB-ESORC','COLF']
ASTRONCOILEF = \
['ERADIAL-2','ETHETA-2','EAXIAL-2','PHASE-2','ERADIAL-3','ETHETA-3','EAXIAL-3','PHASE-3', \
 'ERADIAL-4','ETHETA-4','EAXIAL-4','PHASE-4','ERADIAL-5','ETHETA-5','EAXIAL-5','PHASE-5', \
 'ERADIAL-6','ETHETA-6','EAXIAL-6','PHASE-6','ERADIAL-7','ETHETA-7','EAXIAL-7','PHASE-7', \
 'ERADIAL-8','ETHETA-8','EAXIAL-8','PHASE-8']
ASTRONCOILBF = \
['BT-2','BT-3','BT-4','BT-5','BT-6','BT-7','BT-8', \
 'BRF-2','BRF-3','BRF-4','BRF-5','BRF-6','BRF-7','BRF-8', \
 'PHASEBT-2','PHASEBT-3','PHASEBT-4','PHASEBT-5','PHASEBT-6','PHASEBT-7','PHASEBT-8',]


Ar = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR','AR+','AR2+','AR2*','E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','FZ-AR3S','FR-AR3S','FR-AR+','FZ-AR+','FZ-AR3S','FR-AR3S']
O2 = ['O3','O2','O2+','O','O+','O-','E','S-O3','S-O2+','S-O+','S-O-','SEB-O3','SEB-O+','SEB-O2+','SEB-O-','FR-O+','FZ-O+','FR-O-','FZ-O-']+['O3P3P','O***','S-O3P3P','S-O***','SEB-O3P3P','SEB-O***']
NF3 = ['F']

Ar_Phase = ['S-E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','SRCE-2437','TE','PPOT','FR-E','FZ-E']
O2_Phase = ['S-E','S-O+','S-O-','S-O2+','SEB-O+','SEB-O-','SEB-O2+','TE','PPOT','FR-E','FZ-E']+['S-O3P3P','SEB-O3P3P']

PRCCPAr_PCMC = ['AR^0.35','EB-0.35','ION-TOT0.35']
PRCCPO2_PCMC = ['O^0.35','EB-0.35','ION-TOT0.35']
ESCTAr_PCMC = ['TO BE COMPLETED']
TSHCAr_PCMC = ['AR^0.2S','EB-0.2S','ION-TOT0.2S','AR^0.2T','EB-0.2T','ION-TOT0.2T','AR^4.9U','EB-4.9U','ION-TOT4.9U','AR^9.8V','EB-9.8V','ION-TOT9.8V']					#RM-OUTDATED

EVgeny_PCMC = ['AR^0.1P','EB-0.1P','ION-TOT0.1P','AR^2.1Q','EB-2.1Q','ION-TOT2.1Q']
HYPI_PCMC = ['O^0.2P','EB-0.2P','ION-TOT0.2P','O^4.1Q','EB-4.1Q','ION-TOT4.1Q','O^4.15','EB-4.15','ION-TOT4.15']
HYPII_PCMC = ['O^2.8P','EB-2.8P','ION-TOT2.8P','O^3.5Q','EB-3.5Q','ION-TOT3.5Q']


#Archived variable sets
TSHCOI2019_PCMC = ['AR^0.2S','ION-TOT0.2S','AR^4.4T','ION-TOT4.4T','AR^8.9U','ION-TOT8.9U']
TSHCOI2020_PCMC = ['AR^0.2S','ION-TOT0.2S','AR^4.9T','ION-TOT4.9T','AR^9.8U','ION-TOT9.8U']
MSHC2017Mk0_PCMC = ['AR^0.5S','EB-0.5S','ION-TOT0.5S','AR^1.1B','EB-1.1B','ION-TOT1.1B']
SCCP2018Mk0_PCMC = ['AR^7.7J','ION-TOT7.7J','AR^5.1B','ION-TOT5.1B']
ESCT2018Mk0_PCMC = ['AR^0.3S','EB-0.3S','ION-TOT0.3S']
ParallelPlatePCMC = ['AR^2.67', 'ION-TOT2.67']				#RM-OUTDATED

####################

#Commonly Used Diagnostic Settings:
#### PRCCP ####
#electrodeloc =		[29,44] 					#Reverse [29,62]
#waveformlocs =		[[16,29],[16,44],[16,64],[0,29],[0,44],[0,64]]
#DoFwidth =			R;16,Z;21
#TrendLoc =			H[0,16];R[29,44,64,75]
#thrustloc =		75, 						#stdESCT=76, smlESCT=48/54
#sheathROI =		[34,72]
#sourcewidth =		[0.21]						
#Crop =				R[0.65];Z[1.0,4.0] 

#### PRICP ####	
#electrodeloc = 	[33,33]			#Coil V
#waveformlocs = 	[]
#DoFwidth = 		[]
#TrendLoc = 		H[0];R[36,50]
#thrustloc = 		[79]
#sheathROI = 		[]
#sourcewidth = 		[]
#Crop = 			R[1.0];Z[1.0,9.0]
#Plotmesh = 		'PRCCP'

#### ESCT ####
#electrodeloc = 	[0,5]
#waveformlocs = 	[]
#DoFwidth = 		R;??,Z;??
#TrendLoc =  		H[0];R[??]		#[??] For PROES
#thrustloc = 		[??]
#sheathROI = 		[]
#sourcewidth = 		[??]
#Crop = 			R[0];Z[0,5]

#### TSHC ####
#electrodeloc = 	[5,20]
#waveformlocs = 	[]
#DoFwidth = 		R;20,Z;6
#TrendLoc =  		H[1,22,44,60];R[35,40,44]	#[40] For PROES
#thrustloc = 		[60]
#sheathROI = 		[]
#sourcewidth = 		[5]
#Crop = 			R[0,15];Z[7,14]					#Mk3: R[0,15];Z[10,17.5]

#### SERPENT ####	
#electrodeloc = 	[33,33]			#Coil V
#waveformlocs = 	[]
#DoFwidth = 		[]
#TrendLoc = 		H[0];R[36,50]
#thrustloc = 		[79]
#sheathROI = 		[]
#sourcewidth = 		[]
#Crop = 			R[1.4];Z[1,9]
#Plotmesh = 		False

####################




#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested IEDF/NEDF Variables.
IEDFVariables = PRCCPAr_PCMC				#Requested iprofile_2d variables (no spaces)
NEDFVariables = []							#Requested nprofile_2d variables (no spaces)

#Requested movie1/movie_icp Variables.
IterVariables = ['E','S-E','PPOT','TE']				#Requested Movie_icp (iteration) Variables.		
PhaseVariables = Ar_Phase							#Requested Movie1 (phase) Variables. +['E','AR+']
electrodeloc = [29,44]	#PR[29,44]#THC[35,45]		#Cell location of powered electrode [R,Z].
waveformlocs = []									#Cell locations of additional waveforms [R,Z].

#Requested TECPLOT Variables and plotting locations.
Variables = Phys+Ar + ASTRONCOILEF+ASTRONCOILBF
multivar = []								 #Additional variables plotted ontop of [Variables]
radialprofiles = [10]	#PR[44]#THC[85]#[10] #Radial 1D-Profiles to be plotted (fixed Z-mesh) --
axialprofiles = [10]	#PR[0]#THC[20]		 #Axial 1D-Profiles to be plotted (fixed R-mesh) |
probeloc = [0,44]		#PR[0,44]#THC[20,85] #Cell location For Trend Analysis [R,Z], ([] = min/max)


#Various Diagnostic Settings.
phasecycles = 1.00						#Number of waveform phase cycles to be plotted. (float)
DoFwidth = 10 							#PROES Depth of Field (symmetric about image plane) (cells)
thrustloc = 45			#THC[10]		#Z-axis cell for thrust calculation  (cells)
sheathROI = [34,72]						#Sheath Region of Interest, (Start,End) [cells]				<<< OUTDATED
sourcewidth = [12]						#Source Dimension at ROI, leave empty for auto. [cells]		<<< OUTDATED


#Requested diagnostics and plotting routines.
savefig_convergence = False				#Single-Variables: iter-time axis			Requires movie_icp.pdt
savefig_plot2D = False					#Single-Variables: converged				Requires TECPLOT2D.PDT
#	NEED TO ADD ICOILP-n TOT OPTION WITH ALL COILSETS OVERLAYED

savefig_monoprofiles = False			#Single-Variables; fixed height/radius
savefig_multiprofiles = False			#Multi-Variables; same folder					- NO ASCII OUTPUT
savefig_compareprofiles = True			#Multi-Variables; all folders
savefig_temporalprofiles = False		#Single-Variables; real-time axis
					
savefig_trendphaseaveraged = False		#Converged trends at 'TrendLoc' cell
savefig_trendphaseresolved = False		#Temporal trends at 'TrendLoc' cell				#IN DEVELOPMENT#

savefig_phaseresolve1D = False			#1D Phase Resolved Images
savefig_phaseresolve2D = False			#2D Phase Resolved Images
savefig_sheathdynamics = False			#1D and 2D sheath dynamics images
savefig_PROES =	False					#Simulated PROES Diagnostic

savefig_IEDFangular = False				#2D images of angular IEDF; single folders
savefig_IEDFtrends = False				#1D IEDF trends; all folders
savefig_EEDF = False					#NO PLOTTING ROUTINE							#IN DEVELOPMENT#


#Write processed data to ASCII files.
write_ASCII = True						#All diagnostic output written to ASCII.


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs.
print_Knudsennumber = False				#Print cell averaged Knudsen Number
print_totalpower = False				#Print all requested total powers
print_Reynolds = False					#Print cell averaged sound speed
print_DCbias = False					#Print DC bias at electrodeloc
print_thrust = False					#Print neutral, ion and total thrust
print_sheath = False					#Print sheath width at electrodeloc						

#Image plotting options.
image_extension = '.png'				#Define image extension  ('.png', '.jpg', '.eps')		
image_interp = 'bilinear'				#Define image smoothing  ('none', 'bilinear')
image_cmap = 'plasma'					#Define global colourmap ('jet','plasma','inferno','gnuplot')
# ADD TECPLOT COLOUR SCHEME
# ADD SI/CGS UNIT SCHEMES

image_aspectratio = [10,10]				#Real Size of [X,Y] in cm [Doesn't Rotate - X is always horizontal]
image_radialcrop = []#[0.65]			#Crops 2D images to [R1,R2] in cm
image_axialcrop = []#[1.0,4.0]			#Crops 2D images to [Z1,Z2] in cm
image_cbarlimit = []					#[min,max] colourbar limits

image_plotsymmetry = False#True				#Plot radial symmetry - mirrors across the ISYM axis
image_plotcontours = False#True				#Plot contour lines in 2D images
image_plotoverlay = False				#Plot location(s) of 1D radial/axial profiles onto 2D images
image_plotsheath = False				#Plot sheath extent onto 2D images
image_plotgrid = False					#Plot major/minor gridlines on 1D profiles
image_plotmesh = False#'PRCCP'			#Plot material mesh outlines ('Auto','PRCCP','PRCCPM','ESCT','GEC')
image_numericaxis = False				#### NOT implemented ####

image_rotate = False#True					#Rotate image 90 degrees to the right.			# MAKE [0000-1000, 0-2pi]
image_normalise = False					#Normalise image/profiles to local max
image_logplot = False					#Plot ln(Data), against linear axis.

#Overrides the automatic image labelling.
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = ['Notimplemented']





#============================#




#####TODO#####

### V3.2.0 Release Version ###
#Namelist read-in functions to simplify error analysis. (deal with !!! in .nam)
#Variable Interpolator needs to work with phasedata - Take variables from batch?
#Correct any data 'direction' within readin functions, not within diagnostics.

#Confirm 2DPhaseMovie and PROES radial direction is consistent when rotating:
	#1DPhaseMovie needs no reversal using ExtractRadialProfile
	#2PhaseMovie needs no reversal using ImageExtractor2D
	#PROES requires a reversal using ExtractRadialProfile
#Refactor PROES into 3D array (R,Z,Phase) and perform slice/integration rather than calling plotradial?

#Fix issue with ffmpeg "convert-im6.q16: DistributedPixelCache" and address the ignorance of os.system.

#Remove/simplify the phase-averaged sheath diagnostic?
#Sheathwidth function integrates axially or radially depending on mesh geometry.
#Sheathwidth function has 1D (Scott) and 2D (Greg) capabilities
#SheathWidth function can deal with image rotations.

#IEDF diagnostic capable of comparing between different material surfaces in single image
#IEDF diagnostic saves different material surfaces in different folders
#Add EEDF section and functionalise.

#Thrust diagnostic split into functions performing the same task as before.
#Thrust diagnostic enforces image symmetry, correcting the half-thrust error.

#multivar_profiles diagnostic overhaul - update coding structure and include an ASCII output option

#Add 'probeloc' option to savefig_temporal, where the meshavg is the default but if there are any supplied probeloc cells then those get saved into their own seperate folder.

#Fix Rynolds diagnostic - Current version has issue with NaNs and incorrect averaging
#implement Rynolds number properly (was converted from sound speed diagnostic)
#implement magnetic Rynolds number as well if possible

#implement image_numericaxis, try float(FolderNameTrimmer) as axis.
#implement numerical image_rotate, allow for 000,090,180,270.

#implement tests for key I/O functions
#Implement tests for key mathematical/diagnostic functions (dc self-bias/thrust/soundspeed/knudsen/etc...)
#Implement tests for each diagnostic 
	#AND/OR include a tiny parallel plate folder in github with output figures for comparison



### For Future ###
#include 'garbage ./collection' at the end of each diagnostic.
#Clean up unused functions and ensure homogeneity.
#Split into modules: HELENAIO, HELENADIAGNOSTIC, HELENAPLOTTING
#Update README, include all diagnostics and examples.
#Create developer handbook describing functions.

#Include Andor, OceanOptics, and LeCroy readin functions.














#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Variables and lists for basic processing
Dir = list()					#List of all datafile directories inside folders
Dirlist = list()				#List of all folder directories (not including filenames)
IEDFVariablelist = list()		#List of all variable names in pcmc.prof in header order
Geometrylist = list()			#List containing commonly used geometries [LEGACY: NOT USED]

Globalvarlist = list()			#List of all commonly shared variable names between all folders
Globalnumvars = list()			#Number of commonly shared variables between all folders.

#Variables for mesh_size lists and SI conversion
ISYMlist = list()				#list of ISYM values in folder order in Dirlist
IXZlist = list()				#List of IXZ values in folder order in Dirlist
R_mesh = list()					#List of radial mesh cells for initmesh.out in folder order in Dirlist
Z_mesh = list()					#List of axial mesh cells for initmesh.out in folder order in Dirlist
Raxis = list()					#Radial SI [cm] axis for plotting
Zaxis = list()					#Axial SI [cm] axis for plotting

Depth = list()					#icp.nam Depth input [cm] in folder order in Dirlist
Radius = list()					#icp.nam Radius input [cm] in folder order in Dirlist
Height = list()					#icp.nam Height input [cm] in folder order in Dirlist
dr = list()						#Radial mesh resolution [cm/cell] in folder order in Dirlist
dz = list()						#Axial mesh resolution [cm/cell] in folder order in Dirlist
dy = list()						#Depth mesh resolution [cm/cell] in folder order in Dirlist

#Variables and lists for icp.nam parameters
VRFM,VRFM2   = list(),list()							# Array of reals (In Material Order)
FREQM,FREQM2 = list(),list()							# Array of reals (In Material Order)
FREQC        = list()									# Array of reals (In Material Order)
FREQMAX,FREQMIN  = list(),list()						# Array of reals (In Material Order)
FREQGLOB,FREQALL = list(),list()						# real
IRFPOW = list()											# real
PRESOUT = list()										# real
IMOVIE_FRAMES = list()									# real
NUMMETALS=0; CMETALS,IETRODEM = list(),list()			# int;	Array of strings and ints, respectively
NUMCOILS=0; CCOILS = list()								# int;	Array of strings
IMATSTATS=0; CMATSTATS = list()							# int;	Array of strings
IPCMCSPEC=0; CPCMCSPEC = list()							# int;	Array of strings
IEBINSPCMC=0; EMAXIPCMC=0								# int; 	int

#Lists for icp.dat variables
header_icpdat = list()			#[SpeciesName, Charge, MolecularWeight, StickingCoeff,	- Array of strings
								# Transport, ReturnFrac, ReturnName]
AtomicSpecies = list()			#All species contained within chemistry set				- Array of strings
FluidSpecies  = list() 			#All neutral fluid species (for fluid dynamics) 		- Array of strings
NeutSpecies	= list()			#All neutral and metastable species						- Array of strings
PosSpecies = list()				#All positive ion species								- Array of strings
NegSpecies = list()				#All negative ion species								- Array of strings

#Lists to store raw data
rawdata_2D = list()				#ASCII format TECPLOT2D.pdt data string list			- Variable,Radius,Axis
rawdata_kin = list()			#ASCII format kin.pdt data string list					- Variable,Radius,Axis
rawdata_phasemovie = list()		#ASCII format movie1.pdt data string list				- Variable,Radius,Axis
rawdata_itermovie = list()		#ASCII format movie_icp.pdt data string list 	  		- Variable,Radius,Axis
rawdata_IEDF = list()			#ASCII format iprofile_tec2d.pdt data string list 		- Variable,Radius,Axis
rawdata_mcs = list()			#ASCII format mcs.pdt data string list 			  		- Variable,Radius,Axis

Data = list()					#Data[folder][Variable][Data]						-Data = 2D (R,Z) of reals
DataKin = list()				#Data[folder][Variable][Data]						-Data = 1D (Avg) of reals
DataIEDF = list()				#Data[folder][Variable][Data]						-Data = 2D (R,Z) of reals
DataEEDF = list()				#Data[folder][Variable][Data]						-Data = 2D (R,Z) of Reals
IterMovieData = list()			#ITERMovieData[folder][Timestep][Variable][Data]	-Data = 2D (R,Z) of Reals
PhaseMovieData = list()			#PhaseMovieData[folder][Timestep][Variable][Data]	-Data = 2D (R,Z) of Reals

Moviephaselist = list()			#'CYCL = n'												-int
MovieIterlist = list()			#'ITER = n'												-int
EEDF_TDlist = list()			#'???'													-???

header_itermovie = list()		#Header num rows for movie_icp.pdt						-1D array of ints
header_phasemovie = list()		#Header num rows for movie1.pdt							-1D array of ints
header_IEDFlist = list()		#Header num rows for iprofile_tec2d.pdt					-1D array of ints
header_kinlist = list()			#Header num rows for kin.pdt							-1D array of ints
header_2Dlist = list()			#Header num rows for TECPLOT2D.pdt						-1D array of ints









#====================================================================#
					#WELCOME SPLASH AND INFORMATION#
#====================================================================#

print( '')
print( '--------------------------------------------------------------------')
print( '    __    __   _______  __       _______  __   __      ___          ')
print( '   |  |  |  | |   ____||  |     |   ____||  \ |  |    /   \         ')
print( '   |  |__|  | |  |__   |  |     |  |__   |   \|  |   /  ^  \        ')
print( '   |   __   | |   __|  |  |     |   __|  |  . `  |  /  /_\  \       ')
print( '   |  |  |  | |  |____ |  `----.|  |____ |  |\   | /  _____  \      ')
print( '   |__|  |__| |_______||_______||_______||__| \__|/__/     \__\     ')
print( '                                                              v3.1.4')
print( '--------------------------------------------------------------------')
print( '')
print( 'The following diagnostics were requested:')
print( '-----------------------------------------')
if savefig_plot2D == True:
	print('# 2D Steady-State Image Processing')
if savefig_convergence == True:
	print('# 2D Convergence Movie Processing')
if True in [savefig_phaseresolve2D,savefig_PROES]:
	print('# 2D Phase-Resolved Movie Processing')
if True in [savefig_phaseresolve1D]:
	print('# 1D Phase-Resolved Profile Processing')
if True in [savefig_monoprofiles,savefig_multiprofiles,savefig_compareprofiles,savefig_temporalprofiles]:
	print('# 1D Steady-State Profile Processing')
if True in [print_generaltrends,print_Knudsennumber,print_Reynolds, print_totalpower,print_DCbias,print_thrust]:
	print('# 1D Specific Trend Analysis')
if savefig_trendphaseaveraged == True:
	print('# 1D Steady-State Trend Processing')
if savefig_sheathdynamics == True:
	print('# 1D Phase-Resolved Trend Processing')
if True in [savefig_IEDFangular,savefig_IEDFtrends,savefig_EEDF]:
	print('# Angular Energy Distribution Processing')
print( '-----------------------------------------')
print( '')





#====================================================================#
					#OBTAINING FILE DIRECTORIES#
#====================================================================#

#Obtain system RAM. (and rename enviroment variable)
mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
mem_gib = mem_bytes/(1024.**3)
ext = image_extension

#Define recognized output file data extensions that will be retained in "Dir"
FileExtensions = ['.PDT','.pdt','.nam','.dat','.out']

#Create Directory lists and initialise numfolders to zero.
Dirlist = list() 		#List containing all simulation folder directories relative to HELENA
Dir = list() 			#List containing all output file in each Dirlist folder relative to HELENA
numfolders = 0			#Initiate folder number to zero

#Obtain home directory and contents
HomeDir = list() 		#List of all folders in home.
HomeDirContents = os.listdir( os.path.abspath(".") )
#Determine folders within home directory and add correct 'grammar'.
for i in range(0,len(HomeDirContents)):
	if os.path.isdir(HomeDirContents[i]) == True:
		HomeDir.append('./'+HomeDirContents[i]+'/')
	#endif
#endfor

#Determine number of folders containing accepted file extensions (i.e. simulation folders)
#Extract directories of each sub-folder within home directory
for i in range(0,len(HomeDir)):
	previousnumfolders = numfolders
	CurrentDir = HomeDir[i]
	DirContents = os.listdir(CurrentDir)

	#For each file contained within the subfolders, determine which are datafiles.
	for j in range(0,len(DirContents)):
		Filename = DirContents[j]

		#Save datafiles (with root) to working directory (Dir) and number of datafolders.
		if any([x in Filename for x in FileExtensions]):
			Dir.append(CurrentDir+Filename)
			if (numfolders - previousnumfolders) == 0:
				Dirlist.append(CurrentDir)
				numfolders += 1
			#endif
		else:
			File_Format_Is_Not_Requested = 1
		#endif
	#endfor
#endfor
#Maintain alphanumerical foldername structure (Dirlist) in-sync with dataname structure (Dir)
Dir,Dirlist = sorted(Dir),sorted(Dirlist)

#If no folders detected, end analysis script.
if numfolders == 0:
	print( '-------------------------------------------')
	print( 'No Ouput Files Detected, Aborting Analysis.')
	print( '-------------------------------------------')
	print( '')
	exit()
#endif

#Extract directories for all required data I/O files 
#These directories are relative to HELENA.py directory
icpnam = list(filter(lambda x: 'icp.nam' in x, Dir))
icpdat = list(filter(lambda x: 'icp.dat' in x, Dir))
icpout = list(filter(lambda x: 'icp.out' in x, Dir))
mesh = list(filter(lambda x: 'initmesh.out' in x, Dir))
TEC2D = list(filter(lambda x: 'TECPLOT2D.PDT' in x, Dir))

#Loop over all folders and retrieve mesh sizes and SI sizes.
for l in range(0,numfolders):
	
	#==========##===== INITMESH.OUT READIN =====##==========#
	#==========##===============================##==========#

	#Attempt automated retrieval of mesh sizes.
	try:
		#Identify mesh size from TECPLOT2D file. (Data Array Always Correct Size)
		meshdata = open(TEC2D[l]).readlines()

		#Zone line holds data, split at comma, R&Z values are given by "I=,J=" respectively.
#		R_py3_object = list(filter(lambda x: x.isdigit(), R))
#		Z_py3_object = list(filter(lambda x: x.isdigit(), Z))  
		R = list(filter(lambda x: 'ZONE' in x, meshdata))		#String: 'ZONE I=xxx, J=xxx, F=BLOCK"
		Z = list(filter(lambda x: 'ZONE' in x, meshdata))		#String: 'ZONE I=xxx, J=xxx, F=BLOCK"      
		R = R[0].split(",")[0].strip(' \t\n\r,=ZONE I')			#Split at commas, [0] gives "I=xxx"
		Z = Z[0].split(",")[1].strip(' \t\n\r,=ZONE J') 		#Split at commas, [1] gives "J=xxx"       
		R_mesh.append( int(R) )									#R_mesh (Cells) [int]
		Z_mesh.append( int(Z) )									#Z_mesh (Cells) [int]
	
	#If extraction from TECPLOT2D file fails, attempt to extract from initmesh.out header
	#This is an old method and causes issues with Q-VT meshes and magnified meshes
	except ValueError:
		#Identify mesh size from initmesh.out header:
		meshdata = open(mesh[l]).readline()
		R_mesh.append([int(i) for i in meshdata.split()][1])
		if Magmesh == 1: Z_mesh.append([int(i)+1 for i in meshdata.split()][3])
		elif Magmesh == 2: Z_mesh.append([int(i)+3 for i in meshdata.split()][3])
		elif Magmesh == 3: Z_mesh.append([int(i)+5 for i in meshdata.split()][3])
		#endif

	#If all else fails, request manual input of mesh resolution
	except:
		#If data for current file exists:
		if l <= len(TEC2D)-1:

			#If the initmesh.out file cannot be found, manual input is required.
			print( 'ERR: ICP.NAM GEOMETRY READIN, USING MANUAL MESH CELL INPUT:')
			print( Dirlist[l])
			r_mesh = int(input("DEFINE NUM RADIAL CELLS:"))
			z_mesh = int(input("DEFINE NUM AXIAL CELLS:"))
			print ('')

			R_mesh.append(r_mesh)
			Z_mesh.append(z_mesh)
		#endif
	#endtry

	#Retrieve entire mesh for plotting if requested.	#MESH PLOTTING NOT WORKING#
	if image_plotmesh == True:							#MESH PLOTTING NOT WORKING#
		image_plotmesh = False
		print( 'WARNING: AUTOMESH PLOTTING IS NOT CURRENTLY SUPPORTED')
		print( 'SETTING image_plotmesh = False')
		print( '')
		#Extract mesh data from initmesh.out			#MESH PLOTTING NOT WORKING#
#		mesh = open(mesh[l]).readlines()				#MESH PLOTTING NOT WORKING#
	#endif


	#==========##===== ICP.NAM READIN =====##==========#
	#==========##==========================##==========#

	#Attempt automated retrieval of SI conversion units.
	NamelistData = open(icpnam[l]).readlines()

	#Mesh Geometry Namelist Inputs
	try:
		RADIUS = list(filter(lambda x:'RADIUS=' in x, NamelistData))
		RADIUS = RADIUS[0].split('!!!')[0]
		RADIUS = float(RADIUS.strip(' \t\n\r,=RADIUS'))
		
		RADIUST = list(filter(lambda x:'RADIUST=' in x, NamelistData))
		RADIUST = RADIUST[0].split('!!!')[0]
		RADIUST = float(RADIUST.strip(' \t\n\r,=RADIUST'))

		HEIGHT = list(filter(lambda x:'HEIGHT=' in x, NamelistData))
		HEIGHT = HEIGHT[0].split('!!!')[0]
		HEIGHT = float(HEIGHT.strip(' \t\n\r,=HEIGHT'))
		
		HEIGHTT = list(filter(lambda x:'HEIGHTT=' in x, NamelistData))
		HEIGHTT = HEIGHTT[0].split('!!!')[0]
		HEIGHTT = float(HEIGHTT.strip(' \t\n\r,=HEIGHTT'))

		DEPTH = list(filter(lambda x:'DEPTH=' in x, NamelistData))
		DEPTH = DEPTH[0].split('!!!')[0]
		DEPTH = float(DEPTH.strip(' \t\n\r,=DEPTH'))

		IXZ = list(filter(lambda x:'IXZ=' in x, NamelistData))
		IXZ = IXZ[0].split('!!!')[0]
		IXZ = int(IXZ.strip(' \t\n\r,=IXZ'))
	
		ISYM = list(filter(lambda x:'ISYM=' in x, NamelistData))
		ISYM = ISYM[0].split('!!!')[0]
		ISYM = int(ISYM.strip(' \t\n\r,=ISYM'))
	
		#ISYMlist[l] = 1 if mesh uses radial symmetry, = 0 if not
		if image_plotsymmetry == True: ISYMlist= np.append(ISYMlist,ISYM)
		else: ISYMlist.append(0)
		
		#IXZlist[l] = 1 if mesh uses cartesian coordinates, = 0 if cylindrical
		IXZlist = np.append(IXZlist,IXZ)
		
		#Determine if mesh RADIUS or RADIUST was used, save 'Radius' used for further calculations
		if RADIUS > 0.0: Radius = np.append(Radius, RADIUS)			#[cm]
		elif RADIUST > 0.0: Radius = np.append(Radius, RADIUST)		#[cm]
		if HEIGHT > 0.0: Height = np.append(Height, HEIGHT)			#[cm]
		elif HEIGHTT > 0.0: Height = np.append(Height, HEIGHTT)		#[cm]

		#Determine mesh cell radial (dr), axial (dz), and depth (dy) resolutions 
		dr = np.append(dr, Radius[-1]/(R_mesh[-1]-1))		#[cm/cell]
		dz = np.append(dz, Height[-1]/(Z_mesh[-1]-1))		#[cm/cell]
		dy = np.append(dy, DEPTH)							#[cm/cell] (DEPTH is always 1 cell)
        
	except:
		#If the geometry section cannot be found, manual input is required.
		print('ERR: ICP.NAM GEOMETRY READIN, USING MANUAL MESH SI INPUT:')
		print(Dirlist[l])
		radius = float(input("DEFINE RADIUST [cm]:"))
		height = float(input("DEFINE HEIGHTT [cm]:"))
		depth = float(input("DEFINE DEPTH [cm]:"))
		print('')

		Radius.append(radius)
		Height.append(height)
		Depth.append(depth)
		dr.append(Radius[-1]/(R_mesh[-1]-1))
		dz.append(Height[-1]/(Z_mesh[-1]-1))
	#endtry

	#=====#=====#
	
	#Material Namelist Inputs (frequencies/voltages/powers)
	try:
		NUMMETALS = list(filter(lambda x:'IMETALS' in x,NamelistData))[0].strip(' \t\n\r,=IMETALS')
		NUMMETALS = int(NUMMETALS)+1
		CMETALS = list(filter(lambda x: 'CMETAL=' in x, NamelistData))[0].strip(' \t\n\r,').split()[1:NUMMETALS]
		for i in range(0,len(CMETALS)): CMETALS[i] = str(CMETALS[i].strip(','))
		IETRODEM.append( list(filter(lambda x:'IETRODEM=' in x, NamelistData))[0].split()[1:NUMMETALS])
		for i in range(0,len(IETRODEM[l])): IETRODEM[l][i] = int(IETRODEM[l][i].strip(','))
		##
		NUMCOILS = list(filter(lambda x:'ICOILS' in x,NamelistData))[0].strip(' \t\n\r,=ICOILS')
		NUMCOILS = int(NUMCOILS)
		CCOILS = list(filter(lambda x: 'CCOIL=' in x, NamelistData))[0].strip(' \t\n\r,').split()[1:NUMCOILS]
		for i in range(0,len(CCOILS)): CCOILS[i] = str(CCOILS[i].strip(','))
		##
		VRFM.append( list(filter(lambda x: 'VRFM=' in x, NamelistData))[0].split()[1:NUMMETALS] )
		for i in range(0,len(VRFM[-1])): VRFM[-1][i] = float(VRFM[-1][i].strip(','))
		VRFM2.append( list(filter(lambda x: 'VRFM_2=' in x, NamelistData))[0].split()[1:NUMMETALS] )
		for i in range(0,len(VRFM2[-1])): VRFM2[-1][i] = float(VRFM2[-1][i].strip(','))
		FREQM.append( list(filter(lambda x: 'FREQM=' in x, NamelistData))[0].split()[1:NUMMETALS] )
		for i in range(0,len(FREQM[-1])): FREQM[-1][i] = float(FREQM[-1][i].strip(','))
		FREQM2.append( list(filter(lambda x: 'FREQM_2=' in x, NamelistData))[0].split()[1:NUMMETALS] )
		for i in range(0,len(FREQM2[-1])): FREQM2[-1][i] = float(FREQM2[-1][i].strip(','))
		FREQGLOB.append( list(filter(lambda x:'FREQ=' in x, NamelistData))[0].split() )
		FREQGLOB[l] = float( FREQGLOB[l][0].strip(' \t\n\r,=FREQ') )	#NOTE: Assumes 1 entry for "FREQ="
		##
		FREQC.append( list(filter(lambda x: 'FREQC=' in x, NamelistData))[0].split(',')[0:-1] )
		for i in range(0,len(FREQC[-1])): FREQC[-1][i] = float(FREQC[-1][i].strip(' \t\n\r,=FREQC'))
		##
		IRFPOW.append( list(filter(lambda x:'IRFPOW=' in x, NamelistData))[0].strip(' \t\n\r,=IRFPOW'))
		IRFPOW[l] = float( IRFPOW[l].split('!!!')[0].strip(' \t\n\r,') )
		##
		PRESOUT.append( list(filter(lambda x:'PRESOUT=' in x, NamelistData))[0].strip(' \t\n\r,=PRESOUT'))
		PRESOUT[l] = float( PRESOUT[l].split('!!!')[0].strip(' \t\n\r,') )
	except:
		print( 'ERR: ICP.NAM MATERIAL DEFINITIONS READIN, USING DEFAULT MATERIAL PROPERTIES')
		FREQM.append(13.56E6)		#[Hz]
		FREQM2.append(13.56E6)		#[Hz]
		FREQC.append(13.56E6)		#[Hz]
		FREQGLOB.append(13.56E6)	#[Hz]
		VRFM.append(300.0)			#[Volts]
		VRFM2.append(150.0)			#[Volts]
		IRFPOW.append(100.0)		#[Watts]
		PRESOUT.append(0.85)		#[Torr]
	#endtry

	#No ICP coils are on - Ignore ICP frequencies in FREQALL
	if NUMCOILS == 0: 
		FREQALL = FREQM[l]+FREQM2[l]+FREQGLOB
		FREQALL = [x for x in FREQALL if x > 0]				#Ignore any DC voltages (FREQ = 0)
	#ICP coils are on - include ICP frequencies in FREQALL
	elif NUMCOILS > 0:
		FREQALL = FREQM[l]+FREQM2[l]+FREQC[l]+FREQGLOB
		FREQALL = [x for x in FREQALL if x > 0]				#Ignore any DC voltages (FREQ = 0)
	#endif
	FREQMIN.append( min(FREQALL) )							#Minimum global frequency 	[Hz]
	FREQMAX.append( max(FREQALL) )							#Maximum global frequency	[Hz]
	
	#Debug to check all RF frequencies
	if IDEBUG == 1:	
		print( 'Material Frequencies FREQM', FREQM[l] )
		print( 'Material Frequencies FREQM2', FREQM2[l] )
		print( 'Coil Freqencies FREQC', FREQC[l] )
		print( 'Global Frequency FREQ', FREQGLOB )
	#endif
		
	#=====#=====#

	#Plasma Chemistry Monte-Carlo (PCMC) Namelist Inputs
	try:
		IMATSTATS = list(filter(lambda x:'IMATSTATS' in x,NamelistData))[0].strip(' \t\n\r,=IMATSTATS')
		IMATSTATS = int( IMATSTATS[0].split('!!!')[0] )
		CMATSTATS = list(filter(lambda x: 'CMATSTATS=' in x, NamelistData))[0].strip(' \t\n\r,')
		CMATSTATS = CMATSTATS.split('=')[1].split(',')[0:IMATSTATS]
		##
		IPCMCSPEC = list(filter(lambda x:'IPCMCSPEC' in x,NamelistData))[0].strip(' \t\n\r,=IPCMCSPEC')
		IPCMCSPEC = int( IPCMCSPEC[0].split('!!!')[0])+1		#Add 1 to account for default ION-TOT
		CPCMCSPEC = list(filter(lambda x: 'CPCMCSPEC=' in x, NamelistData))[0].strip(' \t\n\r,')
		CPCMCSPEC = CPCMCSPEC.split('=')[1].split(',')[0:IPCMCSPEC]
		##
		IEBINSPCMC = list(filter(lambda x: 'IEBINSPCMC=' in x, NamelistData))
		IEBINSPCMC = IEBINSPCMC[0].split('!!!')[0]
		IEBINSPCMC = float(IEBINSPCMC.split()[0].strip(' \t\n\r,=IEBINSPCMC'))
		EMAXIPCMC = list(filter(lambda x: 'EMAXIPCMC=' in x, NamelistData))
		EMAXIPCMC = EMAXIPCMC[0].split('!!!')[0]
		EMAXIPCMC = float(EMAXIPCMC.split()[0].strip(' \t\n\r,=EMAXIPCMC'))
	except:
		print( 'ERR: ICP.NAM PCMC READIN, USING DEFAULT PCMC PROPERTIES')
		IEBINSPCMC = 1000
		EMAXIPCMC = 50
	#endtry

	#Phase-Resolved IMOVIE Namelist Inputs
	try:
		imovie_frames = list(filter(lambda x:'IMOVIE_FRAMES=' in x, NamelistData))	#OLD VERSION HAD [0]
		imovie_frames = imovie_frames[0].split('!!!')[0]
		IMOVIE_FRAMES.append( int(imovie_frames.strip(' \t\n\r,=IMOVIE_FRAMES')) )
	except:
		print( 'ERR: ICP.NAM IMOVIE READIN, USING DEFAULT PHASE RESOLVED PROPERTIES')
		IMOVIE_FRAMES.append(180)
	#endtry


	#==========##===== ICP.DAT READIN =====##==========#
	#==========##==========================##==========#

	#Attempt automated retrieval of atomic species
	ChemistryData = open(icpdat[l]).readlines()

	#Plasma chemistry .dat file inputs
	try:
		#Determine end of chemistry set species definition
		for i in range(0,len(ChemistryData)):

			#Atomic Species Defined In Header, read in data line by line from icp.dat
			#len(Header.split()) = 13 for atomic or molecular species definition
			#len(Header.split()) = 8 for material surface interaction definition 
			if len(ChemistryData[i].split()) == 13:
				SpeciesName     = ChemistryData[i].split()[0]
				Charge          = int(ChemistryData[i].split()[2])
				MolecularWeight = float(ChemistryData[i].split()[4])
				StickingCoeff   = float(ChemistryData[i].split()[6])
				TransportBool   = int(ChemistryData[i].split()[8])
				ReturnFrac      = float(ChemistryData[i].split()[10])
				ReturnSpecies   = ChemistryData[i].split()[11]

				#Collect all atomic species (including electrons)
				if SpeciesName not in AtomicSpecies: AtomicSpecies.append(SpeciesName)
				#Seperate species by charge
				if Charge == 0 and SpeciesName not in NeutSpecies: NeutSpecies.append(SpeciesName)
				elif Charge >= 1 and SpeciesName not in PosSpecies:  PosSpecies.append(SpeciesName)
				elif Charge <= -1 and SpeciesName not in NegSpecies: NegSpecies.append(SpeciesName)
				#List of recognized ground-state neutral species for fluid analysis.
				FluidSpecies = ['AR','AR3S','O2','O']	#NOTE :: FLUID SPECIES ARE STILL MANUALLY DEFINED

				#Collect icp.dat header if required for later use
				header_icpdat.append([SpeciesName,Charge,MolecularWeight,StickingCoeff, TransportBool,ReturnFrac,ReturnSpecies])
			#####

			#End of Chemistry Header Denoted By '*', as soon as this is reached, stop reading in.
			elif len(ChemistryData[i].split()) != 13 and len(ChemistryData[i].split()) !=8:
				if ChemistryData[i].split()[0] == '*': break
			#endif
		#endfor
	except:
		print( 'ERR: ICP.COM ATOMIC SPECIES READIN, USING DEFAULT ATOMIC PROPERTIES')
		#List of dafault recognised neutral/metastable atomic sets, add new sets as required.
		ArgonReduced = ['AR','AR+','AR*']
		ArgonFull = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR+','AR2+','AR2*']
		Oxygen = ['O','O+','O-','O*','O2','O2+','O2*']

		AtomicSpecies = ['E']+ArgonReduced+ArgonFull+Oxygen
		NeutSpecies = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR2*','O','O*','O2','O2*']
		PosSpecies = ['AR+','AR2+','O+','O2+']
		NegSpecies = ['E','O-']
		#List of recognized ground-state neutral species for fluid analysis.
		FluidSpecies = ['AR','AR3S','O2','O']
	#endtry 

#==========##========================##==========#
#==========##========================##==========#




#====================================================================#
				#UNPACKING AND ORGANIZATION OF DATA#
#====================================================================#

#VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
#Enumerates requested variables and produces a processlist for plotting.
def VariableEnumerator(Variables,Rawdata,Header):
	processlist = list()
	variablelist = list()

	#For all requested variables, in the requested data header, find which match.
	for j in range(0,len(Variables)):
		for i in range(0,Header):

			#Compare variables and if they match, add to the process list.
			#Default uses [1:-3] slice for the variable string.
			if Variables[j] == Rawdata[i].strip(' ,"\n').replace(' ',''):
				processlist.append(i)
				variablelist.append(Rawdata[i].strip(' ,"\n').replace(' ',''))
				break
			#endif
		#endfor
	#endfor
	return(processlist,variablelist)
#enddef


#Identifies if variable exists in all simulations, rejects if not.
#Allows for the comparison of datasets with different icp.dat files.
#Takes processlist, variablelist, globalcomparisonlist
#Returns processlist and variablelist with largest commonly shared variables.
#proclist,varlist = VariableInterpolator(processlist,Variablelist,Comparisonlist):
def VariableInterpolator(processlist,Variablelist,Comparisonlist):

	#No interpolation needed if variable count is the same for all datasets.
#	if all(map(lambda x: x == Globalnumvars[0], Globalnumvars)) == True:		#Py2.x.x Lambda Method
	if Globalnumvars.count(Globalnumvars[0]) == len(Globalnumvars):				#Py3.x.x Count Method
		return(processlist, Variablelist)
	#endif

	#Identify elements in each variablelist which are not in comparison list.
	interpolation = list()
	for i in range(0,numfolders):
		variablelist = VariableEnumerator(Variables,rawdata_2D[i],header_2Dlist[i])[1]
		inter = set(Comparisonlist).symmetric_difference(variablelist)
		inter = list(inter)
		#Collect list of all variables not present in all folders.
		for i in range(0,len(inter)):
			if inter[i] not in interpolation:
				interpolation.append(inter[i])
			#endif
		#endfor
	#endfor

	#If at least one element not present in all folders, remove from all Proc and Var lists.
	if len(interpolation) != 0:
		for i in range(0,len(interpolation)):
			j = 0
			while j < len(Variablelist):
				#Check for exact string match, e.g to avoid "AR in AR+".
				if interpolation[i] == Variablelist[j]:
					del Variablelist[j]
					del processlist[j]
				else:
					j += 1
				#endif
			#endwhile
		#endfor
	#endif

	return(processlist, Variablelist)
#enddef


#Takes full directory list (Dir) and data filename type (e.g. .png, .txt)
#Returns row-wise list of data and length of datafile.
#rawdata, datalength = ExtractRawData(Dir,'.dat',l)
def ExtractRawData(Dir,NameString,ListIndex=l):
	try:
		DataFileDir = filter(lambda x: NameString in x, Dir)
		DataFileDir = sorted(DataFileDir)
		Rawdata = open(DataFileDir[ListIndex]).readlines()
		nn_data = len(Rawdata)
	except:
		print( 'Unable to extract '+str(NameString))
		exit()
	#endtry
	return(Rawdata,nn_data)
#enddef


#Takes ASCII data in 2/3D format and converts to HELENA friendly structure.
#Requires rawdata(2D/3D), header and variable number and mesh dimensions.
#Allows for an optional offset in 'starting' variable number.
#Returns 2D array of form [Variables,datapoint(R,Z)]
#CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_2D[l],header_2D,numvariables_2D)
def SDFileFormatConvertorHPEM(Rawdata,header,numvariables,offset=0,Zmesh=0,Rmesh=0,Dimension='2D'):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 and Zmesh == 0:
		Rmesh,Zmesh = R_mesh[l],Z_mesh[l]
	#endif

	#Excluding the header, split each row of data and append items to 1D list.
	CurrentFolderData, DataArray1D = list(),list()
	for i in range(header,len(Rawdata)):

		#If end of phasecycle reached, break. (Applicable to 3D Datafiles only)
		if 'CYCL=' in Rawdata[i] or 'ITER=' in Rawdata[i]: 
			break
		else: CurrentRow = Rawdata[i].split()

		#For all elements in the current row, convert to float and save in list.
		for j in range(0,len(CurrentRow)):
			try: DataArray1D.append(float(CurrentRow[j]))
			except: String_Conversion_Error = 1
		#endfor
	#endfor

	#If data is 1D, seperate into 1D chunks using Zmesh as the chunk size
	if Dimension == '1D': 
		#Seperate total 1D array into further 1D sub-arrays with data for each variable.
		for i in range(offset,numvariables):
			numstart = Zmesh*(i)
			numend = Zmesh*(i+1)
			CurrentFolderData.append(list(DataArray1D[numstart:numend]))
		#endfor
		return(CurrentFolderData)

	#If data is 2D, return array in 2D chunks using Zmesh*Rmesh as the chunk size
	elif Dimension == '2D':
		#Seperate total 1D array into 2D array with data for each variable.
		#Offset data by a certain number of variable 'chunks' if requested.
		for i in range(offset,numvariables):
			numstart = int( (Zmesh*Rmesh)*(i) )
			numend = int( (Zmesh*Rmesh)*(i+1) )
			CurrentFolderData.append(list(DataArray1D[numstart:numend]))
		#endfor
		return(CurrentFolderData)
	#endif
#enddef


#Extracts all phase and variable data for the provided folder ID.
#Initial R and Z for CYCL=1 are skipped over and not saved.
#Takes current folder, returns Data[phase][variable][datapoints,R/Z]
#Data,Phaselist = ExtractPhaseData(folder=l,Variables=PhaseVariables)
def ExtractPhaseData(folder=l,Variables=PhaseVariables):
	#Load data from movie_icp file and unpack into 1D array.
	rawdata,filelength = ExtractRawData(Dir,'movie1.pdt',folder)

	#Read through all variables for each file and stop when list ends. 
	#Movie1 has geometry at top, therefore len(header) != len(variables).
	#Only the first encountered geometry is used to define variable zone.
	VariableEndMarker,HeaderEndMarker, = 'GEOMETRY','ZONE'
	variablelist,numvar = list(),0
	for i in range(2,filelength):
		if HeaderEndMarker in str(rawdata[i]):
			header = i+2	#plus 2 to skip to first data line.
			break
		if VariableEndMarker in str(rawdata[i]) and numvar == 0:
			numvar = (i-1-2)	#minus 1 for overshoot, minus 2 for starting at 2.
		if len(rawdata[i]) > 1 and numvar == 0: 
			variablelist.append(str(rawdata[i][:-2].strip(' \t\n\r\"')))
		#endif
	#endfor

	#Enumerate processlist and variablelist, interpolate variables against globalvarlist.
	proclist,varlist = VariableEnumerator(Variables,rawdata,header)
	for i in range(0,len(proclist)): proclist[i] -= 2	#R&Z not included, shift back by two.
	proclist,varlist = VariableInterpolator(proclist,varlist,Comparisonlist)

	#Rough method of obtaining the movie1.pdt cycle locations for data extraction.
	cycleloc = list()
	for i in range(0,len(rawdata)):
		if "CYCL=" in rawdata[i]:
			cycleloc.append(i+1)
		#endif
	#endfor

	#Cycle through all phases for current datafile, appending per cycle.
	#Variables R and Z only saved for first iteration, they are skipped if i == 0.
	FolderData,Phaselist = list(),list()
	for i in range(0,len(cycleloc)-1):	#### -1 IS A HACK TO ALIGN WITH OLD DATA ####
		if i == 0:
			PhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar+2,offset=2)
			FolderData.append(PhaseData[0:numvar])
		else:
			PhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar)
			FolderData.append(PhaseData)
		#endif
		Phaselist.append('CYCL = '+str(i+1))
	#endfor

	return(FolderData,Phaselist,proclist,varlist)
#enddef


#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Three imputs, Data to be written, Filename, 'w'rite or 'a'ppend.
#WriteDataToFile(Image, FolderNameTrimmer(Dirlist[l])+Variablelist[k])
def WriteDataToFile(data,filename,structure='w'):

	#Determine dimensionality of profile.
	if isinstance(data[0], (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, structure)
		for m in range(0,len(data)):
			for n in range(0,len(data[m])):
				datafile.write(str(data[m][n]))
				datafile.write(' ')
			#endfor
			datafile.write('\n')
		#endfor
		datafile.close()

	#Lowest dimention is scalar: ==> 1D array.
	elif isinstance(data, (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, structure)
		for n in range(0,len(data)):
			datafile.write(str(data[n]))
			datafile.write(' ')
		#endfor
		datafile.close()

	return()
#enddef


#Reads 1D or 2D data from textfile in ASCII format, returns data and header.
#Input filename, header length, data dimension and orientation (CSV or RSV).
#Example: OutputData,Header = ReadDataFromFile('/Data.txt', 1, '2D', CSV)
def ReadDataFromFile(Filename,HeaderIdx=0,Dimension='2D',Orientation='CSV'):
	OutputData,Header = list(),list()

	#If data is saved 'Row-wise', use default readin routine.
	if Orientation == 'RSV':
		#Determine dimensionality of profile.
		if Dimension in ['1D','2D']:
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()

			#Extract header and raw data
			for m in range(0,HeaderIdx): Header.append(RawData[m])
			RawData = RawData[HeaderIdx::]

			#Read each row, split it (space delimited) and save.
			for m in range(HeaderIdx,len(RawData)):
				Row = RawData[m].split()
				for n in range(0,len(Row)):
					try: Row[n] = float(Row[n])
					except: Row[n] = str(Row[n])
				#endfor
				OutputData.append(Row)
			#endfor
		#endif

	#=====#

	#If data is saved 'column-wise', transpose the arrays to correct.
	elif Orientation == 'CSV':
		#Determine dimensionality of profile.
		if Dimension in ['1D','2D']:
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()

			#Extract header and raw data
			for m in range(0,HeaderIdx): Header.append(RawData[m])
			RawData = RawData[HeaderIdx::]

			#Enlarge output data array by number of columns
			NumColumns = len(RawData[HeaderIdx+1].split())
			for m in range(0,NumColumns):
				OutputData.append(list())
			#endfor

			#Read each row, split it and save into relevant column of output data.
			for i in range(HeaderIdx,len(RawData)):
				Row = RawData[i].split()
				for j in range(0,len(Row)):
					try: Row[j] = float(Row[j])
					except: Row[j] = str(Row[j])
				#endfor
				for k in range(0,NumColumns):
					OutputData[k].append(Row[k])
				#endfor
			#endfor
		#endif
	#endif

	#Orientation doesn't matter if 0D (scalar data).
	elif Dimension == '0D':
		#Read in 0D data from ASCII formatted file.
		datafile = open(Filename)

		for m in range(0,HeaderIdx): Header.append(RawData[m])
		RawData = datafile.readlines()[HeaderIdx::]
		Row = RawData.split()

		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData,Header)
#enddef


#Creates a new folder if one does not already exist.
#Takes destination dir and namestring, returns new directory.
def CreateNewFolder(Dir,DirString):
	try:
		NewFolderDir = Dir+DirString+'/'
		#os.mkdir(NewFolderDir, 0755);os.mkdir(NewFolderDir, 755);
		os.mkdir(NewFolderDir, 0o755);        
	except:
		a = 1
	#endtry
	return(NewFolderDir)
#enddef


#Takes folder names and returns item after requested underscore index.
#Note, index > 1 will return between two underscores, not the entire string.
def FolderNameTrimmer(DirString,Index=1):
	try:
		for i in range(0,Index):
			underscoreloc = str(DirString[::-1]).index('_')
			cutoff = (len(DirString)-underscoreloc)
			NameString = DirString[cutoff:-1]
			DirString = DirString[:cutoff-1]
		#endfor
	except:
		NameString = str(DirString[2:-1])
	#endtry

	return(NameString)
#enddef


#Takes folder directory and creates a movie from .png images contained within.
def MakeMovie(FolderDir,Output):

	#Break if movies not requested
	if ffmpegMovies == False: return()

	#Obtain current directory and set movie output parameters	
	HomeDir = os.getcwd()
	Output = str(Output)+'.mp4'
	Morph, FPS = 1, 24

	#Ensure folder directory is in BASH-friendly format
	FolderDir = HomeDir+FolderDir[1::]

	#Use ffmpeg to create the movies and save in relevent files.
	os.system("cd "+FolderDir)
	os.system("convert *.png -delay 1 -morph "+str(Morph)+" %05d.morph.jpg > /dev/null")
	os.system("ffmpeg -nostats -loglevel 0 -r "+str(FPS)+" -i %05d.morph.jpg "+Output)
	os.system("rm *.jpg")
	os.system("cd "+HomeDir)

	return()
#enddef


#Runs requested dataconversion script with pre-defined arguments.
#Takes name of convert script, any predefined arguments and newly created files.
#Returns nothing, runs script in each folder, expects to run over all folders.
# *** TECPLOT FILE IS "NPROFILE_TEC2D.PDT" FOR NEUTRAL PROFILE DATA. ***
# *** TECPLOT FILE IS "IPROFILE_TEC2D.PDT FOR ION PROFILE DATA. ***
def AutoConvProf(ConvProfexe,args=[],DirAdditions=[]):
	HomeDir = os.getcwd()
	os.chdir(Dirlist[l])

	#ENTER FILE NAME FOR RAW PLOTTING DATA (simulation pcmc output file, default pcmc.prof)
	pcmcprof = 'pcmc.prof'
	#ENTER TITLE FOR DATA TO BE USED BY TECPLOT (TITLE MUST BE <= 20 CHARACTERS)
	TECPlotTitle = 'title'
	#ENTER 1 TO "APPROVE" WRITING VARIABLES 0 NOT TO "APPROVE"
	WriteVars = '1'
	#Should angular statistics be normalized to f(theta)/solid-angle? (default is f(theta) x d(omega))
	NormaliseAngularStatistics = '1'
	#Should profiles be averaged across 0 degrees?
	AverageAcross0Degrees = '1'
	#[pcmc.prof, title, write, normalise, average] are basic requirements
	InitArgs = [pcmcprof,TECPlotTitle,WriteVars,NormaliseAngularStatistics,AverageAcross0Degrees]
	
	# Write data for each CPCMCSPEC species, for each CMATSTATS material
	# Asks for material first, then for each species seperately, total = IPCMCSPEC + IMATSTATS*IPCMCSPEC
	# EXAMPLE TEXT: SHOULD DATA FOR POSITION  0.25 CM, MATERIAL 5 BE WRITTEN? (1 OR 0)
	SpecArgs = np.ones(IMATSTATS*IPCMCSPEC+IPCMCSPEC, dtype=str).tolist()

	# Should flux(energy) be integrated for angle and flux(angle) be integrated for energy
	IntegrateFluxEnergy = '1'
	IntegrateArgs = [IntegrateFluxEnergy]

	#Concat all pcmc.prof args into a single list of strings
	args = InitArgs+SpecArgs+IntegrateArgs
#	args = [] 						# <<< UNCOMMENT FOR MANUAL ENTRY OF VALUES.
#	print(args)						# UPDATE FEED ARGS IN NAM.READIN AND REMOVE THE ABOVE LINE ONCE FIXED.	
#	print(os.linesep.join(args))	# PRINTS EACH ITEM IN LIST ON SEPERATE LINE
#	#[pcmc.prof,title,1,1,1]+[pcmcmat+pcmcmat*pcmcspecies]+[1] for manual usage

	#=====#=====#

	#Remove old files to avoid any overwrite errors.
	for i in range(0,len(DirAdditions)): os.system('rm -f '+DirAdditions[i])

	#THIS SECTION IS OUT OF DATE, SUBPROCESS.COMMUNICATE REQUIRES "bytes-like objects" NOT STRINGS
	#TypeError: a bytes-like object is required, not 'str'
	#Use predefined arguments if supplied, suppresses output to devnull.
	if len(args) > 0:
		with open(os.devnull, 'w') as fp:
			try:
				subprocess = Popen(ConvProfexe, stdin=PIPE, stdout=fp, encoding='utf8') #noshell=True
				subprocess.communicate(os.linesep.join(args), timeout=15)				#15 second timeout
			except:
				print('')
				print('### Timeout while using conv_prof.exe   ###')
				print('### Check "args" length in AutoConvProf ###')
				print('')
				exit()
		#endwith

	#If no arguments are supplied, run script and allow user inputs.
	elif len(args) == 0:
		os.system(ConvProfexe)
	#endif

	#Update Dir with new filenames, must be supplied manually for now.
	for i in range(0,len(DirAdditions)): Dir.append(Dirlist[l]+DirAdditions[i])

	os.chdir(HomeDir)
	return()
#enddef


#Takes variablenames and checks if variable is radial.
#Returns boolian if variable is radial, used for symmetry options.
#ReverseRadial = IsRadialVariable(Variablelist[l]) 
def IsRadialVariable(variable):

	#List of radially symmetric variable names (i.e. variables that are negative across R=0)
	RadialVariableNames = ['VR-','JR-','FR-','FLUX-R']

	#Return true if supplied variable is within list
	IsRadial = False
	if IsStringInVariable(variable,RadialVariableNames) == True:
		IsRadial = True
	#endif

	return(IsRadial)
#enddef


#Takes array of strings and compares to variable string.
#Returns true if any element of stringarray is in variable.
def IsStringInVariable(variable,stringarray):

	boolian = False
	#Check if each element of string is inside variable.
	for i in range(0,len(stringarray)):
		if stringarray[i] in variable:
			boolian = True
		#endif
	#endfor
	
	return(boolian)
#enddef


#Makeshift way of creating units for each legend entry.
def VariableLabelMaker(variablelist):

	#Define common lists for implicit legend generation.
	Powerlist = ['POW-ALL','POW-TOT','POW-ICP','POW-RF','POW-RF-E']
	Fluxlist = ['FZ-','FR-','EFLUX-R','EFLUX-Z']
	Ionisationlist = ['S-','SEB-']
	Velocitylist = ['VZ-','VR-']
	
	#Define Regular Expression lists for numericised ICP coil set variable names
	RegEx = re.compile('POWICP.')
	POWICPVars = ['POWICP']+[string for string in variablelist if re.match(RegEx, string)]
	
	RegEx = re.compile('ERADIAL.')
	ERADIALVars = ['ERADIAL']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('ETHETA.')
	ETHETAVars = ['ETHETA']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('EAXIAL.')
	EAXIALVars = ['EAXIAL']+[string for string in variablelist if re.match(RegEx, string)]

	RegEx = re.compile('PHASEER.')
	PHASERVars = ['PHASEER']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('PHASE.')
	PHASEVars = ['PHASE']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('PHASEEZ.')
	PHASEZVars = ['PHASEEZ']+[string for string in variablelist if re.match(RegEx, string)]
	
	RegEx = re.compile('BR.')
	BRVars = ['BR']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('BT.')
	BTVars = ['BT']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('BZ.')
	BZVars = ['BZ']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('BRF.')
	BRFVars = ['BRF']+[string for string in variablelist if re.match(RegEx, string)]
	
	RegEx = re.compile('PHASEBR.')
	PHASEBRVars = ['PHASEBR']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('PHASEBT.')
	PHASEBTVars = ['PHASEBT']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('PHASEBZ.')
	PHASEBZVars = ['PHASEBZ']+[string for string in variablelist if re.match(RegEx, string)]
	RegEx = re.compile('J-THETA.')
	JTHETAVars = ['J-THETA']+[string for string in variablelist if re.match(RegEx, string)]
	
	#=====#=====#
	#=====#=====#
	
	Variablelegends = list()
	for i in range(0,len(variablelist)):
		#Explicit Pressure and Species Densities.
		if variablelist[i] in ['PRESSURE']:
			Variable = 'Pressure'
			VariableUnit = '['+str(PressureUnit)+']'			#Default: '[Torr]'
		elif variablelist[i] in ['AR','AR3S']:
			Variable = 'Neutral Ar Density'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] in ['E']:
			Variable = 'Electron Density n$_{e}$'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] in ['AR+']:
			Variable = 'Ar+ Density'
			VariableUnit = '[m$^{-3}$]'

		#Explicit Ionization Rates.
		elif variablelist[i] == 'S-E':
			Variable = 'Bulk e$^-$ Source Rate \n'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-E':
			Variable = 'Secondry e$^-$ Source Rate \n'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'EB-ESORC':
			Variable = 'Secondry e$^-$ Relaxation Rate \n'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'S-AR+':
			Variable = 'Bulk Ar+ Ionization Rate \n'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-AR+':
			Variable = 'Secondry Ar+ Ionization Rate \n'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		
		#Explicit Vibrational States.
		elif variablelist[i] == 'GSH2V1':
			Variable = '1st Vibrational Excited State \n'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] == 'GSH2V4':
			Variable = '4th Vibrational Excited State \n'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] == 'GSH2V14':
			Variable = '14th Vibrational Excited State \n'
			VariableUnit = '[m$^{-3}$]'

		#Explicit Species Temperatures.
		elif variablelist[i] == 'TE':
			Variable = 'Electron Temperature T$_{e}$'
			VariableUnit = '[eV]'
		elif variablelist[i] == 'TG-AVE':
			Variable = 'Neutral Gas Temperature'
			VariableUnit = '[K]'

		#Explicit Species Velocities.
		elif variablelist[i] == 'VZ-NEUTRAL':
			Variable = 'Neutral Axial Velocity'
			VariableUnit = '[ms$^{-1}$]'
		elif variablelist[i] == 'VR-NEUTRAL':
			Variable = 'Neutral Radial Velocity'
			VariableUnit = '[ms$^{-1}$]'
		elif variablelist[i] == 'VZ-ION+':
			Variable = '+Ion Axial Velocity'
			VariableUnit = '[kms$^{-1}$]'
		elif variablelist[i] == 'VR-ION+':
			Variable = '+Ion Radial Velocity'
			VariableUnit = '[kms$^{-1}$]'
		elif variablelist[i] == 'VZ-ION-':
			Variable = '-Ion Axial Velocity'
			VariableUnit = '[kms$^{-1}$]'
		elif variablelist[i] == 'VR-ION-':
			Variable = '-Ion Radial Velocity'
			VariableUnit = '[kms$^{-1}$]'

		#Explicit Species Fluxes.
		elif variablelist[i] == 'EFLUX-Z':
			Variable = 'Electron Axial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif variablelist[i] == 'EFLUX-R':
			Variable = 'Electron Radial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif variablelist[i] == 'FZ-AR+':
			Variable = 'Ar+ Axial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif variablelist[i] == 'FR-AR+':
			Variable = 'Ar+ Radial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif variablelist[i] == 'FZ-AR':
			Variable = 'Ar Axial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif variablelist[i] == 'FR-AR':
			Variable = 'Ar Radial Flux'
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'

		#Explicit Electrodynamic Properties
		elif variablelist[i] in ['PPOT','P-POT']:
			Variable = 'Plasma Potential V$_{p}$'
			VariableUnit = '[V]'
		elif variablelist[i] in ['RHO']:
			Variable = 'Charge Density $\\rho$'
			VariableUnit = '[C cm$^{-3}$]'
			
		elif variablelist[i] in ['EF-TOT']:
			Variable = 'Absolute E-Field Amplitude'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in ERADIALVars+['ER','EAMB-R']:
			Variable = 'Radial E-Field Amplitude $E_{R}$'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in PHASERVars:
			Variable = 'Radial E-Field Phase'
			VariableUnit = '[Radians]'
		elif variablelist[i] in ETHETAVars+['ET','EAMB-T']:
			Variable = 'Azimuthal E-Field Amplitude $E_{\\theta}$'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in PHASEVars:
			Variable = 'Azimuthal E-Field Phase'
			VariableUnit = '[Radians]'
		elif variablelist[i] in EAXIALVars+['EZ','EAMB-Z']:
			Variable = 'Axial E-Field Amplitude $E_{Z}$'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in PHASEZVars:
			Variable = 'Axial E-Field Phase'
			VariableUnit = '[Radians]'

		elif variablelist[i] == 'BRS':
			Variable = 'Radial Static B-field Magnitude $B_{R}$'
			VariableUnit = '[G]'
		elif variablelist[i] == 'BTS':
			Variable = 'Azimuthal Static B-field Magnitude $B_{\\theta}$'
			VariableUnit = '[G]'
		elif variablelist[i] == 'BZS':
			Variable = 'Axial Static B-field Magnitude $B_{Z}$'
			VariableUnit = '[G]'
			
		elif variablelist[i] in BRVars:
			Variable = 'Radial Induced B-field Magnitude $B_{R}$'
			VariableUnit = '[G]'
		elif variablelist[i] in PHASEBRVars:
			Variable = 'Radial Induced B-field Phase'
			VariableUnit = '[Radians]'
		elif variablelist[i] in BTVars+['BTHETA']:
			Variable = 'Azimuthal Induced B-field Magnitude $B_{\\theta}$'
			VariableUnit = '[G]'
		elif variablelist[i] in PHASEBTVars:
			Variable = 'Azimuthal Induced B-field Phase'
			VariableUnit = '[Radians]'
		elif variablelist[i] in BZVars:
			Variable = 'Axial Induced B-field Magnitude $B_{Z}$'
			VariableUnit = '[G]'
		elif variablelist[i] in PHASEBZVars:
			Variable = 'Axial Induced B-field Phase'
			VariableUnit = '[Radians]'
		elif variablelist[i] in BRFVars:
			Variable = 'Total Induced B-field Magnitude $B_{RF}$'
			VariableUnit = '[G]'
			
		elif variablelist[i] in ['JZ-NET']:
			Variable = 'Axial Net Current Density $J_{Z}$'
			VariableUnit = '[mA cm$^{-2}$]'					#Default: '[A cm$^{-2}$]'
		elif variablelist[i] in ['JR-NET']:
			Variable = 'Radial Net Current Density $J_{R}$'
			VariableUnit = '[mA cm$^{-2}$]'					#Default: '[A cm$^{-2}$]'
		elif variablelist[i] in JTHETAVars:
			Variable = 'Azimuthal Net Current Density $J_{\\theta}$'
			VariableUnit = '[mA cm$^{-2}$]'					#Default: '[A cm$^{-2}$]'
		elif variablelist[i] in ['J-TH(MAG)']:
			Variable = 'Azimuthal MCS Electron Current Density $J_{e\\theta}$'
			VariableUnit = '[mA cm$^{-2}$]'					#Default: '[A cm$^{-2}$]'
		elif variablelist[i] in ['J-TH(PHA)']:
			Variable = 'Azimuthal MCS Electron Current Phase $J_{e\\theta}$'
			VariableUnit = '[Radians]'						#Default: '[Radians]'
			
		#Explicit Power Deposition.
		elif variablelist[i] in ['POW-TOT']:
			Variable = 'Total Coupled RF-Power'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] in POWICPVars+['POW-ICP']:		#Note: 	POW-ICP is total icp power, 
			Variable = 'Inductively Coupled RF-Power'		#		POWICP-n is coilset #n icp power
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] in ['POW-RF']:
			Variable = 'RF Power Density'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] in ['POW-RF-E']:
			Variable = 'RF Electron Power Density'
			VariableUnit = '[Wm$^{-3}$]'

		#Explicit Collision Rates.
		elif variablelist[i] == 'COLF':
			Variable = 'Electron Collision Frequency'
			VariableUnit = '[s$^{-1}$]'

		#Implicit Variables.
		elif IsStringInVariable(variablelist[i],Ionisationlist) == True:
			Variable = variablelist[i]
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif IsStringInVariable(variablelist[i],['SRCE-']) == True:
			Variable = variablelist[i]
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif IsStringInVariable(variablelist[i],['T-']) == True:
			Variable = variablelist[i]
			VariableUnit = '[K]'
		elif IsStringInVariable(variablelist[i],Velocitylist) == True:
			Variable = variablelist[i]
			VariableUnit = '[kms$^{-1}$]'
		elif IsStringInVariable(variablelist[i],Fluxlist) == True:
			Variable = variablelist[i]
			VariableUnit = '[m$^{-2}$ s$^{-1}$]'
		elif IsStringInVariable(variablelist[i],['POW-']) == True:
			Variable = variablelist[i]
			VariableUnit = '[Wcm$^{-3}$]'
		elif variablelist[i] in [x.replace('^', '+') for x in AtomicSpecies]:
			Variable = variablelist[i]
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] in AtomicSpecies:
			Variable = variablelist[i]
			VariableUnit = '[m$^{-3}$]'

		#Default if no fitting variable found.
		else:
			Variable = 'Variable'
			VariableUnit = '[Unit]'
		#endif
		Variablelegends.append(Variable+' '+VariableUnit)
	#endfor
	return Variablelegends
#enddef


#Converts units and direction (sign) for input 1D array profiles.
#Takes profile and variable name, returns profile in required SI unit.
#Implicitly calculates for common variables, explicitly for densities.
def VariableUnitConversion(profile,variable):

	#For Pressures, convert to mTorr or Pa as requested, or retain as default Torr.
	if IsStringInVariable(variable,['PRESSURE']) == True:
		for i in range(0,len(profile)):
			if PressureUnit == 'Pa': 		profile[i] = profile[i]*133.333333333	#[Pa]
			elif PressureUnit == 'mTorr': 	profile[i] = profile[i]*1000.0			#[mTorr]
			else: 							profile[i] = profile[i]					#[Torr]
		#endfor
	#endif

	#For ionisation rates, convert from [cm3 s-1] to [m3 s-1]
	if IsStringInVariable(variable,['S-','SEB-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E6
		#endfor
	#endif

	#For fluxes, convert from [cm-2] to [m-2]. (also reverse axial flux)
	if IsStringInVariable(variable,['EFLUX-Z','EFLUX-R','FZ-','FR-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E4
		#endfor
	#endif
	if IsStringInVariable(variable,['EFLUX-Z','FZ-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For velocities, convert from [cms-1] to [ms-1] or [kms-1]. (also reverse axial velocity)
	if IsStringInVariable(variable,['VR-NEUTRAL','VZ-NEUTRAL']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(0.01)			#Neutral [ms-1]
		#endfor
	if IsStringInVariable(variable,['VR-ION+','VZ-ION+','VR-ION-','VZ-ION-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(0.01)*(0.001)	#Ion [kms-1]
		#endfor
	if IsStringInVariable(variable,['VZ-NEUTRAL','VZ-ION+','VZ-ION-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For B-field strengths, convert to Tesla or retain as default Gauss.
	if IsStringInVariable(variable,['BT','BR','BRS','BRF','BTHETA']) == True:
		for i in range(0,len(profile)):
#			if 		Units == 'SI': 	profile[i] = profile[i]/10000			#[T]
#			elif:	Units == 'CGS':	profile[i] = profile[i]					#[G]
			profile[i] = profile[i]											#[G]
		#endfor
	#~~~ AXIAL MAGNETIC FIELD IS NOT REVERSED HERE - RM: NEED TO LOOK INTO THIS... ~~~#
	if IsStringInVariable(variable,['BZ','BZS']) == True:
		for i in range(0,len(profile)):
#			if 		Units == 'SI': 	profile[i] = (profile[i]/10000)#*(-1)	#[T]
#			elif:	Units == 'CGS':	profile[i] = profile[i]					#[G]
			profile[i] = profile[i]#*(-1)										#[G]
		#endfor
	#endif

	#For E-field strengths, convert from [V cm-1] to [V m-1]. (also reverse axial field)
	if IsStringInVariable(variable,['EF-TOT','EAMB-R','EAMB-Z','ETHETA']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]#*100									### [V cm-1] ###
		#endfor
	if IsStringInVariable(variable,['EAMB-Z']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For Current Densities, convert from [A cm-2] to [mA cm-2]. (also reverse axial current)
	if IsStringInVariable(variable,['JZ-NET','JR-NET','J-THETA','J-TH(MAG)']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1000
		#endfor
	if IsStringInVariable(variable,['JZ-NET']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For power densities, convert from [Wcm-3] to [Wm-3].
	if IsStringInVariable(variable,['POW-ALL','POW-TOT','POW-ICP','POW-RF','POW-RF-E']) == True:	
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E6
		#endfor
	#endif

	#For densities, convert from [cm-3] to [m-3]. (AtomicSpecies is defined in icp.nam input)
	if variable in AtomicSpecies or variable in [x.replace('^', '+') for x in AtomicSpecies]:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E6
		#endfor
	#endif

	return(profile)
#enddef


#Applies azimuthal direction (phase) for input 2D data arrays.
#Takes 2D image and variable name, returns image multiplied by sin(phase)
#Only applicable to azimuthal variables, others will be returned unchanged.
def VariableAzimuthalConversion(image,variable):

	#Global toggle to enforce plotting of magnitudes only if requested
	if PlotAzimuthalDirection == False:	return(image)
	
	#=====#=====#
	
	#Determine if supplied profile requires phase conversion
	if variable == 'ETHETA': 		
		phaseprocess,phasevariable = VariableEnumerator(['PHASE'],rawdata_2D[l],header_2Dlist[l])
	elif variable == 'J-THETA':
		phaseprocess,phasevariable = VariableEnumerator(['PHASE'],rawdata_2D[l],header_2Dlist[l])
	elif variable == 'J-TH(MAG)':
		phaseprocess,phasevariable = VariableEnumerator(['J-TH(PHA)'],rawdata_2D[l],header_2Dlist[l])
	else: 
		return(image)
	
	#Extract the appropriate phase data for the supplied variable
	phasemap = ImageExtractor2D(Data[l][phaseprocess[0]],phasevariable[0])

	#Convert from phase [0 --> 2pi] to relative azimuthal direction (0 --> -1 --> +1)
	for i in range(0,len(phasemap)):
		for j in range(0,len(phasemap[i])):
			phasemap[i][j] = np.sin(phasemap[i][j])
		#endfor
	#endfor

	#Apply phasemap to image
	image = image*phasemap
	
	return(image)
#enddef


def ManualPRCCPMesh(Ax=plt.gca()):
	#Plot pocket rocket material dimensions.
	Ax.plot((27*dz[l],27*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((75*dz[l],75*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],27*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((75*dz[l],75*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],75*dz[l]),   ( 0.21, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],75*dz[l]),   (-0.21,-0.21), '-', color='dimgrey', linewidth=2)

	#Alumina Dielectric
	Ax.plot((34*dz[l],74*dz[l]),  ( 0.21,  0.21), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  (-0.21, -0.21), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  ( 0.31,  0.31), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  (-0.31, -0.31), 'c-', linewidth=2)

	#Macor Dielectric
	Ax.plot((34*dz[l],34*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],34*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
	Ax.plot((58*dz[l],58*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((58*dz[l],58*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],58*dz[l]),  ( 47*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],58*dz[l]),  (-47*dr[l],-47*dr[l]), 'b-', linewidth=2)

	#Powered Electrode
	Ax.plot((39*dz[l],49*dz[l]),  ( 25*dr[l], 25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  ( 43*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-43*dr[l],-43*dr[l]), 'r-', linewidth=2)

	#Grounded electrodes
	Ax.plot((34*dz[l],34*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((34*dz[l],34*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
#enddef

#=============#

def ManualPRCCPMMesh(Ax=plt.gca()):
	#Plot pocket rocket material dimensions.
	Ax.plot((27*dz[l],27*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((75*dz[l],75*dz[l]),   (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],27*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((75*dz[l],75*dz[l]),   ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],75*dz[l]),   ( 0.21, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((27*dz[l],75*dz[l]),   (-0.21,-0.21), '-', color='dimgrey', linewidth=2)

	#Alumina Dielectric
	Ax.plot((34*dz[l],74*dz[l]),  ( 0.21,  0.21), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  (-0.21, -0.21), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  ( 0.31,  0.31), 'c-', linewidth=2)
	Ax.plot((34*dz[l],74*dz[l]),  (-0.31, -0.31), 'c-', linewidth=2)

	#Macor Dielectric
	Ax.plot((34*dz[l],34*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],34*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
	Ax.plot((58*dz[l],58*dz[l]),  ( 25*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((58*dz[l],58*dz[l]),  (-25*dr[l],-47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],58*dz[l]),  ( 47*dr[l], 47*dr[l]), 'b-', linewidth=2)
	Ax.plot((34*dz[l],58*dz[l]),  (-47*dr[l],-47*dr[l]), 'b-', linewidth=2)

	#Powered Electrode
	Ax.plot((39*dz[l],49*dz[l]),  ( 25*dr[l], 25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  ( 25*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  (-25*dr[l],-43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-25*dr[l],-25*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  ( 43*dr[l], 43*dr[l]), 'r-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-43*dr[l],-43*dr[l]), 'r-', linewidth=2)

	#Co-axial magnetic ring
	Ax.plot((39*dz[l],49*dz[l]),  ( 48*dr[l], 48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-62*dr[l],-62*dr[l]), 'g-', linewidth=2)

	#Grounded electrodes
	Ax.plot((34*dz[l],34*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((34*dz[l],34*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  (-1.0,-0.21), '-', color='dimgrey', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  ( 1.0, 0.21), '-', color='dimgrey', linewidth=2)
#enddef

#=============#

def ManualHyperionIMesh(Ax=plt.gca()):
	#Plot upstream ICP material dimensions.
	Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68+40)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68-39)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((5*dz[l],51*dz[l]),    ((68+40)*dr[l],(68+40)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((5*dz[l],51*dz[l]),    ((68-39)*dr[l],(68-39)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((51*dz[l],51*dz[l]),   ((68+40)*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((51*dz[l],51*dz[l]),   ((68-39)*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((51*dz[l],111*dz[l]),  ((68+65)*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((51*dz[l],111*dz[l]),  ((68-64)*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((111*dz[l],111*dz[l]), (68*dr[l],(68+65)*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((111*dz[l],111*dz[l]), (68*dr[l],(68-64)*dr[l]), '-', color='lightgrey', linewidth=4)
	#4cm = 74*dz[l], 12cm=111*dz[l]

	#Macor Dielectric
	Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68+40)*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((5*dz[l],5*dz[l]),     (68*dr[l],(68-39)*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((5*dz[l],51*dz[l]),    ((68+40)*dr[l],(68+40)*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((5*dz[l],51*dz[l]),    ((68-39)*dr[l],(68-39)*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((51*dz[l],51*dz[l]),   ((68+40)*dr[l],(68+65)*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((51*dz[l],51*dz[l]),   ((68-39)*dr[l],(68-64)*dr[l]), '-', color='c', linewidth=4)

	#Powered Electrode - 'Metal' Anode
#	Ax.plot((65*dz[l],65*dz[l]),   (2*dr[l],70*dr[l]), '-', color='red', linewidth=4)
#	Ax.plot((65*dz[l],65*dz[l]),   (-2*dr[l],-70*dr[l]), '-', color='red', linewidth=4)

	#Powered ICP Coils - 'Metal'
	Ax.plot((12*dz[l],16*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((12*dz[l],16*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((12*dz[l],12*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((16*dz[l],16*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downstream
	Ax.plot((7*dz[l],11*dz[l]),    (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((7*dz[l],11*dz[l]),    (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((7*dz[l],7*dz[l]),     (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((11*dz[l],11*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

	Ax.plot((22*dz[l],26*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((22*dz[l],26*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((22*dz[l],22*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((26*dz[l],26*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downtream
	Ax.plot((17*dz[l],21*dz[l]),   (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((17*dz[l],21*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((17*dz[l],17*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((21*dz[l],21*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

	Ax.plot((32*dz[l],36*dz[l]),   (114*dr[l],114*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((32*dz[l],36*dz[l]),   (122*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((32*dz[l],32*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((36*dz[l],36*dz[l]),   (114*dr[l],122*dr[l]), '-', color='red', linewidth=5)	#Downstream
	Ax.plot((27*dz[l],31*dz[l]),   (23*dr[l],23*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((27*dz[l],31*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((31*dz[l],31*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((27*dz[l],27*dz[l]),   (23*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Downstream

	#Horseshoe Solenoid
	Ax.plot((44*dz[l],48*dz[l]),   (124*dr[l],124*dr[l]), '-', color='lightgreen', linewidth=5)	#Inboard
	Ax.plot((44*dz[l],48*dz[l]),   (13*dr[l],13*dr[l]), '-', color='lightgreen', linewidth=5)	#Inboard
	Ax.plot((44*dz[l],48*dz[l]),   (132*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Outboard
	Ax.plot((44*dz[l],48*dz[l]),   (5*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)		#Outboard
	Ax.plot((44*dz[l],44*dz[l]),   (124*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Upstream
	Ax.plot((48*dz[l],48*dz[l]),   (124*dr[l],132*dr[l]), '-', color='lightgreen', linewidth=5)	#Downstream
	Ax.plot((44*dz[l],44*dz[l]),   (13*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)	#Upstream
	Ax.plot((48*dz[l],48*dz[l]),   (13*dr[l],5*dr[l]), '-', color='lightgreen', linewidth=5)	#Downstream
#enddef

def ManualHyperionIIMesh(Ax=plt.gca()):
	#Plot upstream ICP material dimensions.
	Ax.plot((1*dz[l],1*dz[l]),       (3*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((1*dz[l],84*dz[l]),      (3*dr[l],3*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((1*dz[l],84*dz[l]),      (109*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)
	Ax.plot((84*dz[l],84*dz[l]),     (3*dr[l],109*dr[l]), '-', color='lightgrey', linewidth=4)

	#Alumina Dielectric
	Ax.plot((2*dz[l],2*dz[l]),     (25*dr[l],87*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((2*dz[l],84*dz[l]),    (25*dr[l],25*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((2*dz[l],84*dz[l]),    (87*dr[l],87*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((14*dz[l],14*dz[l]),   (37*dr[l],75*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((14*dz[l],84*dz[l]),   (37*dr[l],37*dr[l]), '-', color='c', linewidth=4)
	Ax.plot((14*dz[l],84*dz[l]),   (75*dr[l],75*dr[l]), '-', color='c', linewidth=4)

	#Moly Extraction Aperture
	Ax.plot((71*dz[l],84*dz[l]),   (38*dr[l],38*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((71*dz[l],84*dz[l]),   (41*dr[l],41*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((71*dz[l],71*dz[l]),   (38*dr[l],41*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((73*dz[l],73*dz[l]),   (42*dr[l],53*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((74*dz[l],74*dz[l]),   (42*dr[l],53*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((73*dz[l],74*dz[l]),   (53*dr[l],53*dr[l]), '-', color='m', linewidth=4)

	Ax.plot((71*dz[l],84*dz[l]),   (74*dr[l],74*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((71*dz[l],84*dz[l]),   (70*dr[l],70*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((71*dz[l],71*dz[l]),   (70*dr[l],74*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((73*dz[l],73*dz[l]),   (58*dr[l],70*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((74*dz[l],74*dz[l]),   (58*dr[l],70*dr[l]), '-', color='m', linewidth=4)
	Ax.plot((73*dz[l],74*dz[l]),   (58*dr[l],58*dr[l]), '-', color='m', linewidth=4)

	#Powered Electrode - 'Metal' Anode
#	Ax.plot((65*dz[l],65*dz[l]),   (2*dr[l],70*dr[l]), '-', color='red', linewidth=4)
#	Ax.plot((65*dz[l],65*dz[l]),   (-2*dr[l],-70*dr[l]), '-', color='red', linewidth=4)

	#Powered ICP Coils - 'Metal'
	Ax.plot((27*dz[l],27*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((29*dz[l],29*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((27*dz[l],29*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((27*dz[l],29*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)	

	Ax.plot((33*dz[l],33*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((35*dz[l],35*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((33*dz[l],35*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((33*dz[l],35*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)	

	Ax.plot((39*dz[l],39*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((41*dz[l],41*dz[l]),   (9*dr[l],14*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((39*dz[l],41*dz[l]),   (9*dr[l],9*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((39*dz[l],41*dz[l]),   (14*dr[l],14*dr[l]), '-', color='red', linewidth=5)	

	Ax.plot((30*dz[l],30*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((32*dz[l],32*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((30*dz[l],32*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((30*dz[l],32*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

	Ax.plot((36*dz[l],36*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((38*dz[l],38*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((36*dz[l],38*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((36*dz[l],38*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

	Ax.plot((42*dz[l],42*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((44*dz[l],44*dz[l]),   (98*dr[l],103*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((42*dz[l],44*dz[l]),   (98*dr[l],98*dr[l]), '-', color='red', linewidth=5)	
	Ax.plot((42*dz[l],44*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)
	Ax.plot((42*dz[l],44*dz[l]),   (103*dr[l],103*dr[l]), '-', color='red', linewidth=5)

	#Horseshoe Solenoid
	Ax.plot((68*dz[l],68*dz[l]),   (4*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],70*dz[l]),   (4*dr[l],15*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((68*dz[l],70*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],72*dz[l]),   (15*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((72*dz[l],72*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],72*dz[l]),   (22*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)

	Ax.plot((68*dz[l],68*dz[l]),   (4*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],70*dz[l]),   (4*dr[l],15*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((68*dz[l],70*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],72*dz[l]),   (15*dr[l],18*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((72*dz[l],72*dz[l]),   (18*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((70*dz[l],72*dz[l]),   (22*dr[l],22*dr[l]), '-', color='lightgreen', linewidth=5)

	#Iron Magnetic Focus
	Ax.plot((71*dz[l],71*dz[l]),   (43*dr[l],46*dr[l]), '-', color='g', linewidth=6.5)
	Ax.plot((72*dz[l],72*dz[l]),   (43*dr[l],53*dr[l]), '-', color='g', linewidth=6.5)

	Ax.plot((71*dz[l],71*dz[l]),   (65*dr[l],68*dr[l]), '-', color='g', linewidth=6.5)
	Ax.plot((72*dz[l],72*dz[l]),   (58*dr[l],68*dr[l]), '-', color='g', linewidth=6.5)
#enddef

#=============#

def ManualEVgenyMesh(Ax=plt.gca()):
	#Plot upstream ICP material dimensions.
	Ax.plot((2.5*dz[l],2.5*dz[l]), (0*dr[l],20*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((2.5*dz[l],2.5*dz[l]), (-20*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((2.5*dz[l],41*dz[l]),  (20*dr[l],20*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((2.5*dz[l],41*dz[l]),  (-20*dr[l],-20*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((41*dz[l],41*dz[l]),   (20*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((41*dz[l],41*dz[l]),   (-20*dr[l],-90*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((41*dz[l],87*dz[l]),   (90*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((41*dz[l],87*dz[l]),   (-90*dr[l],-90*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((87*dz[l],87*dz[l]),   (0*dr[l],90*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((87*dz[l],87*dz[l]),   (-90*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)

	#Macor Dielectric
	Ax.plot((2.5*dz[l],2.5*dz[l]), (0*dr[l],20*dr[l]), 'c-', linewidth=4)
	Ax.plot((2.5*dz[l],2.5*dz[l]), (-20*dr[l],0*dr[l]), 'c-', linewidth=4)
	Ax.plot((3*dz[l],6*dz[l]),    (20*dr[l],20*dr[l]), 'c-', linewidth=4)
	Ax.plot((3*dz[l],6*dz[l]),    (-20*dr[l],-20*dr[l]), 'c-', linewidth=4)
	Ax.plot((23*dz[l],41*dz[l]),  (20*dr[l],20*dr[l]), 'c-', linewidth=4)
	Ax.plot((23*dz[l],41*dz[l]),  (-20*dr[l],-20*dr[l]), 'c-', linewidth=4)
	Ax.plot((41*dz[l],41*dz[l]),   (20*dr[l],90*dr[l]), 'c-', linewidth=4)
	Ax.plot((41*dz[l],41*dz[l]),   (-20*dr[l],-90*dr[l]), 'c-', linewidth=4)

	#Powered Electrode - LaB6 Cathode
	Ax.plot((6*dz[l],23*dz[l]),   (20*dr[l],20*dr[l]), '-', color='orange', linewidth=5)
	Ax.plot((6*dz[l],23*dz[l]),   (-20*dr[l],-20*dr[l]), '-', color='orange', linewidth=5)

	#Powered Electrode - 'Metal' Anode
	Ax.plot((42*dz[l],86*dz[l]),   (91*dr[l],91*dr[l]), '-', color='red', linewidth=5)
	Ax.plot((42*dz[l],86*dz[l]),   (-91*dr[l],-91*dr[l]), '-', color='red', linewidth=5)

	#Powered ICP Coils - 'Metal'
	Ax.plot((6*dz[l],8*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((6*dz[l],8*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((6*dz[l],8*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((6*dz[l],8*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((6*dz[l],6*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((6*dz[l],6*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((8*dz[l],8*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Downstream
	Ax.plot((8*dz[l],8*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)		#Downstream

	Ax.plot((13.5*dz[l],15.5*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((13.5*dz[l],15.5*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Inboard
	Ax.plot((13.5*dz[l],15.5*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((13.5*dz[l],15.5*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Outboard
	Ax.plot((13.5*dz[l],13.5*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((15.5*dz[l],15.5*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)	#Downtream
	Ax.plot((13.5*dz[l],13.5*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Upstream
	Ax.plot((15.5*dz[l],15.5*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)#Downstream

	Ax.plot((21*dz[l],23*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((21*dz[l],23*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((21*dz[l],23*dz[l]),   (35*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((21*dz[l],23*dz[l]),   (-35*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((21*dz[l],21*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((23*dz[l],23*dz[l]),   (25*dr[l],35*dr[l]), '-', color='red', linewidth=5)		#Downstream
	Ax.plot((21*dz[l],21*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((23*dz[l],23*dz[l]),   (-25*dr[l],-35*dr[l]), '-', color='red', linewidth=5)	#Downstream

	#Solenoid
	Ax.plot((42*dz[l],86*dz[l]),   (93*dr[l],93*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((42*dz[l],86*dz[l]),   (-93*dr[l],-93*dr[l]), '-', color='lightgreen', linewidth=5)
#enddef

def ManualEVgenyMeshOLD(Ax=plt.gca()):
	#Plot upstream ICP material dimensions.
	Ax.plot((0.5*dz[l],0.5*dz[l]), (0*dr[l],10*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((0.5*dz[l],0.5*dz[l]), (-10*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((0.5*dz[l],39*dz[l]),  (10*dr[l],10*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((0.5*dz[l],39*dz[l]),  (-10*dr[l],-10*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((39*dz[l],39*dz[l]),   (10*dr[l], 80*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((39*dz[l],39*dz[l]),   (-10*dr[l], -80*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((39*dz[l],85*dz[l]),   (80*dr[l],80*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((39*dz[l],85*dz[l]),   (-80*dr[l],-80*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((85*dz[l],85*dz[l]),   (0*dr[l],80*dr[l]), '-', color='dimgrey', linewidth=4)
	Ax.plot((85*dz[l],85*dz[l]),   (-80*dr[l],0*dr[l]), '-', color='dimgrey', linewidth=4)

	#Macor Dielectric
	Ax.plot((1*dz[l],4*dz[l]),    (10*dr[l],10*dr[l]), 'c-', linewidth=4)
	Ax.plot((1*dz[l],4*dz[l]),    (-10*dr[l],-10*dr[l]), 'c-', linewidth=4)
	Ax.plot((21*dz[l],39*dz[l]),  (10*dr[l],10*dr[l]), 'c-', linewidth=4)
	Ax.plot((21*dz[l],39*dz[l]),  (-10*dr[l],-10*dr[l]), 'c-', linewidth=4)

	#Powered Electrode - LaB6 Cathode
	Ax.plot((4*dz[l],21*dz[l]),   (10*dr[l],10*dr[l]), '-', color='orange', linewidth=5)
	Ax.plot((4*dz[l],21*dz[l]),   (-10*dr[l],-10*dr[l]), '-', color='orange', linewidth=5)

	#Powered Electrode - 'Metal' Anode
	Ax.plot((40*dz[l],84*dz[l]),   (81*dr[l],81*dr[l]), '-', color='red', linewidth=5)
	Ax.plot((40*dz[l],84*dz[l]),   (-81*dr[l],-81*dr[l]), '-', color='red', linewidth=5)

	#Powered ICP Coils - 'Metal'
	Ax.plot((4*dz[l],6*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((4*dz[l],6*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((4*dz[l],6*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((4*dz[l],6*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((4*dz[l],4*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((4*dz[l],4*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((6*dz[l],6*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Downstream
	Ax.plot((6*dz[l],6*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)		#Downstream

	Ax.plot((11.5*dz[l],13.5*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((11.5*dz[l],13.5*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)#Inboard
	Ax.plot((11.5*dz[l],13.5*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((11.5*dz[l],13.5*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Outboard
	Ax.plot((11.5*dz[l],11.5*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((13.5*dz[l],13.5*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)	#Downtream
	Ax.plot((11.5*dz[l],11.5*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Upstream
	Ax.plot((13.5*dz[l],13.5*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)#Downstream

	Ax.plot((19*dz[l],21*dz[l]),   (15*dr[l],15*dr[l]), '-', color='red', linewidth=5)		#Inboard
	Ax.plot((19*dz[l],21*dz[l]),   (-15*dr[l],-15*dr[l]), '-', color='red', linewidth=5)	#Inboard
	Ax.plot((19*dz[l],21*dz[l]),   (25*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Outboard
	Ax.plot((19*dz[l],21*dz[l]),   (-25*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Outboard
	Ax.plot((19*dz[l],19*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Upstream
	Ax.plot((21*dz[l],21*dz[l]),   (15*dr[l],25*dr[l]), '-', color='red', linewidth=5)		#Downstream
	Ax.plot((19*dz[l],19*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Upstream
	Ax.plot((21*dz[l],21*dz[l]),   (-15*dr[l],-25*dr[l]), '-', color='red', linewidth=5)	#Downstream

	#Solenoid
	Ax.plot((40*dz[l],84*dz[l]),   (83*dr[l],83*dr[l]), '-', color='lightgreen', linewidth=5)
	Ax.plot((40*dz[l],84*dz[l]),   (-83*dr[l],-83*dr[l]), '-', color='lightgreen', linewidth=5)
#enddef

#=============#

def ManualGECMesh(Ax=plt.gca()):	#Greg's GEC overlay code
    thin = 2
    thick = 2
    verythik = 2
    superthick = 3
    
    #Plot upstream ICP material dimensions.
    Ax.plot((56.5*dz[l],56.5*dz[l]), (1.5*dr[l],55.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical right edge
    Ax.plot((34.5*dz[l],34.5*dz[l]), (1.5*dr[l],14.75*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical left edge coil house
    Ax.plot((30.5*dz[l],30.5*dz[l]), (1.5*dr[l],12.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical internal edge coil house
    Ax.plot((40*dz[l],40*dz[l]), (14.75*dr[l],20.75*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical left edge coil rim
    Ax.plot((27.5*dz[l],27.5*dz[l]), (39*dr[l],55.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical left edge electrode
    Ax.plot((27.5*dz[l],27.5*dz[l]), (17.75*dr[l],20.75*dr[l]), '-', color='dimgrey', linewidth=thin)  #vertical right edge coil house rim
    
    Ax.plot((27.5*dz[l],34.5*dz[l]), (55.5*dr[l],55.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal base edge left
    Ax.plot((46*dz[l],56.5*dz[l]), (55.5*dr[l],55.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal base edge right
    Ax.plot((34.5*dz[l],56.5*dz[l]), (1.5*dr[l],1.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal top edge 
    Ax.plot((34.5*dz[l],40*dz[l]), (14.75*dr[l],14.75*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal top edge coil house rim
    Ax.plot((27.5*dz[l],40*dz[l]), (20.75*dr[l],20.75*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal bottom edge coil house rim
    Ax.plot((0.01*dz[l],30.5*dz[l]), (1.5*dr[l],1.5*dr[l]), '-', color='dimgrey', linewidth=thin)  #horizontal top edge internal
    
    	#Dielectric window
    Ax.plot((31.5*dz[l],31.5*dz[l]),  (12.5*dr[l],17.75*dr[l]), 'orange', linewidth=thin)   #Vertical right edege dielectric
    Ax.plot((0.01*dz[l],31.5*dz[l]),  (12.5*dr[l],12.5*dr[l]), 'orange', linewidth=thin)   
    Ax.plot((0.01*dz[l],31.5*dz[l]),  (17.75*dr[l],17.75*dr[l]), 'orange', linewidth=thin)   
    
    	#Powered Electrode
    Ax.plot((0.01*dz[l],27.5*dz[l]),   (39*dr[l],39*dr[l]), 'dimgrey', linewidth=superthick)
    	
    	#Dielectric spacer on Electrode
    #Ax.plot((27.55*dz[l],27.55*dz[l]), (39.5*dr[l],40*dr[l]), 'orange', linewidth=verythik)
    
    	#Gas inlet - 'Metal' A
    Ax.plot((34*dz[l],46*dz[l]), (55.5*dr[l],55.5*dr[l]), 'dimgrey', linewidth=thick)
    	
    
    	#Powered ICP Coils - 'Metal'
    Ax.plot((1.5*dz[l],3.5*dz[l]),   (12.25*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#1st coil bot
    Ax.plot((1.5*dz[l],3.5*dz[l]),   (10*dr[l],10*dr[l]), '-', color='red', linewidth=thick)		#1st coil top
    Ax.plot((1.5*dz[l],1.5*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#1st coil left
    Ax.plot((3.5*dz[l],3.5*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#1st coilright
    
    Ax.plot((7.25*dz[l],9.25*dz[l]),   (12.25*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#2nd coil bot
    Ax.plot((7.25*dz[l],9.25*dz[l]),   (10*dr[l],10*dr[l]), '-', color='red', linewidth=thick)		#2nd coil top
    Ax.plot((7.25*dz[l],7.25*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#2nd coil left
    Ax.plot((9.25*dz[l],9.25*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#2nd coilright
    
    Ax.plot((13.1*dz[l],15*dz[l]),   (12.25*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#3rd coil bot
    Ax.plot((13.1*dz[l],15*dz[l]),   (10*dr[l],10*dr[l]), '-', color='red', linewidth=thick)		#3rd coil top
    Ax.plot((13.1*dz[l],13.1*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#3rd coil left
    Ax.plot((15*dz[l],15*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#3rd coilright
    
    Ax.plot((19*dz[l],20.9*dz[l]),   (12.25*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#4th coil bot
    Ax.plot((19*dz[l],20.9*dz[l]),   (10*dr[l],10*dr[l]), '-', color='red', linewidth=thick)		#4th coil top
    Ax.plot((19*dz[l],19*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#4th coil left
    Ax.plot((20.9*dz[l],20.9*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#4th coilright
    
    Ax.plot((24.75*dz[l],26.5*dz[l]),   (12.25*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#5th coil bot
    Ax.plot((24.75*dz[l],26.5*dz[l]),   (10*dr[l],10*dr[l]), '-', color='red', linewidth=thick)		#5th coil top
    Ax.plot((24.75*dz[l],24.75*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#5th coil left
    Ax.plot((26.5*dz[l],26.5*dz[l]),   (10*dr[l],12.25*dr[l]), '-', color='red', linewidth=thick)		#5th coilright
    
    Ax.plot((0.01*dz[l], 0.2*dz[l]), (57*dr[l],57*dr[l]), color= 'black', linewidth = 14)


#===================##===================#
#===================##===================#

































#====================================================================#
				 	 #READING DATA INTO MEMORY#
#====================================================================#

print('----------------------')
print('Beginning Data Readin:')
print('----------------------')

#Extraction and organization of data from .PDT files.
for l in tqdm(range(0,numfolders)):

	#Load data from TECPLOT2D file and unpack into 1D array.
	rawdata, nn_2D = ExtractRawData(Dir,'TECPLOT2D.PDT',l)
	rawdata_2D.append(rawdata)

	#Read through all variables for each file and stop when list ends.
	Variablelist,HeaderEndMarker = ['Radius','Height'],'ZONE'
	for i in range(2,nn_2D):
		if HeaderEndMarker in str(rawdata_2D[l][i]): break
		else: Variablelist.append(str(rawdata_2D[l][i][:-2].strip(' \t\n\r\"')))
		#endif
	#endfor
	numvariables_2D,header_2D = len(Variablelist),len(Variablelist)+2
	header_2Dlist.append(header_2D)

	#Seperate total 1D data array into sets of data for each variable.
	CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_2D[l],header_2D,numvariables_2D)

	#Save all variables for folder[l] to Data.
	#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
	Data.append(CurrentFolderData)


#===================##===================#
#===================##===================#

#	#Kinetics data readin - NOT CURRENTLY EMPLOYED IN ANY DIAGNOSTICS
#	if True == True:
#
#		#Load data from TECPLOT_KIN file and unpack into 1D array.
#		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
#		rawdata_kin.append(rawdata)
#
#		#Read through all variables for each file and stop when list ends.
#		KinVariablelist,KinHeaderEndMarker = ['T (S)'],'ZONE'
#		for i in range(2,nn_2D):
#			if KinHeaderEndMarker in str(rawdata_kin[l][i]): 
#				I = int(filter(lambda x: x.isdigit(), rawdata_kin[l][i].split(',')[0]))
#				break
#			else: KinVariablelist.append(str(rawdata_kin[l][i][:-2].strip(' \t\n\r\"')))
#			#endif
#		#endfor
#		numvariables_kin,header_kin = len(KinVariablelist),len(KinVariablelist)+2
#		header_kinlist.append(header_kin)
#
#		#Seperate total 1D data array into sets of data for each variable.
#		CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_kin[l],header_kin,numvariables_kin, Zmesh=I,Dimension='1D')
#
#		#Save all variables for folder[l] to Data.
#		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
#		DataKin.append(CurrentFolderData)
#	#endif


#===================##===================#
#===================##===================#

	#IEDF/NEDF file readin.
	if True in [savefig_IEDFangular,savefig_IEDFtrends]:

		#Define arguments and autorun conv_prof.exe if possible.
		#### THIS IS HACKY, WON'T ALWAYS WORK, ARGS LIST NEEDS AUTOMATING ####
		IEDFVarArgs = ['1','1','1','1','1'] 						#Works for 2 species 1 surface.
		ExtraArgs = ['1','1','1','1','1','1','1','1','1','1']#[]	#Hack For Additional Species
		Args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ExtraArgs + ['0','0']
		DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt','iprofile_zones_tec1d.pdt','nprofile_zones_tec1d.pdt']
		#try: AutoConvProf('./conv_prof.exe',Args,DirAdditions)
		#except: print('ConvProf Failure:'+Dirlist[l])
		AutoConvProf('./conv_prof.exe',Args,DirAdditions)

		#Load data from IEDFprofile file and unpack into 1D array.
		rawdata, nn_IEDF = ExtractRawData(Dir,'iprofile_tec2d.pdt',l)
		rawdata_IEDF.append(rawdata)

		#Read through all variables for each file and stop when list ends.
		IEDFVariablelist,HeaderEndMarker = ['Theta [deg]','Energy [eV]'],'ZONE'
		for i in range(2,nn_IEDF):
			#Grab EDFangle(I),EDFbins(J) values from the ZONE line, these outline the datablock size.
			if HeaderEndMarker in str(rawdata_IEDF[l][i]): 
				I = list(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[0]))	#discrete digits
				I = int( ''.join(I) ); EDFangle = I					#Number of EDF angle bins [Integer]
				J = list(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[1]))	#discrete digits
				J = int( ''.join(J) ); EDFbins = J					#Number of EDF energy bins [Integer]
				break
			else: IEDFVariablelist.append(str(rawdata_IEDF[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		numvariables_IEDF,header_IEDF = len(IEDFVariablelist),len(IEDFVariablelist)+2
		header_IEDFlist.append(header_IEDF)

		#Seperate total 1D data array into sets of data for each variable.
		#Data is stored in 2D array of shape: [EDFangle,EDFbins] or [I,J]
		CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_IEDF[l],header_IEDF,numvariables_IEDF,0,I,J)

		#Save all variables for folder[l] to Data.
		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
		DataIEDF.append(CurrentFolderData)
	#endif


#===================##===================#
#===================##===================#

	#EEDF data readin.
	if savefig_EEDF == True:

		#Load data from MCS.PDT file and unpack into 1D array.
		rawdata, nn_mcs = ExtractRawData(Dir,'boltz_tec.pdt',l)
		rawdata_mcs.append(rawdata)

		#Unpack each row of data points into single array of floats.
		#Removing 'spacing' between the floats and ignoring variables above data.
		Energy,Fe = list(),list()
		for i in range(3,len(rawdata_mcs[l])):
			if 'ZONE' in rawdata_mcs[l][i]:
				EEDF_TDlist.append( rawdata_mcs[l][i].split('"')[-2].strip(' ') )
				DataEEDF.append([Energy,Fe])
				Energy,Fe = list(),list()
			#endif
			try:
				Energy.append( float(rawdata_mcs[l][i].split()[0]) )
				Fe.append( float(rawdata_mcs[l][i].split()[1]) )
			except:
				NaN_Value = 1
			#endtry
		#endfor
		a,b = 0,5
		for i in range(a,b):
			plt.plot(DataEEDF[i][0],DataEEDF[i][1], lw=2)
		plt.legend(EEDF_TDlist[a:b])
		plt.xlabel('Energy [eV]')
		plt.ylabel('F(e) [eV-3/2]')
		plt.show()
	#endif


#===================##===================#
#===================##===================#

	if True in [savefig_convergence,savefig_temporalprofiles]:

		#Load data from movie_icp file and unpack into 1D array.
		rawdata,nn_itermovie = ExtractRawData(Dir,'movie_icp.pdt',l)
		rawdata_itermovie.append(rawdata)

		#Read through all variables for each file and stop when list ends. 
		#movie_icp has geometry at top, therefore len(header) != len(variables).
		#Only the first encountered geometry is used to define variable zone.
		VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ITER'
		variablelist,numvar = list(),0
		for i in range(2,nn_itermovie):
			if HeaderEndMarker in str(rawdata[i]): 
				header_iter = i+1		# +1 to skip to data row.
				break
			if VariableEndMarker in str(rawdata[i]) and numvar == 0:
				numvar = (i-3)	# -3 to not include R,Z and remove overflow.
			if len(rawdata[i]) > 1 and numvar == 0: 
				variablelist.append(str(rawdata_itermovie[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		header_itermovie.append(header_iter)

		#Rough method of obtaining the movie_icp iter locations for data extraction.
		Iterloc = list()
		MovieIterlist.append(list())
		for i in range(0,len(rawdata)):
			if "ITER=" in rawdata[i]:
				Iterloc.append(i+1)

				IterStart=rawdata[i].find('ITER')
				MovieIterlist[l].append(rawdata[i][IterStart:IterStart+9])
			#endif
		#endfor
#		print( len(Iterloc))
		#Cycle through all iterations for current datafile, appending per cycle.
		CurrentFolderData,CurrentFolderIterlist = np.array([[]]),np.array([[]])
		for i in range(0,len(Iterloc)):
			if i == 0:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar+2,offset=2)
#				print( 'option1, ITER DATA', len(CurrentIterData[l]))
				CurrentFolderData = np.array([CurrentIterData[0:numvar]])
#				print( 'option1, Current Folder', len(CurrentFolderData[l]), len(CurrentFolderData[l][0]))
			else:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar)
#				print( 'option1, ITER DATA', len(CurrentIterData[l]))                
				CurrentFolderData = np.concatenate((CurrentFolderData, np.array([CurrentIterData])), axis=0)
#				print( 'option2 Current Folder', len(CurrentFolderData[l]))
			#endif
		#endfor
		if l == 0:        
			IterMovieData = np.array([CurrentFolderData])
		else:        
			IterMovieData = np.concatenate((IterMovieData, np.array([CurrentFolderData])), axis=0)
            
#		print('final', IterMovieData, len(IterMovieData), len(IterMovieData[0]), len(IterMovieData[0][0]), len(IterMovieData[0][0][0]),type(IterMovieData))        
	#endif

#===================##===================#
#===================##===================#
#===================##===================#


#Create global list of all variable names and find shortest list.
for l in range(0,numfolders):
	#Alphabetize the Variablelist and keep global alphabetized list.
	tempvarlist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])[1]
	tempvarlist = sort(tempvarlist)
	numvars = len(tempvarlist)

	Globalvarlist.append(tempvarlist)
	Globalnumvars.append(numvars)
#endfor

#Find the folder with the fewest avaliable variables.
val, idx = min((val, idx) for (idx, val) in enumerate(Globalnumvars))
Comparisonlist = Globalvarlist[idx]

#===================##===================#
#===================##===================#


#Empty and delete any non-global data lists.
tempdata,tempdata2 = list(),list()
data_array,templineout = list(),list()
Energy,Fe,rawdata_mcs = list(),list(),list()
Variablelist,variablelist = list(),list()
HomeDir,DirContents = list(),list()
del RADIUS,RADIUST,HEIGHT,HEIGHTT,DEPTH,ISYM,IXZ
del data_array,tempdata,tempdata2,templineout
del Variablelist,variablelist
del Energy,Fe,rawdata_mcs
del HomeDir,DirContents


#Alert user that readin process has ended and continue with selected diagnostics.
if any([savefig_plot2D, savefig_phaseresolve2D, savefig_convergence, savefig_monoprofiles, savefig_multiprofiles, savefig_compareprofiles, savefig_temporalprofiles, savefig_sheathdynamics, savefig_phaseresolve1D, savefig_PROES, savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
	print( '----------------------------------------')
	print( 'Data Readin Complete, Starting Analysis:')
	print( '----------------------------------------')
else:
	print( '------------------')
	print( 'Analysis Complete.')
	print( '------------------')
#endif


#=====================================================================#
#=====================================================================#


























#====================================================================#
				  #COMMONLY USED PLOTTING FUNCTIONS#
#====================================================================#

def Matplotlib_GlobalOptions():
#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html

#	mpl.style.use('classic')								#Resets to classic 1.x.x format
	
	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
	mpl.rcParams['figure.dpi'] = 100						#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
	mpl.rcParams['image.interpolation'] = image_interp		#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
	mpl.rcParams['image.cmap'] = image_cmap 				#Define global colourmap

	#Axis options
	mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'	#View limits coencide with axis ticks
	mpl.rcParams['axes.xmargin'] = 0						#Set default x-axis padding
	mpl.rcParams['axes.ymargin'] = 0						#Set default y-axis padding
	mpl.rcParams['errorbar.capsize'] = 3					#Set error bar end cap width
	mpl.rcParams['font.size'] = 20							#Set global fontsize
	mpl.rcParams['legend.fontsize'] = 'large'				#Set legend fontsize
	mpl.rcParams['figure.titlesize'] = 'medium'				#Set title fontsize

	#Line and Colour options
#	from cycler import cycler								#See below
#	mpl.rcParams['axes.prop_cycle']=cycler(color='bgrcmyk')	#Set default colour names
	mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

	#Maths and Font options
	mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
	mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

	return()
#enddef
Matplotlib_GlobalOptions()	#MUST BE RUN BEFORE ANY DIAGNOSTICS!!!!


#=========================#
#=========================#


def tecplot_cmap():
#	Creates a colourmap closely approximating the Tecplot "modern" map
#
#		Python tutor:	https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
#		RGB Hex Codes:	https://www.rapidtables.com/web/color/RGB_Color.html
#		Colour Picker:	sudo apt-get install gpick 
#
# 	Colour dictionary details the vertices of the linear colour scales for rgb
#		x1		::	Fraction of cbar range (0=min, 1=max)
#		yleft	::  Percent of colour at start of X1 range
#		yright	::	Percent of colour at end of X1 range
#
#	Colour Order (low to high):
#		43006F 27.17% red, 00.00% green, 43.35% blue
#		5BC4CB 35.54% red, 76.56% green, 79.29% blue
#		00C952 00.00% red, 78.51% green, 32.03% blue
#		CDF100 80.07% red, 94.14% green, 00.00% blue
#		AC6313	67.18% red, 38.67% green, 07.42% blue
#		F70000	96.86% red, 00.00% green, 00.00% blue
#
	from matplotlib.colors import ListedColormap, LinearSegmentedColormap

						#x1			#yleft		#yright
	cdict = {'red':   [[0.00,		0.00,		0.26],
			           [0.17,		0.26,		0.36],
			           [0.33,		0.36,		0.00],
			           [0.50,		0.00,		0.80],
			           [0.66,		0.80,		0.67],
			           [0.83,		0.67,		0.97],
			           [1.00,		0.97,		1.00]],
			           
			 'green': [[0.00,		0.00,		0.00],
			           [0.17,		0.00,		0.77],
			           [0.33,		0.77,		0.79],
			           [0.50,		0.79,		0.94],
			           [0.66,		0.94,		0.38],
			           [0.83,		0.38,		0.00],
			           [1.00,		0.00,		0.00]],
			           
			 'blue':  [[0.00,		0.00,		0.43],
			           [0.17,		0.43,		0.79],
			           [0.33,		0.79,		0.32],
			           [0.50,		0.32,		0.00],
			           [0.66,		0.00,		0.07],
			           [0.83,		0.07,		0.00],
			           [1.00,		0.00,		0.00]]}

	cmap = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
	rgba = cmap(np.linspace(0, 1, 256))
	
	return(cmap,cdict)
#enddef           

#Load TECPLOT colourmap
#imshow(Array,cmap=tecplotcmap)
tecplotcmap,tecplotcdict = tecplot_cmap()

#Load IDL colourmap
#Requires std_gamma_II.txt is modules directory
Filename = 'Modules/std_gamma_II.txt'          
map = np.loadtxt(Filename, delimiter=',')
IDL_Gamma_II = col.ListedColormap(map.T, name='IDL_Gamma_II')

#=========================#
#=========================#


def plot_linearmap(cdict):
#	Shows linear rgb colour fractions for cmap colour dictionary
#
#	USAGE:
#		cmap,cdict = tecplot_cmap()
#		plot_linearmap(cdict)

	from matplotlib.colors import ListedColormap, LinearSegmentedColormap

	newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
	rgba = newcmp(np.linspace(0, 1, 256))
	
	fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)

	col = ['r', 'g', 'b']	
	for xx in [0.25, 0.5, 0.75]:
		ax.axvline(xx, color='0.7', linestyle='--')
	#endfor
	for i in range(3):
		ax.plot(np.arange(256)/256, rgba[:, i], color=col[i])
	#endfor
	
	ax.set_xlabel('index')
	ax.set_ylabel('RGB')
	plt.show()
#enddef


#=========================#
#=========================#


def ImageExtractor2D(Data,Variable=[],Rmesh=0,Zmesh=0):
#Returns a 2D array of inputted data with size [R_mesh] x [Z_mesh]
#Can optionally perform variable unit conversion if required.
#Image = ImageExtractor2D(Data,Variable=[]):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Zmesh == 0:
		Rmesh,Zmesh = int(R_mesh[l]),int(Z_mesh[l])
	#endif
	
	#Create empty 2D image of required size.
	numrows = int(len(Data)/Rmesh)
	Image = np.zeros([int(numrows),int(Rmesh)])

	#Reshape data into 2D array for further processing.
	for j in range(0,numrows):
		for i in range(0,Rmesh):
			Start = Rmesh*j
			Row = Zmesh-1-j
			Image[Row,i] = Data[Start+i]
		#endfor
	#endfor

	#Convert units if required.
	Image = VariableUnitConversion(Image,Variable)
	
	#Apply direction to azimuthal magnitudes if required
	Image = VariableAzimuthalConversion(Image,Variable)

	return(Image)
#enddef


#=========================#
#=========================#


def SymmetryConverter(Image,Radial=False):
#Takes 2D image array and produces symmetric image about central axis.
#Returns symmetric 2D image array, allows radially negative values.
#Image = SymmetryConverter(Image,Radial=False)

	#Create new image by reversing and adding itself on the LHS.
	if image_plotsymmetry == True and int(ISYMlist[l]) == 1:
		SymImage = np.zeros([len(Image),2*len(Image[0])])
		if Radial == False:
			for m in range(0,len(Image)):
				SymImage[m] = np.concatenate([Image[m][::-1],Image[m]])
			#endfor
		elif Radial == True:
			for m in range(0,len(Image)):
				SymImage[m] = np.concatenate([-Image[m][::-1],Image[m]])
			#endfor
		#endif
		Image = SymImage
	#endif

	return(Image)
#enddef


#=========================#
#=========================#


def ScaleArray(Array,ScaleFactors):
#Scales a 2D array by the requested X and Y scale factors, can scale asymmetrically.
#Takes a rectilinear 2D array of floats and returns a scaled rectilinear 2D array.
#Employs a linear interpolation scheme when mapping new datapoints.
#ScaledArray = ScaleArray(2DArray,ScaleFactors=[3,3])

	#Naughty module import
	from scipy.ndimage.interpolation import map_coordinates

	#Define old and new array scales based upon supplied factors
	OldScale = [len(Array),len(Array[0])]
	NewScale = [int(OldScale[0]*ScaleFactors[0]),int(OldScale[1]*ScaleFactors[1])]

	#Create new 2D array with required NewScale
	Array,NewDims = np.asarray(Array),list()
	for OriginalLength, NewLength in zip(Array.shape, (NewScale[0],NewScale[1])):
		NewDims.append(np.linspace(0, OriginalLength-1, NewLength))
	#Map old data onto new array employing a linear interpolation
	Coords = np.meshgrid(*NewDims, indexing='ij')
	ScaledArray = map_coordinates(Array, Coords)

	return(ScaledArray)
#enddef


#=========================#
#=========================#

def figure(aspectratio=[],subplots=1,shareX=False):
#Create figure of desired size and with variable axes.
#Returns figure and axes seperately.
#fig,ax = figure(image_aspectratio,1,shareX=False)

	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
	else:
		fig, ax = plt.subplots(subplots, figsize=(10,10), sharex=shareX)
	#endif
	return(fig,ax)
#enddef


#=========================#
#=========================#


def CropImage(ax=plt.gca(),Extent=[],Apply=True,Rotate=True):
#Crops 2D image taking account of image rotation options.
#Takes image axis (assumes default axis), use figure()
#Input Extent format: [ [Rmin,Rmax], [Zmin,Zmax] ] in cm.
#Returns cropping limits in format: [ [R1,R2],[Z1,Z2] ] in cm.
#CropImage(ax[0],Extent=[[R1,R2],[Z1,Z2]],Apply=True,Rotate=True), 

	#Obtain default limits and rotate if needed, doesn't crash if no crop applied.
	#R1,R2 are radial limits of image, Z1,Z2 are axial limits. (non-rotated)
	R1,R2 = ax.get_xlim()[0],ax.get_xlim()[1]
	Z1,Z2 = ax.get_ylim()[0],ax.get_ylim()[1]
	if Rotate == True and image_rotate == True: 
		R1,Z1, R2,Z2 = Z1,R1, Z2,R2
	#endif

	#Set requested cropping limits from function-call or from global.
	if len(Extent) == 2:
		radialcrop,axialcrop = Extent[0],Extent[1]
	else:
		radialcrop = image_radialcrop
		axialcrop = image_axialcrop
	#endif

	#Extract cropping dimentions from image_<input>.
	if len(radialcrop) == 1:
		R1,R2 = -(radialcrop[0]),radialcrop[0]
	elif len(radialcrop) == 2:
		R1,R2 = radialcrop[0],radialcrop[1]
	#endif
	if len(axialcrop) == 1:
		Z1,Z2 = 0,axialcrop[0]
	elif len(axialcrop) == 2:
		Z1,Z2 = axialcrop[0],axialcrop[1]
	#endif

	#Rotate cropping dimentions to match image rotation.
	if Rotate == True:
		if image_rotate == 00:
			Z1,Z2 = Z2,Z1
		elif image_rotate == 90 or image_rotate == True:
			R1,Z1 = Z1,R1
			R2,Z2 = Z2,R2
		elif image_rotate == 180 or image_rotate == False:
			Z1,Z2 = Z1,Z2
		elif image_rotate == 270:
			R1,Z2 = Z2,R1
			R2,Z1 = Z1,R2
		#endif
	#endif

	#Apply cropping dimensions to image.
	if Apply == True:
		ax.set_xlim(R1,R2)
		ax.set_ylim(Z1,Z2)
	#endif

	#Return cropped dimensions in SI units.
	return([[R1,R2],[Z1,Z2]])
#enddef


#=========================#
#=========================#

###
#NB THIS NEEDS A COMPLETE OVERHAUL AND SPLIT INTO SIMPLER FUNCTIONS!
###

def CbarMinMax(Image,PROES=False,Symmetry=image_plotsymmetry):
#Determines min/max colourbar scale within the cropped frame of an image array.
#Takes a 2D image, and returns the min/max value within the cropped region.
#Assumes image symmetry for best results, otherwise R1 is set to zero.
#Works for PROES images too, requires PROES='Axial' or 'Radial'.
#[Minimum,Maximum] = CbarMinMax(Image,PROES=False)

	#Return user defined limits if specified.
	if len(image_cbarlimit) == 2:
		cropmin = image_cbarlimit[0]
		cropmax = image_cbarlimit[1]
		return([cropmin,cropmax])
	#endif

	#Ensure limits are in line with any requested mathematical constraints
#	if image_logplot == True: Image = np.log(Image)
#	if image_normalise == True: Image = Normalise(Image)

	#Extract min/max from cropped region if a region is supplied.
	if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:

		#Import global cell sizes with correct rotation.
		dR,dZ = dr[l],dz[l]
		if image_rotate == True: dR,dZ = dZ,dR
		#endif

		#Convert cropped SI region (CropExtent) to cell region (R1,R2,Z1,Z2).
		#CropExtent is extracted using CropImage function, which applies a rotation.
		CropExtent = CropImage(Apply=False)		#CropExtent applies rotation internally
		R1 = int(CropExtent[0][0]/dR)
		R2 = int(CropExtent[0][1]/dR)
		Z1 = int(CropExtent[1][0]/dZ)
		Z2 = int(CropExtent[1][1]/dZ)
		#endif

		#Re-rotate back so that radial and axial crops align correctly.
		#R1,R2 are radial limits of image, Z1,Z2 are axial limits.
		if image_rotate == True:
			R1,Z1 = Z1,R1
			R2,Z2 = Z2,R2
		#Cell Images have no negative R values, unlike SI (top left cell = 0,0) 
		#To align properly such that R=0 on axis, add R_mesh to both cropping limits.
		if Symmetry == True:
			R1 = int(R_mesh[l] + R1)
			R2 = int(R_mesh[l] + R2)
		#If image has no symmetry applied, replace negative R1 with zero to align axis.
		if Symmetry == False and R1 < 0: 
			R1 = 0
		#endif
		#Reverse axial ROI if required to avoid zero size images, does not affect plotting.
		#HPEM axial zero in top left corner, array must be read with zero in 'bottom left'.
		if Z1 > Z2: 
			Z1,Z2 = Z2,Z1
		#endif

		#Crop the cell region to the desired region
		if PROES == False:
			#Default Orientation, image orientated radially then axially Image[R][Z]
			if Symmetry == False:
				Image = Image[Z1:Z2]
				Image = np.asarray(Image).transpose()
				Image = Image[R1:R2]
				Image = np.asarray(Image).transpose()
			#90 Deg Rotated Orientation, image orientated axially then radially Image[Z][R]
			elif Symmetry == True:
				Image = Image[Z1:Z2]
				Image = np.asarray(Image).transpose()
				Image = Image[R1:R2]
				Image = np.asarray(Image).transpose()
			#endif
		#Crop the cell region for PROES images, they only require y-axis cropping (R,Z).
		elif PROES == 'Axial':
			Image = Image[::-1][Z1:Z2]			#Image[::-1] Accounts for reversed PhaseData order
		elif PROES == 'Radial':
			Image = Image[R1:R2]	
		#endif
	#endif
	
	#=====#=====#
	
	#Flatten image and obtain min/max in region, defaults to full image if no cropping.
	#1D array of 2D image, concats each radial profile in axial order
	flatimage = [item for sublist in Image for item in sublist]
	flatimage = [float(x) for x in flatimage]
	
	#Filtered 1D array, removing any nan's or inf's
	try: flatimage = [x for x in flatimage if not (m.isinf(x) or m.isnan(x))]
	except: print( 'Image Filtering Warning - Cbar may be scaled incorrectly' )
	
	#Take min and max value from the filtered array to use as cbar limits
	cropmin, cropmax = min(flatimage), max(flatimage)
	
	#Return cropped values in list [min,max], as required by colourbar.
	return([cropmin,cropmax])
#enddef


#=========================#
#=========================#


def ImageOptions(fig,ax=plt.gca(),Xlabel='',Ylabel='',Title='',Legend=[],Crop=True,Rotate=True):
#Applies plt.options to current figure based on user input.
#Returns nothing, open figure is required, use figure().
#For best results call immediately before saving/displaying figure.
#ImageOptions(plt.gca(),Xlabel,Ylabel,Title,Legend,Crop=False)

	#Apply user overrides to plots.
	if len(titleoverride) > 0:
		Title = titleoverride
	if len(legendoverride) > 0:
		Legend = legendoverride
	if len(xlabeloverride) > 0:
		Xlabel = xlabeloverride[0]
	if len(ylabeloverride) > 0:
		Ylabel = ylabeloverride[0]
	#endif

	#Set title and legend if one is supplied.
	if len(Title) > 0:
		ax.set_title(Title, fontsize=14, y=1.03)
	if len(Legend) > 0:
		ax.legend(Legend, loc=1, fontsize=16, frameon=False)
	#endif

	#Set labels and ticksize.
	ax.set_xlabel(Xlabel, fontsize=24)
	ax.set_ylabel(Ylabel, fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#Force scientific notation for all axes, accounting for non-scalar x-axes.
	#	RM: THIS SEEMS TO ALSO BE DEPRECIATED - NEED TO REPLACE WITH NEW Matplotlib COMMAND
#	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
#	except: Axes_Contain_Strings = True
	#
#	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')		#Matplotlib v2.x.x TICKFORMAT
#	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')		#Matplotlib v2.x.x TICKFORMAT
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
	#endif

	#Plot mesh outline if requested.	### HACKY ###
	if image_plotmesh == True:
		mesh_auto_plot = 1 #AUTO PLOT MESH NOT implementED!! REQUIRES initmesh.out READER#
	elif image_plotmesh == 'PRCCP' and Crop == True:	
		ManualPRCCPMesh(ax)
	elif image_plotmesh == 'PRCCPM' and Crop == True:	
		ManualPRCCPMMesh(ax)
	elif image_plotmesh == 'GEC' and Crop == True:
		ManualGECMesh(ax)
	elif image_plotmesh == 'EVgeny' and Crop == True:
		ManualEVgenyMesh(ax)
	elif image_plotmesh == 'HyperionI' and Crop == True:
		ManualHyperionIMesh(ax)
	elif image_plotmesh == 'HyperionII' and Crop == True:
		ManualHyperionIIMesh(ax)

	#endif

	#Crop image dimensions, use provided dimensions or default if not provided.
	if isinstance(Crop, (list, np.ndarray) ) == True:
		CropImage(ax,Extent=Crop,Rotate=Rotate)
	elif any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
		if Crop == True:
			CropImage(ax,Rotate=Rotate)
		#endif
	#endif

	#Arrange figure such that labels, legends and titles fit within frame.
	fig.tight_layout()

	return()
#enddef


#=========================#
#=========================#


def Colourbar(ax=plt.gca(),Label='',Ticks=5,Lim=[]):
#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#cbar = Colourbar(ax[0],'Label',5,Lim=[0,1])

	#Set default font and spacing options and modify if required
	Rotation,Labelpad = 270,30
	LabelFontSize,TickFontsize = 24,18
	if '\n' in Label: Labelpad += 25		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)

	#Set number of ticks, label location and define scientific notation.
	cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.locator = ticker.MaxNLocator(nbins=Ticks)
	cbar.ax.yaxis.offsetText.set(size=TickFontsize)
	yticks(fontsize=TickFontsize)

	#Apply colourbar limits if specified.  (lim=[min,max])
	if len(Lim) == 2: im.set_clim(vmin=Lim[0], vmax=Lim[1])

	return(cbar)
#enddef


#=========================#
#=========================#


def InvisibleColourbar(ax='NaN'):
#Creates an invisible colourbar to align subplots without colourbars.
#Takes image axis, returns colourbar axis if further edits are required
#cax = InvisibleColourbar(ax[0])

	#If no axis is supplied, use current open axis
	if ax == 'NaN': ax = plt.gca()

	#Create colourbar axis, ideally should 'find' values of existing cbar! 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)

	#Set new cax to zero size and remove ticks.
	try: cax.set_facecolor('none')				#matplotlib v2.x.x method
	except: cax.set_axis_bgcolor('none')		#matplotlib v1.x.x method
	for axis in ['top','bottom','left','right']:
		cax.spines[axis].set_linewidth(0)
	cax.set_xticks([])
	cax.set_yticks([])

	return(cax)
#enddef



#=========================#
#=========================#


def GenerateAxis(Orientation,Isym=ISYMlist[l],PhaseFrames=range(0,IMOVIE_FRAMES[l])):
#Generates a 1D SI [cm] axis for plotting, includes radial symmetry.
#Takes orientation, symmetry and phasecycle options.
#Returns 1D array in units of [cm] or [omega*t/2pi].
#Raxis=GenerateAxis('Radial',Isym=ISYMlist[l])

	#Create axis list and extract the number of phaseframes
	PhaseResolution = len(PhaseFrames)
	axis = list()
	if Orientation == 'Radial':
		if int(Isym) == 1:
			for i in range(-int(R_mesh[l]),int(R_mesh[l])):
				axis.append(i*dr[l])
		#endfor
		elif int(Isym) != 1:
			for i in range(0,int(R_mesh[l])):
				axis.append(i*dr[l])
			#endfor
		#endif
	elif Orientation == 'Axial':
		for i in range(0,int(Z_mesh[l])):
			axis.append(i*dz[l])
		#endfor
	elif Orientation == 'Phase':
		for i in range(0,int(phasecycles*PhaseResolution)):
			axis.append(  (np.pi*(i*2)/PhaseResolution)/(2*np.pi)  )
		#endfor
	#endif
	return(axis)
#enddef


#=========================#
#=========================#


def Normalise(profile,NormFactor=0):
#Takes 1D or 2D array and returns array Normalised to maximum value.
#If NormFactor is defined, array will be Normalised to this instead.
#Returns Normalised image/profile and the max/min normalization factors.
#NormProfile,Min,Max = Normalise(profile,NormFactor=0)

	NormalisedImage = list()
	
	#determine dimensionality of profile and select normaliztion method.
	if isinstance(profile[0], (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 2D array.
		FlatImage = [item for sublist in profile for item in sublist]
		MaxNormFactor,MinNormFactor = max(FlatImage),min(FlatImage)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0
		#endif

		#Normalise 2D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			NormalisedImage.append( [x/NormFactor for x in profile[i]] )
		#endfor
		profile = NormalisedImage
		return(profile,MaxNormFactor,MinNormFactor)

	#Lowest dimention is still list.
	elif isinstance(profile, (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 1D profile.
		MaxNormFactor,MinNormFactor = max(profile),min(profile)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0

		#Normalise 1D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			profile[i] = profile[i]/NormFactor
		#endfor
	#endif

	return(profile,MinNormFactor,MaxNormFactor)
#enddef


#=========================#
#=========================#


def DataExtent(folder=l,aspectratio=image_aspectratio):
#Takes current image datails and returns extent and rotated aspectratio
#If mesh uses symmetry, will double radius extent centered on zero.
#extent,aspectratio = DataExtent(l)

	#Obtain global variables for current folder.
	Isym = float(ISYMlist[folder])
	radius,height = Radius[folder],Height[folder]
	#Rotated Image: [X,Y] = [Height,Radius]
	if image_rotate == True:
		aspectratio = aspectratio[::-1]
		if Isym == 1: 
			extent= [0,height, -radius, radius]
		elif Isym == 0: 
			extent=[0,height, 0,radius]
		#endif

	#Default mesh orientation: [X,Y] = [Radius,Height]
	elif image_rotate == False:
		if Isym == 1: extent = [-radius,radius, 0,height]
		elif Isym == 0: extent=[0,radius, 0,height]
		#endif
	#endif
	return(extent,aspectratio)
#enddef


#=========================#
#=========================#


def ImagePlotter1D(axis,profile,aspectratio,fig=111,ax=111):
#Create figure and plot a 1D graph with associated image plotting requirements.
#Returns plotted axes and figure if new ones were created.
#Else plots to existing figure and returns the image object.
#ImagePlotter1D(Zaxis,Zprofile,image_aspectratio,fig,ax[0]):

	#Generate new figure if required. {kinda hacky...}
	if fig == 111 and ax == 111:
		fig,ax = figure(aspectratio)
	elif fig == 111:
		fig = figure(aspectratio)
	#endif

	#Apply any required numerical changes to the profile.
	if image_logplot == True:
		profile = np.log(profile)
	if image_normalise == True:
		profile = Normalise(profile)[0]
	#endif

	#Plot profile and return.
	im = ax.plot(axis,profile, lw=2)

	try: return(fig,ax,im)
	except: return()
#enddef


#=========================#
#=========================#


def ImagePlotter2D(Image,extent,aspectratio=image_aspectratio,variable='N/A',fig=111,ax=111):
#Create figure and plot a 2D image with associated image plotting requirements.
#Returns plotted image, axes and figure after applying basic data restructuring.
#fig,ax,im,Image = ImagePlotter2D(Image,extent,image_aspectratio,variablelist[l],fig,ax[0])

	#Generate new figure if required. {kinda hacky...}
	if fig == 111 and ax == 111:
		fig,ax = figure(aspectratio)
	elif fig == 111:
		fig = figure(aspectratio)
	#endif

	#Apply image axis-symmetry, with negative values, if required.
	Radial = IsRadialVariable(variable)  
	Image = SymmetryConverter(Image,Radial)

	#Rotate image if required
	if image_rotate == True:
		Image = np.asarray(Image)
		Image = Image.transpose().tolist()
	#endif

	#Apply any required numerical changes to the image.
	if image_logplot == True:
		Image = np.log(Image)
	elif image_normalise == True:
		Image = Normalise(Image)[0]
	#endif

	#Plot image with or without contour plots, (contour scale = 90% of cbar scale)
	if image_plotcontours == True:
		im = ax.contour(Image,extent=extent,origin="lower")
		im.set_clim(CbarMinMax(Image)[0]*0.90,CbarMinMax(Image)[1]*0.90)
		im = ax.imshow(Image,extent=extent,origin="lower")
	else:
		im = ax.imshow(Image,extent=extent,origin="lower")
	#endif
	return(fig,ax,im,Image)
#enddef


#=========================#
#=========================#


def TrendPlotter(ax=plt.gca(),TrendArray=[],Xaxis=[],Marker='o-',NormFactor=0):
#Creates a 1D image from an array of supplied points.
#Image plotted onto existing axes, figure() should be used.
#NormFactor = 0 will Normalise to maximum of given profile.
#TrendPlotter(ax[0],TrendProfiles,Xaxis,'o-',0)

	#Normalise data to provided normalization factor if required.
	if image_normalise == True:
		TrendArray = Normalise(TrendArray,NormFactor)[0]
	#endif

	#Choose how to plot the trends.
	if image_numericaxis == True:			#!!!!!THIS NEEDS GENERALIZING!!!!!#
		#Plot results against number of cells in each mesh for convergence studies.
		numcells = list()
		for l in range(0,numfolders):
			numcells.append(Z_mesh[l]*R_mesh[l])
		#endfor
		Plot = ax.plot(numcells[::-1],TrendArray[::-1], Marker, lw=2)
	else:
		#Plot results against strings pulled from folder names for batch studies.
		Plot = ax.plot(range(0,numfolders),TrendArray, Marker, lw=2)
		if len(xaxisoverride) > 0:
			ax.set_xticks(np.arange(0,numfolders))
			ax.set_xticklabels(xaxisoverride)
		else:
			ax.set_xticks(np.arange(0,numfolders))
			ax.set_xticklabels(Xaxis)
		#endif
	#endif

	return(Plot)
#enddef


#=========================#
#=========================#


def ExtractRadialProfile(Data,process,variable,Profile,Rmesh='NaN',Isym='NaN'):
#Obtains a radial 1D profile at a requested axial location.
#Returns a 1D array for plotting and performs unit conversion.

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 'NaN' or Isym == 'NaN':
		Rmesh,Isym = R_mesh[l],ISYMlist[l]
	#endif

	#Obtain start location for requested data and perform SI conversion.
	ZStart = int(Rmesh*Profile)
	ZEnd = int(Rmesh*(Profile+1))
	
	#Plot lines for each variable at each requested slice, ignoring ones that fail.
	#If mesh is symmetric, copy the data over and make full plot.
	if int(Isym) == 1:
		Zend = len(Data[process])-ZStart
		Zstart = len(Data[process])-ZEnd
		RProfile = Data[process][Zstart:Zend][::-1]
		#If variable is radially symmetric, add negative to symmetry
		if IsRadialVariable(variable) == True:
			for m in range(0,len(RProfile)):
				RProfile.append(-Data[process][Zstart:Zend][m])
			#endfor
		#If variable is axially symmetric, add positive.
		elif IsRadialVariable(variable) == False:
			for m in range(0,len(RProfile)):
				RProfile.append(Data[process][Zstart:Zend][m])
			#endfor
		#endif
		RProfile = RProfile[::-1]	#Reverse index, negative first then positive values.

	#If the data isn't symmetric, just plot as is.
	elif int(Isym) == 0:
		Zend = len(Data[process])-ZStart
		Zstart = len(Data[process])-ZEnd
		RProfile = Data[process][Zstart:Zend]
		
	#endif

	#Convert units if required and plot.
	RProfile = VariableUnitConversion(RProfile,variable)
	return(RProfile)
#enddef


#=========================#
#=========================#


def ExtractAxialProfile(Data,process,variable,lineout,Rmesh=0,Zmesh=0,Isym=0):
#Obtains an axial 1D profile at a requested radial location.
#Returns a 1D array for plotting and performs unit conversion.

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Zmesh == 0 or Isym == 0:
		Rmesh,Zmesh,ISym = R_mesh[l],Z_mesh[l],ISYMlist[l]
	#endif

	#Pull out Z-data point from each radial line of data and list them.
	Zlineout = list()
	for i in range(0,int(Zmesh)):
		datapoint = Rmesh*i + lineout
		try:
			Zlineout.append(Data[int(process)][int(datapoint)])
		except:
			break
		#endtry
	#endfor

	#Convert units if required
	Zlineout = VariableUnitConversion(Zlineout,variable)
	return(Zlineout)
#enddef


#=========================#
#=========================#


def WaveformLoc(location,origintype):
#Convert cell location for use with WaveformExtractor function.
#Returns mesh location based on input string or [R,Z] list.
#Rcell,Zcell = Waveformloc(electrodeloc,'Phase')

	#if data is from movie1, Zaxis is back to front.
	if len(location) == 2 and origintype == 'Phase':
		RlineoutLoc = (Z_mesh[l]-location[1])-1
		ZlineoutLoc = location[0]
	#if data is from TECPLOT2D, Zaxis is correct orientation.
	elif len(location) == 2 and origintype == '2D':
		RlineoutLoc = location[1]
		ZlineoutLoc = location[0]
	else:
		RlineoutLoc = Z_mesh[l]/2
		ZlineoutLoc = 0
	#endif

	return(RlineoutLoc,ZlineoutLoc)
#enddef


#=========================#
#=========================#


def WaveformExtractor(PhaseData,PPOT,waveformlocation=electrodeloc):
#Takes phasedata for current folder and PPOT process number.
#Returns two arrays: VoltageWaveform at electrode location and average waveform.
#VoltageWaveform,WaveformBias = WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[0])

	#Create required lists and extract electrode location.
	VoltageWaveform,WaveformBias = list(),list()
	RLoc = WaveformLoc(waveformlocation,'Phase')[0]
	ZLoc = WaveformLoc(waveformlocation,'Phase')[1]

	#Create voltage waveform for requested integer number of phasecycles
	for i in range(0,int(phasecycles*len(PhaseData))):
		Index = i % len(PhaseData)		#Use modulo index for additional phase cycles
		VoltageWaveform.append(ExtractAxialProfile(PhaseData[Index],PPOT,'PPOT',ZLoc)[RLoc])
	#endfor

	#Calculate time averaged waveform bias, i.e. waveform symmetry.
	for m in range(0,len(VoltageWaveform)):
		WaveformBias.append(sum(VoltageWaveform)/len(VoltageWaveform))
	#endfor
	
	#Calculate maximum positive and negative waveform amplitudes and compute average Vpp
	PositiveAmp,NegativeAmp = max(VoltageWaveform),min(VoltageWaveform)
	PeakToPeakVoltage = abs(PositiveAmp)+abs(NegativeAmp)
	
	return(VoltageWaveform,WaveformBias,[PositiveAmp,NegativeAmp,PeakToPeakVoltage])
#enddef


#=========================#
#=========================#


def TrendAtGivenLocation(probeloc,process,variable):
#Trend analysis for a given point on a 2D Image.
#Takes global 'probeloc' for safety, process and variable.
#Returns two arrays: One is the X-axis to plot against
#Second is the value of variable at location for all simulations.

	#Refresh lists that change per image.
	R,Z = probeloc[0],probeloc[1]
	Trend = list()
	Xaxis = list()

	#For all simulation folders.
	for l in range(0,numfolders):
	
		#Extract image with given process and variable name.
		Image = ImageExtractor2D(Data[l][process],variable,R_mesh[l],Z_mesh[l])

		#Update X-axis with folder information.
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#TREND ANALYSIS - Value at Given Location Comparison.
		try: Trend.append( Image[Z][R] )
		except: Trend.append(float('NaN'))

		#Display Min/Max value trends to terminal if requested.
		if print_generaltrends == True:
			Location = '('+str(round(R*dr[l],1))+'cm,'+str(round(Z*dz[l],1))+'cm)'
			print( FolderNameTrimmer(Dirlist[l]))
			print( str(variable)+' @ '+Location+':', round(Trend[-1], 5))
		if print_generaltrends == True and l == numfolders-1:
			print( '')
		#endif
	#endfor

	#Normalise to maximum value in each profile if required.
	if image_normalise == True:
		Trend,Min,Max = Normalise(Image)
	#endif

	return(Xaxis,Trend)
#enddef


#=========================#
#=========================#


def MinMaxTrends(lineout,Orientation,process):
#General trend plotting function for use with multiple folders.
#Takes a lineout location and orientation string as input.
#And returns the maximum and minimum values For all folders to be analysed.
#With an associated X-axis composed of the trimmed folder names.

	#Refresh lists that change per profile.
	MaxValueTrend, MinValueTrend = list(), list()
	Xaxis = list()
	p = process

	#For each folder in the directory.
	for l in range(0,numfolders):

		#Create and correct processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])
		processlist,Variablelist = VariableInterpolator(processlist,Variablelist,Comparisonlist)

		#Update X-axis with folder information.
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Obtain radial and axial profiles for further processing.
		if Orientation == 'Radial':
			try: Profile = ExtractRadialProfile(Data[l],processlist[p],Variablelist[p],lineout,R_mesh[l],ISYMlist[l])
			except: Profile = float('NaN')
			#endtry
		elif Orientation == 'Axial':
			try: Profile = ExtractAxialProfile(Data[l],processlist[p],Variablelist[p],lineout,R_mesh[l],Z_mesh[l],ISYMlist[l])
			except: Profile = float('NaN')
			#endtry
		#endif

		#TREND ANALYSIS - Maximum/Minimum Value Comparison.
		try: MaxValueTrend.append(max(Profile))
		except: MaxValueTrend.append(float('NaN'))
		try: MinValueTrend.append(min(Profile))
		except: MinValueTrend.append(float('NaN'))
		#endtry

		#Display Min/Max value trends to terminal if requested.
		if print_generaltrends == True:
			VariableName = VariableLabelMaker(Variablelist)[p]
			print( FolderNameTrimmer(Dirlist[l]))
			print( VariableName+' '+Orientation+'Maximum: ', round(max(MaxValueTrend), 5))
			print( VariableName+' '+Orientation+'Minimum: ', round(min(MinValueTrend), 5))
		if print_generaltrends == True and l == numfolders-1:
			print( '')
		#endif
	#endfor

	#Normalise to maximum value in each profile if required.
	if image_normalise == True:
		MaxValueTrend,MaxMin,MaxMax = Normalise(MaxValueTrend)
		MinValueTrend,MinMin,MinMax = Normalise(MinValueTrend)
	#endif

	return(Xaxis,MaxValueTrend,MinValueTrend)
#enddef


#=========================#
#=========================#

#TREND ANALYSIS - Speed of Sound
def CalcSoundSpeed(NeutralDensity,Pressure,Dimension='2D'):
#Calculates local sound speed via Newton-Laplace equation
#Takes appropriate neutral density [m-3] and pressure [Torr] (0D,1D or 2D)
#Returns same dimensionality array of sound speeds in m/s
#SoundSpeed = CalcSoundSpeed(ArgonDensity,Pressure,Dimension='2D')

	#Initiate required lists and set atomic values
	SoundSpeedArray = list()
	AdiabaticIndex = 5.0/3.0		#Assumes diatomic species (Hydrogen,Helium,Argon, etc...)
	AtomicMass = 39.948*1.66E-27	#[Kg]		#NB HARDCODED FOR ARGON

	#For 0D values:
	if Dimension == '0D':
		ElasticityModulus = AdiabaticIndex*Pressure*133.33					#[Pa] = [kg m-1 s-2]
		MassDensity = NeutralDensity*AtomicMass								#[kg m-3]
		
		#Calculate local sound speed via Newton-Laplace equation
		try: SoundSpeedArray = np.sqrt( ElasticityModulus/MassDensity )		#[m/s]
		except: SoundSpeedArray = np.nan
	#endif
	
	#For 1D arrays:
	if Dimension == '1D':
		for i in range(0,len(NeutralDensity)):
			ElasticityModulus = AdiabaticIndex*Pressure[i]*133.33			#[Pa] = [kg m-1 s-2]
			MassDensity = NeutralDensity[i]*AtomicMass						#[kg m-3]

			#Calculate local sound speed via Newton-Laplace equation
			try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )		#[m/s]
			except: SoundSpeed = np.nan
			SoundSpeedArray.append( SoundSpeed )
		#endfor
	#endif

	#For 2D arrays:
	if Dimension == '2D':
		for i in range(0,len(NeutralDensity)):
			SoundSpeedArray.append(list())
			for j in range(0,len(NeutralDensity[i])):
				ElasticityModulus = AdiabaticIndex*Pressure[i][j]*133.33	#[Pa]
				MassDensity = NeutralDensity[i][j]*AtomicMass				#[kg m-3]

				#Calculate local sound speed via Newton-Laplace equation
				try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )	#[m/s]
				except: SoundSpeed = np.nan
				SoundSpeedArray[i].append( SoundSpeed )
			#endfor
		#endfor
	#endif

	return(SoundSpeedArray)
#enddef


#=========================#
#=========================#

#TREND ANALYSIS - DCbias
def DCbiasMagnitude(PPOTlineout):
#Takes a PPOT profile and calcuates DCbias via difference in voltage drop.
#Can identify DC-bias for parallel plate discharges and dielectric discharges.

	#Identify if radial or axial.
	if len(PPOTlineout) == Z_mesh[l]:
		electrodelocation = WaveformLoc(electrodeloc,'2D')[1]
	elif len(PPOTlineout) in [R_mesh[l],R_mesh[l]*2]:
		electrodelocation = WaveformLoc(electrodeloc,'2D')[0]
	#endif

	#Identify Min/Max Potential magnitudes and location of max potential.
	MinPPOT,MaxPPOT = min(PPOTlineout),max(PPOTlineout)
	MaxIndex = np.argmax(PPOTlineout)

	#Split PPOT profile into each sheath, pre and post max potential
	PreIndex = PPOTlineout[:MaxIndex]
	PostIndex = PPOTlineout[MaxIndex:]


	##=========================================##

	#Metals have flat PPOT profiles, dielectric/plasma have gradients.
	MetalIndices = list([0])
	DielectricIndices = list()
	for i in range(0,len(PPOTlineout)-1):
		if PPOTlineout[i] == PPOTlineout[i+1]:
			MetalIndices.append(i)
		elif PPOTlineout[i] != PPOTlineout[i+1]:
			DielectricIndices.append(i)
		#endif
	#endfor
	MetalIndices.append(len(PPOTlineout)-1)

	#Grounded metal will have a PPOT of zero -- ##INCORRECT IF DC-BIAS == int(0.0)##
	GMetalIndices = list()
	for i in range(0,len(MetalIndices)):
		if PPOTlineout[MetalIndices[i]] == 0:
			GMetalIndices.append(MetalIndices[i])
		#endif
	#endfor

	#Any metal that is not grounded will be powered -- ##BAD ASSUMPTION FOR ALL MESHES##
	PMetalIndices = list()
	for i in range(0,len(MetalIndices)):
		if MetalIndices[i] not in GMetalIndices:
			PMetalIndices.append(MetalIndices[i])
		#endif
	#endfor

	##=========================================##


	#Identify voltage drop through each sheath from max potential.
	try: PreIndexVoltageDrop = MaxPPOT - min(PreIndex)
	except: PreIndexVoltageDrop = MaxPPOT
	#endtry
	try: PostIndexVoltageDrop = MaxPPOT - min(PostIndex)
	except: PostIndexVoltageDrop = MaxPPOT
	#endtry

	#Minimum voltage is not one of the electrodes - "Dielectric Discharge"
	if min(PPOTlineout) not in [ PPOTlineout[0],PPOTlineout[-1] ]:
		try: DCbias = MinPPOT
		except: DCbias = MaxPPOT
		#endtry

	#Minimum voltage is one of the electrodes - "Parallel Plate Discharge"
	else:
		try: DCbias = PPOTlineout[PMetalIndices[0]]
		except: DCbias = PreIndexVoltageDrop - PostIndexVoltageDrop
	#endif

	if IDEBUG == True:
		X1 = range(0,len(PreIndex))
		X2 = range(len(PreIndex),len(PPOTlineout))

		plt.plot(X1,PreIndex, lw=2)
		plt.plot(X2,PostIndex, lw=2)
		plt.plot(np.argmax(PPOTlineout),max(PPOTlineout), 'go',  ms=12)
		for i in range(0,len(GMetalIndices)):
			plt.plot(GMetalIndices[i],PPOTlineout[GMetalIndices[i]], 'ko',  ms=12)
		#endfor
		for i in range(0,len(PMetalIndices)):
			plt.plot(PMetalIndices[i],PPOTlineout[PMetalIndices[i]], 'ro',  ms=12)
		#endfor

		plt.xlabel('Cell Number')
		plt.ylabel('Voltage [V]')
		plt.legend(['PreBulk','PostBulk','Bulk'])
		plt.title('DCBIAS_DEBUG'+str(l+electrodelocation))
		plt.savefig(DirTrends+'DCBIAS_DEBUG'+str(l+electrodelocation)+'.png')
		plt.close('all')
	#endif

	return DCbias
#enddef


#=========================#
#=========================#


#BRINKMANN SHEATH WIDTH CALCULATOR
def CalcSheathExtent(folder=l,Orientation='Axial',Phase='NaN',Ne=list(),Ni=list()):
#Calculates Brinkmann sheath width assuming Child-Langmuir conditions.
#Calculation Methods: 'AbsDensity', 'IntDensity'
#Takes current folder, current axis, movie1 Phase and sheath calc method.
#Returns array of sheath distances from symmetry boundary (or from origin)
#Sx,SxAxis = CalcSheathExtent(folder=l,Phase=moviephaselist[k])

	#Initiate required data storage lists
	NPos,NNeg = list(),list()
	SxAxis,Sx = list(),list()	

	#Import global sheath calculation method and charged particle species names
	SheathMethod=GlobSheathMethod
	if len(SheathIonSpecies) == 0:
		global PosSpecies
		global NegSpecies
	#Force a single sheath species - Legacy Code or for testing purposes
	elif len(SheathIonSpecies) > 0:
		PosSpecies = SheathIonSpecies
		NegSpecies = []
	#endif

	#Return NaN sheath array if no appropriate species are detected
	if len(PosSpecies)+len(NegSpecies) == 0:
		SxAxis = GenerateAxis(Orientation,Isym=ISYMlist[l])
		Sx = np.empty(len(SxAxis))
		[np.nan for x in Sx]
		
		return(Sx,SxAxis)
	#endif

	#Identify charged species and alter names to suit TECPLOT2D nomenclature
	for i in range(0,len(PosSpecies)): PosSpecies[i] = PosSpecies[i] = PosSpecies[i].replace('^','+')
	for i in range(0,len(NegSpecies)): NegSpecies[i] = NegSpecies[i] = NegSpecies[i].replace('^','-')
	if 'E' in NegSpecies: NegSpecies.remove('E')			#Might Cause An Issue With Global!!!
	
	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	#ISSUE WITH THIS REGARDING THE DATA READ-IN
	#PREVIOUS VERSION SENDS Ne and Ni EXPLICITLY INTO FUNCTION, 
	#IDEALLY WOULD HAVE ALL OF THIS RUN ONCE, EXTRACTING THE DATA FROM RAW_PHASEDATA
	#PASSING THE REQUIRED Ne and Neff INTO THE 2ND PART OF THE FUNCTION.
	#IF Ne, Neff don't exist: Run phase-extraction-function
	#Else: run old code below, using global values

	#Obtain current folder ion and electron densities if not already supplied.
	#Default to 2D data format.
#	if Phase == 'NaN' and len(Ne) == 0:
		#Obtain electron density and extract 2D image for further processing.
#		Eproc = VariableEnumerator(['E'],rawdata_2D[folder],header_2Dlist[folder])[0][0]
#		Ne = ImageExtractor2D( Data[folder][Eproc] )

		#Obtain all positive and negative ion densities and extract 2D images for further processing
#		PosSpeciesproc = VariableEnumerator(PosSpecies,rawdata_2D[folder],header_2Dlist[folder])[0]
#		for i in range(0,len(PosSpeciesproc)): 
#			NPos.append( ImageExtractor2D(Data[folder][PosSpeciesproc[i]]) )
		#endfor
#		NegSpeciesproc = VariableEnumerator(NegSpecies,rawdata_2D[folder],header_2Dlist[folder])[0]
#		for i in range(0,len(NegSpeciesproc)): 
#			NNeg.append( ImageExtractor2D(Data[folder][NegSpeciesproc[i]]) )	
		#endfor

	#If phase is supplied, use phase data format.  (Proc=Proc-2 to skip R,Z data in phase data)
#	elif Phase != 'NaN' and len(Ne) == 0:
#		Eproc = VariableEnumerator(['E'],rawdata_phasemovie[folder],header_phasemovie[folder])[0][0]
#		Ne = ImageExtractor2D( PhaseMovieData[folder][Phase][Eproc-2] ) 

		#Obtain all positive and negative ion densities and extract 2D images for further processing
#		PosSpeciesproc=VariableEnumerator(PosSpecies,rawdata_phasemovie[folder],header_phasemovie[folder])[0]
#		for i in range(0,len(PosSpeciesproc)): 
#			NPos.append( ImageExtractor2D(PhaseMovieData[folder][Phase][PosSpeciesproc[i]-2]) )
		#endfor
#		NegSpeciesproc=VariableEnumerator(NegSpecies,rawdata_phasemovie[folder],header_phasemovie[folder])[0]
#		for i in range(0,len(NegSpeciesproc)): 
#			NNeg.append( ImageExtractor2D(PhaseMovieData[folder][Phase][NegSpeciesproc[i]-2]) )
		#endfor

	#If specific electron and ion species densities are supplied, use those
#	elif len(Ne) > 0 or len(Ni) > 0:
#		Ne = ImageExtractor2D( Ne ) 		#Ne[i][j]
#		NPos = [ ImageExtractor2D( Ni ) ]	#Put in array []    (NPos[k][i][j])
#		PosSpeciesproc = ['Ion+']			#Set length to 1
#		NegSpeciesproc = []					#Set length to 0
	#endif

	#Combine 2D images of all positive ion species densities and all negative ion species densitiies
#	NPos = [[sum(x) for x in zip(NPos[0][i],NPos[1][i])] for i in range(len(NPos[0]))]
#	NNeg = [[sum(x) for x in zip(NNeg[0][i],NNeg[1][i])] for i in range(len(NNeg[0]))]
#							HOW TO ZIP ARBITARY NUMBER OF ARRAYS?
#	TotNPos = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
#	for i in range(0,len(TotNPos)):
#		for j in range(0,len(TotNPos[0])):
#			for k in range(0,len(PosSpeciesproc)): TotNPos[i][j] += NPos[k][i][j]
			#endfor
		#endfor
	#endfor
#	TotNNeg = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
#	for i in range(0,len(TotNNeg)):
#		for j in range(0,len(TotNNeg[0])):
#			for k in range(0,len(NegSpeciesproc)): TotNNeg[i][j] += NNeg[k][i][j]
			#endfor
		#endfor
	#endfor

	#Determine effective positive ion density as: Neff = sum(Total NPos)-sum(Total NNeg)
#	Neff = np.zeros( (len(Ne),len(Ne[0])) ).tolist()
#	for i in range(0,len(Neff)):
#		for j in range(0,len(Neff[0])):
#			Neff[i][j] = TotNPos[i][j] - TotNNeg[i][j]
		#endfor
	#endfor

	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	#!!! OLD METHOD !!!
	#Obtain current folder ion and electron densities if not already supplied.
	#Default to 2D data format.
	if Phase == 'NaN' and len(Ne) == 0:
		IONproc = VariableEnumerator(PosSpecies,rawdata_2D[folder],header_2Dlist[folder])[0][0]
		Eproc = VariableEnumerator(['E'],rawdata_2D[folder],header_2Dlist[folder])[0][0]
		Ne,Ni = Data[folder][Eproc], Data[folder][IONproc]
	#If phase is supplied, use phase data format.
	elif Phase != 'NaN' and len(Ne) == 0:
		IONproc = VariableEnumerator(PosSpecies,rawdata_phasemovie[folder],header_phasemovie[folder])[0][0]
		Eproc = VariableEnumerator(['E'],rawdata_phasemovie[folder],header_phasemovie[folder])[0][0]
		IONproc,Eproc = IONproc-2, Eproc-2		#Skip R,Z data inputs in phase data.
		Ne,Ni = PhaseMovieData[folder][Phase][Eproc], PhaseMovieData[folder][Phase][IONproc]
	#endif
	#Extract 2D image for further processing.
	Ne,Neff = ImageExtractor2D(Ne),ImageExtractor2D(Ni)
	#!!! OLD METHOD !!!

	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###
	#Axial sheath array (Sx) is calculated exmploying radial integrations for all axial locations
	#Radial sheath array (Sx) is calculated employing axial integrations for all radial locations
	### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###
	if Orientation == 'Axial':

		#Determine sheath edge through integration of charge density:
		if SheathMethod == 'IntDensity':
			#Sheath extension: integral_(R0->Rwall) ne dR == integral_(Rwall->R0) ni dR
			for i in range(0,len(Neff)):
				#Define wall radius to integrate ions into bulk from.
				for j in range(0,len(Neff[i])):

					#if ion density drops to zero, we've hit a material surface.
					if Neff[i][j] == 0.0 and j == 0:		
						RadialPlasmaExtent = 0
						break
					elif Neff[i][j] == 0.0 and j > 0:
						RadialPlasmaExtent = j-1
						break
					#endif
				#endfor
#				RadialPlasmaExtent = len(Neff[i])	#DEBUG OPTION: Sets RadialPlasmaExtent to max for all Z
				
				#No plasma, all radii are solids, append 'nan' to avoid plotting.
				if RadialPlasmaExtent == 0: 
					Sx.append(np.nan)								#[cm]	
				#If non-zero plasma extent, determine radial cell satisfying Brinkmann Criterion
				elif RadialPlasmaExtent > 0:
					#Refresh sums after every radial profile.
					Neff_sum,Ne_sum = 0.0,0.0
					for j in range(0,RadialPlasmaExtent):
						#Sum density radially for ions and electrons.
						reversed_j = RadialPlasmaExtent-j-1
						Neff_sum += Neff[i][j]		#Sum from R=wall to R=0	[reversed_j] 	####FUDGED####
						Ne_sum += Ne[i][j]			#Sum from R=0 to R=wall [j] 	 		####FUDGED####

						#If ion sum is greater than electron, sheath has begun.
						if Neff_sum/Ne_sum >= 1.0:
							Sx.append(j*dr[l])						#[cm]								
							break
						#If no sheath found within plasma region, append wall location (i.e. Sx=Rwall)
						if j == (RadialPlasmaExtent-1):
							Sx.append((RadialPlasmaExtent+1)*dr[l])	#[cm]
#							Sx.append(np.nan)						#[cm] Nice Plots/Breaks Statistics
						#endif
					#endfor
				#endif
			#endfor

		#==========#

		#Determine sheath edge by 'instantaneous' charge density:
		elif SheathMethod == 'AbsDensity':
			#Sheath extension: ni @R >= ne @R, simplified model.
			for Z in range(0,len(Neff)):
				for R in range(0,len(Neff[Z])):
					#Sheath starts when ion density exceeds electron density.
					if Neff[Z][R]/Ne[Z][R] >= 1.0:
						Sx.append(R*dr[l])
						break
					#If no sheath found, append 'NaN' to avoid plotting.
					if R == (len(Neff[Z])-1):
#						Sx.append(0.0)
						Sx.append(np.nan)
					#endif
				#endfor
			#endfor
		#endif

		#Extract axis to plot sheath against
		SxAxis = GenerateAxis(Orientation,Isym=ISYMlist[l])
	#endif

	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	#NEED TO APPLY RADIAL METHOD!!! FOR NOW THIS IS SET TO ZERO EVERYWHERE#
	if Orientation == 'Radial':
		#Sheath extension: integral_(R0->Rwall) ne dR == integral_(Rwall->R0) ni dR
		for i in range(0,len(Neff)):
			Sx.append(np.nan)
		#endfor
		
		#Extract axis to plot sheath against
		SxAxis = GenerateAxis(Orientation,Isym=ISYMlist[l])
	#endif

	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	#Print sheath characteristics if requested.
	if print_sheath == True:
		#Determine location of interest (loc) for diagnostic
		if Orientation == 'Axial': loc = electrodeloc[0]		#Radial point of interest
		if Orientation == 'Radial': loc = electrodeloc[1]		#Axial point of interest

		#Obtain SheathWidth at electrodeloc
		try: SheathWidth = round(Sx[loc],3)
		except: SheathWidth = 0.0
		print( 'Simulation:', Dirlist[folder])
		print( 'Sheath Location:',SheathWidth*10, 'mm')
		print( 'Sheath Extent:',((sourcewidth[0]*dr[l])-SheathWidth)*10, 'mm')
	#endif

	#Return sheath axis and sheath extent
	return(Sx,SxAxis)
#enddef


#=========================#
#=========================#


def PlotSheathExtent(SxAxis,Sx,ax=plt.gca(),ISymmetry=0):
#PlotSheathExtent(SxAxis,Sx,ax[0],ISYMlist[l])

	#Determine mesh symmetry and create symmetric mesh if required
	if int(ISymmetry) == 1:
		SymSx = list()
		#Symmetric sheath creation
		for i in range(0,len(Sx)): 
			try: SymSx.append(-Sx[i])
			except: SymSx.append(np.nan)
		#endfor
	#endif

	#Plot sheath extent from origin with or without mesh symmetry
	if image_plotsheath == True and int(ISymmetry) == 0:
		ax.plot(SxAxis,Sx, 'w--', lw=2)
	elif image_plotsheath == True and int(ISymmetry) == 1:
		ax.plot(SxAxis,Sx, 'w--', lw=2)
		ax.plot(SxAxis,SymSx, 'w--', lw=2)
	#endif
	
	return()
#enddef


#=========================#
#=========================#

#====================================================================#
#====================================================================#













































#====================================================================#
				  #IMAGE PLOTTING AND DATA ANALYSIS#
#====================================================================#

#====================================================================#
			  #STEADY STATE 2D FIGURES -- CYCLE AVERAGED#
#====================================================================#

#Generate and save image of required variable for given mesh size.
if savefig_plot2D == True:

	for l in range(0,numfolders):
		#Create new folder to keep output plots.
		Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')

		#Create processlist for each folder as required.
		processlist,variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

		#Setting the radial ticks for beauty purposes ~HACKY~
		R = Radius[l]
		SymTicks = [round(-R, 1), round(-R/2, 1), 0, round(R/2, 1), round(R, 1)]
		NoSymTicks = [0, round(R/2, 1), round(R, 1)]

		#Reshape specific part of 1D Data array into 2D image for plotting.
		for k in tqdm(range(0,len(processlist))):

			#Extract full 2D image for further processing.
			Image = ImageExtractor2D(Data[l][processlist[k]],variablelist[k])
			Sx,SxAxis = CalcSheathExtent(folder=l)

			#Generate and rotate figure as requested.
			extent,aspectratio = DataExtent(l)
			fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,variablelist[k])
			PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l])
			
			#Overlay location of 1D profiles if requested, adjusting for image rotation.
			if image_plotoverlay == True:
				for j in range(0,len(radialprofiles)):
					X1,X2 = extent[0],extent[1]
					Y1,Y2 = radialprofiles[j]*dz[l],radialprofiles[j]*dz[l]
					if image_rotate == True: X1,X2,Y1,Y2 = Y1,Y2,X1,X2
					ax.plot((X1,X2),(Y1,Y2),'k--',lw=2)
				#endfor
				for j in range(0,len(axialprofiles)):
					X1,X2 = axialprofiles[j]*dr[l],axialprofiles[j]*dr[l]
					Y1,Y2 = extent[2],extent[3]
					if image_rotate == True: X1,X2,Y1,Y2 = Y1,Y2,X1,X2
					ax.plot((X1,X2),(Y1,Y2),'k--',lw=2)
				#endfor
			#endif

			#Define image beautification variables.
			if image_rotate == True:
				Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
			elif image_rotate == False:
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				plt.gca().invert_yaxis()
			#endif

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = '2D Steady State Plot of '+variablelist[k]+' for \n'+Dirlist[l][2:-1]   #generic titles
			#Add Colourbar (Axis, Label, Bins)
			label,bins = VariableLabelMaker(variablelist),5
			cax = Colourbar(ax,label[k],bins,Lim=CbarMinMax(Image))
			#Finalize image
			ImageOptions(fig,ax,Xlabel,Ylabel,Title)

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(Dir2Dplots, '2Dplots_Data')
				WriteDataToFile(Image, DirWrite+variablelist[k])
				if image_plotsheath == True and k == len(processlist)-1: WriteDataToFile(Sx, DirWrite+'Sx-EXT')
			#endif

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot '+variablelist[k]+ext)
			plt.close('all')
		#endfor
	#endfor

	print('-------------------------------------')
	print('# 2D Steady-State Processing Complete')
	print('-------------------------------------')
#endif



#====================================================================#
			#CONVERGENCE CHECKING MOVIES -- ITERATION BASED#
#====================================================================#

#Plot 2D images at different iterations towards convergence from movie_icp.
if savefig_convergence == True:

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create new folder to keep convergence variable folders in.
		DirConvergence = CreateNewFolder(Dirlist[l],'Convergence/')

		#Create processlist for each folder as required.
		processlist,variablelist = VariableEnumerator(IterVariables,rawdata_itermovie[l],header_itermovie[l])
		#Skip over the R and Z processes as they are not saved properly in iterdata.
		for i in range(0,len(processlist)): processlist[i] = processlist[i]-2

		#Extract saved iteration strings and create list for mean convergence trends
		ConvergenceTrends,IterArray = list(),list()
		for i in range(0,len(MovieIterlist[l])):
			Iter = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))		#List of digits within Iter
			Iter =  int(''.join(Iter))									#Join digits into single integer
			IterArray.append(Iter)										#Append to list of integers
		#endfor

		#for all variables requested by the user.
		for i in tqdm(range(0,len(processlist))):

			#Create new folder to keep output plots.
			DirMovieplots = CreateNewFolder(DirConvergence,variablelist[i]+'_2DConvergence/')
			#Append new list to convergenceTrends for each variable.
			ConvergenceTrends.append(list())

			#Create empty image array based on mesh size and symmetry options.
			try: 
				numrows = len(IterMovieData[l][0][0])/R_mesh[l]
				Image = np.zeros([int(numrows),int(R_mesh[l])])
			except: 
				print( 'No Iteration Data Found For '+Dirlist[l])
				break
			#endtry
			
			#Reshape specific part of 1D Data array into 2D image for plotting.
			if QuickConverge == False:
				for k in range(0,len(MovieIterlist[l])):
					#Extract full 2D image for further processing.
					Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
					Sx,SxAxis = CalcSheathExtent(folder=l)		#NB: CURRENTLY USES TECPLOT2D DATA

					#Take MEAN value of image for general convergence trend.
					ImageMean = sum(Image.flatten())/len(Image.flatten())
					ConvergenceTrends[-1] = np.append(ConvergenceTrends[-1], ImageMean)

					#Generate and rotate figure as requested.
					extent,aspectratio = DataExtent(l)
					fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,variablelist[i])
					PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l])
					
					#Define image axis labels.
					if image_rotate == True:
						Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
					elif image_rotate == False:
						Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
						plt.gca().invert_yaxis()
					#endif

					#Add title, legends, Colourbar (Axis, Label, Bins), etc...
					Title = str(MovieIterlist[l][k])
					label,bins = VariableLabelMaker(variablelist),5
					cax = Colourbar(ax,label[i],bins,Lim=CbarMinMax(Image))
					ImageOptions(fig,ax,Xlabel,Ylabel,Title)

					#Save to seperate folders inside simulation folder.
					#N.B. zfill(3) Asumes Iter never exceeds 999 (i.e. max(IterArray) < 1e4)
					savefig(DirMovieplots+variablelist[i]+'_'+str(IterArray[k]).zfill(3)+ext)
					plt.close('all')
				#endfor

				#Create .mp4 movie from completed images.
				Prefix = FolderNameTrimmer(Dirlist[l])
				MakeMovie(DirMovieplots,Prefix+'_'+variablelist[i])


			#Extract images for convergence trend plotting, but do not plot images.
			if QuickConverge == True:
				for k in range(0,len(MovieIterlist[l])):
					#Extract full 2D image for further processing.
					Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
					#Take MEAN value of image for general convergence trend.
					ImageMean = sum(Image.flatten())/len(Image.flatten())
					ConvergenceTrends[-1] = np.append(ConvergenceTrends[-1], ImageMean)
				#endfor
			#endif
		#endfor

		#=================#

		#Plot a convergence check for all variables in each folder.
		Legend = VariableLabelMaker(variablelist)
		fig, ax = plt.subplots(1, figsize=(10,10))

		#Obtain final iteration values for use as normalisation factors
		FinalIterationValues = list()
		for i in range(0,len(processlist)):
			Image = ImageExtractor2D(IterMovieData[l][-1][processlist[i]],variablelist[i])
			FinalIterationValues.append( sum(Image.flatten())/len(Image.flatten()) )
		#endfor

		#Normalise and plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(ConvergenceTrends)):           
			ConvergenceTrends[i] = Normalise(ConvergenceTrends[i],NormFactor=FinalIterationValues[i])[0]
			ax.plot(IterArray,ConvergenceTrends[i], lw=2)
		#endfor
		#Calculate image Ylims :: also "Limit the limits!" (avoids super zoomed out figures)
		Limits = [min(np.asarray(ConvergenceTrends).flatten()),max(np.asarray(ConvergenceTrends).flatten())]
		if abs(Limits[1]) > 2*Limits[0]: Limits[1] = +4.0
		if abs(Limits[0]) > 2*Limits[1]: Limits[0] = -4.0

		#Image plotting details.
		Title = 'Convergence of '+str(variablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation Iteration','Normalised Mesh-Average Value'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(Limits[0],Limits[1])
		plt.tight_layout()

		#Write data to ASCII files if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirConvergence, 'Convergence_Data')
			SaveString = FolderNameTrimmer(Dirlist[l])+'_ConvergenceData'
			WriteDataToFile(variablelist+['\n'], DirASCII+SaveString, 'w')
			#endif
			for i in range(0,len(ConvergenceTrends)):
				WriteDataToFile('%s \n' %(ConvergenceTrends[i]), DirASCII+SaveString, 'a')
			#endfor
		#endif

		#Print convergence data to terminal if required
		print('')
		print('Percentage Variation At Final Iteration:')
		for i in range(0,len(ConvergenceTrends)):
			ConvergenceFraction = 1-abs( ConvergenceTrends[i][-1]/ConvergenceTrends[i][-2] )
			ConvergencePercentage = round( (ConvergenceFraction*100), 6)
			print( variablelist[i], '\t', ConvergencePercentage, '%')
		#endfor

		#Save figure.
		savefig(DirConvergence+FolderNameTrimmer(Dirlist[l])+'_Convergence'+ext)
		plt.close('all')
	#endfor

	print('------------------------------------')
	print('# 2D Convergence Processing Complete')
	print('------------------------------------')
#endif

#=====================================================================#
#=====================================================================#















































#====================================================================#
				#PROFILE PLOTTING AND DATA ANALYSIS#
#====================================================================#

#====================================================================#
			#AXIAL AND RADIAL PROFILES FROM SINGLE FOLDERS#
#====================================================================#

#Generate and save lineouts of requested variables for given location.
if savefig_monoprofiles == True:

	for l in range(0,numfolders):

		#Create processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

		#Generate SI scale axes for lineout plots and refresh legend.
		Raxis = GenerateAxis('Radial',ISYMlist[l])
		Zaxis = GenerateAxis('Axial',ISYMlist[l])
		Legendlist = list()

		#Generate the radial (horizontal) lineouts for a specific height.
		if len(radialprofiles) > 0:
			#Create folder to keep output plots.
			DirRlineouts = CreateNewFolder(Dirlist[l],'Radial_Profiles/')

			#Loop over all required variables and requested profile locations.
			for i in tqdm(range(0,len(processlist))):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(radialprofiles)):
					#Update legend with location of each lineout.
					if len(Legendlist) < len(radialprofiles):
						Legendlist.append('Z='+str(round((radialprofiles[j])*dz[l], 2))+' cm')
					#endif

					#Plot all requested radial lines on single image per variable.
					Rlineout=ExtractRadialProfile(Data[l],processlist[i],Variablelist[i],radialprofiles[j])
					#Plot lines for each variable at each requested slice.
					ImagePlotter1D(Raxis,Rlineout,image_aspectratio,fig,ax)

					#Write data to ASCII files if requested.
					if write_ASCII == True:
						SaveString = '_R='+str(round((radialprofiles[j])*dz[l], 2))+'cm'
						DirWrite = CreateNewFolder(DirRlineouts, 'Radial_Data')
						WriteDataToFile([Raxis,Rlineout], DirWrite+Variablelist[i]+SaveString)
					#endif
				#endfor

				#Apply image options and axis labels.
				Title = 'Radial Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Radial Distance R [cm]',VariableLabelMaker(Variablelist)[i]
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

				#Save profiles in previously created folder.
				plt.savefig(DirRlineouts+'1D_Radial_'+Variablelist[i]+'_Profiles'+ext)
				plt.close(fig)
			#endfor
			plt.close('all')
		#endif

#===================##===================#

		#Generate the vertical (height) lineouts for a given radius.
		if len(axialprofiles) > 0:
			#Create folder to keep output plots.
			DirZlineouts = CreateNewFolder(Dirlist[l],'Axial_Profiles/')
			Legendlist = list()

			#Collect and plot required data.
			for i in tqdm(range(0,len(processlist))):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(axialprofiles)):
					#Perform SI conversion and save to legend.
					if len(Legendlist) < len(axialprofiles):
						Legendlist.append('R='+str(round(axialprofiles[j]*dr[l], 2))+' cm')
					#endif

					#Plot all requested radial lines on single image per variable.
					Zlineout=ExtractAxialProfile(Data[l],processlist[i],Variablelist[i],axialprofiles[j])
					#Plot lines for each variable at each requested slice.
					ImagePlotter1D(Zaxis,Zlineout[::-1],image_aspectratio,fig,ax)

					#Write data to ASCII files if requested.
					if write_ASCII == True:
						SaveString = '_Z='+str(round((axialprofiles[j])*dr[l], 2))+'cm'
						DirWrite = CreateNewFolder(DirZlineouts, 'Axial_Data')
						WriteDataToFile([Zaxis,Zlineout[::-1]], DirWrite+Variablelist[i]+SaveString)
					#endif
				#endfor

				#Apply image options and axis labels.
				Title = 'Height Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Axial Distance Z [cm]',VariableLabelMaker(Variablelist)[i]
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

				#Save profiles in previously created folder.
				plt.savefig(DirZlineouts+'1D_Axial_'+Variablelist[i]+'_Profiles'+ext)
				plt.close(fig)
			#endfor
			plt.close('all')
		#endif
	#endfor

	print('--------------------------')
	print('# Single Profiles Complete')
	print('--------------------------')
#endif



##====================================================================#
#			#COMPARITIVE PROFILES FROM MULTI-FOLDERS#
##====================================================================#

#Plot comparitive profiles for each variable between folders.
if savefig_compareprofiles == True:

	#Create folder to keep output plots.
	DirComparisons = CreateNewFolder(os.getcwd(),'/1D Comparisons')

	#Generate SI scale axes for lineout plots.
	Raxis = GenerateAxis('Radial',ISYMlist[l])
	Zaxis = GenerateAxis('Axial',ISYMlist[l])

	#Perform radial (horizontal) profile comparisons
	for j in range(0,len(radialprofiles)):

		#Create new folder for each axial or radial slice.
		ProfileFolder = 'Z='+str(round((radialprofiles[j])*dz[l], 2))+'cm'
		DirProfile = CreateNewFolder(DirComparisons,ProfileFolder)

		#For each requested comparison variable.
		for k in tqdm(range(0,len(Variables))):

			#Loop escape if variables that do not exist have been requested.
			if k >= 1 and k > len(Variablelist)-1:
				break
			#endif

			#Create fig of desired size and refresh legendlist.
			fig,ax = figure(image_aspectratio,1)
			Legendlist = list()

			#For each folder in the directory.
			for l in range(0,numfolders):

				#Create processlist for each folder as required.
				processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

				#Correct processlist for folders containing different icp.dat.
				processlist,Variablelist = VariableInterpolator(processlist,Variablelist,Comparisonlist)

				#Update legend with folder information.
				Legendlist.append( FolderNameTrimmer(Dirlist[l]) )
				Ylabels = VariableLabelMaker(Variablelist)

				#Plot all radial profiles for all variables in one folder.
				RProfile = ExtractRadialProfile(Data[l],processlist[k],Variablelist[k], radialprofiles[j],R_mesh[l],ISYMlist[l])

				#Plot radial profile and allow for log y-axis if requested.
				ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)


				#Write data to ASCII files if requested.
				if write_ASCII == True:
					if l == 0:
						WriteFolder = 'Z='+str(round((radialprofiles[j])*dz[l], 2))+'cm_Data'
						DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
						WriteDataToFile(Raxis+['\n'], DirWrite+Variablelist[k], 'w')
					#endif
					WriteDataToFile(RProfile+['\n'], DirWrite+Variablelist[k], 'a')
				#endif

				#Apply image options and axis labels.
				Title = 'Comparison of '+Variablelist[k]+' Profiles at Z='+str(round((radialprofiles[j])*dz[l], 2))+'cm for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel,Legend = 'Radial Distance R [cm]',Ylabels[k],Legendlist
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ Z='+str(round((radialprofiles[j])*dz[l], 2))+'cm profiles'+ext)
			plt.close('all')
		#endfor
	#endfor

##===================##===================#

	#Perform vertical (axial) profile comparisons
	for j in range(0,len(axialprofiles)):

		#Create new folder for each axial or radial slice.
		ProfileFolder = 'R='+str(round((axialprofiles[j])*dr[l], 2))+'cm'
		DirProfile = CreateNewFolder(DirComparisons,ProfileFolder)

		#For each requested comparison variable.
		for k in tqdm(range(0,len(Variables))):

			#Loop escape if variables that do not exist have been requested.
			if k >= 1 and k > len(Variablelist)-1:
				break
			#endif

			#Create fig of desired size and refresh legendlist.
			fig,ax = figure(image_aspectratio,1)
			Legendlist = list()

			#For each folder in the directory.
			for l in range(0,numfolders):

				#Create processlist for each folder as required.
				processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

				#Correct processlist for folders containing different icp.dat.
				processlist,Variablelist = VariableInterpolator(processlist,Variablelist,Comparisonlist)

				#Update legend with folder information.
				Legendlist.append( FolderNameTrimmer(Dirlist[l]) )
				Ylabels = VariableLabelMaker(Variablelist)

				#Obtain axial profile for each folder of the current variable.
				ZProfile = ExtractAxialProfile(Data[l],processlist[k],Variablelist[k], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])

				#Plot axial profile and allow for log y-axis if requested.
				ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)


				#Write data to ASCII files if requested.
				if write_ASCII == True:
					if l == 0:
						WriteFolder = 'R='+str(round((axialprofiles[j])*dr[l], 2))+'cm_Data'
						DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
						WriteDataToFile(Zaxis+['\n'], DirWrite+Variablelist[k], 'w')
					#endif
					WriteDataToFile(ZProfile[::-1]+['\n'], DirWrite+Variablelist[k], 'a')
				#endif

				#Apply image options and axis labels.
				Title = 'Comparison of '+Variablelist[k]+' Profiles at R='+str(round((axialprofiles[j])*dr[l], 2))+'cm for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel,Legend = 'Axial Distance Z [cm]',Ylabels[k],Legendlist
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ R='+str(round((axialprofiles[j])*dr[l], 2))+'cm profiles'+ext)
			plt.close('all')
		#endfor
	#endfor

	print('-------------------------------')
	print('# Comparitive Profiles Complete')
	print('-------------------------------')
#endif



##====================================================================#
#				  #MULTI-PROFILES FROM SAME FOLDER#
##====================================================================#

if savefig_multiprofiles == True:

	#For each folder in turn
	for l in range(0,numfolders):
		#Create global multivar folder.
		Dirlineouts = CreateNewFolder(Dirlist[l],'multivar_Profiles/')

		#Create processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])
		multiprocesslist,multivariablelist = VariableEnumerator(multivar,rawdata_2D[l],header_2Dlist[l])

		#Create variable labels with SI unit conversions if required.
		Ylabel = VariableLabelMaker(Variablelist)
		multiYlabel = VariableLabelMaker(multivariablelist)

		#Generate the vertical (height) lineouts for a given radius.
		if len(axialprofiles) > 0:

			#Generate SI scale axes for lineout plots.
			Zaxis = GenerateAxis('Axial',ISYMlist[l])

			#Perform the plotting for all requested variables.
			for i in tqdm(range(0,len(processlist))):

				#Extract the lineout data from the main data array.
				for j in range(0,len(axialprofiles)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((axialprofiles[j])*dr[l], 2))
					DirZlineouts = CreateNewFolder(Dirlineouts,'R='+Slice+'cm/')	

					#Create legendlist
					Legendlist = list()
					Legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot the initial variable in processlist first.
					ZProfile = ExtractAxialProfile(Data[l],processlist[i],Variablelist[i], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])
					ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):
						#Plot profile for multiplot variables in compareprocesslist.
						ZProfile = ExtractAxialProfile(Data[l],multiprocesslist[m],multivariablelist[m], axialprofiles[j],R_mesh[l],Z_mesh[l],ISYMlist[l])
						ImagePlotter1D(Zaxis,ZProfile[::-1],image_aspectratio,fig,ax)

						#Update legendlist with each variable compared.
						Legendlist.append(VariableLabelMaker(multivariablelist)[m])
					#endfor

					#Apply image options and axis labels.
					Title = str(round((axialprofiles[j])*dr[l], 2))+'cm Height profiles for '+Variablelist[i]+','' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Axial Distance Z [cm]',VariableLabelMaker(Variablelist)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save figures in original folder.
					R = 'R='+str(round((axialprofiles[j])*dr[l], 2))+'_'
					plt.savefig(DirZlineouts+R+Variablelist[i]+'_MultiProfiles'+ext)
					plt.close('all')
				#endfor
			#endfor
		#endif

##===================##===================#

		#Generate the horizontal (Radial) lineouts for a given radius.
		if len(radialprofiles) > 0:
			#Create global multivar folder.
			Dirlineouts = CreateNewFolder(Dirlist[l],'multivar_Profiles/')

			#Generate SI scale axes for lineout plots.
			Raxis = GenerateAxis('Radial',ISYMlist[l])

			#Perform the plotting for all requested variables.
			for i in tqdm(range(0,len(processlist))):

				#Perform the plotting for all requested variables.
				for j in range(0,len(radialprofiles)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((radialprofiles[j])*dz[l], 2))
					DirRlineouts = CreateNewFolder(Dirlineouts,'Z='+Slice+'cm/')

					#Create legendlist
					Legendlist = list()
					Legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot profile for initial variable in processlist.
					RProfile = ExtractRadialProfile(Data[l],processlist[i],Variablelist[i], radialprofiles[j],R_mesh[l],ISYMlist[l])
					ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):
						#Plot profile for multiplot variables in compareprocesslist.
						RProfile = ExtractRadialProfile(Data[l],multiprocesslist[m], multivariablelist[m],radialprofiles[j],R_mesh[l],ISYMlist[l])
						ImagePlotter1D(Raxis,RProfile,image_aspectratio,fig,ax)

						#Update legendlist with each variable compared.
						Legendlist.append(VariableLabelMaker(multivariablelist)[m])
					#endfor


					#Apply image options and axis labels.
					Title = str(round((radialprofiles[j])*dz[l], 2))+'cm Radial Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Radial Distance R [cm]',VariableLabelMaker(Variablelist)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save lines in previously created folder.
					Z = 'Z='+str(round((radialprofiles[j])*dz[l], 2))+'_'
					plt.savefig(DirRlineouts+Z+Variablelist[i]+'_MultiProfiles'+ext)
					plt.close('all')
				#endfor
			#endfor
		#endif
	#endfor

	print('-----------------------------')
	print('# Multiplot Profiles Complete')
	print('-----------------------------')
#endif



##====================================================================#
#			  #ITERMOVIE PROFILES - TEMPORAL ANALYSIS#
##====================================================================#

#Plot 1D profile of itervariables at desired locations
if savefig_temporalprofiles == True:

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create new folder and initiate required lists.
		TemporalTrends,Xaxis = list(),list()
		DirTemporal = CreateNewFolder(Dirlist[l],'Temporal_Profiles/')
		DirMeshAve = CreateNewFolder(DirTemporal,'Mesh_Averaged/')

		#Create processlist for each folder as required.
		processlist,variablelist = VariableEnumerator(Variables,rawdata_itermovie[l],header_itermovie[l])
		#Skip over the R and Z processes as they are not saved properly in iterdata.
		for i in range(0,len(processlist)):
			processlist[i] = processlist[i]-2
		#endfor

		#Create list and x-axis for convergence trend plotting.
		#DtActual is approximate, exact dt per iteration depends upon modules and solvers employed.
		DtActual = (1.0/FREQGLOB[l])*100			# [s]	(~8 microseconds @ 13.56MHz)
		DtActual = DtActual*1000					# [ms]
		for i in range(0,len(MovieIterlist[l])):
			IterDigits = list(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))
			IterDigits = float(''.join(IterDigits))
			IterTime = IterDigits*DtActual
			Xaxis.append(IterTime)
		#endfor

		#for all variables requested by the user.
		for i in tqdm(range(0,len(processlist))):

			#Extract 2D image and take mesh averaged value for iteration trend.
			TemporalProfile = list()
			for k in range(0,len(MovieIterlist[l])):
				Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
				# RM: DEFAULT TEMPORAL PROFILE REPRESENTS THE MESH AVERAGED VALUE
				TemporalProfile.append( sum(Image.flatten())/len(Image.flatten()) )
			#endfor
			TemporalTrends.append(TemporalProfile)

			#Plot each variable against simulation real-time.
			fig, ax = plt.subplots(1, figsize=(10,10))
			ax.plot(Xaxis,TemporalProfile, lw=2)

			#Image plotting details.
			Title = 'Mesh-Averaged Temporal Profile of '+str(variablelist[i])+' for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Simulation time [ms]',VariableLabelMaker(variablelist)[i]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

			#Save figure.
			savefig(DirMeshAve+'Temporal_'+variablelist[i]+ext)
			plt.close('all')

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(DirTemporal, 'Temporal_Data')
				DirWriteMeshAve = CreateNewFolder(DirWrite, 'MeshAveraged_Data')
				WriteDataToFile(Xaxis, DirWriteMeshAve+variablelist[i], 'w')
				WriteDataToFile(['\n']+TemporalProfile, DirWriteMeshAve+variablelist[i], 'a')
			#endif
		#endfor

		#=================#

		#Plot mesh averaged value over 'real-time' in simulation.
		Legend = VariableLabelMaker(variablelist)
		fig, ax = plt.subplots(1, figsize=(14,10))

		#Plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(TemporalTrends)):
			TemporalTrends[i] = Normalise(TemporalTrends[i])[0]
			ax.plot(Xaxis,TemporalTrends[i], lw=2)
		#endfor

		#Image plotting details.
		Title = 'Mesh-Averaged Temporal Profiles of '+str(variablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation time [ms]','Normalised Mesh-Average Value'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(0,1.01+(len(Legend)*0.05))

		#Save figure.
		savefig(DirMeshAve+FolderNameTrimmer(Dirlist[l])+'_Normalised'+ext)
		plt.close('all')
	#endfor
	print('----------------------------')
	print('# Temporal Profiles Complete')
	print('----------------------------')
#endif

#=====================================================================#
#=====================================================================#












































#====================================================================#
				 #GENERAL ENERGY DISTRIBUTION ANALYSIS#
#====================================================================#

#====================================================================#
				#ION-NEUTRAL ANGULAR ENERGY DISTRIBUTIONS#
#====================================================================#

if savefig_IEDFangular == True:

	#For all simulation folders.
	for l in range(0,numfolders):

		#Create new folder for keeping EEDF/IEDF if required.
		DirIEDF = CreateNewFolder(Dirlist[l],'EDFplots')

		#Create processlist for requested EDF species and extract images.
		processlist,variablelist = VariableEnumerator(IEDFVariables,rawdata_IEDF[l],header_IEDFlist[l])
		
		#For all requested variables.
		for i in tqdm(range(0,len(processlist))):

			#Create any required arrays
			EDFprofile = list()
		
			#Extract image from required variable and create required profile lists.
			#Flatten angular distribution across all angles to produce energy distribution.
			EDFImage = ImageExtractor2D(DataIEDF[l][processlist[i]],Rmesh=int(EDFangle),Zmesh=int(EDFbins))
			for j in range(0,len(EDFImage)): EDFprofile.append(sum(EDFImage[j]))
			
			#Smooth kinetic data prior to analysis if requested (Savitzk-Golay filter)
			if PlotKineticFiltering == True:
				WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
				EDFprofile = (savgol_filter(EDFprofile, WindowSize, PolyOrder)).tolist()
				#endfor
			#endif

			#Obtain conversion from energy-bin axis to eV axis and construct energy axis
			deV, eVaxis = (EMAXIPCMC/IEBINSPCMC), list()
			for j in range (0,int(IEBINSPCMC)): eVaxis.append(j*deV)

			#Transpose Image for plotting and reverse both lists to align with other data.
			EDFImage, EDFprofile = EDFImage[::-1].transpose(), EDFprofile[::-1]

			#Determine region of IEDF to plot based on threshold value from array maximum.
			Threshold = EDF_Threshold*max(EDFprofile)
			index_max = np.argmax(EDFprofile) 
			for j in range(int(index_max),len(EDFprofile)): 
				if EDFprofile[j] < Threshold and j != 0: 
					eVlimit = j*deV
 			#break
				elif j == len(EDFprofile)-1:
					eVlimit = EMAXIPCMC
				#endif
			#endfor


			#Plot the angular distribution and EDF of the required species.
			fig,ax = figure([11,9], 2, shareX=True)

			Title = Dirlist[l][2::]+'\n'+variablelist[i]+' Angular Energy Distribution Function'
			Extent=[0,int(EMAXIPCMC), -len(EDFImage)/2,len(EDFImage)/2]

			#Angularly resolved IEDF Figure
			im = ax[0].imshow(EDFImage, extent=Extent, aspect='auto')
			cax = Colourbar(ax[0],variablelist[i]+' EDF($\\theta$)',5)
			Xlabel,Ylabel = '','Angular Dispersion [$\\theta^{\circ}$]'
			ImageCrop = [[0,int(eVlimit)],[-45,45]]					#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Crop=ImageCrop,Rotate=False) 

			#Angle Integrated IEDF figure
			ax[1].plot(eVaxis,EDFprofile, lw=2)
			cax = InvisibleColourbar(ax[1])
			Xlabel,Ylabel = 'Energy [eV]',variablelist[i]+' EDF \n [$\\theta$ Integrated]'
			ImageCrop = [[0,int(eVlimit)],[0,max(EDFprofile)*1.05]]	#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=ImageCrop,Rotate=False)

			plt.savefig(DirIEDF+variablelist[i]+'_EDF'+ext)
			plt.close('all')

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				if i == 0:
					DirASCII = CreateNewFolder(DirIEDF, 'EDF_Data')
					WriteDataToFile(eVaxis, DirASCII+variablelist[i],'w')
				#endif
				WriteDataToFile(EDFImage, DirASCII+variablelist[i]+'_IEDFAngular','w')
				WriteDataToFile(EDFprofile, DirASCII+variablelist[i]+'_IEDFProfile','w')
			#endif
		#endfor
	#endfor
#endif



#====================================================================#
				#ION-NEUTRAL ANGULAR ENERGY ANALYSIS#
#====================================================================#

if savefig_IEDFtrends == True:

	#For all requested IEDF variables
	for i in tqdm(range(0,len(IEDFVariables))):
		#Initiate figure for current variable and any required lists.
		Legendlist,EDFprofiles,EDFEnergyProfile = list(),list(),list()
		Mode_eV,Mean_eV,Median_eV = list(),list(),list()
		Total_eV,Range_eV,FWHM_eV = list(),list(),list()
		GlobRange_eV = [0,0]
		fig,ax = figure()
		
		#Create new global trend folder if it doesn't exist already.
		TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')
		#Create folder for IEDF trends if it doesn't exist already.
		DirIEDFTrends = CreateNewFolder(DirTrends,'IEDF Trends')

		#For all simulation folders.
		for l in range(0,numfolders):
			EDFprofile = list()

			#Create processlist for requested EDF species and extract images.
			processlist,variablelist=VariableEnumerator(IEDFVariables,rawdata_IEDF[l],header_IEDFlist[l])
			Legendlist.append(FolderNameTrimmer(Dirlist[l]))

			#Extract image from required variable and flatten angular distribution profile.
			Image = ImageExtractor2D(DataIEDF[l][processlist[i]],Rmesh=EDFangle,Zmesh=EDFbins)
			for j in range(0,len(Image)): EDFprofile.append(sum(Image[j]))
			#Transpose Image for plotting and reverse both lists due to reading error.
			Image, EDFprofile = Image[::-1].transpose(), EDFprofile[::-1]

			#Obtain conversion from energy-bin axis to eV axis and construct energy axis
			deV, eVaxis = (EMAXIPCMC/IEBINSPCMC), list()
			for j in range (0,int(IEBINSPCMC)): eVaxis.append(j*deV)

			#Plot 1D EDF variable profile on open figure for each simulation folder.
			ax.plot(eVaxis,EDFprofile, lw=2)

			#==========#
			#==========#

			#Perform energy trend diagnostics on current folder (variable 'i') IEDF
			#Smooth kinetic data prior to analysis if requested (Savitzk-Golay filter)
			if KineticFiltering == True:
				WindowSize, PolyOrder = Glob_SavWindow, Glob_SavPolyOrder
				EDFprofile = (savgol_filter(EDFprofile, WindowSize, PolyOrder)).tolist()
				#endfor
			#endif

			#Energy extrema analysis: Returns maximum energy where the fraction is above threshold.
			Threshold = EDF_Threshold*max(EDFprofile)
			ThresholdArray = list()
			for j in range(0,len(EDFprofile)):
				if EDFprofile[j] >= Threshold: ThresholdArray.append(j)
			#endfor
			IndexRange = [min(ThresholdArray),max(ThresholdArray)]
			Range_eV = [IndexRange[0]*deV,IndexRange[1]*deV]
			#Save global extrema between folders for plotting range.
			if Range_eV[0] < GlobRange_eV[0]: GlobRange_eV[0] = Range_eV[0]
			if Range_eV[1] > GlobRange_eV[1]: GlobRange_eV[1] = Range_eV[1]

			#Total energy analysis: Returns total energy contained within EDF profile
			EDFEnergyProfile = list()
			for j in range(0,len(EDFprofile)):
				EDFEnergyProfile.append( EDFprofile[j]*(j*deV) )
			#endfor
			Total_eV.append(sum(EDFEnergyProfile))

			#Average energy analysis: Returns mean/mode/median energies from IEDF.
			#Modal energy taken as most common energy fraction (disregarding EDF_threshold)
			Mode_eV.append( EDFprofile.index(max(EDFprofile))*deV )

			#Mean energy obtained as integrated energy fraction most closely matching total energy,
			#which is averaged over effective energy range determined from EDF_threshold. 
			if GlobMeanCalculation == 'MeanEnergy':
				BinAveragedEnergy = sum(EDFEnergyProfile)/len(EDFEnergyProfile)		#(Total Energy/Num Bins)
				ResidualArray = list()
				#Calculate Residuals using EDF_threshold as upper energy percentile
				for j in range(IndexRange[0],IndexRange[1]):						
					ResidualArray.append(abs(EDFEnergyProfile[j]-BinAveragedEnergy))
				#endfor
				#Capture mean/residual intersections, sort high energy intersection indices first
				NumIntersections = int(len(EDFprofile)/10)
				Intersections = np.argsort(ResidualArray)[:NumIntersections]			
				Intersections = sorted(Intersections,reverse=False)
#				#Single Intersection case - special case, maintained mostly for debugging purposes
#				Intersection = (np.abs(EDFEnergyProfile-BinAveragedEnergy)).argmin()
	
				#Index k defines the IEDF array lower energy edge corresponding to the threshold
				IntersectionRegions,k = list(),0			
				for j in range(0,len(Intersections)-1):
					RegionThreshold = 0.05*len(EDFEnergyProfile)
					#If threshold jump observed, save current intersect region index (k)
					if Intersections[j]+RegionThreshold < Intersections[j+1]:
						IntersectionRegions.append(Intersections[k:j])
						k = j+1
					#Save final intersect region index
					elif j == len(Intersections)-2:
						IntersectionRegions.append(Intersections[k:j])
					#endif
				#endfor
				Intersections = list()
				#If odd number of intersections, likely that low energy one was missed
				if len(IntersectionRegions) % 2 == 1: Intersections.append(1.0)
				#For all intersection regions identify the central index
				for j in range(0,len(IntersectionRegions)):	
					try:	
						RegionCentralIndex = int(len(IntersectionRegions[j])/2.0)
						Intersections.append( IntersectionRegions[j][RegionCentralIndex] )
					except:
						#If no intersections, assume a single zero energy intersection
						if j == 0: Intersections.append(0.0)
					#endtry
				#endfor
				#Extrema represent FWHM of EDF, mean energy lies between extrema
				FWHM_eV.append([min(Intersections)*deV,max(Intersections)*deV])
				MeanEnergyIndex = (max(Intersections)+min(Intersections))/2
			#endif

			#Mean energy obtained as ion energy with population fraction most closely matching IEDF average
			#OLD MEAN ENERGY DEFINITION - MEAN DEFINED BY FRACTION NOT BY ENERGY
			if GlobMeanCalculation == 'MeanFraction':
				BinAveragedFraction = sum(EDFprofile)/len(EDFprofile)
				ResidualArray = list()
				#Calculate Residuals using EDF_threshold as upper energy percentile
				for j in range(IndexRange[0],IndexRange[1]):						
					ResidualArray.append(abs(EDFprofile[j]-BinAveragedFraction))
				#endfor
				#Capture mean/residual intersections, sort high energy intersection indices first
				NumIntersections = int(len(EDFprofile)/8)
				Intersections = np.argsort(ResidualArray)[:NumIntersections]			
				Intersections = sorted(Intersections,reverse=False)	

				#Index k defines the IEDF array lower energy edge corresponding to the threshold
				IntersectionRegions,k = list(),0			
				for j in range(0,len(Intersections)-1):
					RegionThreshold = 0.05*len(EDFprofile)
					#If threshold jump observed, save current intersect region index (k)
					if Intersections[j]+RegionThreshold < Intersections[j+1]:
						IntersectionRegions.append(Intersections[k:j])
						k = j+1
					#Save final intersect region index
					elif j == len(Intersections)-2:
						IntersectionRegions.append(Intersections[k:j])
					#endif
				#endfor

				Intersections = list()
				#If odd number of intersections, likely that low energy one was missed
				if len(IntersectionRegions) % 2 == 1: Intersections.append(1.0)
				#For all intersection regions identify the central index
				for j in range(0,len(IntersectionRegions)):	
					try:	
						RegionCentralIndex = int(len(IntersectionRegions[j])/2.0)
						Intersections.append( IntersectionRegions[j][RegionCentralIndex] )
					except:
						#If no intersections, assume a single zero energy intersection
						if j == 0: Intersections.append(0.0)
					#endtry
				#endfor
				#Extrema represent FWHM of EDF, mean energy lies between extrema
				FWHM_eV.append([min(Intersections)*deV,max(Intersections)*deV])
				MeanEnergyIndex = (max(Intersections)+min(Intersections))/2
			#endif
			Mean_eV.append( MeanEnergyIndex*deV )

#			#Median energy calculated as EDF index representing midpoint of equal integrated energies
#			RisingSum,FallingSum,AbsDiff = 0.0,0.0,list()
#			for j in range(0,len(EDFprofile)):
#				Rising_j, Falling_j = j, (len(EDFprofile)-1-2)-j
#				RisingSum += EDFprofile[Rising_j]*(Rising_j*deV)
#				FallingSum += EDFprofile[Falling_j]*(Falling_j*deV)
#				AbsDiff.append(abs(RisingSum-FallingSum))
#				#endif
#			#endfor
#			MedianIndex = AbsDiff.index(min(AbsDiff))
#			Median_eV.append( MedianIndex*deV ) 				#### NB: MEDIANS' ALL FUCKED UP BRAH! ####

#			#Particle energy variance analysis: Returns FWHM of energy distribution.
#			Take mean and draw line at y = mean 
#			Calculate where y = mean intercepts EDFprofile
#			If only one intercept, first intercept is x = 0
#			Integrate EDFprofile indices between intercepts 
#			Save in 1D array, can be used to get energy spread percentage.

			#==========#
			#==========#

			if IDEBUG == True:
				fig2,ax2 = figure()
				try: BinAveragedValue = BinAveragedEnergy		#MeanEnergy Value
				except: BinAveragedValue = BinAveragedFraction	#MeanFraction Value
				Ylims = [0,max(EDFEnergyProfile)]

				fig2,ax2 = figure()
				ax2.plot(EDFEnergyProfile, 'k-', lw=2)
				ax2.plot(ResidualArray, 'r-', lw=2)
				ax2.plot((0,len(EDFEnergyProfile)),(BinAveragedValue,BinAveragedValue),'b--',lw=2)
				ax2.plot((max(Intersections),max(Intersections)),(Ylims[0],Ylims[1]),'b--',lw=2)
				ax2.plot((min(Intersections),min(Intersections)),(Ylims[0],Ylims[1]),'b--',lw=2)
				ax2.plot((MeanEnergyIndex,MeanEnergyIndex),(Ylims[0],Ylims[1]),'m--',lw=2)
				ax2.legend(['Integrated Energy','abs(ResidualArray)','BinAveragedEnergy/Fraction'])
#				plt.savefig(DirIEDFTrends+'_DebugOutput.png')
			#endif
		#endfor

		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			if i == 0:
				DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
				DirASCIIIEDF = CreateNewFolder(DirASCII,'IEDF_Data')
			#endif
			WriteDataToFile(Legendlist+['\n']+Total_eV, DirASCIIIEDF+variablelist[i]+'_Total', 'w')
			WriteDataToFile(Legendlist+['\n']+Range_eV, DirASCIIIEDF+variablelist[i]+'_Range', 'w')
			WriteDataToFile(Legendlist+['\n']+Mode_eV, DirASCIIIEDF+variablelist[i]+'_Mode', 'w')
			WriteDataToFile(Legendlist+['\n']+Mean_eV, DirASCIIIEDF+variablelist[i]+'_Mean', 'w')
		#endif

		##IEDF PROFILES##
		#===============#

		#Apply image options to IEDF plot generated in the above loop.
		Title = Dirlist[l][2::]+'\n'+variablelist[i]+' Angular Energy Distribution Function Profiles'
		Xlabel,Ylabel = 'Energy [eV]',variablelist[i]+' EDF [$\\theta$ Integrated]'
		ImageCrop = [ [0,GlobRange_eV[1]], [] ]		#[[X1,X2],[Y1,Y2]]
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=ImageCrop,Rotate=False)

		plt.savefig(DirIEDFTrends+variablelist[i]+'_EDFprofiles'+ext)
		plt.close('all')

		##ENERGY ANALYSIS##
		#=================#

		#Plot IEDF average energies with respect to simulation folder names.
		fig,ax = figure()
		TrendPlotter(ax,Mean_eV,Legendlist,NormFactor=0)
		TrendPlotter(ax,Mode_eV,Legendlist,NormFactor=0)
##		TrendPlotter(ax,Range_eV[0],Legendlist,NormFactor=0)
##		TrendPlotter(ax,Range_eV[1],Legendlist,NormFactor=0)

		Title = Dirlist[l][2::]+'\n'+'Average '+variablelist[i]+' Energies'
		Legend = ['EDF Mean Energy','EDF Mode Energy','EDF Min Energy','EDF Max Energy']
		Xlabel,Ylabel = 'Varied Property',variablelist[i]+' Energy [eV]'
		ImageCrop = [[],[0,max(Mean_eV+Mode_eV)*1.15]]		#[[X1,X2],[Y1,Y2]]
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=ImageCrop,Rotate=False)

		plt.savefig(DirIEDFTrends+variablelist[i]+'_AverageEnergies'+ext)
		plt.close('all')
	#endfor
#endif


if any([savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
	print('--------------------------------')
	print('# EEDF/IEDF Processing Complete.')
	print('--------------------------------')
#endif

#=====================================================================#
#=====================================================================#










































#====================================================================#
				 #GENERAL TREND PLOTTING ANALYSIS#
#====================================================================#

#====================================================================#
				#COMPARATIVE TRENDS -- MULTI-FOLDER#
#====================================================================#

if savefig_trendphaseaveraged == True or print_generaltrends == True:

	#Create trend folder for outputs - Assumes that scanned variable is within trimmed foler name
	TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
	TrendVariable = ''.join(TrendVariable)												#Single string of chars
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#For each requested comparison variable.
	for k in tqdm(range(0,len(Variables))):

		#Create processlist for largest output, only compare variables shared between all folders.
		processlist,Variablelist = VariableEnumerator(Variables,max(rawdata_2D),max(header_2Dlist))
		processlist,Variablelist = VariableInterpolator(processlist,Variablelist,Comparisonlist)

		#Create Y-axis legend for each variable to be plotted.
		YaxisLegend = VariableLabelMaker(Variablelist)

		#Loop escape if variables that do not exist have been requested.
		if k >= 1 and k > len(Variablelist)-1:
			break
		#endif

		#Create fig of desired size and refresh legendlist.
		fig,ax = figure(image_aspectratio,1)
		Legendlist = list()


		##AXIAL TRENDS##
		#===============#

		#Perform trend analysis on requested axial profiles.
		for j in range(0,len(axialprofiles)):

			#Create folder for axial trends if needed.
			DirAxialTrends = CreateNewFolder(DirTrends,'Axial Trends')

			#Take Trend at Given Location or Default to Min/Max Trends.
			if len(probeloc) == 2:
				#Append requested position to the legendlist.
				R,Z = probeloc[0],probeloc[1]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			elif len(probeloc) == 1:
				#Append requested position to the legendlist.
				R,Z = probeloc[0],axialprofiles[j]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			else:
				#Obtain min/max trend values for requested profile over all folders.
				Xaxis,MaxTrend,MinTrend = MinMaxTrends(axialprofiles[j],'Axial',k)
				Trend = MaxTrend
				#Append the radial position to the legendlist.
				Legendlist.append( 'R='+str(round((axialprofiles[j]*dr[l]), 2))+'cm' )
			#endif

			#Plot trends for each variable over all folders, applying image options.
			TrendPlotter(ax,Trend,Xaxis,NormFactor=0)
			Title='Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

			#Write data to ASCII format datafile if requested.
			if write_ASCII == True:
				if j == 0:
					DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
					DirASCIIAxial = CreateNewFolder(DirASCII,'Axial_Data')
					WriteDataToFile(Xaxis+['\n'], DirASCIIAxial+Variablelist[k]+'_Trends', 'w')
				#endif
				WriteDataToFile(Trend+['\n'], DirASCIIAxial+Variablelist[k]+'_Trends', 'a')
			#endif

		#Save one image per variable with data from all simulations.
		if len(axialprofiles) > 0:
			plt.savefig(DirAxialTrends+'Axial Trends in '+Variablelist[k]+ext)
			plt.close('all')
		#endif


		##RADIAL TRENDS##
		#===============#

		#Create fig of desired size and refresh legendlist.
		fig,ax = figure(image_aspectratio,1)
		Legendlist = list()

		#Perform trend analysis on requested radial profiles.
		for j in range(0,len(radialprofiles)):

			#Create folder for axial trends if needed.
			DirRadialTrends = CreateNewFolder(DirTrends,'Radial Trends')

			#Take Trend at Given Location or Default to Min/Max Trends.
			if len(probeloc) == 2:
				#Append requested position to the legendlist.
				R,Z = probeloc[0],probeloc[1]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			elif len(probeloc) == 1:
				#Append requested position to the legendlist.
				R,Z = radialprofiles[j],probeloc[0],
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			else:
				#Obtain min/max trend values for requested profile over all folders.
				Xaxis,MaxTrend,MinTrend = MinMaxTrends(radialprofiles[j],'Radial',k)
				Trend = MaxTrend
				#Append the axial position to the legendlist.
				Legendlist.append( 'Z='+str(round((radialprofiles[j]*dz[l]), 2))+'cm' )
			#endif

			#Plot trends for each variable over all folders, applying image options.
			TrendPlotter(ax,Trend,Xaxis,NormFactor=0)
			Title='Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

			#Write data to ASCII format datafile if requested.
			if write_ASCII == True:
				if j == 0:
					DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
					DirASCIIRadial = CreateNewFolder(DirASCII,'Radial_Data')
					WriteDataToFile(Xaxis+['\n'], DirASCIIRadial+Variablelist[k]+'_Trends', 'w')
				#endif
				WriteDataToFile(Trend+['\n'], DirASCIIRadial+Variablelist[k]+'_Trends', 'a')
			#endif

		#Save one image per variable with data from all simulations.
		if len(radialprofiles) > 0:
			plt.savefig(DirRadialTrends+'Radial Trends in '+Variablelist[k]+ext)
			plt.close('all')
		#endif
	#endfor
#endif




#=====================================================================#
						# DC-BIAS CALCULATOR #
#=====================================================================#

if savefig_trendphaseaveraged == True or print_DCbias == True:

	#Create Trend folder to keep output plots.
	TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
	TrendVariable = ''.join(TrendVariable)												#Single string of chars
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initiate lists required for storing data.
	Xaxis = list()
	DCbias = list()

	#For all folders.
	for l in range(0,numfolders):

		#Create processlist for each folder as required.
		Process,Variable = VariableEnumerator(['P-POT'],rawdata_2D[l],header_2Dlist[l])

		#Update X-axis with folder information.
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Locate powered electrode for bias extraction.
		Rlineoutloc = WaveformLoc(electrodeloc,'2D')[0]
		Zlineoutloc = WaveformLoc(electrodeloc,'2D')[1]

		#Obtain radial and axial profiles for further processing.
		try: Rlineout = ExtractRadialProfile(Data[l],Process[0],Variable[0],Rlineoutloc,R_mesh[l],  ISYMlist[l])
		except: Rlineout = float('NaN')
		#endtry
		try: Zlineout = ExtractAxialProfile(Data[l],Process[0],Variable[0],Zlineoutloc,R_mesh[l],Z_mesh[l],ISYMlist[l])
		except: Zlineout = float('NaN')
		#endtry

		#Obtain DCbias on axis and across the centre radius of the mesh.
		AxialDCbias = DCbiasMagnitude(Zlineout[::-1])
		RadialDCbias = DCbiasMagnitude(Rlineout)

		#Choose axial or radial DCbias based on user input, else autoselect most probable.
		if DCbiasaxis == 'Radial':
			DCbias.append(RadialDCbias)
		elif DCbiasaxis == 'Axial':
			DCbias.append(AxialDCbias)
		elif DCbiasaxis == 'Auto':
			#Compare Axial and Radial DCbias, if same pick Axial, if not pick the largest.
			if AxialDCbias != RadialDCbias:
				if abs(AxialDCbias) > abs(RadialDCbias):
					DCbias.append(AxialDCbias)
				else:
					DCbias.append(RadialDCbias)
				#endif
			else:
				DCbias.append(AxialDCbias)
			#endif
		#endif

		#Display DCbias to terminal if requested.
		if print_DCbias == True:
			print(Dirlist[l])
			print('DC Bias:',round(DCbias[l],5),'V')
		#endif
	#endfor

	#Write data to ASCII format datafile if requested.
	if write_ASCII == True:
		DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
		DCASCII = [Xaxis,DCbias]
		WriteDataToFile(DCASCII, DirASCII+'DCbias_Trends')
	#endif

	#Plot and beautify the DCbias, applying normalization if requested.
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(ax,DCbias,Xaxis,NormFactor=0)

	#Apply image options and axis labels.
	Title = 'Trend in DCbias with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','DC bias [V]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirTrends+'Powered Electrode DCbias'+ext)
	plt.close('all')
#endif



#====================================================================#
					#POWER DEPOSITED DIAGNOSTIC#
#====================================================================#

if savefig_trendphaseaveraged == True or print_totalpower == True:

	#Create Trend folder to keep output plots.
	TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
	TrendVariable = ''.join(TrendVariable)												#Single string of chars
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Create required lists.
	RequestedPowers,DepositedPowerList = list(),list()
	Xaxis,Powers = list(),list()

	#Update X-axis with folder information.
	for l in range(0,numfolders): Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

	#Create list of requested power variables and ensure they also exist in all compared folders
	for i in range(0,len(Variables)):
		if 'POW-' in Variables[i] and Variables[i] in Comparisonlist:
			RequestedPowers.append(Variables[i])
		#endif
	#endfor

	#For each different power deposition mechanism requested.
	for k in range(0,len(RequestedPowers)):
		#For all folders.
		for l in range(0,numfolders):

			#Create extract data for the neutral flux and neutral velocity.
			processlist,Variablelist = VariableEnumerator(RequestedPowers,rawdata_2D[l],header_2Dlist[l])
			
			#Extract full 2D power density image. [W/m3]
			PowerDensity = ImageExtractor2D(Data[l][processlist[k]])
			PowerDensity = VariableUnitConversion(PowerDensity,Variablelist[k])

			#=====#=====#

			#Cylindrical integration of power per unit volume ==> total coupled power.
			if IXZlist[l] == 0:
			
				Power = 0												#[W]
				#For each axial slice
				for j in range(0,Z_mesh[l]):
					#For each radial slice
					for i in range(0,R_mesh[l]-1):
						#Calculate radial plane volume of a ring at radius [i], correcting for central r=0.
						InnerArea = np.pi*( (i*(dr[l]/100))**2 )		#[m^2]
						OuterArea = np.pi*( ((i+1)*(dr[l]/100))**2 )	#[m^2]
						RingVolume = (OuterArea-InnerArea)*(dz[l]/100)	#[m^3]

						#Calculate Power by multiplying power density for ring[i] by volume of ring[i]
						Power += PowerDensity[j][i]*RingVolume 			#[W]
					#endfor
				#endfor
				DepositedPowerList.append(Power)

				#Display power to terminal if requested.
				if print_totalpower == True:
					print(Dirlist[l])
					print(RequestedPowers[k]+' Deposited:',round(Power,4),'W')
				#endif
				
			#=====#=====#
				
			#Cartesian integration of power per unit volume ==> total coupled power.	
			elif IXZlist[l] == 1:
			
				Power = 0												#[W]
				#For each axial slice
				for j in range(0,Z_mesh[l]):
					#For each radial slice
					for i in range(0,R_mesh[l]-1):
					
						#Calculate cell area, all cells have the same area in a cartesian grid
						CellArea = (dr[l]/100.)*(dz[l]/100.)			#[m^2]
						CellVolume = CellArea*(dy[l]/100.)				#[m^3]

						#Calculate Power by multiplying power density for ring[i] by volume of ring[i]
						Power += PowerDensity[j][i]*CellVolume 			#[W]
					#endfor
				#endfor
				DepositedPowerList.append(Power)

				#Display power to terminal if requested.
				if print_totalpower == True:
					print(Dirlist[l])
					print(RequestedPowers[k]+' Deposited:',round(Power,4),'W')
				#endif
			#endif	
		#endfor
		
		#==========#==========#

		#Plot and beautify each requested power deposition seperately.
		fig,ax = figure(image_aspectratio,1)
		Powers.append( DepositedPowerList[k*numfolders:(k+1)*numfolders] )
		TrendPlotter(ax,Powers[k],Xaxis,NormFactor=0)

		#Apply image options and axis labels.
		Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

		plt.savefig(DirTrends+RequestedPowers[k]+' Deposition Trends'+ext)
		plt.close('all')
	#endfor

	#Write data to ASCII format datafile if requested.
	if write_ASCII == True:
		DirASCII, TotalPowerASCII = CreateNewFolder(DirTrends,'Trend_Data'), [Xaxis]
		for k in range(0,len(RequestedPowers)): TotalPowerASCII.append(Powers[k])
		WriteDataToFile(TotalPowerASCII, DirASCII+'RFPower_Trends')
	#endif

	#Plot a comparison of all power depositions requested.
	fig,ax = figure(image_aspectratio,1)
	for k in range(0,len(RequestedPowers)):TrendPlotter(ax,Powers[k],Xaxis,NormFactor=0)

	#Apply image options and axis labels.
	Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

	plt.savefig(DirTrends+'Power Deposition Comparison'+ext)
	plt.close('all')
#endif



#====================================================================#
				  	#ION/NEUTRAL THRUST ANALYSIS#
#====================================================================#

#ABORT DIAGNOSTIC UNLESS ARGON IS SUPPLIED, WILL FIX LATER!!!
if 'AR3S' in list(set(FluidSpecies).intersection(Variables)):
	if savefig_trendphaseaveraged == True or print_thrust == True:

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		NeutralThrustlist,IonThrustlist,Thrustlist = list(),list(),list()
		NeutralIsplist,IonIsplist,ThrustIsplist = list(),list(),list()
		Xaxis = list()
	
		#Extract Positive, Negative and Neutral species names (Excluding electrons)
		NeutralSpecies = list(set(FluidSpecies).intersection(Variables))
		PositiveIons = PosSpecies
		NegativeIons = NegSpecies[:-1]

		#For all folders.
		for l in range(0,numfolders):

			#Update X-axis with folder information.
			Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

			#Extract data required for Thrust calculations, discharge plane (Z) = thrustloc.
			processlist,variablelist = VariableEnumerator(['VZ-NEUTRAL'],rawdata_2D[l],header_2Dlist[l])
			NeutralVelocity = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
			processlist,variablelist = VariableEnumerator(['VZ-ION+'],rawdata_2D[l],header_2Dlist[l])
			IonVelocity = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
			processlist,variablelist = VariableEnumerator(['FZ-AR3S'],rawdata_2D[l],header_2Dlist[l])
			NeutralAxialFlux = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
			processlist,variablelist = VariableEnumerator(['FZ-AR+'],rawdata_2D[l],header_2Dlist[l])
			IonAxialFlux = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
			processlist,variablelist = VariableEnumerator(['TG-AVE'],rawdata_2D[l],header_2Dlist[l])
			NeutGasTemp = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
			processlist,variablelist = VariableEnumerator(['PRESSURE'],rawdata_2D[l],header_2Dlist[l])
			try: 
				Pressure = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc)
				PressureDown = ExtractRadialProfile(Data[l],processlist[0],variablelist[0],thrustloc+1)
			except: 
				Pressure = np.zeros(R_mesh[l]*2)
				PressureDown = np.zeros(R_mesh[l]*2)
			#endtry

			#Convert pressure to Torr if required (delta pressure in thrust calculations expect Torr)
			if PressureUnit == 'Pa':
				for i in range(0,len(Pressure)): Pressure[i] = Pressure[i]/133.33
				for i in range(0,len(PressureDown)): PressureDown[i] = PressureDown[i]/133.33
			elif PressureUnit == 'mTorr':
				for i in range(0,len(Pressure)): Pressure[i] = Pressure[i]/1000.0
				for i in range(0,len(PressureDown)): PressureDown[i] = PressureDown[i]/1000.0
			#endif

			#Define which gas is used and calculate neutral mass per atom.
			NeutralIsp,IonIsp = list(),list()
			Argon,Xenon = 39.948,131.29			 #amu
			NeutralMass = Argon*1.67E-27		 #Kg

			#Choose which method to solve for thrust: 'ThermalVelocity','AxialMomentum'
			if GlobThrustMethod == 'ThermalVelocity':
				#Technique assumes cylindrical geometry, cartesian geometry will be overestimated.
				#Integrates neutral momentum loss rate based on neutral gas temperature.
				#Assumes angularly symmetric temperature and Maxwellian velocity distribution.
				NeutralThrust = 0
				for i in range(0,R_mesh[l]):
					#Calculate radial plane area of a ring at radius [i], correcting for central r=0.
					Circumference = 2*np.pi*(i*(dr[l]/100))		#m
					CellArea = Circumference*(dr[l]/100)		#m^2
					if CellArea == 0:
						CellArea = np.pi*(dr[l]/100)**2			#m^2
					#endif  

					#Calculate most probable neutral velocity based on temperature
					MeanVelocity = np.sqrt( (2*1.38E-23*NeutGasTemp[i])/(NeutralMass) )  	#m/s

					#If current cell is gas phase (Pressure > 0.0), calculate thrust
					if Pressure[i] > 0.0:
						#Calculate Neutral mass flow rate and integrate thrust via F = (dm/dt)Ve.
						NeutralMassFlowRate = NeutralAxialFlux[i]*NeutralMass*CellArea	#Kg/s
						NeutralExitVelocity = NeutralVelocity[i]						#m/s
						NeutralThrust += NeutralMassFlowRate * NeutralExitVelocity 		#N
						if NeutralExitVelocity > 0:
							NeutralIsp.append(NeutralExitVelocity)
						#endif
					#endif
				#endfor

				#Add neutral thrust and Isp to arrays (dummy variables not calculated)
				NeutralThrustlist.append( round(NeutralThrust*1000,5) )		#mN
				Thrustlist.append( round( NeutralThrust*1000,5) )			#mN
				NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81			#s
				Thrust,ThrustIsp = NeutralThrust,NeutralIsp					#N,s
				IonThrust,IonIsp = 1E-30,1E-30								#'Not Calculated'
				DiffForce = 1E-30											#'Not Calculated'
			#endif

			#====================#

			elif GlobThrustMethod == 'AxialMomentum':
				#Technique assumes cylindrical geometry, cartesian geometry will be overestimated.
				#Integrates ion/neutral momentum loss rate and differental pressure for concentric rings.
				#Assumes pressure differential, ion/neutral flux equal for all angles at given radii.

				#CellArea increases from central R=0.
				#Ensure pressure index aligns with radial index for correct cell area.
				#ONLY WORKS WHEN SYMMETRY OPTION IS ON, NEED A MORE ROBUST METHOD!
				Pressure,PressureDown = Pressure[0:R_mesh[l]][::-1],PressureDown[0:R_mesh[l]][::-1]
				NeutralVelocity = NeutralVelocity[0:R_mesh[l]][::-1]
				NeutralAxialFlux = NeutralAxialFlux[0:R_mesh[l]][::-1]
				IonVelocity = IonVelocity[0:R_mesh[l]][::-1]
				IonAxialFlux = IonAxialFlux[0:R_mesh[l]][::-1]
				#ONLY WORKS WHEN SYMMETRY OPTION IS ON, NEED A MORE ROBUST METHOD!

				DiffForce,NeutralThrust,IonThrust = 0,0,0
				for i in range(0,R_mesh[l]):
					#Calculate radial plane area of a ring at radius [i], correcting for central r=0.
					Circumference = 2*np.pi*(i*(dr[l]/100))		#m
					CellArea = Circumference*(dr[l]/100)		#m^2
					if CellArea == 0:
						CellArea = np.pi*((dr[l]/100)**2)		#m^2
					#endif

					#Calculate differential pressure between thrustloc-(thrustloc+1)
					if Pressure[i] > 0.0:
						DiffPressure = (Pressure[i]-PRESOUT[l])*133.33			#N/m^2
#						DiffPressure = (Pressure[i]-PressureDown[i])*133.33		#N/m^2
						DiffForce += DiffPressure*CellArea						#N
					else:
						DiffForce += 0.0
					#endif

					#Calculate Neutral mass flow rate and integrate thrust via F = (dm/dt)Ve.
					NeutralMassFlowRate = NeutralAxialFlux[i]*NeutralMass*CellArea	#Kg/s
					NeutralExitVelocity = NeutralVelocity[i]						#m/s
					NeutralThrust += NeutralMassFlowRate * NeutralExitVelocity 		#N
					if NeutralExitVelocity > 0:
						NeutralIsp.append(NeutralExitVelocity)
					#endif

					#Calculate Ion mass flow rate and integrate thrust via F = (dm/dt)Ve.
					IonMassFlowRate = IonAxialFlux[i]*NeutralMass*CellArea	#Kg/s
					IonExitVelocity = IonVelocity[i]*1000					#m/s
					IonThrust += IonMassFlowRate * IonExitVelocity 			#N
					if IonExitVelocity > 0: 
						IonIsp.append(IonExitVelocity)
					#endif
				#endfor
				if len(IonIsp) == 0: IonIsp.append(np.nan)
				if len(NeutralIsp) == 0: NeutralIsp.append(np.nan)

				#Add total thrust and calculate Isp of each component
				Thrust = DiffForce + NeutralThrust + IonThrust				#N
				NeutralFraction = NeutralThrust/(Thrust-DiffForce)			#Ignore dP/dz
				IonFraction = IonThrust/(Thrust-DiffForce)					#Ignore dP/dz

				IonIsp = (sum(IonIsp)/len(IonIsp))/9.81						#s
				NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81			#s
				ThrustIsp = NeutralFraction*NeutralIsp+IonFraction*IonIsp 	#s

				NeutralThrustlist.append( round(NeutralThrust*1000,5) )		#mN
				IonThrustlist.append( round(IonThrust*1000,5) )				#mN
				Thrustlist.append( round(Thrust*1000,5) )					#mN
				NeutralIsplist.append( round(NeutralIsp,5) )				#s
				IonIsplist.append( round(IonIsp,5) )						#s
				ThrustIsplist.append( round(ThrustIsp,5) )					#s
			#endif

			#====================#

			#Display thrust to terminal if requested.
			if print_thrust == True:
				print(Dirlist[l], '@ Z=',round(thrustloc*dz[l],2),'cm')
				print('NeutralThrust', round(NeutralThrust*1000,2), 'mN @ ', round(NeutralIsp,2),'s')
				print('IonThrust:', round(IonThrust*1000,4), 'mN @ ', round(IonIsp,2),'s')
				print('D-Pressure:', round(DiffForce*1000,4), 'mN')
				print('Thrust:',round(Thrust*1000,4),'mN @ ', round(ThrustIsp,2),'s')
				print('')
			#endif
		#endfor

		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
			WriteDataToFile(Xaxis+['\n'], DirASCII+'Thrust_Trends','w')
			WriteDataToFile(Thrustlist+['\n'], DirASCII+'Thrust_Trends','a')
			WriteDataToFile(ThrustIsplist, DirASCII+'Thrust_Trends','a')
		#endif

		#=====#=====#

		#Plot total thrust and ion/neutral components.
		fig,ax1 = figure(image_aspectratio,1)
		TrendPlotter(ax1,Thrustlist,Xaxis,Marker='ko-',NormFactor=0)
		TrendPlotter(ax1,NeutralThrustlist,Xaxis,Marker='r^-',NormFactor=0)
#		TrendPlotter(ax1,IonThrustlist,Xaxis,Marker='bs-',NormFactor=0)

		#Apply image options and save figure.
		Title='Thrust at Z='+str(round(thrustloc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Thrust F$_{T}$ [mN]'
		ax1.legend(['Total Thrust','Neutral Component','Ion Component'], fontsize=18, frameon=False)
		ImageOptions(fig,ax1,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirTrends+'Thrust Trends'+ext)
		plt.close('all')

		#=====#=====#

		#Plot Specific Impulse for total thrust and ion/neutral components.
		fig,ax1 = figure(image_aspectratio,1)
		TrendPlotter(ax1,ThrustIsplist,Xaxis,Marker='ko-',NormFactor=0)
		TrendPlotter(ax1,NeutralIsplist,Xaxis,Marker='r^-',NormFactor=0)
#		TrendPlotter(ax1,IonIsplist,Xaxis,Marker='bs-',NormFactor=0)

		#Apply image options and save figure.
		Title = 'Specific Impulse at Z='+str(round(thrustloc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Specific Impulse I$_{sp}$ [s]'
		ax1.legend(['Total I$_{sp}$','Neutral Component','Ion Component'], fontsize=18, frameon=False)
		ImageOptions(fig,ax1,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirTrends+'Isp Trends'+ext)
		plt.close('all')
	#endif
#endif



#====================================================================#
			 		#PHASE-AVERAGED SHEATH TRENDS#
#====================================================================#

if savefig_trendphaseaveraged == True or print_sheath == True:
#NB: 	This diagnostic is very out of date, particularily the sheathROI treatment...
#		Could be easily worked into the savefig_temporaltrends diagnostic or simply re-written
#		Might be worth considering an overhaul... but it works for now.

	#Create Trend folder to keep output plots.
	TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
	TrendVariable = ''.join(TrendVariable)												#Single string of chars
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initialize any required lists.
	Xaxis,SxLocExtent,SxMaxExtent = list(),list(),list()	

	#Obtain sheathROI and sourcewidth automatically if none are supplied.
	if len(sheathROI) != 2:
		#image_radialcrop Convert to Cells 
		#image_axialcrop Convert to Cells
		#Use axialcrop or radialcrop to set automatic ROI!
		Start,End = 34,72				#AUTOMATIC ROUTINE REQUIRED#
		sheathROI = [Start,End]			#AUTOMATIC ROUTINE REQUIRED#
	#endif
	if len(sourcewidth) == 0:
		#Take Variable that is zero in metals (Density?)
		#Take Axial/Radial slice depending on sheath direction.
		#Find Cell distance from zero to 'wall' at electrodeloc.
		#Convert to SI [cm], set to automatic width.
		sourcewidth = [0.21]			#AUTOMATIC ROUTINE REQUIRED#
	#endif

	SxMeanExtent,SxMeanExtentArray = list(),list()
	#For all selected simulations, obtain Xaxis, sheath value and save to array.
	for l in range(0,numfolders):
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Obtain sheath thickness array for current folder 
		Sx,SxAxis = CalcSheathExtent(folder=l)

		#Calculate mean sheath extent across ROI. On failure provide null point for sheath thickness.
		try:
			SxMeanExtentArray = list()
			for i in range(sheathROI[0],sheathROI[1]):	SxMeanExtentArray.append(Sx[i])
			SxMeanExtent.append(sum(SxMeanExtentArray)/len(SxMeanExtentArray))
		except:
			SxMeanExtent.append( np.nan )
		#endtry

		#Extract maximum sheath thickness from within region of interest
		try: SxMaxExtent.append( ((sourcewidth[0]*dr[l])-max(Sx[sheathROI[0]:sheathROI[1]]))*10 )
		except: SxMaxExtent.append( np.nan )

		#Extract sheath width adjacent to powered electrode
		#loc = electrodeloc[0]		#Radial
		loc = electrodeloc[1] 		#Axial
		try: SxLocExtent.append( ((sourcewidth[0]*dr[l])-Sx[loc])*10 )
		except:	SxLocExtent.append( np.nan )
	#endfor

	#===============================#

	#Write trend data to ASCII format datafile if requested.
	if write_ASCII == True:
		DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
		WriteDataToFile(Xaxis, DirASCII+'Sx-Avg_Trends','w')
		WriteDataToFile('\n', DirASCII+'Sx-Avg_Trends','w')
		WriteDataToFile(SxMaxExtent, DirASCII+'Sx-Avg_Trends','a')
	#endif

	#Generate figure and plot trends.	
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(ax,SxMeanExtent,Xaxis,NormFactor=0)

	#Apply image options and axis labels.
	Title = 'Maximum Sheath Extension With Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Extension [mm]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

	plt.savefig(DirTrends+'Sheath Extension (Phase-Averaged)'+ext)
	plt.close('all')
#endif



#====================================================================#
				  		#KNUDSEN NUMBER ANALYSIS#
#====================================================================#

#Only perform on bulk fluid dynamics relevent species.
if bool(set(FluidSpecies).intersection(Variables)) == True:
	if savefig_trendphaseaveraged == True or print_Knudsennumber == True:

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		KnudsenAverage,Xaxis = list(),list()

		#For all folders - calculate Knudsen number for each cell of TECPLOT2D data
		for l in range(0,numfolders):

			#Using effective radius of argon in this calculation.
			Dimentionality = 2*(Radius[l]/100)		#[m]
			CrossSection = np.pi*((7.1E-11)**2)		#[m^2]

			#Extract data for the neutral flux and neutral velocity.
			processlist,Variablelist = VariableEnumerator(FluidSpecies,rawdata_2D[l],header_2Dlist[l])
			Sx,SxAxis = CalcSheathExtent(folder=l)

			#Update X-axis with folder information.
			Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

			#Create empty image array based on mesh size and symmetry options.
			numrows = len(Data[l][0])/R_mesh[l]
			Image = np.zeros([Z_mesh[l],R_mesh[l]])

			#Produce Knudsen number 2D image using density image.
			for j in range(0,Z_mesh[l]):
				for i in range(0,R_mesh[l]):
					Start = R_mesh[l]*j
					Row = Z_mesh[l]-1-j

					LocalDensity = (Data[l][processlist[0]][Start+i])*1E6					#[m-3]
					try:
						KnudsenNumber = (1/(LocalDensity*CrossSection*Dimentionality))		#[-]
					except:
						KnudsenNumber = 0													#[-]
					#endtry
					Image[Row,i] = KnudsenNumber
				#endfor
			#endfor

			#Display average Knudsen number to terminal if requested.
			KnudsenAverage.append( sum(Image)/(len(Image[0])*len(Image)) )
			if print_Knudsennumber == True:
				print(Dirlist[l])
				print('Average Knudsen Number:', KnudsenAverage[l])
			#endif

			#Create new folder to keep 2D output plots.
			Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')
			#Write image data to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = CreateNewFolder(Dir2Dplots,'2Dplots_Data')
				WriteDataToFile(Image, DirASCII+'Kn','w')
			#endif

			#Label and save the 2D Plots.
			extent,aspectratio = DataExtent(l)
			fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio)
			PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l])

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = 'Knudsen Number Image for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
			cax = Colourbar(ax,'Knudsen Number $K_{n}$',5,Lim=CbarMinMax(Image))
			ImageOptions(fig,ax,Xlabel,Ylabel,Title)

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot Kn'+ext)
			plt.close('all')
		#endfor


		#Write trend data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
			WriteDataToFile(Xaxis, DirASCII+'Kn_Trends','w')
			WriteDataToFile('\n', DirASCII+'Kn_Trends','w')
			WriteDataToFile(KnudsenAverage, DirASCII+'Kn_Trends','a')
		#endif

		#Plot a comparison of all average Knudsen numbers.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(ax,KnudsenAverage,Xaxis,NormFactor=0)

		#Image plotting details.
		Title = 'Average Knudsen Number with Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Average Knudsen Number $K_{n}$'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		#Save figure.
		plt.savefig(DirTrends+'KnudsenNumber_Comparison'+ext)
		plt.close('all')
	#endif
#endif



#====================================================================#
				#REYNOLDS NUMBER / SOUND SPEED ANALYSIS#
#====================================================================#

#Only perform on bulk fluid dynamics relevent species.
if bool(set(FluidSpecies).intersection(Variables)) == True:
	if savefig_trendphaseaveraged == True or print_Reynolds == True:

		#Create Trend folder to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		AverageSoundSpeed,Xaxis = list(),list()
		NeutralDensities = list()

		#For all folders.
		for l in range(0,numfolders):

			#Extract spatially resolved pressure and neutral densities.
			processlist,variablelist = VariableEnumerator(['PRESSURE'],rawdata_2D[l],header_2Dlist[l])
			Pressure = ImageExtractor2D(Data[l][processlist[0]],variablelist[0])
			Sx,SxAxis = CalcSheathExtent(folder=l)
			
			#If only single neutral species - extract that density
			processlist,Variablelist = VariableEnumerator(FluidSpecies,rawdata_2D[l],header_2Dlist[l])
			if len(processlist) == 1: 
				NeutralDensity = ImageExtractor2D(Data[l][processlist[0]],variablelist[0])
			#If there are multiple neutral species, combine them to get total neutral density
			elif len(processlist) > 1:
				for i in range(0,len(processlist)):
					NeutralDensities.append( ImageExtractor2D(Data[l][processlist[i]],variablelist[i]) )
				#endfor

				#Create empty neutral density array based on mesh size and symmetry options.
				numrows = len(Data[l][0])/R_mesh[l]
				NeutralDensity = np.zeros([Z_mesh[l],R_mesh[l]])

				#Combine all neutral densities to get total neutral density - if required.
				for i in range(0,len(NeutralDensities)):
					for j in range(0,len(NeutralDensities[i])):
						for k in range(0,len(NeutralDensities[i][j])):
							NeutralDensity[j][k] += NeutralDensities[i][j][k]
						#endfor
					#endfor
				#endfor
			#endif

			#Update X-axis with folder information.
			Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

			#Calculate 2D sound speed image using neutral density and pressure
			Image = CalcSoundSpeed(NeutralDensity,Pressure,Dimension='2D')

			#Display mesh-averaged sound speed to terminal if requested.
			AverageSoundSpeed.append( sum(Image)/(len(Image[0])*len(Image)) )
			if print_Reynolds == True:
				print(Dirlist[l])
				print('Average Sound Speed:', AverageSoundSpeed[l])
			#endif

			#Create new folder to keep 2D output plots.
			Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')
			#Write image data to ASCII format datafile if requested.
			if write_ASCII == True:
				DirASCII = CreateNewFolder(Dir2Dplots,'2Dplots_Data')
				WriteDataToFile(Image, DirASCII+'Cs','w')
			#endif

			#Label and save the 2D Plots.
			extent,aspectratio = DataExtent(l)
			fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio)
			PlotSheathExtent(SxAxis,Sx,ax,ISYMlist[l])

			#Image plotting details, invert Y-axis to fit 1D profiles.
			#ERROR WITH IMAGE LIMIT - LIKELY DUE TO NANS - #Lim=CbarMinMax(Image)
			Title = 'Sound Speed Image for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
			cax = Colourbar(ax,'Sound Speed $C_{s}$ [m/s]',5,Lim=[]) 
			ImageOptions(fig,ax,Xlabel,Ylabel,Title)

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot Cs'+ext)
			plt.close('all')
		#endfor
		

		#Write trend data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
			WriteDataToFile(Xaxis, DirASCII+'Cs_Trends','w')
			WriteDataToFile('\n', DirASCII+'Cs_Trends','w')
			WriteDataToFile(AverageSoundSpeed, DirASCII+'Cs_Trends','a')
		#endif

		#Plot a comparison of all average Knudsen numbers.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(ax,AverageSoundSpeed,Xaxis,NormFactor=0)

		#Image plotting details.
		Title = 'Average Sound Speed with Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Average Sound Speed'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

		#Save figure.
		plt.savefig(DirTrends+'SoundSpeed_Comparison'+ext)
		plt.close('all')
	#endif
#endif

#===============================#
#===============================#

if any([savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_Reynolds, print_totalpower, print_DCbias, print_thrust, print_sheath]) == True:
	print('---------------------------')
	print('# Trend Processing Complete')
	print('---------------------------')
#endif

#=====================================================================#
#=====================================================================#

















































#====================================================================#
				#PHASE RESOLVED DIAGNOSTICS (REQ MOVIE1)#
#====================================================================#

#====================================================================#
						#1D PHASE RESOLVED MOVIES#
#====================================================================#

#Plot Phase-Resolved profiles with electrode voltage and requested variables.
if savefig_phaseresolve1D == True:

	#Tnitiate any required lists.
	VoltageWaveforms,WaveformBiases,VariedValuelist = list(),list(),list()

	#for all folders.
	for l in range(0,numfolders):

		#Create global folders to keep output plots and collect graph title.
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'1DPhase')
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get 'E','AR+','PPOT')
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables)
		SxData,SxPhase,Sxproc,Sxvar = ExtractPhaseData(folder=l,Variables=['E','AR+'])
		PPOT = ExtractPhaseData(folder=l,Variables=['PPOT'])[2][0]

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)
		Raxis = GenerateAxis('Radial',ISYMlist[l])
		Zaxis = GenerateAxis('Axial',ISYMlist[l])

		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)):
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

		plt.savefig(DirPhaseResolved+VariedValuelist[l]+' Waveform'+ext)
		plt.close('all')

		#Write waveform data in ASCII format if required.
		if write_ASCII == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'1DPhase_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPhase+'VoltageWaveforms')
		#endif

		#==============#

		#for all requested variables.
		for i in tqdm(range(0,len(proclist))):

			#Create new folder to keep specific plots.
			DirMovieplots = CreateNewFolder(DirPhaseResolved,varlist[i]+'_1DPhaseResolved')

			#Refresh lineout lists between variables.
			Lineouts,ProfileOrientation = list(),list()

			#Concatinate all requested lineouts together, keeping seperate orientation.
			for m in range(0,len(radialprofiles)):
				Lineouts.append(radialprofiles[m])
				ProfileOrientation.append('Radial')
			#endfor
			for m in range(0,len(axialprofiles)):
				Lineouts.append(axialprofiles[m])
				ProfileOrientation.append('Axial')
			#endfor

			#For all requested lineouts and orientations.
			for k in range(0,len(Lineouts)):

				#Refresh required lists.
				VariableMax,VariableMin = list(),list()

				#Create folders to keep output plots for each variable.
				if ProfileOrientation[k] == 'Axial':
					NameString= varlist[i]+'_'+str(round(Lineouts[k]*dr[l],2))+'cm[R]'
				if ProfileOrientation[k] == 'Radial':
					NameString= varlist[i]+'_'+str(round(Lineouts[k]*dz[l],2))+'cm[Z]'
				if savefig_phaseresolve1D == True:
					Dir1DProfiles = CreateNewFolder(DirMovieplots,NameString)
				#endif

				#Collect Normalization data for plotting.
				for j in range(0,len(Phaselist)):
					#Record local maximum and minimum for each phase.
					if ProfileOrientation[k] == 'Axial': 
						Profile = ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],Lineouts[k])
					elif ProfileOrientation[k] == 'Radial': 
						Profile = ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],Lineouts[k])
					#endif
					Profile,Minimum,Maximum = Normalise(Profile)
					VariableMax.append(Maximum)
					VariableMin.append(Minimum)
				#endfor
				#Find global maximum and minimum for all phases.
				VariableMax = max(VariableMax)
				VariableMin = min(VariableMin)

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Phaselist)):
					Phase = int( round(Phaseaxis[j]*360.0,3) )		#[Deg]

					if ProfileOrientation[k] == 'Axial':
						ZlineoutLoc,axis = Lineouts[k],Zaxis
						PhaseResolvedProfile = ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],ISYMlist[l])[::-1]
						ProfileString = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm \n'
						Xlabel = 'Axial Distance Z [cm]'
					elif ProfileOrientation[k] == 'Radial':
						RlineoutLoc,axis = Lineouts[k],Raxis
						PhaseResolvedProfile = ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],RlineoutLoc,R_mesh[l],ISYMlist[l])
						ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm \n'
						Xlabel = 'Radial Distance R [cm]'
					#endif

					#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
					fig,ax = figure(image_aspectratio,2)
					Ylabel = VariableLabelMaker(varlist)
					fig.suptitle('Phase-Resolved '+varlist[i]+' for '+VariedValuelist[l]+ProfileString+str(Phaselist[j]), y=0.97, fontsize=16)

					#Plot profile and apply image options.
					ax[0].plot(axis, PhaseResolvedProfile, lw=2)
					ImageOptions(fig,ax[0],Xlabel,Ylabel[i],Crop=False)
					ax[0].set_ylim(VariableMin,VariableMax*1.02)

					#Plot waveform and apply image options.
					ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
					ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
					Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
					ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)

					#Clean up image and save with relevent filename.
					#NOTE:	zfill assumes phase < 999 degrees	(i.e. < 1e4)
					fig.tight_layout()
					plt.subplots_adjust(top=0.90)
					plt.savefig(Dir1DProfiles+NameString+'_'+str(Phase).zfill(3)+ext)
					plt.close('all')

					#Write Phase data in ASCII format if required.
					if write_ASCII == True:
						DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'1DPhase_Data')
						DirASCIIPhaseloc = CreateNewFolder(DirASCIIPhase,ProfileString[3:-2])
						Cycle = str( Phaselist[j].replace(" ", "") )
						SaveString = DirASCIIPhaseloc+varlist[i]+'_'+Cycle
						WriteDataToFile(PhaseResolvedProfile[::-1], SaveString)
					#endif
				#endfor

				#Create .mp4 movie from completed images.
				Prefix = FolderNameTrimmer(Dirlist[l])+'_'+NameString
				MakeMovie(Dir1DProfiles,Prefix)
			#endfor
		#endfor
	#endfor

	print('---------------------------------------')
	print('# 1D Phase-Resolved Processing Complete')
	print('---------------------------------------')
#endif




#====================================================================#
					#2D PHASE RESOLVED MOVIES#
#====================================================================#

#Plot 2D images over all saved phase cycles with included wavevform guide.
if savefig_phaseresolve2D == True:

	#Initialize required lists.
	VoltageWaveforms,WaveformBiases,VariedValuelist = list(),list(),list()

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create global folder to keep output plots and collect graph title.
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'2DPhase/')
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get 'E','AR+','PPOT')
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables)
		SxData,SxPhase,Sxproc,Sxvar = ExtractPhaseData(folder=l,Variables=['E','AR+'])
		PPOT = ExtractPhaseData(folder=l,Variables=['PPOT'])[2][0]

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)
		Raxis = GenerateAxis('Radial',ISYMlist[l])
		Zaxis = GenerateAxis('Axial',ISYMlist[l])

		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)): 
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

		plt.savefig(DirPhaseResolved+VariedValuelist[l]+' Waveform'+ext)
		plt.close('all')

		#Write PROES data in ASCII format if required.
		if write_ASCII == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'2DPhase_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPhase+'VoltageWaveforms')
		#endif

		#===============#


		#for all variables requested by the user.
		for i in tqdm(range(0,len(proclist))):

			#Create new folder to keep specific plots.
			DirMovieplots = CreateNewFolder(DirPhaseResolved,varlist[i]+'_2DPhaseResolved/')

			#Obtain maximum and minimum values of current variable over all phases.
			MinLim,MaxLim = list(),list()
			for j in range(0,len(Phaselist)):
				Image = ImageExtractor2D(PhaseData[j][proclist[i]],varlist[i])
				MinLim.append( CbarMinMax(Image,Symmetry=False)[0] )
				MaxLim.append( CbarMinMax(Image,Symmetry=False)[1] )
			#endfor
			Limits = [min(MinLim),max(MaxLim)]

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for j in range(0,len(Phaselist)):
				Phase = int( round(Phaseaxis[j]*360.0,3) )		#[Deg]

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(PhaseData[j][proclist[i]],varlist[i])
				#Extract Ni and Ne variables for sheath processing.
				if image_plotsheath == True:
					Ne = SxData[j][Sxproc[Sxvar.index('E')]]
					Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
					Sx,SxAxis = CalcSheathExtent(folder=l,Phase=j,Ne=Ne,Ni=Ni)
				#endif

				#Obtain image extent and axis labels based on image symmetry and rotation.
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				if image_rotate == True: Xlabel,Ylabel = Ylabel,Xlabel
				extent,aspectratio = DataExtent(l)

				#Create figure and axes, plot image on top and waveform underneath.
				fig,ax = figure(aspectratio,2)
				Title = 'Phase-Resolved '+varlist[i]+'\n'+str(Phaselist[j])
				fig.suptitle(Title, y=0.97, fontsize=18)

				#Plot 2D image, applying image options and cropping as required.
				fig,ax[0],im,Image = ImagePlotter2D(Image,extent,aspectratio,varlist[i],fig,ax[0])
				PlotSheathExtent(SxAxis,Sx,ax[0],ISYMlist[l])
				
				#Add Colourbar (Axis, Label, Bins)
				ImageOptions(fig,ax[0],Xlabel,Ylabel,Crop=True)
				Ylabel = VariableLabelMaker(varlist)
				cax = Colourbar(ax[0],Ylabel[i],5,Lim=Limits)

				#Plot waveform and apply image options.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
				ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
				#Add Invisible Colourbar to sync X-axis
				InvisibleColourbar(ax[0])
				
				#Clean up image and save with relevent filename.
				#NOTE:	zfill assumes phase < 999 degrees	(i.e. < 1e4)
				fig.tight_layout()
				plt.subplots_adjust(top=0.90)
				savefig(DirMovieplots+varlist[i]+'_'+str(Phase).zfill(3)+ext)
				plt.close('all')


				#Write Phase data in ASCII format if required.
				if write_ASCII == True:
					DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'2DPhase_Data')
					DirASCIIPhaseVar = CreateNewFolder(DirASCIIPhase,varlist[i])
					Cycle = str( Phaselist[j].replace(" ", "") )
					SaveString = DirASCIIPhaseVar+varlist[i]+'_'+Cycle
					WriteDataToFile(Image, SaveString)
				#endif
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			MakeMovie(DirMovieplots,Prefix+'_'+varlist[i])
		#endfor
	#endfor

	print('---------------------------------------')
	print('# 2D Phase-Resolved Processing Complete')
	print('---------------------------------------')
#endif




#====================================================================#
				#PHASE TRENDS & SHEATH DYNAMICS#
#====================================================================#

#Process phase resolved data from multiple folders to extract trends.
if savefig_sheathdynamics == True:

	#Initiate arrays between folders
	VoltageWaveforms,WaveformBiases = list(),list()
	SxMaxExtTrend,SxMeanExtTrend= list(),list()
	SxMaxVelTrend,SxMeanVelTrend = list(),list()
	SxDynRangeTrend = list()
	VariedValuelist = list()

	#Read data from each simulation folder individually, saves on RAM.
	for l in tqdm(range(0,numfolders)):

		#Create global folders to keep output plots.
		TrendVariable = list(filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0])))	#List of discrete chars
		TrendVariable = ''.join(TrendVariable)												#Single string of chars
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'2DPhase/')
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')
		DirSheath = CreateNewFolder(DirTrends,'Sheath Trends')

		#Create processlist for each folder as required. (Always get 'E','AR+','PPOT')
		SxData,SxPhase,Sxproc,Sxvar = ExtractPhaseData(folder=l,Variables=PhaseVariables+['E','AR+','PPOT'])
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Extract waveforms from desired electrode locations.
		PPOT = Sxproc[Sxvar.index('PPOT')]
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(SxData,PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(SxData,PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(SxData,PPOT)

		### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###
		#Axial sheath array (Sx) is calculated using radial integrations.
		#Radial sheath array (Sx) is calculated using axial integrations.
		Orientation = 'Axial'
		if Orientation == 'Axial': loc = electrodeloc[1]		#Axial depth where sheath plotted
		elif Orientation == 'Radial': loc = electrodeloc[0]		#Radial depth where sheath plotted
		Phaseaxis = GenerateAxis('Phase',Isym=ISYMlist[l])

		#=============#

		SxLoc = list()
		#For all phases, calculate sheath width and record sheath width at electrodeloc
		for k in range(0,len(SxPhase)):
			#Extract Ni and Ne variables for sheath processing.
			Ne = SxData[k][Sxproc[Sxvar.index('E')]]
			Ni = SxData[k][Sxproc[Sxvar.index('AR+')]]

			#calculate sheath width employing 'E' and 'AR+'
			Sx,SxAxis = CalcSheathExtent(folder=l,Phase=j,Ne=Ne,Ni=Ni)
			for j in range(0,len(Sx)): 
				Sx[j] = ((sourcewidth[0]*dr[l])-Sx[j])*10	#Convert to mm
			#endfor
			
			SxLoc.append(Sx[loc])
		#endfor

		#Determine phase-averaged sheath proprties, removing 'nans' unless all data is 'nan'
		SxLocNoNaN = [x for x in SxLoc if np.isnan(x) == False]
		if len(SxLocNoNaN) > 0:
			SxDynRangeTrend.append(max(SxLocNoNaN)-min(SxLocNoNaN))		#Dynamic Range
			SxMeanExtTrend.append(sum(SxLocNoNaN)/len(SxLocNoNaN))		#Mean Extension
			SxMaxExtTrend.append(max(SxLocNoNaN))						#Max Extension
		else:
			SxDynRangeTrend.append( np.nan )
			SxMeanExtTrend.append( np.nan )
			SxMaxExtTrend.append( np.nan )
		#endif

		#Calculate phase-averaged (mean) sheath velocity.
		#Assumes one sheath collapse and one sheath expansion per rf-cycle.
		RFPeriod = 1.0/FREQMIN[l]										#[s]
		MeanSheathExtent = (sum(SxLoc)/len(SxLoc))/1E6					#[km]
		MeanSheathVelocity = (2*MeanSheathExtent)/RFPeriod				#[km/s]
		SxMeanVelTrend.append( MeanSheathVelocity )						#[km/s]		

		#Calculate maximum instantaneous sheath velocity.
		#Assumes sheath collapse velocity > sheath expansion velocity
		Collapsed,CollapsedPhase = min(SxLoc),SxLoc.index(min(SxLoc))
		Extended,ExtendedPhase = max(SxLoc),SxLoc.index(max(SxLoc))

		SheathExtension = (Extended-Collapsed)/1000.0  					#[m]
		PhaseResolution = 1.0/(FREQMIN[l]*len(Phaseaxis))				#[s]
		SheathTime = (ExtendedPhase-CollapsedPhase)*PhaseResolution		#[s]
		try:
			MaxSheathVelocity = SheathExtension/SheathTime				#[m/s]
			SxMaxVelTrend.append(MaxSheathVelocity/1000.0)				#[km/s]
		except:
			SxMaxVelTrend.append(0.0)									#[km/s]
		#endtry

		#=============#

		#Scale sheath extension by required number of phasecycles.
		ScaledSxLoc = list()
		for n in range(0,int(phasecycles*len(SxLoc))):
			Index = n % len(SxLoc)		#Modulo index for multiple phasecycles
			ScaledSxLoc.append(SxLoc[Index])
		#endfor
		SxLoc=ScaledSxLoc

		#Print results to terminal if requested.
		if print_sheath == True:
			print(Dirlist[l][2:-1])
			print('Sheath Extension:',round(SheathExtension*1E3,2),'[mm]')
			print('Sheath Collapse Time:',round(SheathTime*1E9,2),'[ns]')
			print('Mean Sheath Velocity:',round(MeanSheathVelocity/1E3,2),'[km/s]')
			print('Max Sheath Velocity:',round(MaxSheathVelocity/1E3,2),'[km/s]')
			print('')
		#endif

		#=============#

		#Plot phase-resolved sheath extension [ax0] and voltage waveform [ax1] for current folder
		fig,ax = figure(image_aspectratio,2,shareX=True)
		ax[0].plot(Phaseaxis,SxLoc, lw=2)		
		Title = 'Phase-Resolved Sheath Extension for '+VariedValuelist[l]
		Ylabel = 'Sheath Extension [mm]'
		ImageOptions(fig,ax[0],Title=Title,Ylabel=Ylabel,Crop=False)
		
		ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
		ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)

		plt.savefig(DirPhaseResolved+VariedValuelist[l]+' SheathDynamics'+ext)
		plt.close('all')

		#Write phase-resolved sheath dynamics to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirPhaseResolved,'2DPhase_Data')
			DirASCIISheath = CreateNewFolder(DirASCII,'SheathDynamics')
			WriteDataToFile(Phaseaxis+['\n']+SxLoc, DirASCIISheath+VariedValuelist[l]+'SheathDynamics')
		#endif
	#endfor

	#Write sheath trends to ASCII format datafile if requested.
	if write_ASCII == True:
		DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
		DirASCIISheath = CreateNewFolder(DirASCII,'Sheath_Data')
		WriteDataToFile(VariedValuelist+['\n']+SxMaxExtTrend, DirASCIISheath+'MaxExtent_Trends')
		WriteDataToFile(VariedValuelist+['\n']+SxMeanExtTrend, DirASCIISheath+'MeanExtent_Trends')
		WriteDataToFile(VariedValuelist+['\n']+SxDynRangeTrend, DirASCIISheath+'DynamicRange_Trends')
		WriteDataToFile(VariedValuelist+['\n']+SxMaxVelTrend, DirASCIISheath+'MaxVelocity_Trends')
		WriteDataToFile(VariedValuelist+['\n']+SxMeanVelTrend, DirASCIISheath+'MeanVelocity_Trends')
	#endif


	#Plot maximum sheath extension trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMaxExtTrend,VariedValuelist,NormFactor=0)
	Title='Maximum Sheath Extension W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Max Sheath Extension [mm]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Max Sheath Extension Trends'+ext)
	plt.close('all')

	#Plot mean sheath extension trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMeanExtTrend,VariedValuelist,NormFactor=0)
	Title='Mean Sheath Extension W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Mean Sheath Extension [mm]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Mean Sheath Extension Trends'+ext)
	plt.close('all')

	#Plot sheath dynamic range (max-min extension) trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxDynRangeTrend,VariedValuelist,NormFactor=0)
	Title='Sheath Dynamic Range W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Dynamic Range [mm]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Sheath Dynamic Range Trends'+ext)
	plt.close('all')

	#Plot maximum sheath velocity trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMaxVelTrend,VariedValuelist,NormFactor=0)
	Title='Maximum Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Max Sheath Velocity Trends'+ext)
	plt.close('all')

	#Plot mean sheath velocity trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMeanVelTrend,VariedValuelist,NormFactor=0)
	Title='Phase-averaged Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Mean Sheath Velocity Trends'+ext)
	plt.close('all')


	print('----------------------------------------')
	print('# Phase-Resolved Trend Analysis Complete')
	print('----------------------------------------')
#endif



#====================================================================#
					#SIMULATED PROES DIAGNOSTIC#
#====================================================================#

#Plot Phase-Resolved profiles with electrode voltage and requested variables.
if savefig_PROES == True:
	VariedValuelist = list()

	#for all folders.
	for l in range(0,numfolders):

		#Create global folders to keep output plots and collect graph title.
		VoltageWaveforms,WaveformBiases = list(),list()
		DirPROES = CreateNewFolder(Dirlist[l],'PROES')

		#Update X-axis with folder information.
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get 'E','AR+','PPOT')
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables)
		SxData,SxPhase,Sxproc,Sxvar = ExtractPhaseData(folder=l,Variables=['E','AR+'])
		PPOT = ExtractPhaseData(folder=l,Variables=['PPOT'])[2][0]

		#Generate SI scale axes for lineout plots.
		Phaseaxis = GenerateAxis('Phase',ISYMlist[l],Phaselist)		#[omega*t/2pi]
		Raxis = GenerateAxis('Radial',ISYMlist[l])					#[cm]
		Zaxis = GenerateAxis('Axial',ISYMlist[l])					#[cm]

		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias,ElectrodeVpp = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)): 
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['Waveform Vpp: '+str(round(ElectrodeVpp[2],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))

		plt.savefig(DirPROES+VariedValuelist[l]+' Waveform'+ext)
		plt.close('all')

		#Write PROES data in ASCII format if required.
		if write_ASCII == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPROES = CreateNewFolder(DirPROES,'PROES_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPROES+'VoltageWaveforms')
		#endif

		#==============#


		#for all requested variables.
		for i in tqdm(range(0,len(proclist))):
			#Refresh lineout lists between variables.
			Lineouts,ProfileOrientation = list(),list()

			#Concatinate all requested lineouts together, keeping seperate orientation.
			for m in range(0,len(radialprofiles)):
				Lineouts.append(radialprofiles[m])
				ProfileOrientation.append('Radial')
			#endfor
			for m in range(0,len(axialprofiles)):
				Lineouts.append(axialprofiles[m])
				ProfileOrientation.append('Axial')
			#endfor

			#For all requested profiles and orientations.
			for k in range(0,len(Lineouts)):
				#Refresh required lists.
				VariableMax,VariableMin = list(),list()
				PhaseSx = list()
				PROES = list()

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Phaselist)):
					#Refresh lists between each phasecycle.
					IntegratedDoFArray,DoFArrays = list(),list()

					#Collect each profile for stitching into a PROES image if required.
					if DoFwidth > 0:
						#Determine range of lineouts within the depth of field.
						DOFRegion = [(Lineouts[k]-DoFwidth),(Lineouts[k]+DoFwidth)]
						#If DOF extends beyond mesh, alert user and abort diagnostic.
						if any(DOFRegion) < 0:
							print('----------------------------------')
							print('Depth-of-Field Extends Beyond Mesh')
							print('----------------------------------')
							break
						#endif

						#Collect profiles from DOF region and transpose to allow easy integration.
						for LineoutLoc in range(DOFRegion[0],DOFRegion[1]):
							if ProfileOrientation[k] == 'Radial': DoFArrays.append(ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							elif ProfileOrientation[k] == 'Axial': DoFArrays.append(ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							#endif
						#endfor
						DoFArrays = np.asarray(DoFArrays).transpose().tolist()

						#Integrate DoF profiles spatially, obtaining PROES profile for phase 'j'
						for n in range(0,len(DoFArrays)):
							IntegratedDoFArray.append( sum(DoFArrays[n])/(DoFwidth*2+1) )
						#endif
						PROES.append(IntegratedDoFArray)

					#If no DoF then simply collect lineout from required location.
					elif DoFwidth == 0:
						LineoutLoc = Lineouts[k]
						if ProfileOrientation[k] == 'Radial':
							PROES.append(ExtractRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc))
						elif ProfileOrientation[k] == 'Axial':
							PROES.append(ExtractAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
						#endif
					#endif

					if image_plotsheath == True:
						#Extract Ni and Ne variables and perform sheath processing.
						Ne = SxData[j][Sxproc[Sxvar.index('E')]]
						Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
						#Extract the full spatial sheath extent for phase 'j' 
						#and save sheath extent at PROES profile location
						Sx,SxAxis = CalcSheathExtent(folder=l,Orientation=ProfileOrientation[k],Phase=j,Ne=Ne,Ni=Ni)
						PhaseSx.append(Sx[Lineouts[k]])
					#endif
				#endfor
				
				#Scale PROES image and Sx arrays by required number of phasecycles.
				ScaledPROES,ScaledPhaseSx,ScaledPhaseSxAxis = list(),list(),list()
				for n in range(0,int(phasecycles*len(PROES))):
					Index = n % len(PROES)		#Modulo index for multiple phasecycles
					ScaledPROES.append(PROES[Index])
					if image_plotsheath == True:
						ScaledPhaseSx.append(PhaseSx[Index])
					#endif
				#endfor
				try: PROES,PhaseSx = ScaledPROES,ScaledPhaseSx
				except: PROES = ScaledPROES

				#Create figure and rotate PROES such that phaseaxis aligns with waveform.
				fig,ax = figure(image_aspectratio,2,shareX=True)
				PROES = ndimage.rotate(PROES, 90)
				#Choose correct axial or radial distance axis and create associated folder.
				x1,x2 = Phaseaxis[0],Phaseaxis[-1]
				if ProfileOrientation[k] == 'Axial':
					ProfileString = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm'
					NameString = varlist[i]+'_'+ProfileString[2::]
					Ylabel = 'Axial Distance Z [cm]'
					Crop = [image_axialcrop[::-1],image_radialcrop] #Reversed accounting for rotation.
					y1,y2 = Zaxis[-1],Zaxis[0]						#Reversed accounting for top origin.
				elif ProfileOrientation[k] == 'Radial' and int(ISYMlist[l]) == 1:
					ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+ProfileString[2::]
					Ylabel = 'Radial Distance R [cm]'
					Crop = [image_radialcrop,image_axialcrop]
					y1,y2 = Raxis[-1],-Raxis[-1]
				elif ProfileOrientation[k] == 'Radial' and int(ISYMlist[l]) == 0:
					ProfileString = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+ProfileString[2::]
					Ylabel = 'Radial Distance R [cm]'
					Crop = [image_radialcrop,image_axialcrop]
					y1,y2 = Raxis[-1],0
				#endif
				DirPROESloc = CreateNewFolder(DirPROES,ProfileString[3::])

				#Create PROES image along line of sight with phase-locked waveform.
				fig.suptitle( 'Simulated '+varlist[i]+' PROES for '+VariedValuelist[l]+ProfileString+'\n DoF = '+str(round(((2*DoFwidth)+1)*dz[l],2))+' cm', y=0.95, fontsize=18)
				im = ax[0].contour(PROES,extent=[x1,x2,y1,y2],origin='lower')
				im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='bottom',aspect='auto')
				PlotSheathExtent(Phaseaxis,PhaseSx,ax[0],ISYMlist[l])
				
				#Beautify Image
				ImageOptions(fig,ax[0],Xlabel='',Ylabel=Ylabel,Crop=Crop)
				ax[0].set_xticks([])
				ax[0].set_xlim(x1,x2)
				#Add Colourbar (Axis, Label, Bins)
				label = VariableLabelMaker(varlist)
				Colourbar(ax[0],label[i],5,Lim=CbarMinMax(PROES,ProfileOrientation[k]))

				#Plot Voltage Waveform.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
				ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
				ax[1].set_xlim(x1,x2)
				#Add Invisible Colourbar to sync X-axis
				InvisibleColourbar(ax[1])

				#Cleanup layout and save images.
				fig.tight_layout()
				plt.subplots_adjust(top=0.85)
				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' PROES'+ext)
				plt.close('all')

				#Write PROES data in ASCII format if required.
				if write_ASCII == True:
					DirASCIIPROES = CreateNewFolder(DirPROES,'PROES_Data')
					DirASCIIPROESloc = CreateNewFolder(DirASCIIPROES,ProfileString[3::])
					WriteDataToFile(PROES, DirASCIIPROESloc+varlist[i]+'_PROES')
				#endif

				#==============##==============#
				#==============##==============#
				
				#Spatially Collapse 2D PROES image along the line of sight.
				PROES,TemporalPROES = PROES.transpose().tolist(),list()
				for m in range(0,len(PROES) ):
					TemporalPROES.append( (sum(PROES[m][::]))/(len(PROES)/2) )
				#endfor

				#Plot Spatially Collapsed PROES with phase axis.
				fig,ax = figure(image_aspectratio,2,shareX=True)
				ax[0].plot(Phaseaxis,TemporalPROES, lw=2)
				Title = 'Spatially Integrated '+varlist[i]+' for '+VariedValuelist[l]+ProfileString+'\n DoF = '+str(round(((2*DoFwidth)+1)*dz[l],2))+' cm'
				Ylabel = 'Spatially Integrated '+varlist[i]
				ImageOptions(fig,ax[0],Title=Title,Ylabel=Ylabel,Crop=False)
#				ax[0].set_xlim(x1,x2)

				#Plot Waveform onto Temporally collapsed PROES.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
#				ax[1].set_xlim(x1,x2)

				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' TemporalPROES'+ext)
				plt.close('all')
			
				#=====#=====#
				
				#					!!!NOTE!!!
				#RM: 	SpatialPROES array length mismatch with R,Z Axes
				#		try/except Hack no longer working as both axes are incorrect
				#		Not sure why data axes is different length, something wrong with sum?
				
				#Temporally collapse 2D PROES image through defined phase fraction.
#				SpatialPROES = list()
#				for m in range(0,len(PROES)):
#					SpatialPROES.append( (sum(PROES[m][::])) ) 		#Needs 'maths'-ing
				#endfor

				#Plot Temporally Collapsed PROES with required axis.
#				fig,ax = figure(image_aspectratio,2,shareX=True)
#				try: ax[0].plot(Raxis,SpatialPROES, lw=2)		### HACKY ###
#				except: ax[0].plot(Zaxis,SpatialPROES, lw=2)	### HACKY ###
#				Xlabel = 'Phase [$\omega$t/2$\pi$]'
#				Ylabel = 'Temporally Integrated '+varlist[i]
#				ImageOptions(fig,ax[0],Xlabel,Ylabel,Crop=False)

				#Plot Waveform onto Spatially collapsed PROES.
#				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
#				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
#				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
#				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)

#				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' SpatialPROES'+ext)
#				plt.close('all')

				#==============##==============#

			if l == numfolders and k == len(Lineouts):
				print('-------------------------------')
				print('# PROES Image analysis Complete')
				print('-------------------------------')
			#endif
		#endfor
	#endfor
#endif

#===============================#

if any([savefig_sheathdynamics, savefig_phaseresolve1D, savefig_phaseresolve2D, savefig_PROES]) == True:
	print('----------------------------------')
	print('# Phase Resolved Profiles Complete')
	print('----------------------------------')
#endif

#=====================================================================#
#=====================================================================#



















































#====================================================================#
							  #CODE DUMP#
#====================================================================#
#
##===============================#
##     Paper Trend Locations     #
##===============================#
#
##SDoyle2017a: 
##dz(5.50/118), dr(2.55/102) height=[24,43], Trend=[19]
#
##SDoyle2018a: 
##PROES, (Z=14.2, 21.0, 31.0) radialprofiles = [29,44,64]
##Dielectric locations: [16,29],[16,44],[16,64]
#
##===============================#
## Archived Switchboard Settings #
##===============================#
#
#### PP-SCCP ####
##electrodeloc = 	[0,3]
##waveformlocs = 	[]
##DoFwidth = 		R;??,Z;??
##TrendLoc =  		H[0];R[]
##thrustloc = 		[]
##sheathROI = 		[]
##sourcewidth = 		[]
##Crop = 			R[];Z[]
#
##### TSHC-2017 ####
##electrodeloc = 	[0,15]
##waveformlocs = 	[]
##DoFwidth = 		R;5,Z;10
##TrendLoc =  		H[0,20];R[30,60,90]
##thrustloc = 		[]
##sheathROI = 		[]
##sourcewidth = 		[]
##Crop = 			R[0.0,1.0];Z[0.5,2.5]
#
#### TSHC-OI Mk3 ###
##electrodeloc = 	[58,15]
##waveformlocs = 	[]
##DoFwidth = 		R;??,Z;??
##TrendLoc =  		H[0,23,45];R[46,55,64]			#R,Z = 0.2cm/cell,0.1cm/cell
##thrustloc = 		[]
##sheathROI = 		[]
##sourcewidth = 		[]
##Crop = 			R[0,12];Z[4,7]
#
#### HYPERION-I Mk1 ###
##electrodeloc = 	[51,14]HYPI OR [12,28]HYPII			#Upper(Positive) ICP coil
##waveformlocs = 	[[51,24][51,34]]					#Middle ICP Coil, Lower(Negative) coil
##DoFwidth = 		R;??,Z;??
##TrendLoc =  		H[0];R[50]HYPI OR H[56];R[50]HYPII	#R,Z = 0.2cm/cell,0.1cm/cell
##thrustloc = 		[]
##sheathROI = 		[]
##sourcewidth = 		[]
##Crop = 			R[];Z[]
#
#### EVgeny Mk1 ###
##electrodeloc = 	[31,14]							#Middle ICP Coil
##waveformlocs = 	[[31,6],[31,23],[20,6]]			#[UpstreamCoil],[DownstreamCoil],[DielectricSurface]
##DoFwidth = 		R;??,Z;??
##TrendLoc =  		H[0,21,41];R[]					#R,Z = 0.2cm/cell,0.1cm/cell
##thrustloc = 		[]
##sheathROI = 		[45,85]							#Downstream
##sourcewidth = 		[90]							#Downstream
##Crop = 			R[];Z[]
##Plotmesh = 		'EVgeny'
#
##===============================#
##             Notes             #
##===============================#
#
## Disabled the following warning message regarding scalar assignment of 'arr' axis.
## /home/sjd549/.local/lib/python2.7/site-packages/numpy/ma/core.py:6385
#
#
##====================================================================#
#				  	#GRAPHICAL USER INTERFACE#
##====================================================================#
#
#
#use_GUI = False
#if use_GUI == True:
#	try:
#		# Python2
#		import Tkinter as tk
#		from ttk import *
#	except ImportError:
#		# Python3
#		import tkinter as tk
#		from ttk import *
#	#endtry
#
#	#Create switchboard directory for GUI.
#	Switchboard = {}
#
##=============#
#
#	#Toggles button font between green/red and outputs corresponding true/false.
#	def toggle(string,btn):
#		#to get the present state of the toggle button
#		if btnlist[btn].config('fg')[-1] == 'green':
#			btnlist[btn].config(fg='red')
#		else:
#			btnlist[btn].config(fg='green')
#		#endif
#
#		#Add name of button and true/false value to switchboard dictionary.
#		global Switchboard
#		if btnlist[btn].config('fg')[-1] == 'green':
#			Switchboard[string] = 'True'
#		elif btnlist[btn].config('fg')[-1] == 'red':
#			Switchboard[string] = 'False'
#		#endif
#	#enddef
#
##=============#
#
#	#Initialize and configure root window.
#	root = tk.Tk()
#	root.title("Hpem ELectronic ENgine Analysis v0.8.4")
#	root.minsize(width=800, height=600)
#	root.maxsize(width=4*1080, height=1080)
#
#	#Buttonlist as displayed in GUI.
#	btnlist = ['2D Images','2D Converge','2D Phase']
#
#	#Add toggle buttons.  (lambda returns a reference to a nameless function)
#	btnlist[0] = tk.Button(text=btnlist[0], width=12, fg='red')
#	btnlist[0]["command"] = lambda: toggle('savefig_plot2D',0)
#	btnlist[0].grid(row=1, column=0)
#
#	btnlist[1] = tk.Button(text=btnlist[1], width=12, fg='red')
#	btnlist[1]["command"] = lambda: toggle('savefig_convergence',1)
#	btnlist[1].grid(row=1, column=1)
#
#	btnlist[2] = tk.Button(text=btnlist[2], width=12, fg='red')
#	btnlist[2]["command"] = lambda: toggle('savefig_phaseresolve',2)
#	btnlist[2].grid(row=1, column=2)
#
#	#Add run button, disables GUI and progresses to the main program.
#	Run = tk.Button(text='Run Analysis', height=4, width=18, fg='black')
#	Run["command"] = root.destroy
#	Run.grid(row=3,column=3,padx=0,pady=100)
#
#	root.mainloop()
#	print Switchboard
##endif
#
##=====================================================================#
##=====================================================================#
#
#
#
#
#
#
#
#
