#!/usr/bin/env python

#################################
#		Point of Contact		#
#								#
#	   Dr. Scott J. Doyle		#
#	   University of York		#
#	   York Plasma Institute	#
#	   1&2 Genesis Building		#
#	   North Yorkshire, UK		#
#	  Scott.Doyle@Physics.org	#
#								#
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

	print ''
	print 'First time use requires installation of additional python modules'
	print 'Please type your password when prompted to allow installation:'
	print ''
	try:
		os.system('sudo apt-get install python-pip')
		os.system('sudo apt-get install python-matplotlib')
		os.system('sudo apt-get install python-numpy')
		os.system('sudo apt-get install python-scipy')
		os.system('sudo apt-get install ffmpeg')
		os.system('pip install tqdm')
	except:
		print ''
		print 'Error installing required packages'
		print 'Please attempt manual installation'
		print ''
	#endtry
	print ''
	print ''
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

#Enforce matplotlib to avoid instancing undisplayed windows
#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib			#matplotlib.use('Agg')

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *
 

#====================================================================#
				  		 #LOW LEVEL INPUTS#
#====================================================================#

#Various debug and streamlining options.
Magmesh = 1							#initmesh.exe magnification factor. (almost obsolete)
ffmpegMovies = True					#If False: Suppresses ffmpeg routines, saves RAM.
DebugMode = False					#Produces debug outputs for relevent diagnostics.
QuickConverge = False				#Supresses 2D Convergence images in savefig_convergence

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#List of recognized data extensions for file readin
FileExtensions = ['.PDT','.pdt','.nam','.dat','.out']

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

#Define units for particular variables
PressureUnit = 'Torr'#'Pa'							#'Torr','mTorr','Pa'
BFieldUnit	=  'Gauss'						#'Gauss','Tesla'


####################

#Commonly used variable sets.
Phys = ['P-POT','TE','EF-TOT','EAMB-Z','EAMB-R','RHO','BR','BRS','BZ','BZS','BT','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','EFLUX-R','EFLUX-Z','JZ-NET','JR-NET','TG-AVE','PRESSURE','POW-RF','POW-RF-E','POW-ICP','EB-ESORC','COLF']
Ar = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR','AR+','AR2+','AR2*','E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','FZ-AR3S','FR-AR3S','FR-AR+','FZ-AR+','FZ-AR3S','FR-AR3S']+Phys
O2 = ['O3','O2','O2+','O','O+','O-','E','S-O3','S-O2+','S-O+','S-O-','SEB-O3','SEB-O+','SEB-O2+','SEB-O-','FR-O+','FZ-O+','FR-O-','FZ-O-']+['O3P3P','O***','S-O3P3P','S-O***','SEB-O3P3P','SEB-O***']+Phys

Ar_Phase = ['S-E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','SRCE-2437','TE','PPOT','FR-E','FZ-E']
O2_Phase = ['S-E','S-O+','S-O-','S-O2+','SEB-O+','SEB-O-','SEB-O2+','TE','PPOT','FR-E','FZ-E']+['S-O3P3P','SEB-O3P3P']

PRCCPAr_PCMC = ['AR^0.35','EB-0.35','ION-TOT0.35']
PRCCPO2_PCMC = ['O^0.35','EB-0.35','ION-TOT0.35']
ESCTAr_PCMC = ['TO BE COMPLETED']
TSHCMk4a_PCMC = ['AR^1.1T','EB-1.1T','ION-TOT1.1T','AR^5.4U','EB-5.4U','ION-TOT5.4U','AR^9.8V','EB-9.8V','ION-TOT9.8V']

EVgeny_PCMC = ['AR^0.1P','EB-0.1P','ION-TOT0.1P','AR^2.1Q','EB-2.1Q','ION-TOT2.1Q']
HYPI_PCMC = ['O^0.2P','EB-0.2P','ION-TOT0.2P','O^4.1Q','EB-4.1Q','ION-TOT4.1Q','O^4.15','EB-4.15','ION-TOT4.15']
HYPII_PCMC = ['O^2.8P','EB-2.8P','ION-TOT2.8P','O^3.5Q','EB-3.5Q','ION-TOT3.5Q']


#Archived variable sets
TSHCOI2019_PCMC = ['AR^0.2S','ION-TOT0.2S','AR^4.4T','ION-TOT4.4T','AR^8.9U','ION-TOT8.9U']
TSHCOI2020_PCMC = ['AR^0.2S','ION-TOT0.2S','AR^4.9T','ION-TOT4.9T','AR^9.8U','ION-TOT9.8U']
MSHC2017Mk0_PCMC = ['AR^0.5S','EB-0.5S','ION-TOT0.5S','AR^1.1B','EB-1.1B','ION-TOT1.1B']
MSHC2021Mk4_PCMC = ['AR^0.2S','EB-0.2S','ION-TOT0.2S','AR^0.2T','EB-0.2T','ION-TOT0.2T']
SCCP2018Mk0_PCMC = ['AR^7.7J','ION-TOT7.7J','AR^5.1B','ION-TOT5.1B']
ESCT2018Mk0_PCMC = ['AR^0.3S','EB-0.3S','ION-TOT0.3S']

####################



#Commonly Used Diagnostic Settings:
#### PRCCP ####
#electrodeloc =		[29,44] 					#Reverse [29,62]
#waveformlocs =		[16,29],[16,44],[16,64],[0,29],[0,44],[0,64]]
#DOFWidth =			R;16,Z;21
#TrendLoc =			H[0];R[29,44,64,75]			#R[44] for PROES
#ThrustLoc =		75, 						#stdESCT=76, smlESCT=48/54
#SheathROI =		[34,72]
#SourceWidth =		[0.21]						
#image_radialcrop = [0.65]						#[R1,R2] in cm
#image_axialcrop = 	[1.0,4.0]					#[Z1,Z2] in cm

#### PRuICP ####	
#electrodeloc = 	[33,33]			#Coil V
#waveformlocs = 	[]
#DOFWidth = 		[]
#TrendLoc = 		H[0];R[36,50]
#ThrustLoc = 		[79]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[0.0,1.0];Z[1.5,10]
#Plotmesh = 		'PRCCP'

#### TSHC-Ar ####
#electrodeloc = 	[5,20]
#waveformlocs = 	[]
#DOFWidth = 		R;20,Z;6
#TrendLoc =  		H[1,22,44,60];R[35,40,44]		#R[39] for PROES
#ThrustLoc = 		[45]
#SheathROI = 		[]
#SourceWidth = 		[5]
#Crop = 			R[0,14];Z[7,13]					#Mk3: R[0,15];Z[10,17.5]

#### TSHC-Mk4PPa ####
#electrodeloc = 	[5,30]
#waveformlocs = 	[]
#DOFWidth = 		R;20,Z;5
#TrendLoc =  		H[0,24,48];R[34,39,44,46]		#R[40] for PROES
#ThrustLoc = 		[45]
#SheathROI = 		[]
#SourceWidth = 		[5]
#Crop = 			R[0,14];Z[5,10]					#R[2,14] to avoid on-axis craziness

#### SERPENT ####	
#electrodeloc = 	[33,33]			#Coil V
#waveformlocs = 	[]
#DOFWidth = 		R;60,Z;5
#TrendLoc = 		H[0];R[36,50]
#ThrustLoc = 		[79]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[1.4];Z[1,9]
#Plotmesh = 		False

####################




#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested IEDF/NEDF Variables.
IEDFVariables = TSHCMk4a_PCMC				#Requested iprofile_2d variables (no spaces)
NEDFVariables = []							#Requested nprofile_2d variables (no spaces)

#Requested movie1/movie_icp Variables.
IterVariables = ['E','S-E','PPOT','TE']				#Requested Movie_icp (iteration) Variables.		
PhaseVariables = Ar_Phase							#Requested Movie1 (phase) Variables. +['E','AR+']
electrodeloc = [30,15]								#Cell location of powered electrode [R,Z].
waveformlocs = []									#Cell locations of additional waveforms [R,Z].

#Requested TECPLOT Variables and plotting locations.
Variables = Ar
MultiVar = []							#Additional variables plotted ontop of [Variables]
radialineouts = []#[85]		 			#Radial 1D-Profiles to be plotted (fixed Z-mesh) --
heightlineouts = [33]#[33] 	#[66]?		#Axial 1D-Profiles to be plotted (fixed R-mesh) |
TrendLocation = [] 						#Cell location For Trend Analysis [R,Z], ([] = min/max)


#Various Diagnostic Settings.
phasecycles = 1.00						#Number of waveform phase cycles to be plotted. (float)
DoFWidth = 10 				#20?		#PROES Depth of Field (symmetric about image plane) (cells)
ThrustLoc = 45							#Z-axis cell for thrust calculation  (cells)
SheathROI = [34,72]						#Sheath Region of Interest, (Start,End) [cells]
SourceWidth = [12]						#Source Dimension at ROI, leave empty for auto. [cells]
EDF_Threshold = 0.01					#Maximum Recognised EEDF/IEDF energy fraction (Plot all: 0.0)


#Requested diagnostics and plotting routines.
savefig_convergence = False				#Requires movie_icp.pdt
savefig_plot2D = True					#Requires TECPLOT2D.PDT

savefig_monoprofiles = False			#Single-Variables; fixed height/radius
savefig_multiprofiles = False			#Multi-Variables; same folder
savefig_comparelineouts = False			#Multi-Variables; all folders
savefig_trendphaseaveraged = False		#Single-Variables; fixed cell location (or max/min)
savefig_trendphaseresolved = False		#Single-Variables; Phase-resolved data.
savefig_pulseprofiles = False			#Single-Variables; plotted against real-time axis

savefig_phaseresolve1D = False			#1D Phase Resolved Images
savefig_phaseresolve2D = False			#2D Phase Resolved Images
savefig_PROES = False					#Simulated PROES Diagnostic

savefig_IEDFangular = False				#2D images of angular IEDF; Single simulation directory
savefig_IEDFtrends = False				#1D IEDF trends; 			All simulation directories
savefig_EEDF = False					#NO PLOTTING ROUTINE		#IN DEVELOPMENT#

#Write processed data to ASCII files.
write_ASCII = True						#All diagnostic output written to ASCII.


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs.
print_Knudsennumber = False				#Print cell averaged Knudsen Number
print_soundspeed = False				#Print cell averaged sound speed
print_totalpower = False				#Print all requested total powers
print_DCbias = False					#Print DC bias at electrodeloc
print_thrust = False					#Print neutral, ion and total thrust
print_sheath = False					#Print sheath width at electrodeloc


#Image plotting options.
image_extension = '.png'				#Extensions ('.png', '.jpg', '.eps')
image_aspectratio = [10,10]				#[x,y] in cm [Doesn't rotate dynamically]
image_radialcrop = []#[0,12.5]			#[R1,R2] in cm								CROPPING IS TOTALLY FUCKED
image_axialcrop = []#[5,10]				#[Z1,Z2] in cm								CROPPING IS TOTALLY FUCKED
image_cbarlimit = []					#[min,max] colourbar limits	

image_plotsymmetry = False#True			#Toggle radial symmetry
image_numericaxis = False				#### NOT IMPLIMENTED ####
image_contourplot = False#True			#Toggle contour Lines in images
image_1Doverlay = False					#Overlay location(s) of radialineouts/heightlineouts
image_plotgrid = False					#Plot major/minor gridlines on profiles
image_plotmesh = False					#Plot material mesh outlines ('Auto','PRCCP','ESCT')
image_rotate = False#True				#Rotate image 90 degrees to the right.

image_normalize = False					#Normalize image/profiles to local max
image_logplot = False					#Plot ln(Data), against linear axis.
image_sheath = False#True				#Plot sheath width onto 2D images.


#Overrides the automatic image labelling.
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = ['NotImplimented']





#============================#




#V1.1.0 Release Version To Do list:
#Fix sound speed diagnostic - Current version has issue with NaNs and incorrect averaging

#Clarified SheathWidth function Axial/Radial definition, also corrected this in the phase-resolved sheath trends diagnostic.

#Corrected 1DPhaseMovie, 2DPhaseMovie and PROES, however the radial direction is not consistent. 
#1DPhaseMovie needs no reversal using plotradialprofile
#2PhaseMovie needs no reversal using ImageExtractor2D
#PROES requires a reversal using plotradialprofile
#Refactor PROES into 3D array (R,Z,Phase) and perform slice/integration rather than calling plotradial?

#Sheathwidth function integrates axially or radially depending on mesh geometry.
#Sheathwidth function has 1D (Scott) and 2D (Greg) capabilities
#SheathWidth function can deal with image rotations.

#Thrust diagnostic split into functions performing the same task as before.
#Thrust diagnostic enforces image symmetry, correcting the half-thrust error.

#Automate AutoConvProfData() and provide a timeout setting. 
#Automate IEDFVarArgs input for varying PCMC settings.

#IEDF diagnostic capable of comparing between different material surfaces in single image
#IEDF diagnostic saves different material surfaces in different folders

#Variable Interpolator needs to work with phasedata - Take variables from batch?

#Fix issue with ffmpeg "convert-im6.q16: DistributedPixelCache" and address the ignorance of os.system.














#####TODO#####

#For Future:
#include 'garbage ./collection' at the end of each diagnostic.
#Correct data 'direction' in readin functions, not diagnostics.
#Include Andor, OO and LeCroy readin functions.
#Namelist read-in functions to simplify error analysis. (deal with !!! in .nam)

#Impliment image_numericaxis, try float(FolderNameTrimmer) as axis.
#Impliment numerical image_rotate, allow for 000,090,180,270.
#Functionalise PROES images, 3D PROES cube model?
#Remove/simplify the phase-averaged sheath diagnostic?
#Add EEDF section and Functionalise.

#Clean up unused functions and ensure homogeneity.
#Update README, include all diagnostics and examples.
#Create developer handbook describing functions.
#Python 3.x compatable.





#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()					#List of all datafile directories inside folders
Dirlist = list()				#List of all folder directories (not including filenames)
IEDFVariablelist = list()		#List of all variable names in pcmc.prof in header order
Geometrylist = list()			#List containing commonly used geometries [LEGACY: NOT USED]

Globalvarlist = list()			#List of all commonly shared variable names between all folders
Globalnumvars = list()			#Number of commonly shared variables between all folders.

#Create mesh_size lists and SI conversion
Isymlist = list()				#Boolian list of ISYM values in folder order in in Dirlist
R_mesh = list()					#List of radial mesh cells for initmesh.out in folder order in Dirlist
Z_mesh = list()					#List of axial mesh cells for initmesh.out in folder order in Dirlist
Raxis = list()					#Radial SI [cm] axis for plotting
Zaxis = list()					#Axial SI [cm] axis for plotting

Depth = list()					#icp.nam Depth input [cm] in folder order in Dirlist
Radius = list()					#icp.nam Radius input [cm] in folder order in Dirlist
Height = list()					#icp.nam Height input [cm] in folder order in Dirlist
dr = list()						#Radial mesh resolution [cm/cell] in folder order in Dirlist
dz = list()						#Axial mesh resolution [cm/cell] in folder order in Dirlist

#Lists for icp.nam variables
VRFM,VRFM2 = list(),list()
FREQM,FREQM2 = list(),list()
FREQC        = list()
FREQGLOB,IRFPOW = list(),list()
MAXFREQ,MINFREQ = list(),list()
PRESOUT = list()
IETRODEM = list()
IMOVIE_FRAMES = list()

#Lists for icp.dat variables
header_icpdat = list()			#[SpeciesName, Charge, MolecularWeight, StickingCoeff,
								# Transport, ReturnFrac, ReturnName]
AtomicSpecies = list()			#All species contained within chemistry set
FluidSpecies  = list() 			#All 'bulk' fluid species (for fluid dynamics analysis)
NeutSpecies	= list()			#All neutral and metastable species
PosSpecies = list()				#All positive ion species
NegSpecies = list()				#All negative ion species

#Lists to store raw data
rawdata_2D = list()				#ASCII format TECPLOT2D data string list 		  - Variable,Radius,Axis
rawdata_kin = list()			#ASCII format kin.pdt data string list 			  - Variable,Radius,Axis
rawdata_phasemovie = list()		#ASCII format movie1.pdt data string list 		  - Variable,Radius,Axis
rawdata_itermovie = list()		#ASCII format movie_icp.pdt data string list 	  - Variable,Radius,Axis
rawdata_IEDF = list()			#ASCII format iprofile_tec2d.pdt data string list - Variable,Radius,Axis
rawdata_mcs = list()			#ASCII format mcs.pdt data string list 			  - Variable,Radius,Axis

Data = list()					#Data[folder][Variable][Datapoint]						-2D Data
DataKin = list()				#Data[folder][Variable][Datapoint]						-1D Data
DataIEDF = list()				#Data[folder][Variable][Datapoint]						-2D Data
DataEEDF = list()				#Data[folder][Variable][Datapoint]						-2D Data
IterMovieData = list()			#ITERMovieData[folder][timestep][variable][datapoints]	-3D Data
PhaseMovieData = list()			#PhaseMovieData[folder][timestep][variable][datapoints]	-3D Data

Moviephaselist = list()			#'CYCL = n'
MovieIterlist = list()			#'ITER = n'
EEDF_TDlist = list()			#'???'

header_itermovie = list()
header_phasemovie = list()
header_IEDFlist = list()
header_kinlist = list()
header_2Dlist = list()









#====================================================================#
					#WELCOME TEXT AND INFORMATION#
#====================================================================#

print ''
print '--------------------------------------------------------------------'
print '    __    __   _______  __       _______  __   __      ___          '
print '   |  |  |  | |   ____||  |     |   ____||  \ |  |    /   \         '
print '   |  |__|  | |  |__   |  |     |  |__   |   \|  |   /  ^  \        '
print '   |   __   | |   __|  |  |     |   __|  |  . `  |  /  /_\  \       '
print '   |  |  |  | |  |____ |  `----.|  |____ |  |\   | /  _____  \      '
print '   |__|  |__| |_______||_______||_______||__| \__|/__/     \__\     '
print '                                                              v1.0.4'
print '--------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if savefig_plot2D == True:
	print'# 2D Steady-State Image Processing'
if savefig_convergence == True:
	print'# 2D Convergence Movie Processing'
if True in [savefig_phaseresolve2D,savefig_PROES]:
	print'# 2D Phase-Resolved Movie Processing'
if True in [savefig_phaseresolve1D]:
	print'# 1D Phase-Resolved Profile Processing'
if True in [savefig_monoprofiles,savefig_multiprofiles,savefig_comparelineouts,savefig_pulseprofiles]:
	print'# 1D Steady-State Profile Processing'
if True in [print_generaltrends,print_Knudsennumber,print_soundspeed, print_totalpower,print_DCbias,print_thrust]:
	print'# 1D Specific Trend Analysis'
if savefig_trendphaseaveraged == True:
	print'# 1D Steady-State Trend Processing'
if savefig_trendphaseresolved == True:
	print'# 1D Phase-Resolved Trend Processing'
if True in [savefig_IEDFangular,savefig_IEDFtrends,savefig_EEDF]:
	print'# Angular Energy Distribution Processing'
print '-----------------------------------------'
print ''





#====================================================================#
					#OBTAINING FILE DIRECTORIES#
#====================================================================#

#Obtain system RAM. (and rename enviroment variable)
mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
mem_gib = mem_bytes/(1024.**3)
ext = image_extension

#Create Directory lists and initialise numfolders to zero.
Dirlist = list() 		#List of all folders
Dir = list() 			#Directory of files within folders
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

#Determine number of folders containing accepted file extensions (i.e. data folders)
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
	print '-------------------------------------------'
	print 'No Ouput Files Detected, Aborting Analysis.'
	print '-------------------------------------------'
	print ''
	exit()
#endif


#Begin the retrieval of geometry from mesh and input files.
icpnam = filter(lambda x: 'icp.nam' in x, Dir)
icpdat = filter(lambda x: 'icp.dat' in x, Dir)
icpout = filter(lambda x: 'icp.out' in x, Dir)
mesh = filter(lambda x: 'initmesh.out' in x, Dir)
TEC2D = filter(lambda x: 'TECPLOT2D.PDT' in x, Dir)

#Loop over all folders and retrieve mesh sizes and SI sizes.
for l in range(0,numfolders):
	
	#==========##===== INITMESH.OUT READIN =====##==========#
	#==========##===============================##==========#

	#Attempt automated retrieval of mesh sizes.
	try:
		#Identify mesh size from TECPLOT2D file. (Data Array Always Correct Size)
		meshdata = open(TEC2D[l]).readlines()

		#Zone line holds data, split at comma, R&Z values are given by "I=,J=" respectively.
		R = filter(lambda x: 'ZONE' in x, meshdata)[0].split(",")[0].strip(' \t\n\r,=ZONE I')
		Z = filter(lambda x: 'ZONE' in x, meshdata)[0].split(",")[1].strip(' \t\n\r,=ZONE J')
		R_mesh.append( int(filter(lambda x: x.isdigit(), R)) )
		Z_mesh.append( int(filter(lambda x: x.isdigit(), Z)) )

	except ValueError:
		#Identify mesh size from initmesh.out file. (Issues with Q-VT and Magmesh)
		meshdata = open(mesh[l]).readline()
		R_mesh.append([int(i) for i in meshdata.split()][1])
		if Magmesh == 1: Z_mesh.append([int(i)+1 for i in meshdata.split()][3])
		elif Magmesh == 2: Z_mesh.append([int(i)+3 for i in meshdata.split()][3])
		elif Magmesh == 3: Z_mesh.append([int(i)+5 for i in meshdata.split()][3])
		#endif

	except:
		#If data for current file exists, ask for manual input.
		if l <= len(TEC2D)-1:

			#If the initmesh.out file cannot be found, manual input is required.
			print '#======================================================================#'
			print 'INITMESH GEOMETRY READIN ERROR, PLEASE MANUALLY DEFINE MESH GEOMETRY FOR'
			print '#======================================================================#'
			print Dirlist[l]
			r_mesh = int(raw_input("DEFINE NUM RADIAL CELLS:"))
			z_mesh = int(raw_input("DEFINE NUM AXIAL CELLS:"))
			print ''

			R_mesh.append(r_mesh)
			Z_mesh.append(z_mesh)
		#endif
	#endtry

	#Retrieve entire mesh for plotting if requested.	#MESH PLOTTING NOT WORKING#
	if image_plotmesh == True:							#MESH PLOTTING NOT WORKING#
		print '#================================================#'
		print 'Mesh Outline Plotting Does Not Currently Function.'
		print '#================================================#'
		print ''
		#Extract mesh data from initmesh.out			#MESH PLOTTING NOT WORKING#
		mesh = open(mesh[l]).readlines()				#MESH PLOTTING NOT WORKING#
	#endif


	#==========##===== ICP.NAM READIN =====##==========#
	#==========##==========================##==========#

	#Attempt automated retrieval of SI conversion units.
	NamelistData = open(icpnam[l]).readlines()

	#Mesh Geometry Namelist Inputs
	try:
		RADIUS = float(filter(lambda x:'RADIUS=' in x, NamelistData)[0].strip(' \t\n\r,=RADIUS'))
		RADIUST = float(filter(lambda x:'RADIUST=' in x, NamelistData)[0].strip(' \t\n\r,=RADIUST'))
		HEIGHT = float(filter(lambda x:'HEIGHT=' in x, NamelistData)[0].strip(' \t\n\r,=HEIGHT'))
		HEIGHTT = float(filter(lambda x:'HEIGHTT=' in x, NamelistData)[0].strip(' \t\n\r,=HEIGHTT'))
		DEPTH = float(filter(lambda x:'DEPTH=' in x, NamelistData)[0].strip(' \t\n\r,=DEPTH'))
		SYM = float(filter(lambda x:'ISYM=' in x, NamelistData)[0].strip(' \t\n\r,=ISYM'))
		if image_plotsymmetry == True: Isymlist.append(SYM)
		else: Isymlist.append(0)
		if RADIUS > 0.0: Radius.append(RADIUS)
		elif RADIUST > 0.0: Radius.append(RADIUST)
		if HEIGHT > 0.0: Height.append(HEIGHT)
		elif HEIGHTT > 0.0: Height.append(HEIGHTT)
		Depth.append(DEPTH)
		#endif
		dr.append(Radius[-1]/(R_mesh[-1]-1))
		dz.append(Height[-1]/(Z_mesh[-1]-1))
	except:
		#If the geometry section cannot be found, manual input is required.
		print '#====================================================================#'
		print 'ICP.NAM GEOMETRY READIN ERROR, PLEASE MANUALLY DEFINE MESH SI SIZE FOR'
		print '#====================================================================#'
		print Dirlist[l]
		radius = float(raw_input("DEFINE RADIUST [cm]:"))
		height = float(raw_input("DEFINE HEIGHTT [cm]:"))
		depth = float(raw_input("DEFINE DEPTH [cm]:"))
		print ''

		Radius.append(radius)
		Height.append(height)
		Depth.append(depth)
		dr.append(Radius[-1]/(R_mesh[-1]-1))
		dz.append(Height[-1]/(Z_mesh[-1]-1))
	#endtry

	#Material Namelist Inputs (frequencies/voltages/powers)   [FREQGLOB ONLY READS 10 CHARACTERS]
	try:
		NUMMETALS = int(filter(lambda x: x.isdigit(),filter(lambda x:'IMETALS' in x,NamelistData)[0]))+1
		CMETALS = filter(lambda x: 'CMETAL=' in x, NamelistData)[0].split()[1:NUMMETALS]
		VRFM.append(filter(lambda x: 'VRFM=' in x, NamelistData)[0].split()[1:NUMMETALS])
		VRFM2.append(filter(lambda x: 'VRFM_2=' in x, NamelistData)[0].split()[1:NUMMETALS])
		FREQM.append(filter(lambda x: 'FREQM=' in x, NamelistData)[0].split()[1:NUMMETALS])
		FREQM2.append(filter(lambda x: 'FREQM_2=' in x, NamelistData)[0].split()[1:NUMMETALS])
		FREQC.append(filter(lambda x: 'FREQC=' in x, NamelistData)[0].split()[1:NUMMETALS])
		FREQGLOB.append(float(filter(lambda x:'FREQ=' in x, NamelistData)[0].strip(' \t\n\r,=FREQ')[0:10]))
		IRFPOW.append(float(filter(lambda x:'IRFPOW=' in x, NamelistData)[0].strip(' \t\n\r,=IRFPOW')))
		IETRODEM.append(filter(lambda x:'IETRODEM=' in x, NamelistData)[0].split()[1:NUMMETALS])
		for i in range(0,len(IETRODEM[l])): IETRODEM[l][i] = int(IETRODEM[l][i].strip(','))
		PRESOUT.append(  float(filter(lambda x:'PRESOUT=' in x, NamelistData)[0].strip(' \t\n\r,=PRESOUT')))
	except:
		print '#==========================================================================#'
		print 'ICP.NAM MATERIAL DEFINITIONS READIN ERROR, USING DEFAULT MATERIAL PROPERTIES'
		print '#===========================================================================#'
		FREQM.append(13.56E6)
		FREQM2.append(13.56E6)
		FREQC.append(13.56E6)
		FREQGLOB.append(13.56E6)
		VRFM.append(300.0)
		VRFM2.append(150.0)
		IRFPOW.append(100.0)
		PRESOUT.append(0.85)
	#endtry

	#Plasma Chemistry Monte-Carlo (PCMC) Namelist Inputs
	try:
		IEBINSPCMC = float(filter(lambda x: 'IEBINSPCMC=' in x, NamelistData)[0].split()[0].strip(' \t\n\r,=IEBINSPCMC'))
		EMAXIPCMC = float(filter(lambda x: 'EMAXIPCMC=' in x, NamelistData)[0].split()[0].strip(' \t\n\r,=EMAXIPCMC '))
	except:
		print '#======================================================#'
		print 'ICP.NAM PCMC READIN ERROR, USING DEFAULT PCMC PROPERTIES'
		print '#======================================================#'
		IEBINSPCMC = 250
		EMAXIPCMC = 160
	#endtry

	#Phase-Resolved IMOVIE Namelist Inputs
	try:
		IMOVIE_FRAMES.append(int(filter(lambda x:'IMOVIE_FRAMES=' in x, NamelistData)[0].strip(' \t\n\r,=IMOVIE_FRAMES')))
	except:
		print '#==================================================================#'
		print 'ICP.NAM IMOVIE READIN ERROR, USING DEFAULT PHASE RESOLVED PROPERTIES'
		print '#==================================================================#'
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
				FluidSpecies = ['AR','AR3S','O2','O']	#FLUID SPECIES ARE STILL MANUALLY DEFINED

				#Collect icp.dat header if required for later use
				header_icpdat.append([SpeciesName,Charge,MolecularWeight,StickingCoeff, TransportBool,ReturnFrac,ReturnSpecies])
			#####

			#End of Chemistry Header Denoted By '*', as soon as this is reached, stop reading in.
			elif len(ChemistryData[i].split()) != 13 and len(ChemistryData[i].split()) !=8:
				if ChemistryData[i].split()[0] == '*': break
			#endif
		#endfor
	except:
		print '#==================================================================#'
		print 'ICP.COM ATOMIC SPECIES READIN ERROR, USING DEFAULT ATOMIC PROPERTIES'
		print '#==================================================================#'
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


	#clean up variables and assign required types.
	try:
#		for i in range(0,len(CMETALS[l])): CMETALS[l][i] = CMETALS[i].strip(',\'') #!!!BROKEN!!!
		VRFM[l] = float( VRFM[l][IETRODEM[l].index(1)].strip(',') )
		VRFM2[l] = float( VRFM2[l][IETRODEM[l].index(1)].strip(',') )
		FREQM[l] = float( FREQM[l][IETRODEM[l].index(1)].strip(',') )
		FREQM2[l] = float( FREQM2[l][IETRODEM[l].index(1)].strip(',') )
		try: FREQC[l] = float( FREQMC[l][IETRODEM[l].index(1)].strip(',') )
		except: ICP_Material_Not_Found=1

		MINFREQ.append( min([FREQM[l],FREQM2[l],FREQC[l],FREQGLOB[l]]) )
		MAXFREQ.append( max([FREQM[l],FREQM2[l],FREQC[l],FREQGLOB[l]]) )
	except:
		Material_Property_Conversion_Error=1
	#endtry
#endfor

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

	#Return default if atomic physics is the same in all datasets.
	if all(map(lambda x: x == Globalnumvars[0], Globalnumvars)) == True:
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
		print 'Unable to extract '+str(NameString)
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
			numstart = (Zmesh*Rmesh)*(i)
			numend = (Zmesh*Rmesh)*(i+1)
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
		os.mkdir(NewFolderDir, 0755);
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
def AutoConvProfData(Convertexe,args=[],DirAdditions=[]):
	HomeDir = os.getcwd()
	os.chdir(Dirlist[l])

	#Remove old files to avoid any overwrite errors.
	for i in range(0,len(DirAdditions)): os.system('rm -f '+DirAdditions[i])

	#Use predefined arguments if supplied, suppresses output to devnull.
	if len(args) > 0:
		with open(os.devnull, 'w') as fp:
			subprocess = Popen(Convertexe, stdin=PIPE, stdout=fp) #noshell=True
			subprocess.communicate(os.linesep.join(args))
		#endwith

	#If no arguments supplied, run script and allow user inputs.
	elif len(args) == 0:
		os.system(Convertexe)
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

	Radial = False
	if IsStringInVariable(variable,['VR-','JR-','FR-','FLUX-R']) == True:
		Radial = True
	#endif

	return(Radial)
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

	Variablelegends = list()
	for i in range(0,len(variablelist)):
		#Explicit Species Densities.
		if variablelist[i] == 'E':
			Variable = 'Electron Density'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] in ['AR','AR3S']:
			Variable = 'Neutral Ar Density'
			VariableUnit = '[m$^{-3}$]'
		elif variablelist[i] == 'AR+':
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

		#Explicit Species Temperatures.
		elif variablelist[i] == 'TE':
			Variable = 'Electron Temperature'
			VariableUnit = '[eV]'
		elif variablelist[i] == 'TG-AVE':
			Variable = 'Neutral Gas Temperature'
			VariableUnit = '[K]'

		#Explicit Species Velocities and Resulting Pressure.
		elif variablelist[i] == 'PRESSURE':
			Variable = 'Pressure'
			VariableUnit = '['+str(PressureUnit)+']'			#Default: '[Torr]'
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
		elif variablelist[i] == 'P-POT':
			Variable = 'Plasma Potential'
			VariableUnit = '[V]'
		elif variablelist[i] == 'PPOT':
			Variable = 'Plasma Potential'
			VariableUnit = '[V]'
		elif variablelist[i] == 'RHO':
			Variable = 'Charge Density'
			VariableUnit = '[C cm$^{-3}$]'
		elif variablelist[i] == 'EF-TOT':
			Variable = 'E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in ['ER','EAMB-R']:
			Variable = 'Radial E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] in ['EZ','EAMB-Z']:
			Variable = 'Axial E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] == 'BT':
			Variable = 'B-field Strength'
			VariableUnit = '['+str(BFieldUnit)+']'			#Default: '[G]'
		elif variablelist[i] == 'BR':
			Variable = 'Radial B-field Strength'
			VariableUnit = '['+str(BFieldUnit)+']'			#Default: '[G]'
		elif variablelist[i] == 'BRS':
			Variable = 'Radial B-field Strength'
			VariableUnit = '['+str(BFieldUnit)+']'			#Default: '[G]'
		elif variablelist[i] == 'BZ':
			Variable = 'Axial B-field Strength'
			VariableUnit = '['+str(BFieldUnit)+']'			#Default: '[G]'
		elif variablelist[i] == 'BZS':
			Variable = 'Axial B-field Strength'
			VariableUnit = '['+str(BFieldUnit)+']'			#Default: '[G]'
		elif variablelist[i] == 'JZ-NET':
			Variable = 'Axial Current Density'
			VariableUnit = '[mA cm$^{-2}$]'
		elif variablelist[i] == 'JR-NET':
			Variable = 'Radial Current Density'
			VariableUnit = '[mA cm$^{-2}$]'

		#Explicit Power Deposition.
		elif variablelist[i] == 'POW-TOT':
			Variable = 'Capacitively Coupled RF-Power'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] == 'POW-ICP':
			Variable = 'Inductively Coupled RF-Power'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] == 'POW-RF':
			Variable = 'RF-Power Deposited'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] == 'POW-RF-E':
			Variable = 'RF-Power Deposited by e$^-$'
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
	if IsStringInVariable(variable,['PRESSURE']) == True and PressureUnit == 'Pa':
		for i in range(0,len(profile)):
			profile[i] = profile[i]*133.333333333
		#endfor
	if IsStringInVariable(variable,['PRESSURE']) == True and PressureUnit == 'mTorr':
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1000.0
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

	#For B-field strengths, convert to Tesla or retain as default Gauss. (also reverse axial field)
	if IsStringInVariable(variable,['BT','BRS']) == True and BFieldUnit == 'Tesla':
		for i in range(0,len(profile)):
			profile[i] = profile[i]/10000
		#endfor
	if IsStringInVariable(variable,['BZS']) == True and BFieldUnit == 'Tesla':
		for i in range(0,len(profile)):
			profile[i] = (profile[i]/10000)*(-1)
		#endfor
	if IsStringInVariable(variable,['BZS']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For E-field strengths, convert from [V cm-1] to [V m-1]. (also reverse axial field)
	if IsStringInVariable(variable,['EF-TOT','EAMB-R','EAMB-Z']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]#*100	### [V cm-1] ###
		#endfor
	if IsStringInVariable(variable,['EAMB-Z']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(-1)
		#endfor
	#endif

	#For Current Densities, convert from [A cm-2] to [mA cm-2]. (also reverse axial current)
	if IsStringInVariable(variable,['JZ-NET','JR-NET']) == True:
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

	#Co-axial magnetic ring
	Ax.plot((39*dz[l],49*dz[l]),  ( 48*dr[l], 48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],39*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  ( 48*dr[l], 62*dr[l]), 'g-', linewidth=2)
	Ax.plot((49*dz[l],49*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-48*dr[l]), 'g-', linewidth=2)
	Ax.plot((39*dz[l],49*dz[l]),  (-48*dr[l],-62*dr[l]), 'g-', linewidth=2)

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

#===================##===================#
#===================##===================#

































#====================================================================#
				 	 #READING DATA INTO MEMORY#
#====================================================================#

print'-----------------------'
print'Beginning Data Read-in.'
print'-----------------------'

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

	#Kinetics data readin - NOT CURRENTLY EMPLOYED IN ANY DIAGNOSTICS
	if True == True:

		#Load data from TECPLOT_KIN file and unpack into 1D array.
		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
		rawdata_kin.append(rawdata)

		#Read through all variables for each file and stop when list ends.
		KinVariablelist,KinHeaderEndMarker = ['T (S)'],'ZONE'
		for i in range(2,nn_2D):
			if KinHeaderEndMarker in str(rawdata_kin[l][i]): 
				I = int(filter(lambda x: x.isdigit(), rawdata_kin[l][i].split(',')[0]))
				break
			else: KinVariablelist.append(str(rawdata_kin[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		numvariables_kin,header_kin = len(KinVariablelist),len(KinVariablelist)+2
		header_kinlist.append(header_kin)

		#Seperate total 1D data array into sets of data for each variable.
		CurrentFolderData = SDFileFormatConvertorHPEM(rawdata_kin[l],header_kin,numvariables_kin, Zmesh=I,Dimension='1D')

		#Save all variables for folder[l] to Data.
		#Data is now 3D array of form [folder,variable,datapoint(R,Z)]
		DataKin.append(CurrentFolderData)
	#endif


#===================##===================#
#===================##===================#

	#IEDF/NEDF file readin.
	if True in [savefig_IEDFangular,savefig_IEDFtrends]:

		#Define arguments and autorun conv_prof.exe if possible.
		#### THIS IS HACKY, WON'T ALWAYS WORK, ARGS LIST NEEDS AUTOMATING ####
		IEDFVarArgs = ['1','1','1','1','1'] 	#Works for 2 species 1 surface.
		ExtraArgs = ['1','1','1','1','1','1','1','1','1','1']#[]	#Hack For Additional Species
		Args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ExtraArgs + ['0','0']
		DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt','iprofile_zones_tec1d.pdt','nprofile_zones_tec1d.pdt']
		try: AutoConvProfData('./conv_prof.exe',Args,DirAdditions)
		except: print 'ConvProf Failure:'+Dirlist[l]

		#Load data from IEDFprofile file and unpack into 1D array.
		rawdata, nn_IEDF = ExtractRawData(Dir,'iprofile_tec2d.pdt',l)
		rawdata_IEDF.append(rawdata)

		#Read through all variables for each file and stop when list ends.
		IEDFVariablelist,HeaderEndMarker = ['Theta [deg]','Energy [eV]'],'ZONE'
		for i in range(2,nn_IEDF):
			#Grab EDFangle(I),EDFbins(J) values from the ZONE line, these outline the datablock size.
			if HeaderEndMarker in str(rawdata_IEDF[l][i]): 
				I = int(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[0]))
				J = int(filter(lambda x: x.isdigit(), rawdata_IEDF[l][i].split(',')[1]))
				EDFangle, EDFbins = I,J
				break
			else: IEDFVariablelist.append(str(rawdata_IEDF[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		numvariables_IEDF,header_IEDF = len(IEDFVariablelist),len(IEDFVariablelist)+2
		header_IEDFlist.append(header_IEDF)

		#Seperate total 1D data array into sets of data for each variable.
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

	if True in [savefig_convergence,savefig_pulseprofiles]:

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

		#Cycle through all iterations for current datafile, appending per cycle.
		CurrentFolderData,CurrentFolderIterlist = list(),list()
		for i in range(0,len(Iterloc)):
			if i == 0:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar+2,offset=2)
				CurrentFolderData.append(CurrentIterData[0:numvar])
			else:
				CurrentIterData = SDFileFormatConvertorHPEM(rawdata,Iterloc[i],numvar)
				CurrentFolderData.append(CurrentIterData)
			#endif
		#endfor
		IterMovieData.append(CurrentFolderData)
	#endif


#===================##===================#
#===================##===================#

	Batch=False
	if True in [savefig_phaseresolve2D,savefig_phaseresolve1D,savefig_PROES] and Batch==True:

		#Load data from movie_icp file and unpack into 1D array.
		rawdata,nn_phasemovie = ExtractRawData(Dir,'movie1.pdt',l)
		rawdata_phasemovie.append(rawdata)

		#Read through all variables for each file and stop when list ends. 
		#Movie1 has geometry at top, therefore len(header) != len(variables).
		#Only the first encountered geometry is used to define variable zone.
		VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ZONE'
		variablelist,numvar = list(),0
		for i in range(2,nn_phasemovie):
			if HeaderEndMarker in str(rawdata_phasemovie[l][i]): 
				header_phase = i+2		# +2 to skip to data row.
				break
			if VariableEndMarker in str(rawdata_phasemovie[l][i]) and numvar == 0:
				numvar = (i-3)	# -3 to not include R,Z and remove overflow.
			if len(rawdata_phasemovie[l][i]) > 1 and numvar == 0: 
				variablelist.append(str(rawdata_phasemovie[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		header_phasemovie.append(header_phase)

		#Rough method of obtaining the movie1.pdt cycle locations for data extraction.
		cycleloc = list()
		for i in range(0,len(rawdata_phasemovie[l])):
			if "CYCL=" in rawdata_phasemovie[l][i]:
				cycleloc.append(i+1)
			#endif
		#endfor

		#Cycle through all phases for current datafile, appending per cycle.
		CurrentFolderData,CurrentFolderPhaselist = list(),list()
		for i in range(0,len(cycleloc)):
			if i == 0:
				CurrentPhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar+2,offset=2)
				CurrentFolderData.append(CurrentPhaseData[0:numvar])
			else:
				CurrentPhaseData = SDFileFormatConvertorHPEM(rawdata,cycleloc[i],numvar)
				CurrentFolderData.append(CurrentPhaseData)
			#endif
			CurrentFolderPhaselist.append('CYCL = '+str(i+1))
		#endfor
		Moviephaselist.append(CurrentFolderPhaselist)
		PhaseMovieData.append(CurrentFolderData)
	#endif
#endfor


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
del RADIUS,RADIUST,HEIGHT,HEIGHTT,DEPTH,SYM
del data_array,tempdata,tempdata2,templineout
del Variablelist,variablelist
del Energy,Fe,rawdata_mcs
del HomeDir,DirContents


#Alert user that readin process has ended and continue with selected diagnostics.
if any([savefig_plot2D, savefig_phaseresolve2D, savefig_convergence, savefig_monoprofiles, savefig_multiprofiles, savefig_comparelineouts, savefig_pulseprofiles, savefig_trendphaseresolved, savefig_phaseresolve1D, savefig_PROES, savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, savefig_IEDFangular, savefig_IEDFtrends, savefig_EEDF]) == True:
	print '----------------------------------------'
	print 'Data Readin Complete, Starting Analysis:'
	print '----------------------------------------'
else:
	print '------------------'
	print 'Analysis Complete.'
	print '------------------'
#endif


#=====================================================================#
#=====================================================================#



























#====================================================================#
				  #COMMONLY USED PLOTTING FUNCTIONS#
#====================================================================#

#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html
def Matplotlib_GlobalOptions():

#	mpl.style.use('classic')								#Resets to classic 1.x.x format
	
	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
	mpl.rcParams['figure.dpi'] = 100						#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
	mpl.rcParams['image.interpolation'] = 'bilinear'		#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
	mpl.rcParams['image.cmap'] = 'plasma'					#Select global colourmap 
	#'jet','plasma','gnuplot'

	#Axis options
	mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'	#View limits coencide with axis ticks
	mpl.rcParams['axes.xmargin'] = 0						#Set default x-axis padding
	mpl.rcParams['axes.ymargin'] = 0						#Set default y-axis padding
	mpl.rcParams['errorbar.capsize'] = 3					#Set error bar end cap width
	mpl.rcParams['font.size'] = 12							#Set global fontsize
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



#Returns a 2D array of inputted data with size [R_mesh] x [Z_mesh]
#Can optionally perform variable unit conversion if required.
#Image = ImageExtractor2D(Data,Variable=[]):
def ImageExtractor2D(Data,Variable=[],Rmesh=0,Zmesh=0):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Zmesh == 0:
		Rmesh,Zmesh = R_mesh[l],Z_mesh[l]
	#endif

	#Create empty 2D image of required size.
	numrows = len(Data)/Rmesh
	Image = np.zeros([numrows,Rmesh])

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

	return(Image)
#enddef



#=========================#
#=========================#



#Takes 2D image array and produces symmetric image about central axis.
#Returns symmetric 2D image array, allows radially negative values.
#Image = SymmetryConverter(Image,Radial=False)
def SymmetryConverter(Image,Radial=False):

	#Create new image by reversing and adding itself on the LHS.
	if image_plotsymmetry == True and Isymlist[l] == 1:
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



#Scales a 2D array by the requested X and Y scale factors, can scale asymmetrically.
#Takes a rectilinear 2D array of floats and returns a scaled rectilinear 2D array.
#Employs a linear interpolation scheme when mapping new datapoints.
#ScaledArray = ScaleArray(2DArray,ScaleFactors=[3,3])
def ScaleArray(Array,ScaleFactors):
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



#Create figure of desired size and with variable axes.
#Returns figure and axes seperately.
#fig,ax = figure(image_aspectratio,1,shareX=False)
def figure(aspectratio=[],subplots=1,shareX=False):
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
	else:
		fig, ax = plt.subplots(subplots, figsize=(10,10), sharex=shareX)
	#endif
	return(fig,ax)
#enddef



#=========================#
#=========================#



#Crops 2D image taking account of image rotation options.
#Takes image axis (assumes default axis), use figure()
#Input Extent format: [ [Rmin,Rmax], [Zmin,Zmax] ] in cm.
#Returns cropping limits in format: [ [R1,R2],[Z1,Z2] ] in cm.
#CropImage(ax[0],Extent=[[R1,R2],[Z1,Z2]],Apply=True,Rotate=True), 
def CropImage(ax=plt.gca(),Extent=[],Apply=True,Rotate=True):

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



#Provides a new colourbar scale for cropped images.
#Takes a 2D image, and returns the min/max value within the cropped region.
#Assumes image symmetry for best results, otherwise R1 is set to zero.
#Works for PROES images too, requires PROES='Axial' or 'Radial'.
#[Minimum,Maximum] = CbarMinMax(Image,PROES=False)
def CbarMinMax(Image,PROES=False,Symmetry=image_plotsymmetry):

	#Return user defined limits if specified.
	if len(image_cbarlimit) == 2:
		cropmin = image_cbarlimit[0]
		cropmax = image_cbarlimit[1]
		return([cropmin,cropmax])
	#endif

	#Ensure limits are in line with any requested mathematical constraints
#	if image_logplot == True: Image = np.log(Image)
#	if image_normalize == True: Image = Normalize(Image)

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
			R1 = R_mesh[l] + R1
			R2 = R_mesh[l] + R2
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
			Image = Image[::-1][Z1:Z2]			#Account for reversed PhaseData
		elif PROES == 'Radial':
			Image = Image[R1:R2]	
		#endif
	#endif

	#Flatten image and obtain min/max in region, defaults to full image if no cropping.
	try:	
		flatimage = [item for sublist in Image for item in sublist]
		cropmin,cropmax = min(flatimage),max(flatimage)
	except:	
		print 'IMAGE CROPPING OUTSIDE MESH BOUNDARIES: CHECK IMAGE_RADIALCROP,IMAGE_AXIALCROP'
		exit()
	#endtry

	#Remove any nan's or infs from image limits (potentially caused by normalisation or logging)
	if np.isinf(cropmin) == True and cropmax < 0.0: cropmin = cropmax		#Cropmin > Cropmax
	elif np.isinf(cropmin) == True: cropmin = 0.0
	if np.isinf(cropmax) == True: cropmax = 0.0	

	#Return cropped values in list [min,max], as required by colourbar.
	return([cropmin,cropmax])
#enddef



#=========================#
#=========================#



#Applies plt.options to current figure based on user input.
#Returns nothing, open figure is required, use figure().
#For best results call immediately before saving/displaying figure.
#ImageOptions(plt.gca(),Xlabel,Ylabel,Title,Legend,Crop=False)
def ImageOptions(fig,ax=plt.gca(),Xlabel='',Ylabel='',Title='',Legend=[],Crop=True,Rotate=True):

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
	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
	except: Axes_Contain_Strings = True
#	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')	#Old tickformat.
#	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')	#Old tickformat.
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
	#endif

	#Plot mesh outline if requested.	### HACKY ###
	if image_plotmesh == True:
		mesh_auto_plot = 1 #AUTO PLOT MESH#
		#NOT IMPLIMENTED!! REQUIRES initmesh.out READER#
	elif image_plotmesh == 'PRCCP' and Crop == True:	
		ManualPRCCPMesh(ax)
	elif image_plotmesh == 'EVgeny' and Crop == True:
		ManualEVgenyMesh(ax)
	elif image_plotmesh == 'HyperionI' and Crop == True:
		ManualHyperionIMesh(ax)
	elif image_plotmesh == 'HyperionII' and Crop == True:
		ManualHyperionIIMesh(ax)
	#endif

	#Crop image dimensions:	If Crop is a list, use cropping dimensions provided to this function
	if isinstance(Crop, (list, np.ndarray) ) == True:
		CropImage(ax,Extent=Crop,Rotate=Rotate)
	#Else if crop is a boolian, use default cropping dimensions from switchboard
	elif Crop == True:
		if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
			CropImage(ax,Rotate=Rotate)
		#endif
	#endif

	#Arrange figure such that labels, legends and titles fit within frame.
	fig.tight_layout()

	return()
#enddef



#=========================#
#=========================#



#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#cbar = Colourbar(ax[0],'Label',5,Lim=[0,1])
def Colourbar(ax=plt.gca(),Label='',Ticks=5,Lim=[]):

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



#Creates an invisible colourbar to align subplots without colourbars.
#Takes image axis, returns colourbar axis if further edits are required
#cax = InvisibleColourbar(ax[0])
def InvisibleColourbar(ax='NaN'):
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



#Generates a 1D SI [cm] axis for plotting, includes radial symmetry.
#Takes orientation, symmetry and phasecycle options.
#Returns 1D array in units of [cm] or [omega*t/2pi].
#Raxis=GenerateAxis('Radial',Isym=Isymlist[l])
def GenerateAxis(Orientation,Isym=Isymlist[l],PhaseFrames=range(0,IMOVIE_FRAMES[l])):

	#Create axis list and extract the number of phaseframes
	PhaseResolution = len(PhaseFrames)
	axis = list()

	if Orientation == 'Radial':
		if Isym == 1:
			for i in range(-R_mesh[l],R_mesh[l]):
				axis.append(i*dr[l])
		#endfor
		elif Isym != 1:
			for i in range(0,R_mesh[l]):
				axis.append(i*dr[l])
			#endfor
		#endif
	elif Orientation == 'Axial':
		for i in range(0,Z_mesh[l]):
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



#Takes 1D or 2D array and returns array normalized to maximum value.
#If NormFactor is defined, array will be normalized to this instead.
#Returns normalized image/profile and the max/min normalization factors.
#NormProfile,Min,Max = Normalize(profile,NormFactor=0)
def Normalize(profile,NormFactor=0):
	NormalizedImage = list()

	#determine dimensionality of profile and select normaliztion method.
	if isinstance(profile[0], (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 2D array.
		FlatImage = [item for sublist in profile for item in sublist]
		MaxNormFactor,MinNormFactor = max(FlatImage),min(FlatImage)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0
		#endif

		#Normalize 2D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			NormalizedImage.append( [x/NormFactor for x in profile[i]] )
		#endfor
		profile = NormalizedImage
		return(profile,MaxNormFactor,MinNormFactor)

	#Lowest dimention is still list.
	elif isinstance(profile, (list, np.ndarray) ) == True:

		#Obtain max and min normalization factors for 1D profile.
		MaxNormFactor,MinNormFactor = max(profile),min(profile)

		#Fix for division by zero and other infinity related things...
		if 'inf' in str(MaxNormFactor) or MaxNormFactor == 0.0: MaxNormFactor = 1.0
		if 'inf' in str(MinNormFactor) or MinNormFactor == 0.0: MinNormFactor = 0.0

		#Normalize 1D array to local maximum.
		if NormFactor == 0: NormFactor = MaxNormFactor
		for i in range(0,len(profile)):
			profile[i] = profile[i]/NormFactor
		#endfor
	#endif

	return(profile,MinNormFactor,MaxNormFactor)
#enddef



#=========================#
#=========================#



#Takes current image datails and returns extent and rotated aspectratio
#If mesh uses symmetry, will double radius extent centered on zero.
#extent,aspectratio = DataExtent(l)
def DataExtent(folder=l,aspectratio=image_aspectratio,rotate=image_rotate):

	#Obtain global variables for current folder.
	Isym = Isymlist[folder]
	radius,height = Radius[folder],Height[folder]

	#Rotated Image: [X,Y] = [Height,Radius]
	if image_rotate == True:
		aspectratio = aspectratio[::-1]
		if Isym == 1: extent=[0,height, -radius,radius]
		elif Isym == 0: extent=[0,height, 0,radius]
		#endif

	#Default mesh orientation: [X,Y] = [Radius,Height]
	elif rotate == False:
		if Isym == 1: extent = [-radius,radius, 0,height]
		elif Isym == 0: extent=[0,radius, 0,height]
		#endif
	#endif

	return(extent,aspectratio)
#enddef



#=========================#
#=========================#



#Create figure and plot a 1D graph with associated image plotting requirements.
#Returns plotted axes and figure if new ones were created.
#Else plots to existing figure and returns the image object.
#ImagePlotter1D(Zlineout,Zaxis,image_aspectratio,fig,ax[0]):
def ImagePlotter1D(profile,axis,aspectratio,fig=111,ax=111):

	#Generate new figure if required. {kinda hacky...}
	if fig == 111 and ax == 111:
		fig,ax = figure(aspectratio)
	elif fig == 111:
		fig = figure(aspectratio)
	#endif

	#Apply any required numerical changes to the profile.
	if image_logplot == True:
		profile = np.log(profile)
	if image_normalize == True:
		profile = Normalize(profile)[0]
	#endif

	#Plot profile and return.
	im = ax.plot(axis,profile, lw=2)

	try: return(fig,ax,im)
	except: return()
#enddef



#=========================#
#=========================#



#Create figure and plot a 2D image with associated image plotting requirements.
#Returns plotted image, axes and figure after applying basic data restructuring.
#fig,ax,im,Image = ImagePlotter2D(Image,extent,image_aspectratio,variablelist[l],fig,ax[0])
def ImagePlotter2D(Image,extent,aspectratio=image_aspectratio,variable='N/A',fig=111,ax=111,Origin="lower"):

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
	elif image_normalize == True:
		Image = Normalize(Image)[0]
	#endif

	#Plot image with or without contour plots, (contour scale = 90% of cbar scale)
	if image_contourplot == True:
		im = ax.contour(Image,extent=extent,origin=Origin)
		im.set_clim(CbarMinMax(Image)[0]*0.90,CbarMinMax(Image)[1]*0.90)
		im = ax.imshow(Image,extent=extent,origin=Origin)
	else:
		im = ax.imshow(Image,extent=extent,origin=Origin)
	#endif
	return(fig,ax,im,Image)
#enddef



#=========================#
#=========================#



#Creates a 1D image from an array of supplied points.
#Image plotted onto existing axes, figure() should be used.
#NormFactor = 0 will normalize to maximum of given profile.
#TrendPlotter(ax[0],TrendProfiles,Xaxis,'o-',0)
def TrendPlotter(ax=plt.gca(),TrendArray=[],Xaxis=[],Marker='o-',NormFactor=0):

	#Normalize data to provided normalization factor if required.
	if image_normalize == True:
		TrendArray = Normalize(TrendArray,NormFactor)[0]
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



#Obtains a radial 1D profile at a requested axial location.
#Returns a 1D array for plotting and performs unit conversion.
def PlotRadialProfile(Data,process,variable,lineout,Rmesh=0,Isym=0):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Isym == 0:
		Rmesh,Isym = R_mesh[l],Isymlist[l]
	#endif

	#Obtain start location for requested data and perform SI conversion.
	ZStart = Rmesh*lineout
	ZEnd = Rmesh*(lineout+1)

	#Plot lines for each variable at each requested slice, ignoring ones that fail.
	#If mesh is symmetric, copy the data over and make full plot.
	if Isym == 1:
		Zend = len(Data[process])-ZStart
		Zstart = len(Data[process])-ZEnd
		Rlineout = Data[process][Zstart:Zend][::-1]
		#If variable is radially symmetric, add negative to symmetry
		if IsRadialVariable(variable) == True:
			for m in range(0,len(Rlineout)):
				Rlineout.append(-Data[process][Zstart:Zend][m])
			#endfor
		#If variable is axially symmetric, add positive.
		elif IsRadialVariable(variable) == False:
			for m in range(0,len(Rlineout)):
				Rlineout.append(Data[process][Zstart:Zend][m])
			#endfor
		#endif
		Rlineout = Rlineout[::-1]	#Reverse index, negative first then positive values.		

	#If the data isn't symmetric, just plot as is.
	elif Isym == 0:
		Zend = len(Data[process])-ZStart
		Zstart = len(Data[process])-ZEnd
		Rlineout = Data[process][Zstart:Zend]
	#endif

	#Convert units if required and plot.
	Rlineout = VariableUnitConversion(Rlineout,variable)
	return(Rlineout)
#enddef



#=========================#
#=========================#



#Obtains an axial 1D profile at a requested radial location.
#Returns a 1D array for plotting and performs unit conversion.
def PlotAxialProfile(Data,process,variable,lineout,Rmesh=0,Zmesh=0,Isym=0):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Zmesh == 0 or Isym == 0:
		Rmesh,Zmesh,ISym = R_mesh[l],Z_mesh[l],Isymlist[l]
	#endif

	#Pull out Z-data point from each radial line of data and list them.
	Zlineout = list()
	for i in range(0,Zmesh):
		datapoint = Rmesh*i + lineout
		try:
			Zlineout.append(Data[process][datapoint])
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



#Convert cell location for use with WaveformExtractor function.
#Returns mesh location based on input string or [R,Z] list.
#Rcell,Zcell = Waveformloc(electrodeloc,'Phase')
def WaveformLoc(location,origintype):

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



#Takes phasedata for current folder and PPOT process number.
#Returns two arrays: VoltageWaveform at electrode location and average waveform.
#VoltageWaveform,WaveformBias = WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[0])
def WaveformExtractor(PhaseData,PPOT,waveformlocation=electrodeloc):

	#Create required lists and extract electrode location.
	VoltageWaveform,WaveformBias = list(),list()
	RLoc = WaveformLoc(waveformlocation,'Phase')[0]
	ZLoc = WaveformLoc(waveformlocation,'Phase')[1]

	#Create voltage waveform for requested integer number of phasecycles
	for i in range(0,int(phasecycles*len(PhaseData))):
		Index = i % len(PhaseData)		#Use modulo index for additional phase cycles
		VoltageWaveform.append(PlotAxialProfile(PhaseData[Index],PPOT,'PPOT',ZLoc)[RLoc])
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



#Trend analysis for a given point on a 2D Image.
#Takes global 'TrendLocation' for safety, process and variable.
#Returns two arrays: One is the X-axis to plot against
#Second is the value of variable at location for all simulations.
def TrendAtGivenLocation(TrendLocation,process,variable):

	#Refresh lists that change per image.
	R,Z = TrendLocation[0],TrendLocation[1]
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
			print FolderNameTrimmer(Dirlist[l])
			print str(variable)+' @ '+Location+':', round(Trend[-1], 5)
		if print_generaltrends == True and l == numfolders-1:
			print ''
		#endif
	#endfor

	#Normalize to maximum value in each profile if required.
	if image_normalize == True:
		Trend,Min,Max = Normalize(Image)
	#endif

	return(Xaxis,Trend)
#enddef



#=========================#
#=========================#



#General trend plotting function for use with multiple folders.
#Takes a lineout location and orientation string as input.
#And returns the maximum and minimum values For all folders to be analysed.
#With an associated X-axis composed of the trimmed folder names.
def MinMaxTrends(lineout,Orientation,process):

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
			try: Profile = PlotRadialProfile(Data[l],processlist[p],Variablelist[p],lineout,R_mesh[l],Isymlist[l])
			except: Profile = float('NaN')
			#endtry
		elif Orientation == 'Axial':
			try: Profile = PlotAxialProfile(Data[l],processlist[p],Variablelist[p],lineout,R_mesh[l],Z_mesh[l],Isymlist[l])
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
			print FolderNameTrimmer(Dirlist[l])
			print VariableName+' '+Orientation+'Maximum: ', round(max(MaxValueTrend), 5)
			print VariableName+' '+Orientation+'Minimum: ', round(min(MinValueTrend), 5)
		if print_generaltrends == True and l == numfolders-1:
			print ''
		#endif
	#endfor

	#Normalize to maximum value in each profile if required.
	if image_normalize == True:
		MaxValueTrend,MaxMin,MaxMax = Normalize(MaxValueTrend)
		MinValueTrend,MinMin,MinMax = Normalize(MinValueTrend)
	#endif

	return(Xaxis,MaxValueTrend,MinValueTrend)
#enddef



#=========================#
#=========================#


#TREND ANALYSIS - Speed of Sound
#Calculates local sound speed via Newton-Laplace equation
#Takes appropriate neutral density [m-3] and pressure [Torr] (0D,1D or 2D)
#Returns same dimensionality array of sound speeds in m/s
#SoundSpeed = LocalSoundSpeed(ArgonDensity,Pressure,Dimension='2D')
def LocalSoundSpeed(NeutralDensity,Pressure,Dimension='2D'):
	#Initiate required lists and set atomic values
	SoundSpeedArray = list()
	AdiabaticIndex = 5.0/3.0		#For Argon/Helium
	AtomicMass = 39.948*1.66E-27	#Kg

	#For 0D values:
	if Dimension == '0D':
		ElasticityModulus = AdiabaticIndex*Pressure*133.33
		MassDensity = NeutralDensity*AtomicMass
		try: SoundSpeedArray = np.sqrt( ElasticityModulus/MassDensity )
		except: SoundSpeedArray = np.nan
	#endif

	#For 1D arrays:
	if Dimension == '1D':
		for i in range(0,len(NeutralDensity)):
			ElasticityModulus = AdiabaticIndex*Pressure[i]*133.33
			MassDensity = NeutralDensity[i]*AtomicMass

			#Calculate local sound speed via Newton-Laplace equation
			try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )
			except: SoundSpeed = np.nan
			SoundSpeedArray.append( SoundSpeed )
		#endfor
	#endif

	#For 2D arrays:
	if Dimension == '2D':
		for i in range(0,len(NeutralDensity)):
			SoundSpeedArray.append(list())
			for j in range(0,len(NeutralDensity[i])):
				ElasticityModulus = AdiabaticIndex*Pressure[i][j]*133.33
				MassDensity = NeutralDensity[i][j]*AtomicMass

				#Calculate local sound speed via Newton-Laplace equation
				try: SoundSpeed = np.sqrt( ElasticityModulus/MassDensity )
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
#Takes a PPOT profile and calcuates DCbias via difference in voltage drop.
#Can identify DC-bias for parallel plate discharges and dielectric discharges.
def DCbiasMagnitude(PPOTlineout):

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

	if DebugMode == True:
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



#Calculates Brinkmann sheath width assuming Child-Langmuir conditions.
#Calculation Methods: 'AbsDensity', 'IntDensity'
#Takes current folder, current axis, movie1 Phase and sheath calc method.
#Returns array of sheath distances from origin and can plot this if requested.
#Sx = SheathExtent(folder=l,Phase=moviephaselist[k])
def SheathExtent(folder=l,ax='NaN',Orientation='Axial',Phase='NaN',Ne=list(),Ni=list()):
	#Return null array if sheath plotting is not required:
	if image_sheath == False: return([np.nan],[np.nan])

	#Initiate required data storage lists
	NPos,NNeg = list(),list()
	Sx,SymSx = list(),list()	

	#Import global sheath calculation method and charged particle species names
	SheathMethod=GlobSheathMethod
	if len(SheathIonSpecies) == 0:
		global PosSpecies
		global NegSpecies
	#Force single sheath species - Legacy Code or for testing purposes
	elif len(SheathIonSpecies) > 0:
		PosSpecies = SheathIonSpecies
		NegSpecies = []
	#endif

	#Identify charged species and alter names to suit TECPLOT2D nomenclature
	for i in range(0,len(PosSpecies)): PosSpecies[i] = PosSpecies[i] = PosSpecies[i].replace('^','+')
	for i in range(0,len(NegSpecies)): NegSpecies[i] = NegSpecies[i] = NegSpecies[i].replace('^','-')
	if 'E' in NegSpecies: NegSpecies.remove('E')			#Might Cause An Issue With Global!!!


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
						anti_j = RadialPlasmaExtent-j-1
						Neff_sum += Neff[i][j]		#Sum from R=wall to R=0	[anti_j] 	####FUDGED####
						Ne_sum += Ne[i][j]			#Sum from R=0 to R=wall [j] 	 	####FUDGED####

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

		#Create symmetric sheath boundary
		for i in range(0,len(Sx)): 
			try: SymSx.append(-Sx[i])
			except: SymSx.append(np.nan)
		#endfor
	#endif

	#=======#	#=======#	#=======#
	#=======#	#=======#	#=======#

	#NEED TO APPLY RADIAL METHOD!!! FOR NOW THIS IS SET TO ZERO EVERYWHERE#
	if Orientation == 'Radial':
		#Sheath extension: integral_(R0->Rwall) ne dR == integral_(Rwall->R0) ni dR
		for i in range(0,len(Neff)):
			Sx.append(np.nan)
		#endfor

		#Create symmetric sheath boundary
		for i in range(0,len(Sx)): 
			try: SymSx.append(-Sx[i])
			except: SymSx.append(np.nan)
		#endfor
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
		print 'Simulation:', Dirlist[folder]
		print 'Sheath Location:',SheathWidth*10, 'mm'
		print 'Sheath Extent:',((SourceWidth[0]*dr[l])-SheathWidth)*10, 'mm'
	#endif

	#THIS SHOULD PROBABLY BE A SEPERATE FUNCTION
	if ax != 'NaN':
		#Generate axis and plot sheath characteristics if requested.
		Axis=GenerateAxis(Orientation,Isym=Isymlist[folder])
		if image_sheath == True:
			ax.plot(Axis,Sx, 'w--', lw=2)
			ax.plot(Axis,SymSx, 'w--', lw=2)
		#endif
	#endif

	#Return sheath expansion
	return(Sx,SymSx)
#enddef


#=========================#
#=========================#
































#====================================================================#
				  #IMAGE PLOTTING AND DATA ANALYSIS#
#====================================================================#

#====================================================================#
				#STEADY STATE MOVIES -- TIME AVERAGED#
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

			#Generate and rotate figure as requested.
			extent,aspectratio = DataExtent(l)
			fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,variablelist[k])
			#Add sheath thickness to figure if requested.
			if image_sheath == True: Sx = SheathExtent(folder=l,ax=ax)[0]

			#Overlay location of 1D profiles if requested, adjusting for image rotation.
			if image_1Doverlay == True:
				for j in range(0,len(radialineouts)):
					X1,X2 = extent[0],extent[1]
					Y1,Y2 = radialineouts[j]*dz[l],radialineouts[j]*dz[l]
					if image_rotate == True: X1,X2,Y1,Y2 = Y1,Y2,X1,X2
					ax.plot((X1,X2),(Y1,Y2),'k--',lw=2)
				#endfor
				for j in range(0,len(heightlineouts)):
					X1,X2 = heightlineouts[j]*dr[l],heightlineouts[j]*dr[l]
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
			Title = '2D Steady State Plot of '+variablelist[k]+' for \n'+Dirlist[l][2:-1]
			#Add Colourbar (Axis, Label, Bins)
			label,bins = VariableLabelMaker(variablelist),5
			cax = Colourbar(ax,label[k],bins,Lim=CbarMinMax(Image))
			#Finalize image
			ImageOptions(fig,ax,Xlabel,Ylabel,Title)

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(Dir2Dplots, '2Dplots_Data')
				WriteDataToFile(Image, DirWrite+variablelist[k])
				if image_sheath == True and k == len(processlist)-1: WriteDataToFile(Sx, DirWrite+'Sx-EXT')
			#endif

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot '+variablelist[k]+ext)
			plt.close('all')
		#endfor
	#endfor

	print'-------------------------------------'
	print'# 2D Steady-State Processing Complete'
	print'-------------------------------------'
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
		for i in range(0,len(processlist)):
			processlist[i] = processlist[i]-2
		#endfor

		#Create list and x-axis for convergence trend plotting.
		ConvergenceTrends,Xaxis = list(),list()
		for i in range(0,len(MovieIterlist[l])):
			Iter = filter(lambda x: x.isdigit(), MovieIterlist[l][i])
			Xaxis.append(float(Iter))
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
				Image = np.zeros([numrows,R_mesh[l]])
			except: 
				print 'No Iteration Data Found For '+Dirlist[l]
				break
			#endtry
			
			#Reshape specific part of 1D Data array into 2D image for plotting.
			if QuickConverge == False:
				for k in range(0,len(MovieIterlist[l])):

					#Extract full 2D image for further processing.
					Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
					#Take Max value of image for general convergence trend.
					ConvergenceTrends[-1].append( sum(Image.flatten())/len(Image.flatten()) )

					
					#Generate and rotate figure as requested.
					extent,aspectratio = DataExtent(l)
					fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,variablelist[i])
					#Add sheath thickness to figure if requested. (CURRENTLY USES TECPLOT2D DATA)
					Sx = SheathExtent(folder=l,ax=ax)[0]

					#Define image axis labels.
					if image_rotate == True:
						Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
					elif image_rotate == False:
						Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
						plt.gca().invert_yaxis()
					#endif

					#Image plotting details.
					Title = str(MovieIterlist[l][k])
					#Add Colourbar (Axis, Label, Bins)
					label,bins = VariableLabelMaker(variablelist),5
					cax = Colourbar(ax,label[i],bins,Lim=CbarMinMax(Image))
					#Finalize image
					ImageOptions(fig,ax,Xlabel,Ylabel,Title)

					#Save to seperate folders inside simulation folder.
					num1,num2,num3 = k % 10, k/10 % 10, k/100 % 10
					Number = str(num3)+str(num2)+str(num1)
					savefig(DirMovieplots+variablelist[i]+'_'+Number+ext)
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
					#Take Max value of image for general convergence trend.
					ConvergenceTrends[-1].append( sum(Image.flatten())/len(Image.flatten()) )
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

		#Normalize and plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(ConvergenceTrends)):
			ConvergenceTrends[i] = Normalize(ConvergenceTrends[i],NormFactor=FinalIterationValues[i])[0]
			ax.plot(Xaxis,ConvergenceTrends[i], lw=2)
		#endfor
		Limits = [min(np.asarray(ConvergenceTrends).flatten()),max(np.asarray(ConvergenceTrends).flatten())]
		if Limits[1] > 3*Limits[0]: Limits[1] = 3			#Limit the limits! (avoids early messy figures)

		#Image plotting details.
		Title = 'Convergence of '+str(variablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation Iteration','Normalized Mesh-Average Value'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(Limits[0],Limits[1])

		#Write data to ASCII files if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirConvergence, 'Convergence_Data')
			SaveString = FolderNameTrimmer(Dirlist[l])+'_ConvergenceData'
			WriteDataToFile(variablelist+['\n'], DirASCII+SaveString, 'w')
			#endif
			for i in range(0,len(ConvergenceTrends)):
				WriteDataToFile(ConvergenceTrends[i]+['\n'], DirASCII+SaveString, 'a')
			#endfor
		#endif

		#Print convergence data to terminal if required
		print ''
		print 'Percentage Variation At Last Iteration:'
		for i in range(0,len(ConvergenceTrends)):
			ConvergenceFraction = 1-abs( ConvergenceTrends[i][-1]/ConvergenceTrends[i][-2] )
			ConvergencePercentage = round( (ConvergenceFraction*100), 6)
			print variablelist[i], ConvergencePercentage, '%'
		#endfor

		#Save figure.
		savefig(DirConvergence+FolderNameTrimmer(Dirlist[l])+'_Convergence'+ext)
		plt.close('all')
	#endfor

	print'------------------------------------'
	print'# 2D Convergence Processing Complete'
	print'------------------------------------'
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
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])
		Legendlist = list()

		#Generate the radial (horizontal) lineouts for a specific height.
		if len(radialineouts) > 0:
			#Create folder to keep output plots.
			DirRlineouts = CreateNewFolder(Dirlist[l],'Radial_Profiles/')

			#Loop over all required variables and requested profile locations.
			for i in tqdm(range(0,len(processlist))):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(radialineouts)):
					#Update legend with location of each lineout.
					if len(Legendlist) < len(radialineouts):
						Legendlist.append('Z='+str(round((radialineouts[j])*dz[l], 2))+' cm')
					#endif

					#Plot all requested radial lines on single image per variable.
					Rlineout=PlotRadialProfile(Data[l],processlist[i],Variablelist[i],radialineouts[j])
					#Plot lines for each variable at each requested slice.
					ImagePlotter1D(Rlineout,Raxis,image_aspectratio,fig,ax)

					#Write data to ASCII files if requested.
					if write_ASCII == True:
						SaveString = '_R='+str(round((radialineouts[j])*dz[l], 2))+'cm'
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
		if len(heightlineouts) > 0:
			#Create folder to keep output plots.
			DirZlineouts = CreateNewFolder(Dirlist[l],'Axial_Profiles/')
			Legendlist = list()

			#Collect and plot required data.
			for i in tqdm(range(0,len(processlist))):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(heightlineouts)):
					#Perform SI conversion and save to legend.
					if len(Legendlist) < len(heightlineouts):
						Legendlist.append('R='+str(round(heightlineouts[j]*dr[l], 2))+' cm')
					#endif

					#Plot all requested radial lines on single image per variable.
					Zlineout=PlotAxialProfile(Data[l],processlist[i],Variablelist[i],heightlineouts[j])
					#Plot lines for each variable at each requested slice.
					ImagePlotter1D(Zlineout[::-1],Zaxis,image_aspectratio,fig,ax)

					#Write data to ASCII files if requested.
					if write_ASCII == True:
						SaveString = '_Z='+str(round((heightlineouts[j])*dr[l], 2))+'cm'
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

	print'--------------------------'
	print'# Single Profiles Complete'
	print'--------------------------'
#endif



#====================================================================#
			#COMPARITIVE PROFILES FROM MULTI-FOLDERS#
#====================================================================#



#Plot comparitive profiles for each variable between folders.
#Perform horizontal profile comparisons
if savefig_comparelineouts == True:

	#Create folder to keep output plots.
	DirComparisons = CreateNewFolder(os.getcwd(),'/1D Comparisons')

	#Generate SI scale axes for lineout plots.
	Raxis = GenerateAxis('Radial',Isymlist[l])
	Zaxis = GenerateAxis('Axial',Isymlist[l])

	#Perform radial profile comparisons
	for j in range(0,len(radialineouts)):

		#Create new folder for each axial or radial slice.
		ProfileFolder = 'Z='+str(round((radialineouts[j])*dz[l], 2))+'cm'
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
				Rlineout = PlotRadialProfile(Data[l],processlist[k],Variablelist[k],radialineouts[j],R_mesh[l],Isymlist[l])

				#Plot radial profile and allow for log y-axis if requested.
				ImagePlotter1D(Rlineout,Raxis,image_aspectratio,fig,ax)


				#Write data to ASCII files if requested.
				if write_ASCII == True:
					if l == 0:
						WriteFolder = 'Z='+str(round((radialineouts[j])*dz[l], 2))+'cm_Data'
						DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
						WriteDataToFile(Raxis+['\n'], DirWrite+Variablelist[k], 'w')
					#endif
					WriteDataToFile(Rlineout+['\n'], DirWrite+Variablelist[k], 'a')
				#endif

				#Apply image options and axis labels.
				Title = 'Comparison of '+Variablelist[k]+' Profiles at Z='+str(round((radialineouts[j])*dz[l], 2))+'cm for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel,Legend = 'Radial Distance R [cm]',Ylabels[k],Legendlist
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ Z='+str(round((radialineouts[j])*dz[l], 2))+'cm profiles'+ext)
			plt.close('all')
		#endfor
	#endfor

#===================##===================#

	#Perform vertical lineout comparisons
	for j in range(0,len(heightlineouts)):

		#Create new folder for each axial or radial slice.
		ProfileFolder = 'R='+str(round((heightlineouts[j])*dr[l], 2))+'cm'
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
				Zlineout = PlotAxialProfile(Data[l],processlist[k],Variablelist[k],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])

				#Plot axial profile and allow for log y-axis if requested.
				ImagePlotter1D(Zlineout[::-1],Zaxis,image_aspectratio,fig,ax)


				#Write data to ASCII files if requested.
				if write_ASCII == True:
					if l == 0:
						WriteFolder = 'R='+str(round((heightlineouts[j])*dr[l], 2))+'cm_Data'
						DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
						WriteDataToFile(Zaxis+['\n'], DirWrite+Variablelist[k], 'w')
					#endif
					WriteDataToFile(Zlineout[::-1]+['\n'], DirWrite+Variablelist[k], 'a')
				#endif

				#Apply image options and axis labels.
				Title = 'Comparison of '+Variablelist[k]+' Profiles at R='+str(round((heightlineouts[j])*dr[l], 2))+'cm for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel,Legend = 'Axial Distance Z [cm]',Ylabels[k],Legendlist
				ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ R='+str(round((heightlineouts[j])*dr[l], 2))+'cm profiles'+ext)
			plt.close('all')
		#endfor
	#endfor

	print'-------------------------------'
	print'# Comparitive Profiles Complete'
	print'-------------------------------'
#endif



#====================================================================#
				  #MULTI-PROFILES FROM SAME FOLDER#
#====================================================================#



if savefig_multiprofiles == True:

	#For each folder in turn
	for l in range(0,numfolders):
		#Create global multivar folder.
		Dirlineouts = CreateNewFolder(Dirlist[l],'MultiVar_Profiles/')

		#Create processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])
		multiprocesslist,multiVariablelist = VariableEnumerator(MultiVar,rawdata_2D[l],header_2Dlist[l])

		#Create variable labels with SI unit conversions if required.
		Ylabel = VariableLabelMaker(Variablelist)
		multiYlabel = VariableLabelMaker(multiVariablelist)

		#Generate the vertical (height) lineouts for a given radius.
		if len(heightlineouts) > 0:

			#Generate SI scale axes for lineout plots.
			Zaxis = GenerateAxis('Axial',Isymlist[l])

			#Perform the plotting for all requested variables.
			for i in tqdm(range(0,len(processlist))):

				#Extract the lineout data from the main data array.
				for j in range(0,len(heightlineouts)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((heightlineouts[j])*dr[l], 2))
					DirZlineouts = CreateNewFolder(Dirlineouts,'R='+Slice+'cm/')

					#Create legendlist
					Legendlist = list()
					Legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot the initial variable in processlist first.
					Zlineout = PlotAxialProfile(Data[l],processlist[i],Variablelist[i],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])
					ImagePlotter1D(Zlineout[::-1],Zaxis,image_aspectratio,fig,ax)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):

						#Plot profile for multiplot variables in compareprocesslist.
						Zlineout = PlotAxialProfile(Data[l],multiprocesslist[m],multiVariablelist[m],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])
						ImagePlotter1D(Zlineout[::-1],Zaxis,image_aspectratio,fig,ax)

						#Update legendlist with each variable compared.
						Legendlist.append(VariableLabelMaker(multiVariablelist)[m])
					#endfor


					#Apply image options and axis labels.
					Title = str(round((heightlineouts[j])*dr[l], 2))+'cm Height profiles for '+Variablelist[i]+','' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Axial Distance Z [cm]',VariableLabelMaker(Variablelist)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save figures in original folder.
					R = 'R='+str(round((heightlineouts[j])*dr[l], 2))+'_'
					plt.savefig(DirZlineouts+R+Variablelist[i]+'_MultiProfiles'+ext)
					plt.close('all')
				#endfor
			#endfor
		#endif

#===================##===================#

		#Generate the horizontal (Radial) lineouts for a given radius.
		if len(radialineouts) > 0:
			#Create global multivar folder.
			Dirlineouts = CreateNewFolder(Dirlist[l],'MultiVar_Profiles/')

			#Generate SI scale axes for lineout plots.
			Raxis = GenerateAxis('Radial',Isymlist[l])

			#Perform the plotting for all requested variables.
			for i in tqdm(range(0,len(processlist))):

				#Perform the plotting for all requested variables.
				for j in range(0,len(radialineouts)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((radialineouts[j])*dz[l], 2))
					DirRlineouts = CreateNewFolder(Dirlineouts,'Z='+Slice+'cm/')

					#Create legendlist
					Legendlist = list()
					Legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot profile for initial variable in processlist.
					Rlineout = PlotRadialProfile(Data[l],processlist[i],Variablelist[i],radialineouts[j],R_mesh[l],Isymlist[l])
					ImagePlotter1D(Rlineout,Raxis,image_aspectratio,fig,ax)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):

						#Plot profile for multiplot variables in compareprocesslist.
						Rlineout = PlotRadialProfile(Data[l],multiprocesslist[m],multiVariablelist[m],radialineouts[j],R_mesh[l],Isymlist[l])
						ImagePlotter1D(Rlineout,Raxis,image_aspectratio,fig,ax)

						#Update legendlist with each variable compared.
						Legendlist.append(VariableLabelMaker(multiVariablelist)[m])
					#endfor


					#Apply image options and axis labels.
					Title = str(round((radialineouts[j])*dz[l], 2))+'cm Radial Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]
					Xlabel,Ylabel = 'Radial Distance R [cm]',VariableLabelMaker(Variablelist)[i]
					ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

					#Save lines in previously created folder.
					Z = 'Z='+str(round((radialineouts[j])*dz[l], 2))+'_'
					plt.savefig(DirRlineouts+Z+Variablelist[i]+'_MultiProfiles'+ext)
					plt.close('all')
				#endfor
			#endfor
		#endif
	#endfor

	print'-----------------------------'
	print'# Multiplot Profiles Complete'
	print'-----------------------------'
#endif



#====================================================================#
			  #ITERMOVIE PROFILES - PULSE ANALYSIS#
#====================================================================#



#Plot 1D profile of itervariables at desired locations
if savefig_pulseprofiles == True:

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create new folder and initiate required lists.
		PulseTrends,Xaxis = list(),list()
		DirPulse = CreateNewFolder(Dirlist[l],'Pulse_Profiles/')
		DirMeshAve = CreateNewFolder(DirPulse,'Mesh_Averaged/')

		#Create processlist for each folder as required.
		processlist,variablelist = VariableEnumerator(Variables,rawdata_itermovie[l],header_itermovie[l])
		#Skip over the R and Z processes as they are not saved properly in iterdata.
		for i in range(0,len(processlist)):
			processlist[i] = processlist[i]-2
		#endfor

		#Create list and x-axis for convergence trend plotting.
		#DtActual is approximate, exact dt per iteration depends upon modules and solvers employed.
		DtActual = (1.0/FREQM[l])*100				#s	(~8 microseconds @ 13.56MHz)
		DtActual = DtActual*1000					#ms
		for i in range(0,len(MovieIterlist[l])):
			Xaxis.append( float(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))*DtActual )
		#endfor

		#for all variables requested by the user.
		for i in tqdm(range(0,len(processlist))):

			#Extract 2D image and take mesh averaged value for iteration trend.
			PulseProfile = list()
			for k in range(0,len(MovieIterlist[l])):
				#for further processing.
				Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
				PulseProfile.append( sum(Image.flatten())/len(Image.flatten()) )
			#endfor
			PulseTrends.append(PulseProfile)

			#Plot each variable against simulation real-time.
			fig, ax = plt.subplots(1, figsize=(10,10))
			ax.plot(Xaxis,PulseProfile, lw=2)

			#Image plotting details.
			Title = 'Simulation Time Profile of '+str(variablelist[i])+' for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Simulation time [ms]',VariableLabelMaker(variablelist)[i]
			ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

			#Save figure.
			savefig(DirMeshAve+FolderNameTrimmer(Dirlist[l])+'_'+variablelist[i]+ext)
			plt.close('all')

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(DirPulse, 'Pulse_Data')
				DirWriteMeshAve = CreateNewFolder(DirWrite, 'MeshAveraged_Data')
				WriteDataToFile(Xaxis, DirWriteMeshAve+variablelist[i], 'w')
				WriteDataToFile(['\n']+PulseProfile, DirWriteMeshAve+variablelist[i], 'a')
			#endif
		#endfor

		#=================#

		#Plot mesh averaged value over 'real-time' in simulation.
		Legend = VariableLabelMaker(variablelist)
		fig, ax = plt.subplots(1, figsize=(10,10))

		#Plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(PulseTrends)):
			PulseTrends[i] = Normalize(PulseTrends[i])[0]
			ax.plot(Xaxis,PulseTrends[i], lw=2)
		#endfor

		#Image plotting details.
		Title = 'Simulation Time Profile of '+str(variablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation time [ms]','Normalized Mesh-Average Value'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(0,1.01+(len(Legend)*0.05))

		#Save figure.
		savefig(DirMeshAve+'Normalized_'+FolderNameTrimmer(Dirlist[l])+ext)
		plt.close('all')
	#endfor
	print'-------------------------'
	print'# Pulse Profiles Complete'
	print'-------------------------'
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
		DirEDF = CreateNewFolder(Dirlist[l],'EDFplots')

		#Create processlist for requested EDF species and extract images.
		processlist,variablelist = VariableEnumerator(IEDFVariables,rawdata_IEDF[l],header_IEDFlist[l])
		
		#For all requested variables.
		for i in tqdm(range(0,len(processlist))):
			EDFprofile = list()

			#Extract image from required variable and create required profile lists.
			#Flatten angular distribution across all angles to produce energy distribution.
			Image = ImageExtractor2D(DataIEDF[l][processlist[i]],Rmesh=EDFangle,Zmesh=EDFbins)
			for j in range(0,len(Image)): EDFprofile.append(sum(Image[j]))

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
			Image, EDFprofile = Image[::-1].transpose(), EDFprofile[::-1]

			#Determine region of IEDF to plot based on threshold value from array maximum.
			Threshold = EDF_Threshold*max(EDFprofile)
			for j in range(EDFprofile.index(max(EDFprofile)),len(EDFprofile)): 
				if EDFprofile[j] < Threshold and j != 0: 
					eVlimit = j*deV
 					break
				elif j == len(EDFprofile)-1:
					eVlimit = EMAXIPCMC
				#endif
			#endfor


			#Plot the angular distribution and EDF of the required species.
			fig,ax = figure([11,9], 2, shareX=True)

			Title = Dirlist[l][2::]+'\n'+variablelist[i]+' Angular Energy Distribution Function'
			Extent=[0,EMAXIPCMC, -len(Image)/2,len(Image)/2]

			#Angularly resolved IEDF Figure
			im = ax[0].imshow(Image, extent=Extent, aspect='auto')
			cax = Colourbar(ax[0],variablelist[i]+' EDF($\\theta$)',5)
			Xlabel,Ylabel = '','Angular Dispersion [$\\theta^{\circ}$]'
			ImageCrop = [[0,eVlimit],[-45,45]]					#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax[0],Xlabel,Ylabel,Title,Crop=ImageCrop,Rotate=False) 

			#Angle Integrated IEDF figure
			ax[1].plot(eVaxis,EDFprofile, lw=2)
			cax = InvisibleColourbar(ax[1])
			Xlabel,Ylabel = 'Energy [eV]',variablelist[i]+' EDF \n [$\\theta$ Integrated]'
			ImageCrop = [[0,eVlimit],[0,max(EDFprofile)*1.05]]	#[[X1,X2],[Y1,Y2]]
			ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=ImageCrop,Rotate=False)

			plt.savefig(DirEDF+variablelist[i]+'_EDF'+ext)
			plt.close('all')

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				if i == 0:
					DirASCII = CreateNewFolder(DirEDF, 'EDF_Data')
#					WriteDataToFile(eVaxis, DirASCII+variablelist[i],'w')
				#endif
				WriteDataToFile(Image, DirASCII+variablelist[i]+'_IEDFAngular','w')
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
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
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

			#Mean energy obtained as integrated energy fraction most closely matching total energy
			#averaged over effective energy range determined from EDF_threshold. 
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
				#Single Intersection case
#				Intersection = (np.abs(EDFEnergyProfile-BinAveragedEnergy)).argmin()
	
				#k defines region lower energy edge index
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

				#k defines region lower energy edge index
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

			#Median energy calculated as EDF index representing midpoint of equal integrated energies
#			RisingSum,FallingSum,AbsDiff = 0.0,0.0,list()
#			for j in range(0,len(EDFprofile)):
#				Rising_j, Falling_j = j, (len(EDFprofile)-1-2)-j
#				RisingSum += EDFprofile[Rising_j]*(Rising_j*deV)
#				FallingSum += EDFprofile[Falling_j]*(Falling_j*deV)
#				AbsDiff.append(abs(RisingSum-FallingSum))
				#endif
			#endfor
#			MedianIndex = AbsDiff.index(min(AbsDiff))
#			Median_eV.append( MedianIndex*deV ) 				#### MEDIANS' FUCKED UP BRAH! ####

			#Particle energy variance analysis: Returns FWHM of energy distribution.
#			Take mean and draw line at y = mean 
#			Calculate where y = mean intercepts EDFprofile
#			If only one intercept, first intercept is x = 0
#			Integrate EDFprofile indices between intercepts 
#			Save in 1D array, can be used to get energy spread percentage.

			#==========#
			#==========#

			if DebugMode == True:
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
		ImageCrop = [ [0,GlobRange_eV[1]],[] ]		#[[X1,X2],[Y1,Y2]]
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legendlist,Crop=ImageCrop,Rotate=False)

		plt.savefig(DirIEDFTrends+variablelist[i]+'_EDFprofiles'+ext)
		plt.close('all')


		##ENERGY ANALYSIS##
		#=================#

		#Plot IEDF average energies with respect to simulation folder names.
		fig,ax = figure()
		TrendPlotter(ax,Mean_eV,Legendlist,NormFactor=0)
		TrendPlotter(ax,Mode_eV,Legendlist,NormFactor=0)
#		TrendPlotter(ax,Range_eV[0],Legendlist,NormFactor=0)
#		TrendPlotter(ax,Range_eV[1],Legendlist,NormFactor=0)

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
	print'--------------------------------'
	print'# EEDF/IEDF Processing Complete.'
	print'--------------------------------'
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

	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
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
		for j in range(0,len(heightlineouts)):

			#Create folder for axial trends if needed.
			DirAxialTrends = CreateNewFolder(DirTrends,'Axial Trends')

			#Take Trend at Given Location or Default to Min/Max Trends.
			if len(TrendLocation) == 2:
				#Append requested position to the legendlist.
				R,Z = TrendLocation[0],TrendLocation[1]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			elif len(TrendLocation) == 1:
				#Append requested position to the legendlist.
				R,Z = TrendLocation[0],heightlineouts[j]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			else:
				#Obtain min/max trend values for requested profile over all folders.
				Xaxis,MaxTrend,MinTrend = MinMaxTrends(heightlineouts[j],'Axial',k)
				Trend = MaxTrend
				#Append the radial position to the legendlist.
				Legendlist.append( 'R='+str(round((heightlineouts[j]*dr[l]), 2))+'cm' )
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
		if len(heightlineouts) > 0:
			plt.savefig(DirAxialTrends+'Axial Trends in '+Variablelist[k]+ext)
			plt.close('all')
		#endif


		##RADIAL TRENDS##
		#===============#

		#Create fig of desired size and refresh legendlist.
		fig,ax = figure(image_aspectratio,1)
		Legendlist = list()

		#Perform trend analysis on requested radial profiles.
		for j in range(0,len(radialineouts)):

			#Create folder for axial trends if needed.
			DirRadialTrends = CreateNewFolder(DirTrends,'Radial Trends')

			#Take Trend at Given Location or Default to Min/Max Trends.
			if len(TrendLocation) == 2:
				#Append requested position to the legendlist.
				R,Z = TrendLocation[0],TrendLocation[1]
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			elif len(TrendLocation) == 1:
				#Append requested position to the legendlist.
				R,Z = radialineouts[j],TrendLocation[0],
				Location = '(R'+str(round(R*dr[l],1))+'cm, Z'+str(round(Z*dz[l],1))+'cm)'
				Legendlist.append(Location)
				#Take trend at given location if specified.
				Xaxis,Trend = TrendAtGivenLocation([R,Z],processlist[k],Variablelist[k])

			else:
				#Obtain min/max trend values for requested profile over all folders.
				Xaxis,MaxTrend,MinTrend = MinMaxTrends(radialineouts[j],'Radial',k)
				Trend = MaxTrend
				#Append the axial position to the legendlist.
				Legendlist.append( 'Z='+str(round((radialineouts[j]*dz[l]), 2))+'cm' )
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
		if len(radialineouts) > 0:
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
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
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
		try: Rlineout = PlotRadialProfile(Data[l],Process[0],Variable[0],Rlineoutloc,R_mesh[l],  Isymlist[l])
		except: Rlineout = float('NaN')
		#endtry
		try: Zlineout = PlotAxialProfile(Data[l],Process[0],Variable[0],Zlineoutloc,R_mesh[l],Z_mesh[l],Isymlist[l])
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
			print Dirlist[l]
			print 'DC Bias:',round(DCbias[l],5),'V'
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
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Create required lists.
	RequestedPowers,DepositedPowerList = list(),list()
	Xaxis,Powers = list(),list()

	#Update X-axis with folder information.
	for l in range(0,numfolders): Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

	#Identify which power densities have been requested.
	for i in range(0,len(Variables)):
		#Avoids POW-ICP unless coils are used
		if len(FREQC[l]) > 0 and Variables[i] in ['POW-ALL','POW-TOT','POW-ICP','POW-RF','POW-RF-E']:
			RequestedPowers.append(Variables[i])	
		#All powers here should be saved in all simulations (Including DC ones)		
		elif Variables[i] in ['POW-ALL','POW-TOT','POW-RF','POW-RF-E']:
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

			Power = 0
			#Cylindrical integration of power per unit volume ==> total coupled power.
			for j in range(0,Z_mesh[l]):
				#For each radial slice
				for i in range(0,R_mesh[l]-1):
					#Calculate radial plane volume of a ring at radius [i], correcting for central r=0.
					InnerArea = np.pi*( (i*(dr[l]/100))**2 )		#m^2
					OuterArea = np.pi*( ((i+1)*(dr[l]/100))**2 )	#m^2
					RingVolume = (OuterArea-InnerArea)*(dz[l]/100)	#m^3

					#Calculate Power by multiplying power density for ring[i] by volume of ring[i]
					Power += PowerDensity[j][i]*RingVolume 			#W
				#endfor
			#endfor
			DepositedPowerList.append(Power)

			#Display power to terminal if requested.
			if print_totalpower == True:
				print Dirlist[l]
				print RequestedPowers[k]+' Deposited:',round(Power,4),'W'
			#endif
		#endfor

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
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
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

			#Extract data required for Thrust calculations, discharge plane (Z) = ThrustLoc.
			processlist,variablelist = VariableEnumerator(['VZ-NEUTRAL'],rawdata_2D[l],header_2Dlist[l])
			NeutralVelocity = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
			processlist,variablelist = VariableEnumerator(['VZ-ION+'],rawdata_2D[l],header_2Dlist[l])
			IonVelocity = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
			processlist,variablelist = VariableEnumerator(['FZ-AR3S'],rawdata_2D[l],header_2Dlist[l])
			NeutralAxialFlux = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
			processlist,variablelist = VariableEnumerator(['FZ-AR+'],rawdata_2D[l],header_2Dlist[l])
			IonAxialFlux = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
			processlist,variablelist = VariableEnumerator(['TG-AVE'],rawdata_2D[l],header_2Dlist[l])
			NeutGasTemp = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
			processlist,variablelist = VariableEnumerator(['PRESSURE'],rawdata_2D[l],header_2Dlist[l])
			try: 
				Pressure = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc)
				PressureDown = PlotRadialProfile(Data[l],processlist[0],variablelist[0],ThrustLoc+1)
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

					#Calculate differential pressure between ThrustLoc-(ThrustLoc+1)
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
				print Dirlist[l], '@ Z=',round(ThrustLoc*dz[l],2),'cm'
				print 'NeutralThrust', round(NeutralThrust*1000,2), 'mN @ ', round(NeutralIsp,2),'s'
				print 'IonThrust:', round(IonThrust*1000,4), 'mN @ ', round(IonIsp,2),'s'
				print 'D-Pressure:', round(DiffForce*1000,4), 'mN'
				print 'Thrust:',round(Thrust*1000,4),'mN @ ', round(ThrustIsp,2),'s'
				print ''
			#endif
		#endfor

		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
			WriteDataToFile(Xaxis+['\n'], DirASCII+'Thrust_Trends','w')
			WriteDataToFile(Thrustlist+['\n'], DirASCII+'Thrust_Trends','a')
			WriteDataToFile(ThrustIsplist, DirASCII+'Thrust_Trends','a')
		#endif


		#Plot total thrust and ion/neutral components.
		fig,ax1 = figure(image_aspectratio,1)
		TrendPlotter(ax1,Thrustlist,Xaxis,Marker='ko-',NormFactor=0)
		TrendPlotter(ax1,NeutralThrustlist,Xaxis,Marker='r^-',NormFactor=0)
#		TrendPlotter(ax1,IonThrustlist,Xaxis,Marker='bs-',NormFactor=0)

		#Apply image options and save figure.
		Title='Thrust at Z='+str(round(ThrustLoc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Thrust F$_{T}$ [mN]'
		ax1.legend(['Total Thrust','Neutral Component','Ion Component'], fontsize=18, frameon=False)
		ImageOptions(fig,ax1,Xlabel,Ylabel,Title,Crop=False)

		plt.savefig(DirTrends+'Thrust Trends'+ext)
		plt.close('all')


		#Plot Specific Impulse for total thrust and ion/neutral components.
		fig,ax1 = figure(image_aspectratio,1)
		TrendPlotter(ax1,ThrustIsplist,Xaxis,Marker='ko-',NormFactor=0)
		TrendPlotter(ax1,NeutralIsplist,Xaxis,Marker='r^-',NormFactor=0)
#		TrendPlotter(ax1,IonIsplist,Xaxis,Marker='bs-',NormFactor=0)

		#Apply image options and save figure.
		Title = 'Specific Impulse at Z='+str(round(ThrustLoc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
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

	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initialize any required lists.
	Xaxis,SxLocExtent,SxMaxExtent = list(),list(),list()	

	#Obtain SheathROI and SourceWidth automatically if none are supplied.
	if len(SheathROI) != 2:
		#image_radialcrop Convert to Cells 
		#image_axialcrop Convert to Cells
		#Use axialcrop or radialcrop to set automatic ROI!
		Start,End = 34,72			#AUTOMATIC ROUTINE REQUIRED#
		SheathROI = [Start,End]		#AUTOMATIC ROUTINE REQUIRED#
	#endif
	if len(SourceWidth) == 0:
		#Take Variable that is zero in metals (Density?)
		#Take Axial/Radial slice depending on sheath direction.
		#Find Cell distance from zero to 'wall' at electrodeloc.
		#Convert to SI [cm], set to automatic width.
		SourceWidth = [0.21]			#AUTOMATIC ROUTINE REQUIRED#
	#endif

	SxMeanExtent,SxMeanExtentArray = list(),list()
	#For all selected simulations, obtain Xaxis, sheath value and save to array.
	for l in range(0,numfolders):
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Obtain sheath thickness array for current folder 
		Sx = SheathExtent(folder=l)[0]

		#Calculate mean sheath extent across ROI. On failure provide null point for sheath thickness.
		try:
			SxMeanExtentArray = list()
			for i in range(SheathROI[0],SheathROI[1]):	SxMeanExtentArray.append(Sx[i])
			SxMeanExtent.append(sum(SxMeanExtentArray)/len(SxMeanExtentArray))
		except:
			SxMeanExtent.append( np.nan )
		#endtry

		#Extract maximum sheath thickness from within region of interest
		try: SxMaxExtent.append( ((SourceWidth[0]*dr[l])-max(Sx[SheathROI[0]:SheathROI[1]]))*10 )
		except: SxMaxExtent.append( np.nan )

		#Extract sheath width adjacent to powered electrode
		#loc = electrodeloc[0]		#Radial
		loc = electrodeloc[1] 		#Axial
		try: SxLocExtent.append( ((SourceWidth[0]*dr[l])-Sx[loc])*10 )
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
	TrendPlotter(ax,SxMaxExtent,Xaxis,NormFactor=0)

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
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		KnudsenAverage,Xaxis = list(),list()

		#For all folders.
		for l in range(0,numfolders):

			#Using effective radius of argon in this calculation.
			Dimentionality = 2*(Radius[l]/100)		#meters
			CrossSection = np.pi*((7.1E-11)**2)		#meters

			#Extract data for the neutral flux and neutral velocity.
			processlist,Variablelist = VariableEnumerator(FluidSpecies,rawdata_2D[l],header_2Dlist[l])

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

					LocalDensity = (Data[l][processlist[0]][Start+i])*1E6
					try:
						KnudsenNumber = (1/(LocalDensity*CrossSection*Dimentionality))
					except:
						KnudsenNumber = 0
					#endtry
					Image[Row,i] = KnudsenNumber
				#endfor
			#endfor

			#Display average Knudsen number to terminal if requested.
			KnudsenAverage.append( sum(Image)/(len(Image[0])*len(Image)) )
			if print_Knudsennumber == True:
				print Dirlist[l]
				print 'Average Knudsen Number:', KnudsenAverage[l]
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
			#Add sheath thickness to figure if requested.
			Sx = SheathExtent(folder=l,ax=ax)[0]

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
				  #LOCAL/GLOBAL SOUND SPEED ANALYSIS#
#====================================================================#

PERFORMSOUNDSPEED = False
#Only perform on bulk fluid dynamics relevent species.
if bool(set(FluidSpecies).intersection(Variables)) == True and PERFORMSOUNDSPEED == True:
	if savefig_trendphaseaveraged == True or print_soundspeed == True:

		#Create Trend folder to keep output plots.
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		AverageSoundSpeed,Xaxis = list(),list()
		NeutralDensities = list()

		#For all folders.
		for l in range(0,numfolders):

			#Extract spatially resolved pressure and neutral densities.
			processlist,variablelist = VariableEnumerator(['PRESSURE'],rawdata_2D[l],header_2Dlist[l])
			Pressure = ImageExtractor2D(Data[l][processlist[0]],variablelist[0])
			#If single neutral species - extract density
			processlist,Variablelist = VariableEnumerator(FluidSpecies,rawdata_2D[l],header_2Dlist[l])
			if len(processlist) == 1: 
				NeutralDensity = ImageExtractor2D(Data[l][processlist[0]],variablelist[0])
			#If multiple neutral species, combine them to get total neutral density
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
			Image = LocalSoundSpeed(NeutralDensity,Pressure,Dimension='2D')

			#Display mesh-averaged sound speed to terminal if requested.
			AverageSoundSpeed.append( sum(Image)/(len(Image[0])*len(Image)) )
			if print_soundspeed == True:
				print Dirlist[l]
				print 'Average Sound Speed:', AverageSoundSpeed[l]
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
			#Add sheath thickness to figure if requested.
			Sx = SheathExtent(folder=l,ax=ax)[0]

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

if any([savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_soundspeed, print_totalpower, print_DCbias, print_thrust, print_sheath]) == True:
	print'---------------------------'
	print'# Trend Processing Complete'
	print'---------------------------'
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
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Phaselist)
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])

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
			Lineouts,LineoutsOrientation = list(),list()

			#Concatinate all requested lineouts together, keeping seperate orientation.
			for m in range(0,len(radialineouts)):
				Lineouts.append(radialineouts[m])
				LineoutsOrientation.append('Radial')
			for m in range(0,len(heightlineouts)):
				Lineouts.append(heightlineouts[m])
				LineoutsOrientation.append('Axial')

			#For all requested lineouts and orientations.
			for k in range(0,len(Lineouts)):

				#Refresh required lists.
				VariableMax,VariableMin = list(),list()

				#Create folders to keep output plots for each variable.
				if LineoutsOrientation[k] == 'Axial':
					NameString= varlist[i]+'_'+str(round(Lineouts[k]*dr[l],2))+'cm[R]'
				if LineoutsOrientation[k] == 'Radial':
					NameString= varlist[i]+'_'+str(round(Lineouts[k]*dz[l],2))+'cm[Z]'
				if savefig_phaseresolve1D == True:
					Dir1DProfiles = CreateNewFolder(DirMovieplots,NameString)
				#endif

				#Collect Normalization data for plotting.
				for j in range(0,len(Phaselist)):
					#Record local maximum and minimum for each phase.
					if LineoutsOrientation[k] == 'Axial': Profile = PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],Lineouts[k])
					elif LineoutsOrientation[k] == 'Radial': Profile = PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],Lineouts[k])
					#endif
					Profile,Minimum,Maximum = Normalize(Profile)
					VariableMax.append(Maximum)
					VariableMin.append(Minimum)
				#endfor
				#Find global maximum and minimum for all phases.
				VariableMax = max(VariableMax)
				VariableMin = min(VariableMin)

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Phaselist)):

					if LineoutsOrientation[k] == 'Axial':
						ZlineoutLoc,axis = Lineouts[k],Zaxis
						PhaseResolvedProfile = PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])[::-1]
						lineoutstring = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm \n'
						Xlabel = 'Axial Distance Z [cm]'
					elif LineoutsOrientation[k] == 'Radial':
						RlineoutLoc,axis = Lineouts[k],Raxis
						PhaseResolvedProfile = PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],RlineoutLoc,R_mesh[l],Isymlist[l])
						lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm \n'
						Xlabel = 'Radial Distance R [cm]'
					#endif

					#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
					fig,ax = figure(image_aspectratio,2)
					Ylabel = VariableLabelMaker(varlist)
					fig.suptitle('Phase-Resolved '+varlist[i]+' for '+VariedValuelist[l]+lineoutstring+str(Phaselist[j]), y=0.97, fontsize=16)

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
					fig.tight_layout()
					plt.subplots_adjust(top=0.90)
					num1,num2,num3 = j % 10, j/10 % 10, j/100 % 10
					Number = str(num3)+str(num2)+str(num1)
					plt.savefig(Dir1DProfiles+NameString+'_'+Number+ext)
					plt.close('all')

					#Write Phase data in ASCII format if required.
					if write_ASCII == True:
						DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'1DPhase_Data')
						DirASCIIPhaseloc = CreateNewFolder(DirASCIIPhase,lineoutstring[3:-2])
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

	print'---------------------------------------'
	print'# 1D Phase-Resolved Processing Complete'
	print'---------------------------------------'
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
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Phaselist)
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])

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

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(PhaseData[j][proclist[i]],varlist[i])
				#Extract Ni and Ne variables for sheath processing.
				if image_sheath == True:
					Ne = SxData[j][Sxproc[Sxvar.index('E')]]
					Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
				else: Ne,Ni = np.nan, np.nan
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
				fig,ax[0],im,Image = ImagePlotter2D(Image,extent,aspectratio,varlist[i],fig,ax[0],"upper")
				Sx = SheathExtent(folder=l,ax=ax[0],Phase=j,Ne=Ne,Ni=Ni)[0]
				ImageOptions(fig,ax[0],Xlabel,Ylabel,Crop=True)
				#Add Colourbar (Axis, Label, Bins)
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

				#Cleanup layout and save images.
				fig.tight_layout()
				plt.subplots_adjust(top=0.90)
				num1,num2,num3 = j % 10, j/10 % 10, j/100 % 10
				Number = str(num3)+str(num2)+str(num1)
				savefig(DirMovieplots+varlist[i]+'_'+Number+ext)
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

	print'---------------------------------------'
	print'# 2D Phase-Resolved Processing Complete'
	print'---------------------------------------'
#endif




#====================================================================#
				#PHASE TRENDS & SHEATH DYNAMICS#
#====================================================================#

#Process phase resolved data from multiple folders to extract trends.
if savefig_trendphaseresolved == True:

	#Initiate arrays between folders
	VoltageWaveforms,WaveformBiases = list(),list()
	SxMaxExtTrend,SxMeanExtTrend= list(),list()
	SxMaxVelTrend,SxMeanVelTrend = list(),list()
	SxDynRangeTrend = list()
	VariedValuelist = list()

	#Read data from each simulation folder individually, saves on RAM.
	for l in tqdm(range(0,numfolders)):

		#Create global folders to keep output plots.
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
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
		Phaseaxis = GenerateAxis('Phase',Isym=Isymlist[l])

		#=============#

		SxLoc = list()
		#For all phases, calculate sheath width and record sheath width at electrodeloc
		for k in range(0,len(SxPhase)):
			#Extract Ni and Ne variables for sheath processing.
			Ne = SxData[k][Sxproc[Sxvar.index('E')]]
			Ni = SxData[k][Sxproc[Sxvar.index('AR+')]]

			#calculate sheath width employing 'E' and 'AR+'
			Sx = SheathExtent(folder=l,Phase=k,Ne=Ne,Ni=Ni)[0]
			for j in range(0,len(Sx)): 
				Sx[j] = ((SourceWidth[0]*dr[l])-Sx[j])*10	#Convert to mm
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
		RFPeriod = 1.0/MINFREQ[l]										#[s]
		MeanSheathExtent = (sum(SxLoc)/len(SxLoc))/1E6					#[km]
		MeanSheathVelocity = (2*MeanSheathExtent)/RFPeriod				#[km/s]
		SxMeanVelTrend.append( MeanSheathVelocity )						#[km/s]		

		#Calculate maximum instantaneous sheath velocity.
		#Assumes sheath collapse velocity > sheath expansion velocity
		Collapsed,CollapsedPhase = min(SxLoc),SxLoc.index(min(SxLoc))
		Extended,ExtendedPhase = max(SxLoc),SxLoc.index(max(SxLoc))

		SheathExtension = (Extended-Collapsed)/1000.0  					#[m]
		PhaseResolution = 1.0/(MINFREQ[l]*len(Phaseaxis))				#[s]
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
			print Dirlist[l][2:-1]
			print 'Sheath Extension:',round(SheathExtension*1E3,2),'[mm]'
			print 'Sheath Collapse Time:',round(SheathTime*1E9,2),'[ns]'
			print 'Mean Sheath Velocity:',round(MeanSheathVelocity/1E3,2),'[km/s]'
			print 'Max Sheath Velocity:',round(MaxSheathVelocity/1E3,2),'[km/s]'
			print ''
		#endif

		#=============#

		#Plot phase-resolved sheath extension for current folder
		fig,ax = figure(image_aspectratio,2,shareX=True)
		ax[0].plot(Phaseaxis,SxLoc, lw=2)
		Ylabel = 'Sheath Extension [mm]'
		ImageOptions(fig,ax[0],Ylabel=Ylabel,Crop=False)
		ax[0].set_xticks([])

		#Plot Waveform onto Temporally collapsed PROES.
		ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
		ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
		ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.5))

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


	print'----------------------------------------'
	print'# Phase-Resolved Trend Analysis Complete'
	print'----------------------------------------'
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

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Phaselist)
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])

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
			Lineouts,LineoutsOrientation = list(),list()

			#Concatinate all requested lineouts together, keeping seperate orientation.
			for m in range(0,len(radialineouts)):
				Lineouts.append(radialineouts[m])
				LineoutsOrientation.append('Radial')
			for m in range(0,len(heightlineouts)):
				Lineouts.append(heightlineouts[m])
				LineoutsOrientation.append('Axial')

			#For all requested lineouts and orientations.
			for k in range(0,len(Lineouts)):
				#Refresh required lists.
				VariableMax,VariableMin = list(),list()
				PhaseSx,SymPhaseSx = list(),list()
				PROES = list()

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Phaselist)):
					#Refresh lists between each phasecycle.
					IntegratedDoFArray,DoFArrays = list(),list()

					#Collect each profile for stitching into a PROES image if required.
					if DoFWidth > 0:
						#Determine range of lineouts within the depth of field.
						DOFRegion = [(Lineouts[k]-DoFWidth),(Lineouts[k]+DoFWidth)]
						#If DOF extends beyond mesh, alert user and abort diagnostic.
						if any(DOFRegion) < 0:
							print '----------------------------------'
							print 'Depth-of-Field Extends Beyond Mesh'
							print '----------------------------------'
							break
						#endif

						#Collect profiles from DOF region and transpose to allow easy integration.
						for LineoutLoc in range(DOFRegion[0],DOFRegion[1]):
							if LineoutsOrientation[k] == 'Radial': DoFArrays.append(PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							elif LineoutsOrientation[k] == 'Axial': DoFArrays.append(PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							#endif
						#endfor
						DoFArrays = np.asarray(DoFArrays).transpose().tolist()

						#Integrate DoF profiles spatially, obtaining PROES profile for phase 'j'
						for n in range(0,len(DoFArrays)):
							IntegratedDoFArray.append( sum(DoFArrays[n])/(DoFWidth*2+1) )
						#endif
						PROES.append(IntegratedDoFArray)

					#If no DoF then simply collect lineout from required location.
					elif DoFWidth == 0:
						LineoutLoc = Lineouts[k]
						if LineoutsOrientation[k] == 'Radial':
							PROES.append(PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc))
						elif LineoutsOrientation[k] == 'Axial':
							PROES.append(PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
						#endif
					#endif 

					if image_sheath == True:
						#Extract Ni and Ne variables and perform sheath processing.
						Ne = SxData[j][Sxproc[Sxvar.index('E')]]
						Ni = SxData[j][Sxproc[Sxvar.index('AR+')]]
						#Extract the full spatial sheath extent for phase 'j', save location for PROES.
						Sx = SheathExtent(folder=l,Orientation=LineoutsOrientation[k],Phase=j,Ne=Ne,Ni=Ni)
						PhaseSx.append(Sx[0][Lineouts[k]])
						SymPhaseSx.append(Sx[1][Lineouts[k]])
					#endif
				#endfor

				#Scale PROES image and Sx arrays by required number of phasecycles.
				ScaledPROES,ScaledPhaseSx,ScaledSymPhaseSx = list(),list(),list()
				for n in range(0,int(phasecycles*len(PROES))):
					Index = n % len(PROES)		#Modulo index for multiple phasecycles
					ScaledPROES.append(PROES[Index])
					if image_sheath == True:
						ScaledPhaseSx.append(PhaseSx[Index])
						ScaledSymPhaseSx.append(SymPhaseSx[Index])
					#endif
				#endfor
				try: PROES,PhaseSx,SymPhaseSx = ScaledPROES,ScaledPhaseSx,ScaledSymPhaseSx
				except: PROES = ScaledPROES

				#Create figure and rotate PROES such that phaseaxis aligns with waveform.
				fig,ax = figure(image_aspectratio,2,shareX=True)
				PROES = ndimage.rotate(PROES, 90)

				#Choose correct axial or radial axis and set image dimensions and labels
				#IF PROES image plane is axial:
				if LineoutsOrientation[k] == 'Axial':
					lineoutstring = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm'
					NameString = varlist[i]+'_'+lineoutstring[2::]
					Ylabel = 'Axial Distance Z [cm]'
					x1,x2 = Phaseaxis[0],Phaseaxis[-1]
					y1,y2 = Zaxis[-1],Zaxis[0] 						#Y1,Y2 Reversed 	- see below
					Crop = [image_radialcrop,image_axialcrop[::-1]]	#AxialCrop Reversed - see below
					origins = ['upper','top']						#HPEM uses upper origin (Z=0.0cm at top of mesh)

				#IF PROES image plane is radial, with symmetry:
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 1:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					Crop = [image_radialcrop,image_axialcrop]		#Crop is not reversed !!this is how PhD was set!!
					x1,x2 = Phaseaxis[0],Phaseaxis[-1]
					y1,y2 = Raxis[-1],-Raxis[-1]
					origins = ['lower','bottom']					#HPEM has R=0.0 on-axis, increasing outwards

				#IF PROES image plane is radial, without symmetry:
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 0:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					x1,x2 = Phaseaxis[0],Phaseaxis[-1]
					y1,y2 = Raxis[-1],0
					Crop = [image_axialcrop,image_radialcrop]		#Crop is reversed as Y-axis represents R for PROES
					origins = ['lower','bottom']					#HPEM has R=0.0 on-axis, increasing outwards
				#endif
				DirPROESloc = CreateNewFolder(DirPROES,lineoutstring[3::])

				#Create PROES image along line of sight with phase-locked waveform.
				fig.suptitle( 'Simulated '+varlist[i]+' PROES for '+VariedValuelist[l]+lineoutstring+'\n DoF = '+str(round(((2*DoFWidth)+1)*dz[l],2))+' cm', y=0.95, fontsize=18)
				if image_contourplot == True: im = ax[0].contour(PROES,extent=[x1,x2,y1,y2],origin='lower')
				im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='bottom',aspect='auto')
				if image_sheath == True:
					ax[0].plot(Phaseaxis,PhaseSx,'w--',lw=1.0)
					ax[0].plot(Phaseaxis,SymPhaseSx,'w--',lw=1.0)
				#endif
				ImageOptions(fig,ax[0],Xlabel='',Ylabel=Ylabel,Crop=Crop)
				ax[0].set_xticks([])
				ax[0].set_xlim(x1,x2)
				#Add Colourbar (Axis, Label, Bins)
				label = VariableLabelMaker(varlist)
				Colourbar(ax[0],label[i],5,Lim=CbarMinMax(PROES,LineoutsOrientation[k]))

				#Plot Waveform.
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
					DirASCIIPROESloc = CreateNewFolder(DirASCIIPROES,lineoutstring[3::])
					WriteDataToFile(PROES, DirASCIIPROESloc+varlist[i]+'_PROES')
				#endif



				#==============##==============#

				#Temporally collapse 2D PROES image through defined phase fraction.
				SpatialPROES = list()
				for m in range(0,len(PROES)):
					SpatialPROES.append( (sum(PROES[m][::])) ) 	#NEEDS 'MATHS-ING'
				#endfor

				#Spatially Collapse 2D PROES image along the line of sight.
				PROES,TemporalPROES = PROES.transpose().tolist(),list()
				for m in range(0,len(PROES) ):
					TemporalPROES.append( (sum(PROES[m][::]))/(len(PROES)/2) )
				#endfor

				#Plot Spatially Collapsed PROES with phase axis.
				fig,ax = figure(image_aspectratio,2,shareX=True)
				ax[0].plot(Phaseaxis,TemporalPROES, lw=2)
				Ylabel = 'Spatially Integrated '+varlist[i]
				ImageOptions(fig,ax[0],Ylabel=Ylabel,Crop=False)
				ax[0].set_xticks([])
				ax[0].set_xlim(x1,x2)

				#Plot Waveform onto Temporally collapsed PROES.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(fig,ax[1],Xlabel,Ylabel,Crop=False)
				ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.25))
				ax[1].set_xlim(x1,x2)

				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' TemporalPROES'+ext)
				plt.close('all')

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

				#plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' SpatialPROES'+ext)
#				plt.close('all')

				#==============##==============#

			#endfor
			if l == numfolders and k == len(Lineouts):
				print'-------------------------------'
				print'# PROES Image analysis Complete'
				print'-------------------------------'
			#endif
		#endfor
	#endfor
#endif

#===============================#

if any([savefig_trendphaseresolved, savefig_phaseresolve1D, savefig_phaseresolve2D, savefig_PROES]) == True:
	print'----------------------------------'
	print'# Phase Resolved Profiles Complete'
	print'----------------------------------'
#endif

#=====================================================================#
#=====================================================================#
























































#====================================================================#
							  #CODE DUMP#
#====================================================================#

#===============================#
#     Paper Trend Locations     #
#===============================#

#SDoyle2017a: 
#dz(5.50/118), dr(2.55/102) height=[24,43], Trend=[19]

#SDoyle2018a: 
#PROES, (Z=14.2, 21.0, 31.0) radialineouts = [29,44,64]
#Dielectric locations: [16,29],[16,44],[16,64]

#===============================#
# Archived Switchboard Settings #
#===============================#

#### PP-SCCP ####
#electrodeloc = 	[0,3]
#waveformlocs = 	[]
#DOFWidth = 		R;??,Z;??
#TrendLoc =  		H[0];R[]
#ThrustLoc = 		[]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[];Z[]

#### TSHC-2017 ####
#electrodeloc = 	[0,15]
#waveformlocs = 	[]
#DOFWidth = 		R;5,Z;10
#TrendLoc =  		H[0,20];R[30,60,90]
#ThrustLoc = 		[]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[0.0,1.0];Z[0.5,2.5]

### TSHC-OI Mk3 ###
#electrodeloc = 	[58,15]
#waveformlocs = 	[]
#DOFWidth = 		R;??,Z;??
#TrendLoc =  		H[0,23,45];R[46,55,64]			#R,Z = 0.2cm/cell,0.1cm/cell
#ThrustLoc = 		[]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[0,12];Z[4,7]

### HYPERION-I Mk1 ###
#electrodeloc = 	[51,14]HYPI OR [12,28]HYPII			#Upper(Positive) ICP coil
#waveformlocs = 	[[51,24][51,34]]					#Middle ICP Coil, Lower(Negative) coil
#DOFWidth = 		R;??,Z;??
#TrendLoc =  		H[0];R[50]HYPI OR H[56];R[50]HYPII	#R,Z = 0.2cm/cell,0.1cm/cell
#ThrustLoc = 		[]
#SheathROI = 		[]
#SourceWidth = 		[]
#Crop = 			R[];Z[]

### EVgeny Mk1 ###
#electrodeloc = 	[31,14]							#Middle ICP Coil
#waveformlocs = 	[[31,6],[31,23],[20,6]]			#[UpstreamCoil],[DownstreamCoil],[DielectricSurface]
#DOFWidth = 		R;??,Z;??
#TrendLoc =  		H[0,21,41];R[]					#R,Z = 0.2cm/cell,0.1cm/cell
#ThrustLoc = 		[]
#SheathROI = 		[45,85]							#Downstream
#SourceWidth = 		[90]							#Downstream
#Crop = 			R[];Z[]
#Plotmesh = 		'EVgeny'

#===============================#
#             Notes             #
#===============================#

# Disabled the following warning message regarding scalar assignment of 'arr' axis.
# /home/sjd549/.local/lib/python2.7/site-packages/numpy/ma/core.py:6385


#====================================================================#
				  	#GRAPHICAL USER INTERFACE#
#====================================================================#


use_GUI = False
if use_GUI == True:
	try:
		# Python2
		import Tkinter as tk
		from ttk import *
	except ImportError:
		# Python3
		import tkinter as tk
		from ttk import *
	#endtry

	#Create switchboard directory for GUI.
	Switchboard = {}

#=============#

	#Toggles button font between green/red and outputs corresponding true/false.
	def toggle(string,btn):
		#to get the present state of the toggle button
		if btnlist[btn].config('fg')[-1] == 'green':
			btnlist[btn].config(fg='red')
		else:
			btnlist[btn].config(fg='green')
		#endif

		#Add name of button and true/false value to switchboard dictionary.
		global Switchboard
		if btnlist[btn].config('fg')[-1] == 'green':
			Switchboard[string] = 'True'
		elif btnlist[btn].config('fg')[-1] == 'red':
			Switchboard[string] = 'False'
		#endif
	#enddef

#=============#

	#Initialize and configure root window.
	root = tk.Tk()
	root.title("Hpem ELectronic ENgine Analysis v0.8.4")
	root.minsize(width=800, height=600)
	root.maxsize(width=4*1080, height=1080)

	#Buttonlist as displayed in GUI.
	btnlist = ['2D Images','2D Converge','2D Phase']

	#Add toggle buttons.  (lambda returns a reference to a nameless function)
	btnlist[0] = tk.Button(text=btnlist[0], width=12, fg='red')
	btnlist[0]["command"] = lambda: toggle('savefig_plot2D',0)
	btnlist[0].grid(row=1, column=0)

	btnlist[1] = tk.Button(text=btnlist[1], width=12, fg='red')
	btnlist[1]["command"] = lambda: toggle('savefig_convergence',1)
	btnlist[1].grid(row=1, column=1)

	btnlist[2] = tk.Button(text=btnlist[2], width=12, fg='red')
	btnlist[2]["command"] = lambda: toggle('savefig_phaseresolve',2)
	btnlist[2].grid(row=1, column=2)

	#Add run button, disables GUI and progresses to the main program.
	Run = tk.Button(text='Run Analysis', height=4, width=18, fg='black')
	Run["command"] = root.destroy
	Run.grid(row=3,column=3,padx=0,pady=100)

	root.mainloop()
	print Switchboard
#endif

#=====================================================================#
#=====================================================================#








