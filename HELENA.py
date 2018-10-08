#!/usr/bin/env python

#################################
#		Point of Contact		#
#								#
#	   Mr. Scott J. Doyle		#
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

#Import core modules
import matplotlib.cm as cm
import numpy as np
import scipy as sp
import math as m
import os, sys
import os.path

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *



#====================================================================#
				  		#DEFAULT PARAMETERS#
#====================================================================#

#Various debug and streamlining options.
Magmesh = 1							#initmesh.exe magnification factor. (almost obsolete)
DisableMovie = False				#Suppresses ffmpeg routines, saves RAM.
DebugMode = False					#Produces debug outputs for relevent diagnostics.

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images

#List of recognized data extensions for file readin
FileExtensions = ['.PDT','.pdt','.nam']


#Calculation Methods:
GlobSheathMethod = 'AbsDensity'		#Set Global Sheath Calculation Method.
#Choices: ('AbsDensity','IntDensity')
GlobThrustMethod = 'AxialMomentum'	#Set Global Thrust Calculation Method. 
#Choices:('ThermalVelocity','AxialMomentum')
DCbiasaxis = 'Auto'					#Direction to calculate dc bias over.
#Choices:('Axial','Radial','Auto')



#List of recognised neutral/metastable atomic density sets, add new sets as required.
ArgonReduced = ['AR','AR+','AR*']
ArgonFull = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR+','AR2+','AR2*']
Oxygen = ['O','O+','O-','O*','O2','O2+','O2*']
AtomicSet = ['E']+ArgonReduced+ArgonFull+Oxygen

#List of recognized ground-state neutral species for fluid analysis.
NeutSpecies = ['AR','AR3S','O2']

#Commonly used variable sets.
Phys = ['P-POT','TE','EF-TOT','EAMB-Z','EAMB-R','RHO','BT','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','EFLUX-R','EFLUX-Z','TG-AVE','PRESSURE','POW-RF','POW-RF-E','EB-ESORC']
Ar = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR','AR+','AR2+','AR2*','E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','FZ-AR3S','FR-AR3S','FR-AR+','FZ-AR+','FZ-AR3S','FR-AR3S']+Phys
O2 = ['O2','O2+','O','O+','O-','E','FR-O-','FZ-O-']+Phys

Ar_Phase = ['S-E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','SRCE-2437','TE','PPOT','FR-E','FZ-E','SEB-AR4P','SEB-AR+']

ESCT_PCMC = ['AR^0.3S','EB-0.3S','ION-TOT0.3S']
MSHC_PCMC = ['AR^0.5S','EB-0.5S','ION-TOT0.5S','AR^1.1B','EB-1.1B','ION-TOT1.1B']
SCCP_PCMC = ['AR^7.7J','ION-TOT7.7J','AR^5.1B','ION-TOT5.1B']
PR_PCMC = ['AR^0.35','EB-0.35','ION-TOT0.35']





        
        
        



#Commonly Used Diagnostic Settings
#electrodeloc	#YPR [29,44],[16,44] #SPR [0,107] 	#MSHC [0,12]
#waveformlocs 	#YPR [[16,29],[16,44],[16,64]]
#DOFWidth		#YPR R;16,Z;41   					#MSHC R;5,Z;10
#TrendLoc		#YPR H[0];R[29,44,64] 				#MSHC H[0,20];R[20]
#ThrustLoc		#YPR=74, stdESCT=76, smlESCT=48/54
#SheathROI		#YPR=[34,72]
#SourceWidth	#YPR=R[0.21]						#MSHC R[]

#Commonly Used Image Settings
#Crop YPR R[0.6];Z[1,4.5]   
#Crop #MSHC R[0.0,1.0];Z[0.5,2.5]

#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested IEDF/NEDF Variables.
IEDFVariables = PR_PCMC		#Requested iprofile_2d variables (no spaces)
NEDFVariables = []			#Requested nprofile_2d variables (no spaces)

#Requested movie1/movie_icp Variables.
IterVariables = ['E','S-E','PPOT','TE']		#Requested Movie_icp (iteration) Variables.		
PhaseVariables = Ar_Phase				#Requested Movie1 (phase) Variables. +['E','AR+']
electrodeloc = [29,44]						#Cell location of powered electrode [R,Z].
waveformlocs = [[16,29],[16,44],[16,64]]	#Cell locations of additional waveforms [R,Z].

#Various Diagnostic Settings.
phasecycles = 2							#Number of waveform phase cycles to be plotted. [number]
DoFWidth = 41							#PROES Depth of Field (symmetric on image plane) [cells]
ThrustLoc = 74							#Z-axis cell for thrust calculation  [cells]
SheathROI = [34,72]						#Sheath Region of Interest, (Start,End) [cells]
SourceWidth = [16]						#Source Dimension at ROI, leave empty for auto. [cells]

#Requested TECPLOT Variables and plotting locations.
Variables = Ar
MultiVar = []							#Additional variables plotted ontop of [Variables]
radialineouts = [29,44,64,74] 						#Radial 1D-Profiles to be plotted (fixed Z-mesh) --
heightlineouts = [0]						#Axial 1D-Profiles to be plotted (fixed R-mesh) |
TrendLocation = [] 						#Cell location For Trend Analysis [R,Z], ([] = min/max)


#Requested diagnostics and plotting routines.
savefig_convergence = False				#Requires movie_icp.pdt
savefig_plot2D = False					#Requires TECPLOT2D.PDT

savefig_monoprofiles = False			#Single-Variables; fixed height/radius
savefig_multiprofiles = False			#Multi-Variables; same folder
savefig_comparelineouts = False			#Multi-Variables; all folders
savefig_trendphaseaveraged = False		#Single-Variables; fixed cell location (or max/min)
savefig_trendphaseresolved = False		#Single-Variables; Phase-resolved data.
savefig_pulseprofiles = False			#Single-Variables; plotted against real-time axis

savefig_phaseresolve1D = False			#1D Phase Resolved Images
savefig_phaseresolve2D = False			#2D Phase Resolved Images
savefig_PROES = False					#Simulated PROES Diagnostic

savefig_IEDFangular = False				#2D images of angular IEDF; single folders.
savefig_IEDFtrends = False				#1D IEDF trends; all folders.
savefig_EEDF = False					#NO PLOTTING ROUTINE		#IN DEVELOPMENT#

#Write processed data to ASCII files.
write_ASCII = True						#All diagnostic output written to ASCII.


#Steady-State diagnostics terminal output toggles.
print_generaltrends = False				#Verbose Min/Max Trend Outputs.
print_Knudsennumber = False				#Print cell averaged Knudsen Number
print_totalpower = False				#Print all requested total powers
print_DCbias = False					#Print DC bias at electrodeloc
print_thrust = False					#Print neutral, ion and total thrust
print_sheath = False					#Print sheath width at electrodeloc


#Image plotting options.
image_extension = '.png'				#Extensions ('.png', '.jpg', '.eps')
image_aspectratio = [10,10]				#[x,y] in cm [Doesn't rotate dynamically]
image_radialcrop = [0.6]				#[R1,R2] in cm
image_axialcrop = [1,4]					#[Z1,Z2] in cm
image_cbarlimit = []					#[min,max] colourbar limits	

image_plotsymmetry = True				#Toggle radial symmetry
image_numericaxis = False				#### NOT IMPLIMENTED ####
image_contourplot = True				#Toggle contour Lines in images
image_plotgrid = False					#Plot major/minor gridlines on profiles
image_plotmesh = 'PR'					#### NOT IMPLIMENTED ####	('Auto','PR')
image_rotate = True						#Rotate image 90 degrees to the right.

image_normalize = False					#Normalize image/profiles to local max
image_logplot = False					#Plot ln(Data), against linear axis.
image_sheath = True						#Plot sheath width onto 2D images.


#Overrides the automatic image labelling.
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = ['NotImplimented']





#============================#












#####TODO#####

#For V 0.11.n:
#FIX CBARMINMAX FUNCTION!!! PROES NEEDS CROPPING AND RADIAL MINMAX IS WRONG.
#ADD if DOFWIDTH < LINEOUT LOCATION SKIP AND WARNING IN PROES
#rotate data at read-in and remove confusing [::-1] from diagnostics.
#Thrust diagnostic uses Rlineout function which employs symmetry <-- need to force symmetry!
#Functionalize thrust calculation with options for neutral/ion/pressure diff.
#Impliment image_numericaxis, try float(FolderNameTrimmer) as axis.
#Impliment numerical image_rotate, allow for 000,090,180,270.
#Introduce Dirlist creating function, using os.module (remove findtools)
#IEDF ASCII routine needs to save the IEDF energy axis. (current assumes 1ev - 250ev)
#Remove the need for ['E','AR+'] as default in phasevariables - use SheathData list?
#Variable Interpolator needs to work with phasedata - Take variables from batch?
#Convert EnumerateVariable Function to use header, not full rawdata.
#Use Headerlists to store the actual headers
#header_2Dlist.append(rawdata[0:header_2D])


#For V 1.0.0:
#SheathWidth function needs to be able to work axially and radially
#SheathWidth function needs to be able to deal with image rotations
#SheathWidth function requires automatic ROI calculation.
#Remove/simplify the phase-averaged sheath diagnostic. 
#Functionalise PROES images, Radial/Axial PROES profile collectors?
#Complete IEDF/NEDF section and Functionalise
#Add EEDF section and Functionalise.
#Clean up unused functions and ensure homogeneity.


#For Future:
#introduce seaborn into the program en-masse.
#introduce 'garbage collection' at the end of each diagnostic.
#Update README, include all diagnostics and examples.
#Create developer handbook describing functions.
#Python 3.x compatable.





#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()
Dirlist = list()
IEDFVariablelist = list()
Geometrylist = list()

Globalvarlist = list()
Globalnumvars = list()

#Create mesh_size lists and SI conversion
Isymlist = list()
R_mesh = list()
Z_mesh = list()
Raxis = list()
Zaxis = list()

Depth = list()
Radius = list()
Height = list()
dr = list()
dz = list()

VRFM,VRFM2 = list(),list()
FREQM,FREQM2 = list(),list()
FREQICP,IRFPOW = list(),list()
MAXFREQ,MINFREQ = list(),list()
IETRODEM = list()

#Create lists to store data
rawdata_2D = list()
rawdata_kin = list()
rawdata_phasemovie = list()
rawdata_itermovie = list()
rawdata_IEDF = list()
rawdata_mcs = list()

Data = list()					#Data[folder][Variable][Datapoint]
DataIEDF = list()				#Data[folder][Variable][Datapoint]
DataEEDF = list()				#Data[folder][Variable][Datapoint]
IterMovieData = list()			#ITERMovieData[folder][timestep][variable][datapoints]
PhaseMovieData = list()			#PhaseMovieData[folder][timestep][variable][datapoints]

Moviephaselist = list()			#'CYCL = n'
MovieIterlist = list()			#'ITER = n'
EEDF_TDlist = list()			#'???'

header_itermovie = list()
header_phasemovie = list()
header_IEDFlist = list()
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
print '                                                            v0.12.4 '
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
if True in [print_generaltrends,print_Knudsennumber,print_totalpower,print_DCbias,print_thrust]:
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
numfolders = 0

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


#Begin the retrieval of geometry from mesh and input files.
icpnam = filter(lambda x: 'icp.nam' in x, Dir)
icpout = filter(lambda x: 'icp.out' in x, Dir)
mesh = filter(lambda x: 'initmesh.out' in x, Dir)
TEC2D = filter(lambda x: 'TECPLOT2D.PDT' in x, Dir)

#Loop over all folders and retrieve mesh sizes and SI sizes.
for l in range(0,numfolders):

	#Attempt automated retrieval of mesh sizes.
	try:
		#Identify mesh size from TECPLOT2D file.
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
			print '#========================================================#'
			print 'Cannot extract mesh, please manually define mesh geometry.'
			print '#========================================================#'
			r_mesh = int(raw_input("Please Define R_mesh: "))
			z_mesh = int(raw_input("Please Define Z_mesh: "))
			print ''

			R_mesh.append(r_mesh)
			Z_mesh.append(z_mesh)
		#endif
	#endtry


				#MESH PLOTTING NOT WORKING#
#################################################################
	#Retrieve entire mesh for plotting if requested.
	if image_plotmesh == True:
		print '#================================================#'
		print 'Mesh Outline Plotting Does Not Currently Function.'
		print '#================================================#'
		print ''
		#Extract mesh data from initmesh.out
		mesh = open(mesh[l]).readlines()
	#endif

	#Inform Mesh Size for inputdeck plotting purposes.
	print Dirlist[l]
	print 'R_mesh:', R_mesh[-1], 'Z_mesh:', Z_mesh[-1]
	print ''
#################################################################


	#Attempt automated retrieval of SI conversion units.
	SImeshdata = open(icpnam[l]).readlines()

	#Retrieve useful input variables from icp.nam.
	try:
		NUMPHASE = int(filter(lambda x: x.isdigit(),filter(lambda x:'IMOVIE_FRAMES' in x,SImeshdata)[0]))
		NUMMETALS = int(filter(lambda x: x.isdigit(),filter(lambda x:'IMETALS' in x,SImeshdata)[0]))+1
		MATERIALS = filter(lambda x: 'CMETAL=' in x, SImeshdata)[0].split()[1:NUMMETALS]
	except:
		print 'ICP.NAM READIN ERROR, IGNORING MESH MATERIAL TYPES'
	#endtry

	#Input frequencies/voltages/powers   [FREQICP ONLY READS 10 CHARACTERS]
	try:
		VRFM.append(filter(lambda x: 'VRFM=' in x, SImeshdata)[0].split()[1:NUMMETALS])
		VRFM2.append(filter(lambda x: 'VRFM_2=' in x, SImeshdata)[0].split()[1:NUMMETALS])
		FREQM.append(filter(lambda x: 'FREQM=' in x, SImeshdata)[0].split()[1:NUMMETALS])
		FREQM2.append(filter(lambda x: 'FREQM_2=' in x, SImeshdata)[0].split()[1:NUMMETALS])
		FREQICP.append(float(filter(lambda x:'FREQ=' in x, SImeshdata)[0].strip(' \t\n\r,=FREQ')[0:10]))
		IRFPOW.append(float(filter(lambda x:'IRFPOW=' in x, SImeshdata)[0].strip(' \t\n\r,=IRFPOW')))
		IETRODEM.append(filter(lambda x:'IETRODEM=' in x, SImeshdata)[0].split()[1:NUMMETALS])
		for i in range(0,len(IETRODEM[l])): IETRODEM[l][i] = int(IETRODEM[l][i].strip(','))
	except:
		print 'ICP.NAM READIN ERROR, USING DEFAULT MATERIAL PROPERTIES'
		FREQM.append(13.56E6)
		FREQM2.append(13.56E6)
		FREQICP.append(13.56E6)
		VRFM.append(240.0)
		VRFM2.append(240.0)
		IRFPOW.append(100.0)
	#endtry

	#SI Conversion unit extraction.
	try:
		RADIUS = float(filter(lambda x:'RADIUS=' in x, SImeshdata)[0].strip(' \t\n\r,=RADIUS'))
		RADIUST = float(filter(lambda x:'RADIUST=' in x, SImeshdata)[0].strip(' \t\n\r,=RADIUST'))
		HEIGHT = float(filter(lambda x:'HEIGHT=' in x, SImeshdata)[0].strip(' \t\n\r,=HEIGHT'))
		HEIGHTT = float(filter(lambda x:'HEIGHTT=' in x, SImeshdata)[0].strip(' \t\n\r,=HEIGHTT'))
		DEPTH = float(filter(lambda x:'DEPTH=' in x, SImeshdata)[0].strip(' \t\n\r,=DEPTH'))
		SYM = float(filter(lambda x:'ISYM=' in x, SImeshdata)[0].strip(' \t\n\r,=ISYM'))
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
		print '#=================================================#'
		print 'icp.nam not found, please manually define variables'
		print '#=================================================#'
		radius = float(raw_input("Please Define SI radius: "))
		height = float(raw_input("Please Define SI height: "))
		depth = float(raw_input("Please Define SI depth: "))
		print ''

		Radius.append(radius)
		Height.append(height)
		Depth.append(depth)
		dr.append(Radius[-1]/(R_mesh[-1]-1))
		dz.append(Height[-1]/(Z_mesh[-1]-1))
	#endtry

	#clean up variables and assign required types.
	try:
#		for i in range(0,len(MATERIALS[l])): MATERIALS[l][i] = MATERIALS[i].strip(',\'') <-Broken
		VRFM[l] = float( VRFM[l][IETRODEM[l].index(1)].strip(',') )
		VRFM2[l] = float( VRFM2[l][IETRODEM[l].index(1)].strip(',') )
		FREQM[l] = float( FREQM[l][IETRODEM[l].index(1)].strip(',') )
		FREQM2[l] = float( FREQM2[l][IETRODEM[l].index(1)].strip(',') )

		MINFREQ.append( min([FREQM[l],FREQM2[l],FREQICP[l]]) )
		MAXFREQ.append( max([FREQM[l],FREQM2[l],FREQICP[l]]) )
	except:
		Material_Property_Conversion_Error=1
	#endtry
#endfor






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


#Takes directory list and data filename type (e.g. .png, .txt)
#Returns datalist of contents and length of datalist.
#rawdata, datalength = ExtractRawData(Dir,'.dat',l)
def ExtractRawData(Dirlist,NameString,ListIndex=l):
	try:
		DataFileDir = filter(lambda x: NameString in x, Dirlist)
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
def SDFileFormatConvertorHPEM(Rawdata,header,numvariables,offset=0,Zmesh=0,Rmesh=0):

	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Zmesh == 0:
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

	#Seperate total 1D array into 2D array with data for each variable.
	#Offset data by a certain number of variable 'chunks' if requested.
	for i in range(offset,numvariables):
		numstart = (Zmesh*Rmesh)*(i)
		numend = (Zmesh*Rmesh)*(i+1)
		CurrentFolderData.append(list(DataArray1D[numstart:numend]))
	#endfor

	return(CurrentFolderData)
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


#Reads 1D or 2D data from textfile in ASCII format.
#One input, filename string, returns data array.
def ReadDataFromFile(Filename,Dimension='1D'):
	datafile = open(Filename)
	OutputData = list()

	#Determine dimensionality of profile.
	if Dimension == '2D':
		#Read in 2D data from ASCII formatted file.	
		RawData = datafile.readlines()
		for m in range(0,len(RawData)):
			Row = RawData[m].split()
			for n in range(0,len(Row)):
				#Convert to float if possible.
				try: Row[n] = float(Row[n])
				except: Row[n] = Row[n]
			#endfor
			OutputData.append(Row)
		#endfor

	#Lowest dimention is scalar: ==> 1D array.
	elif Dimension == '1D':
		#Read in 1D data from ASCII formatted file.
		Row = datafile.readline().split()
		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData)
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
def Automovie(FolderDir,Output):

	#Break if movies not requested
	if DisableMovie == True: return()

	#Correct file extention on name.
	HomeDir = os.getcwd()
	Output = Output+'.mp4'
	Morph, FPS = 1, 24

	#Use ffmpeg to create the movies and save in relevent files.
	os.chdir(FolderDir)
	os.system("convert *.png -delay 1 -morph "+str(Morph)+" %05d.morph.jpg > /dev/null")
	os.system("ffmpeg -nostats -loglevel 0 -r "+str(FPS)+" -i %05d.morph.jpg "+Output)
	os.system("rm *.jpg")
	os.chdir(HomeDir)
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
	Ionizationlist = ['S-','SEB-']
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
			Variable = 'Bulk e$^-$ Source Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-E':
			Variable = 'Secondry e$^-$ Source Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'EB-ESORC':
			Variable = 'Secondry e$^-$ Relaxation Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'S-AR+':
			Variable = 'Bulk Ar+ Ionization Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-AR+':
			Variable = 'Secondry Ar+ Ionization Rate'
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
			VariableUnit = '[Torr]'
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
			VariableUnit = '[G]'
		elif variablelist[i] == 'JZ-NET':
			Variable = 'Axial Current Density'
			VariableUnit = '[mA cm$^{-2}$]'
		elif variablelist[i] == 'JR-NET':
			Variable = 'Radial Current Density'
			VariableUnit = '[mA cm$^{-2}$]'

		#Explicit Power Deposition.
		elif variablelist[i] == 'POW-TOT':
			Variable = 'Total RF-Power Deposited'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] == 'POW-RF':
			Variable = 'RF-Power Deposited'
			VariableUnit = '[Wm$^{-3}$]'
		elif variablelist[i] == 'POW-RF-E':
			Variable = 'RF-Power Deposited by e$^-$'
			VariableUnit = '[Wm$^{-3}$]'


		#Implicit Variables.
		elif IsStringInVariable(variablelist[i],Ionizationlist) == True:
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
		elif variablelist[i] in AtomicSet:
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

	#For ionization rates, convert from [cm3 s-1] to [m3 s-1]
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

	#For densities, convert from [cm-3] to [m-3]. (AtomicSet is defined in default parameters)
	if variable in AtomicSet:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E6
		#endfor
	#endif

	return(profile)
#enddef


def ManualPRMesh(Ax=plt.gca()):
	#Plot pocket rocket material dimensions.
	Ax.plot((1.32,1.32),   (-1.0,-0.21), 'w-', linewidth=2)
	Ax.plot((3.7,3.7),     (-1.0,-0.21), 'w-', linewidth=2)
	Ax.plot((1.32,1.32),   ( 1.0, 0.21), 'w-', linewidth=2)
	Ax.plot((3.7,3.7),     ( 1.0, 0.21), 'w-', linewidth=2)
	Ax.plot((1.32,3.7),    ( 0.21, 0.21), 'w-', linewidth=2)
	Ax.plot((1.32,3.7),    (-0.21,-0.21), 'w-', linewidth=2)

	#Alumina Dielectric
	Ax.plot((34.2*dz[l],73.8*dz[l]),  ( 0.21,  0.21), 'c-', linewidth=2)
	Ax.plot((34.2*dz[l],73.8*dz[l]),  (-0.21, -0.21), 'c-', linewidth=2)
	Ax.plot((34.2*dz[l],73.8*dz[l]),  ( 0.31,  0.31), 'c-', linewidth=2)
	Ax.plot((34.2*dz[l],73.8*dz[l]),  (-0.31, -0.31), 'c-', linewidth=2)

	#Powered Electrode
	Ax.plot((40*dz[l],50*dz[l]),  ( 0.31, 0.31), 'r-', linewidth=2)
	Ax.plot((40*dz[l],50*dz[l]),  (-0.31,-0.31), 'r-', linewidth=2)
	Ax.plot((40*dz[l],40*dz[l]),  ( 0.31, 0.60), 'r-', linewidth=2)
	Ax.plot((40*dz[l],40*dz[l]),  (-0.31,-0.60), 'r-', linewidth=2)
	Ax.plot((50*dz[l],50*dz[l]),  ( 0.31, 0.60), 'r-', linewidth=2)
	Ax.plot((50*dz[l],50*dz[l]),  (-0.31,-0.60), 'r-', linewidth=2)

	#Grounded electrodes
	Ax.plot((34*dz[l],34*dz[l]),  (-1.0,-0.21), 'w-', linewidth=2)
	Ax.plot((34*dz[l],34*dz[l]),  ( 1.0, 0.21), 'w-', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  (-1.0,-0.21), 'w-', linewidth=2)
	Ax.plot((74*dz[l],74*dz[l]),  ( 1.0, 0.21), 'w-', linewidth=2)
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

	#Kinetics data readin - NOT CURRENTLY USED
	if True == False:

		#Load data from TECPLOT_KIN file and unpack into 1D array.
		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
		rawdata_kin.append(rawdata)
	#endif


#===================##===================#
#===================##===================#

	#IEDF/NEDF file readin.
	if True in [savefig_IEDFangular,savefig_IEDFtrends]:

		#Define arguments and autorun conv_prof.exe if possible.
		IEDFVarArgs = ['1','1','1','1','1'] #### THIS IS HACKY, WON'T ALWAYS WORK ####
		args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ['0','0']
		DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt','iprofile_zones_tec1d.pdt','nprofile_zones_tec1d.pdt']
		try: AutoConvProfData('./conv_prof.exe',args,DirAdditions)
		except: print Dirlist[l]

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
	mpl.rcParams['image.cmap'] = 'jet'						#Select global colourmap 
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
	if image_rotate == True: 
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
	if image_rotate == 00:
		Z1,Z2 = Z1,Z2
	elif image_rotate == 90 or image_rotate == True:
		R1,Z1 = Z1,R1
		R2,Z2 = Z2,R2
	elif image_rotate == 180 or image_rotate == False:
		Z1,Z2 = Z2,Z1
	elif image_rotate == 270:
		R1,Z2 = Z2,R1
		R2,Z1 = Z1,R2
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
#Works for PROES images too, requires PROES='Axial' or 'Radial'.
#[Minimum,Maximum] = CbarMinMax(Image,PROES=False)
def CbarMinMax(Image,PROES=False):

	#Return user defined limits if specified.
	if len(image_cbarlimit) == 2:
		cropmin = image_cbarlimit[0]
		cropmax = image_cbarlimit[1]
		return([cropmin,cropmax])
	#endif

	#Ensure limits are in line with any requested mathematical constraints
	if image_logplot == True: Image = np.log(Image)
	if image_normalize == True: Image = normalize(Image)

	#Modify image to cropped region if a region is supplied.
	if any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:

		#Import global cell sizes and apply rotation for the maths stage.
		dR,dZ = dr[l],dz[l]
		if image_rotate == True: dR,dZ = dZ,dR
		#endif

		#Convert cropped SI region (CropExtent) to cell region (R1,R2,Z1,Z2).
		CropExtent = CropImage(Apply=False)		#CropExtent applies rotation internally
		R1 = int(CropExtent[0][0]/dR)
		R2 = int(CropExtent[0][1]/dR)
		Z1 = int(CropExtent[1][0]/dZ)
		Z2 = int(CropExtent[1][1]/dZ)
		#endif

		#Re-rotate back so that the image is cut in the correct order.
		#R1,R2 are radial limits of image, Z1,Z2 are axial limits.
		if image_rotate == True:
			R1,Z1 = Z1,R1
			R2,Z2 = Z2,R2
		#Replace negative R1 with zero, Images here have no symmetry.
		if R1 < 0: R1 = 0
		#endif

		#Crop the cell region to the desired region, axial first, then radial.
		if PROES == False:
			Image = Image[Z1:Z2]
			Image = np.asarray(Image).transpose()
			Image = Image[R1:R2]
			Image = np.asarray(Image).transpose()
		#Crop the cell region for PROES images, they only require one cropped axis.
		elif PROES == 'Axial':
			Image = Image[::-1][Z1:Z2]			#Reverse image, origin at top.
		elif PROES == 'Radial':
			Image = np.asarray(Image).transpose()
			Image = Image[R1:R2]
			Image = np.asarray(Image).transpose()
		#endif
	#endif

	#Flatten image and obtain min/max in region, use full image if no cropping.
	flatimage = [item for sublist in Image for item in sublist]
	cropmin,cropmax = min(flatimage),max(flatimage)

	#Return cropped values as list [min,max], as required by colourbar.
	return([cropmin,cropmax])
#enddef



#=========================#
#=========================#



#Applies plt.options to current figure based on user input.
#Returns nothing, current image is required, use figure().
#ImageOptions(plt.gca(),Xlabel,Ylabel,Title,Legend,Crop=False)
def ImageOptions(ax=plt.gca(),Xlabel='',Ylabel='',Title='',Legend=[],Crop=True):

	#Apply user overrides to plots. ## Experimental ##
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
	elif image_plotmesh == 'PR' and Crop == True:	
		ManualPRMesh(ax)
	#endif

	#Crop image dimensions, use provided dimensions or default if not provided.
	if isinstance(Crop, (list, np.ndarray) ) == True:
		CropImage(ax,Crop)
	elif any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
		if Crop == True:
			CropImage(ax)
		#endif
	#endif

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

	#Colourbar plotting details
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)
	#Set number of ticks, label location and scientific notation.
	tick_locator = ticker.MaxNLocator(nbins=Ticks)
	cbar.locator = tick_locator
	cbar.set_label(Label, rotation=270,labelpad=30,fontsize=24)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.update_ticks()
	#Size of font
	cbar.ax.yaxis.offsetText.set(size=18)
	yticks(fontsize=18)

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
def GenerateAxis(Orientation,Isym=Isymlist[l],phasepoints=range(0,180)):
	
	#Extract number of phase datapoints and create axis list.
	phasepoints = len(phasepoints)
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
		for i in range(0,phasecycles*phasepoints):
			axis.append(  (np.pi*(i*2)/phasepoints)/(2*np.pi)  )
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
def DataExtent(folder=l,aspectratio=image_aspectratio):

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
	elif image_rotate == False:
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
def ImagePlotter2D(Image,extent,aspectratio=image_aspectratio,variable='N/A',fig=111,ax=111):

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
	if image_numericaxis == True:
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

	#WHY Rmesh=R_mesh[l],Zmesh=Z_mesh[l] NOT WORKING FOR MULTIPLE FOLDERS???
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

	#Obtain applied voltage waveform, Refresh list between folders if needed.
	for j in range(0,phasecycles):
		for i in range(0,len(PhaseData)):
			VoltageWaveform.append(PlotAxialProfile(PhaseData[i],PPOT,'PPOT',ZLoc)[RLoc])
		#endfor
	#endfor

	#Calculate time averaged waveform bias, i.e. waveform symmetry.
	for m in range(0,len(VoltageWaveform)):
		WaveformBias.append(sum(VoltageWaveform)/len(VoltageWaveform))
	#endfor
	
	return(VoltageWaveform,WaveformBias)
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



#Calculates sheath width assuming Child-Langmuir conditions.
#Calculation Methods: 'AbsDensity', 'IntDensity'
#Takes current folder, current axis, movie1 Phase and sheath calc method.
#Returns array of sheath distances from origin and can plot this if requested.
#Sx = SheathThickness(folder=l,Phase=moviephaselist[k])
def SheathThickness(folder=l,ax=plt.gca(),Phase='NaN',Ne=list(),Ni=list()):
	#Initiate required lists and set sheath method.
	SheathMethod=GlobSheathMethod
	Sx,SymSx = list(),list()	

	#Obtain current folder ion and electron densities if not already supplied.
	#Default to 2D data format.
	if Phase == 'NaN' and len(Ne) == 0:
		IONproc = VariableEnumerator(['AR+'],rawdata_2D[folder],header_2Dlist[folder])[0][0]
		Eproc = VariableEnumerator(['E'],rawdata_2D[folder],header_2Dlist[folder])[0][0]
		Ne,Ni = Data[folder][Eproc], Data[folder][IONproc]
	#If phase is supplied, use phase data format.
	elif Phase != 'NaN' and len(Ne) == 0:
		IONproc = VariableEnumerator(['AR+'],rawdata_phasemovie[folder],header_phasemovie[folder])[0][0]
		Eproc = VariableEnumerator(['E'],rawdata_phasemovie[folder],header_phasemovie[folder])[0][0]
		IONproc,Eproc = IONproc-2, Eproc-2		#Skip R,Z data inputs in phase data.
		Ne,Ni = PhaseMovieData[folder][Phase][Eproc], PhaseMovieData[folder][Phase][IONproc]
	#endif
	#Extract 2D image for further processing.
	Ne,Ni = ImageExtractor2D(Ne),ImageExtractor2D(Ni)

	#=======#

	### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###
	Orientation = 'Axial'
	#Determine electrode location.
	if Orientation == 'Radial': loc = electrodeloc[0]
	elif Orientation == 'Axial': loc = electrodeloc[1]
	### CURRENTLY ONLY AXIAL METHOD IS EMPLOYED ###

	#Determine sheath edge through integration of charge.
	if SheathMethod == 'IntDensity':
		#Sheath extension: integral_(0-R) ne dR == integral_(0-R) ni dR (Gibson 2015)
		for j in range(0,len(Ni)):
			#Refresh sums after every radial profile.
			Ni_sum,Ne_sum = 0.0,0.0
			for i in range(0,len(Ni[j])):
				#Sum density outward radially for ions and electrons.
				anti_i = len(Ni[j])-i-1
				Ni_sum += Ni[j][i]
				Ne_sum += Ne[j][i]

				#If ion sum is greater than electron, sheath has begun.
				if Ni_sum/Ne_sum >= 1.0: 
					Sx.append(i*dr[l])
					break
				#If no sheath found, append 'NaN' to avoid plotting.
				if i == (len(Ni[j])-1):
					Sx.append('NaN')
				#endif
			#endfor
		#endfor

	#Determine sheath edge by 'instantaneous' charge
	elif SheathMethod == 'AbsDensity':
		#Sheath extension: ni @R >= ne @R, simplified model.
		for j in range(0,len(Ni)):
			for i in range(0,len(Ni[j])):
				#Sheath starts when ion density exceeds electron density.
				if Ni[j][i]/Ne[j][i] >= 1.0:		###SUPPRESS OUTPUT?###
					Sx.append(i*dr[l])
					break
				#If no sheath found, append 'NaN' to avoid plotting.
				if i == (len(Ni[j])-1):
					Sx.append('NaN')
				#endif
			#endfor
		#endfor
	#endif

	#=======#

	#Print sheath characteristics if requested.
	if print_sheath == True:
		#Obtain SheathWidth at electrodeloc
		try: SheathWidth = round(Sx[loc],3)
		except: SheathWidth = 0.0
		print 'Simulation:', Dirlist[folder]
		print 'Sheath Location:',SheathWidth*10, 'mm'
		print 'Sheath Extent:',((SourceWidth[0]*dr[l])-SheathWidth)*10, 'mm'
	#endif

	#=======#

	#THIS SHOULD PROBABLY BE A SEPERATE FUNCTION
	#THAT USES SHEATHTHICKNESS AS AN INPUT

	#Generate Axis and Remove NaN's from data, 
	Axis=GenerateAxis(Orientation,Isym=Isymlist[folder])
	for i in range(0,len(Sx)): 
		if Sx[i] == 'NaN': Sx[i],Axis[i] = np.nan,np.nan
		#endif
	#endfor

	#Create symmetric sheath boundary
	for i in range(0,len(Sx)): SymSx.append(-Sx[i])
	#endfor

	#Plot sheath characteristics if requested.
	if image_sheath == True:
		ax.plot(Axis,Sx, 'w--', lw=2)
		ax.plot(Axis,SymSx, 'w--', lw=2)
	#endif

	#Return sheath expansion
	return(Sx)
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
			Sx = SheathThickness(folder=l,ax=ax)

			#Define image beautification variables.
			if image_rotate == True:
				Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
			elif image_rotate == False:
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				plt.gca().invert_yaxis()
			#endif

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = '2D Steady State Plot of '+variablelist[k]+' for \n'+Dirlist[l][2:-1]
			ImageOptions(ax,Xlabel,Ylabel,Title)

			#Add Colourbar (Axis, Label, Bins)
			label,bins = VariableLabelMaker(variablelist),5
			cax = Colourbar(ax,label[k],bins,Lim=CbarMinMax(Image))

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(Dir2Dplots, '2Dplots_Data')
				WriteDataToFile(Image, DirWrite+variablelist[k])
				if k == len(processlist)-1: WriteDataToFile(Sx, DirWrite+'Sx-EXT')
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
			Xaxis.append(filter(lambda x: x.isdigit(), MovieIterlist[l][i]))
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
			for k in range(0,len(MovieIterlist[l])):

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(IterMovieData[l][k][processlist[i]],variablelist[i])
				#Take Max value of image for general convergence trend.
				ConvergenceTrends[-1].append( sum(Image.flatten())/len(Image.flatten()) )

				#Generate and rotate figure as requested.
				extent,aspectratio = DataExtent(l)
				fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio,variablelist[i])
				#Add sheath thickness to figure if requested.
				Sx = SheathThickness(folder=l,Ax=ax)

				#Define image axis labels.
				if image_rotate == True:
					Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
				elif image_rotate == False:
					Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
					plt.gca().invert_yaxis()
				#endif

				#Image plotting details.
				Title = str(MovieIterlist[l][k])
				ImageOptions(ax,Xlabel,Ylabel,Title)

				#Add Colourbar (Axis, Label, Bins)
				label,bins = VariableLabelMaker(variablelist),5
				cax = Colourbar(ax,label[i],bins,Lim=CbarMinMax(Image))

				#Save to seperate folders inside simulation folder.
				num1,num2,num3 = k % 10, k/10 % 10, k/100 % 10
				Number = str(num3)+str(num2)+str(num1)
				savefig(DirMovieplots+variablelist[i]+'_'+Number+ext)
				plt.close('all')
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+variablelist[i])
		#endfor


		#=================#


		#Plot a convergence check for all variables in each folder.
		Legend = VariableLabelMaker(variablelist)
		fig, ax = plt.subplots(1, figsize=(10,10))

		#Normalize and plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(ConvergenceTrends)):
			ConvergenceTrends[i] = Normalize(ConvergenceTrends[i])[0]
			ax.plot(Xaxis,ConvergenceTrends[i], lw=2)
		#endfor

		#Image plotting details.
		Title = 'Convergence of '+str(variablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation Iteration','Normalized Mesh-Average Value'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(0,1.01+(len(Legend)*0.05))

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
				ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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
				ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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
				ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
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
				ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
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
					ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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
					ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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

		#Create processlist for each folder as required.
		processlist,variablelist = VariableEnumerator(Variables,rawdata_itermovie[l],header_itermovie[l])
		#Skip over the R and Z processes as they are not saved properly in iterdata.
		for i in range(0,len(processlist)):
			processlist[i] = processlist[i]-2
		#endfor

		#Create list and x-axis for convergence trend plotting.
		#DtActual is approximate, exact dt per 'iteration' is not known.
		DtActual = (1/FREQM[l])*100					#S	(~8 microseconds @ 13.56MHz)
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
			Xlabel,Ylabel = 'Simulation time [S]',VariableLabelMaker(variablelist)[i]
			ImageOptions(ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

			#Save figure.
			savefig(DirPulse+FolderNameTrimmer(Dirlist[l])+'_'+variablelist[i]+ext)
			plt.close('all')

			#Write data to ASCII files if requested.
			if write_ASCII == True:
				DirWrite = CreateNewFolder(DirPulse, 'Pulse_Data')
				WriteDataToFile(Xaxis, DirWrite+variablelist[i], 'w')
				WriteDataToFile(['\n']+PulseProfile, DirWrite+variablelist[i], 'a')
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
		Xlabel,Ylabel = 'Simulation time [S]','Normalized Mesh-Average Value'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(0,1.01+(len(Legend)*0.05))

		#Save figure.
		savefig(DirPulse+'Normalized_'+FolderNameTrimmer(Dirlist[l])+ext)
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
			Image = ImageExtractor2D(DataIEDF[l][processlist[i]],Rmesh=EDFangle,Zmesh=EDFbins)

			#Flatten angular distribution across all angles to produce energy distribution.
			for j in range(0,len(Image)): EDFprofile.append(sum(Image[j]))

			#Transpose Image for plotting and reverse both lists to align with other data.
			Image, EDFprofile = Image[::-1].transpose(), EDFprofile[::-1]

			#Plot the angular distribution and EDF of the required species.
			fig,ax = figure([11,9], 2, shareX=True)

			Title = Dirlist[l][2::]+'\n'+variablelist[i]+' Angular Energy Distribution Function'
			Extent=[0,len(Image[0]), -len(Image)/2,len(Image)/2]
			fig.suptitle(Title, y=0.995, fontsize=16)

			#Angular Figure
			im = ax[0].imshow(Image,extent=Extent)
			ImageOptions(ax[0],Ylabel='Angular Dispersion [$\\theta^{\circ}$]', Crop=False)						
			cax = Colourbar(ax[0],variablelist[i]+' EDF($\\theta$)',5)

			#Integrated IEDF figure
			ax[1].plot(EDFprofile, lw=2)
			Xlabel,Ylabel = 'Energy [eV]',variablelist[i]+' EDF \n [$\\theta$ Integrated]'
			ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

			plt.tight_layout()
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
		Legendlist,EDFprofiles = list(),list()
		Mode_eV,Mean_eV,Max_eV = list(),list(),list()
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

			#Plot current variable profile to figure for each simulation folder.
			ax.plot(EDFprofile, lw=2)


			#Perform a trend analysis on current folder variable i IEDF
			#Average energy analysis: Returns mean/mode energies from IEDF.
			mean = (max(EDFprofile)+min(EDFprofile))/2
			meanindex = (np.abs(EDFprofile[2::]-mean)).argmin()
			Mean_eV.append( EDFprofile.index(EDFprofile[meanindex]) ) 
			Mode_eV.append( EDFprofile.index(max(EDFprofile)) )

			#Maximum energy analysis: Returns maximum energy below a set threshold.
			threshold = 0.01*max(EDFprofile)
			for j in range(0,len(EDFprofile)):
				if EDFprofile[j] >= threshold: tempMax_eV = j
			#endfor
			Max_eV.append(tempMax_eV)

			#Particle energy variance analysis: Returns FWHM of energy distribution.
			#Take mean and draw line at y = mean 
			#Calculate where y = mean intercepts EDFprofile
			#If only one intercept, first intercept is x = 0
			#Integrate EDFprofile indices between intercepts 
			#Save in 1D array, can be used to get energy spread percentage.
		#endfor


		#Write data to ASCII format datafile if requested.
		if write_ASCII == True:
			if i == 0:
				DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
				DirASCIIIEDF = CreateNewFolder(DirASCII,'IEDF_Data')
			#endif
			WriteDataToFile(Legendlist+['\n']+Mode_eV, DirASCIIIEDF+variablelist[i]+'_Mode', 'w')
			WriteDataToFile(Legendlist+['\n']+Mean_eV, DirASCIIIEDF+variablelist[i]+'_Mean', 'w')
			WriteDataToFile(Legendlist+['\n']+Max_eV, DirASCIIIEDF+variablelist[i]+'_Max', 'w')
		#endif

		##IEDF PROFILES##
		#===============#

		#Apply image options to IEDF plot generated in the above loop.
		Title = Dirlist[l][2::]+'\n'+variablelist[i]+' Angular Energy Distribution Function Profiles'
		Xlabel,Ylabel = 'Energy [eV]',variablelist[i]+' EDF [$\\theta$ Integrated]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)
#		ax.set_xlim(0,100)

		plt.savefig(DirIEDFTrends+variablelist[i]+'_EDFprofiles'+ext)
		plt.close('all')


		#ENERGY ANALYSIS#
		#===============#

		#Plot average energy analysis profiles against simulation folder names.
		fig,ax = figure()
		TrendPlotter(ax,Mean_eV,Legendlist,NormFactor=0)
		TrendPlotter(ax,Mode_eV,Legendlist,NormFactor=0)
#		TrendPlotter(ax,Max_eV,Legendlist,NormFactor=0)

		Title = Dirlist[l][2::]+'\n'+'Average '+variablelist[i]+' Energies'
		Legend = ['EDF Mean Energy','EDF Mode Energy','EDF Max Energy']
		Xlabel,Ylabel = 'Varied Property',variablelist[i]+' Energy [eV]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)

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
			ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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
			ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

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
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

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
		if Variables[i] in ['POW-ALL','POW-TOT','POW-ICP','POW-RF','POW-RF-E']:
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
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

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
	ImageOptions(ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

	plt.savefig(DirTrends+'Power Deposition Comparison'+ext)
	plt.close('all')
#endif



#====================================================================#
				  	#ION/NEUTRAL THRUST ANALYSIS#
#====================================================================#


if savefig_trendphaseaveraged == True or print_thrust == True:

	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initiate lists required for storing data.
	NeutralThrustlist,IonThrustlist,Thrustlist = list(),list(),list()
	Xaxis = list()

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

			#Add neutral thrust and Isp to arrays (save dummy variables not calculated)
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
			#Correct extracted profiles to agree with this direction.
			Pressure,PressureDown = Pressure[0:64][::-1],PressureDown[0:64][::-1]
			NeutralVelocity,NeutralAxialFlux = NeutralVelocity[0:64][::-1],NeutralAxialFlux[0:64][::-1]
			IonVelocity,IonAxialFlux = IonVelocity[0:64][::-1],IonAxialFlux[0:64][::-1]
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
				#Ensure pressure index aligns with radial index for correct cell area.
				if Pressure[i] > 0.0:
#					DiffPressure = (Pressure[i]-0.85)*133.33				#N/m^2
					DiffPressure = (Pressure[i]-PressureDown[i])*133.33		#N/m^2
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
			if len(IonIsp) == 0: IonIsp.append(1E-30)
			if len(NeutralIsp) == 0: NeutralIsp.append(1E-30)

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
		WriteDataToFile(Thrustlist, DirASCII+'Thrust_Trends','a')
	#endif

	#Plot requested thrusts, neutral and ion, to 1st and 2nd y-axes.
	fig,ax1 = figure(image_aspectratio,1)
#	ax2 = ax1.twinx()
	P1 = TrendPlotter(ax1,Thrustlist,Xaxis,Marker='ko-',NormFactor=0)
	P2 = TrendPlotter(ax1,NeutralThrustlist,Xaxis,Marker='r^-',NormFactor=0)
#	P3 = TrendPlotter(ax2,IonThrustlist,Xaxis,Marker='bs-',NormFactor=0)
	Pn = P1+P2#+P3

	#Apply image options and save figure.
	Title = 'Thrust at Z='+str(round(ThrustLoc*dz[0],2))+'cm with varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Net Thrust [mN]'
	ImageOptions(ax1,Xlabel,Ylabel,Title,Crop=False)
#	ImageOptions(ax2,Ylabel='Ion Thrust [mN]',Crop=False)
	ax1.legend(Pn, ['Total Thrust','Neutral Component','Ion Component'], fontsize=18, frameon=False)

	plt.savefig(DirTrends+'Thrust Trends'+ext)
	plt.close('all')
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
	#loc = electrodeloc[0]		#Radial
	loc = electrodeloc[1] 		#Axial

	#For all selected simulations, obtain Xaxis, sheath value and save to array.
	for l in range(0,numfolders):
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Obtain sheath thickness array for current folder 
		Sx = SheathThickness(folder=l)

		#Extract maximum sheath thickness from region of interest and width at electrodeloc [Loc].
		#If this fails, provide null point for sheath thickness.
		try: 
			SxMaxExtent.append( ((SourceWidth[0]*dr[l])-max(Sx[SheathROI[0]:SheathROI[1]]))*10 )
			SxLocExtent.append( ((SourceWidth[0]*dr[l])-Sx[loc])*10 )
		except:
			SxMaxExtent.append( 'NaN' )
			SxLocExtent.append( 'NaN' )
		#endtry
	#endfor

	#===============================#

	#Generate figure and plot trends.	
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(ax,SxLocExtent,Xaxis,NormFactor=0)

	#Apply image options and axis labels.
	Title = 'Maximum Sheath Extension With Varying '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Extension [mm]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Legend=[],Crop=False)

	plt.savefig(DirTrends+'Sheath Extension (Phase-Averaged)'+ext)
	plt.close('all')
#endif



#====================================================================#
				  		#KNUDSEN NUMBER ANALYSIS#
#====================================================================#


#Only perform if a neutralspecies is included within the atomic set.
if bool(set(NeutSpecies).intersection(Variables)) == True:
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

			#Create extract data for the neutral flux and neutral velocity.
			processlist,Variablelist = VariableEnumerator(NeutSpecies,rawdata_2D[l],header_2Dlist[l])

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

			#Create new folder to keep 2D output plots.
			Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')

			#Display average Knudsen number to terminal if requested.
			KnudsenAverage.append( sum(Image)/(len(Image[0])*len(Image)) )
			if print_Knudsennumber == True:
				print Dirlist[l]
				print 'Average Knudsen Number:', KnudsenAverage[l]
			#endif

			#Label and save the 2D Plots.
			extent,aspectratio = DataExtent(l)
			fig,ax,im,Image = ImagePlotter2D(Image,extent,aspectratio)
			#Add sheath thickness to figure if requested.
			Sx = SheathThickness(folder=l,ax=ax)

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = 'Knudsen Number Image for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
			ImageOptions(ax,Xlabel,Ylabel,Title)
			plt.gca().invert_yaxis()

			#Add Colourbar (Axis, Label, Bins)
			label,bins = 'Knudsen Number',5
			cax = Colourbar(ax,label,bins,Lim=CbarMinMax(Image))

			#Save Figure
			plt.savefig(Dir2Dplots+'KnudsenNumber'+ext)
			plt.close('all')
		#endfor

		#Plot a comparison of all average Knudsen numbers.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(ax,KnudsenAverage,Xaxis,NormFactor=0)

		#Image plotting details.
		Title = 'Average Knudsen Number with Changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Average Knudsen Number'
		ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

		#Save figure.
		plt.savefig(DirTrends+'KnudsenNumber Comparison'+ext)
		plt.close('all')
	#endif
#endif

#===============================#

if any([savefig_trendphaseaveraged, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, print_sheath]) == True:
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
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'1DPhase/')
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get PPOT)
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables+['E','AR+'])
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
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)):
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['rf self-bias: '+str(round(ElectrodeBias[0],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)

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
			DirMovieplots = CreateNewFolder(DirPhaseResolved,varlist[i]+'_1DPhaseResolved/')

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
						PhaseResolvedlineout = PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])
						lineoutstring = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm \n'
						Xlabel = 'Axial Distance Z [cm]'
					elif LineoutsOrientation[k] == 'Radial':
						RlineoutLoc,axis = Lineouts[k],Raxis
						PhaseResolvedlineout = PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],RlineoutLoc,R_mesh[l],Isymlist[l])
						lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm \n'
						Xlabel = 'Radial Distance R [cm]'
					#endif

					#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
					fig,ax = figure(image_aspectratio,2)
					Ylabel = VariableLabelMaker(varlist)
					fig.suptitle('Phase-Resolved '+varlist[i]+' for '+VariedValuelist[l]+lineoutstring+str(Phaselist[j]), y=0.97, fontsize=16)

					#Plot profile and apply image options.
					ax[0].plot(axis, PhaseResolvedlineout[::-1], lw=2)
					ImageOptions(ax[0],Xlabel,Ylabel[i],Crop=False)
					ax[0].set_ylim(VariableMin,VariableMax*1.02)

					#Plot waveform and apply image options.
					ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
					ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
					Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
					ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

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
						WriteDataToFile(PhaseResolvedlineout[::-1], SaveString)
					#endif
				#endfor

				#Create .mp4 movie from completed images.
				Prefix = FolderNameTrimmer(Dirlist[l])+'_'+NameString
				Automovie(Dir1DProfiles,Prefix)
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

		#Extract Data[Phase][Variable][R,Z] and update trend axis.
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables+['E','AR+'])
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
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)): 
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['rf self-bias: '+str(round(ElectrodeBias[0],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)

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
				MinLim.append( CbarMinMax(Image)[0] )
				MaxLim.append( CbarMinMax(Image)[1] )
			#endfor
			Limits = [min(MinLim),max(MaxLim)]

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for j in range(0,len(Phaselist)):

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(PhaseData[j][proclist[i]],varlist[i])
				#Extract Ni and Ne variables for sheath processing.
				Ne = PhaseData[j][proclist[varlist.index('E')]]
				Ni = PhaseData[j][proclist[varlist.index('AR+')]]

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
				SheathThickness(folder=l,Ax=ax[0],Phase=j,Ne=Ne,Ni=Ni)
				ImageOptions(ax[0],Xlabel,Ylabel,Crop=True)
				#Add Colourbar (Axis, Label, Bins)
				Ylabel = VariableLabelMaker(varlist)
				cax = Colourbar(ax[0],Ylabel[i],5,Lim=Limits)

				#Plot waveform and apply image options.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)
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
					Cycle = str( Phaselist[j].replace(" ", "") )
					SaveString = DirASCIIPhase+varlist[i]+'_'+Cycle
					WriteDataToFile(Image, SaveString)
				#endif
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+varlist[i])
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

		#Create processlist for each folder as required. (Always get PPOT)
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables+['E','AR+'])
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Extract waveforms from desired electrode locations.
		PPOT = proclist[varlist.index('PPOT')]
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseData,PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseData,PPOT)

		#Select axial or radial electrode location (loc) and create axis.
		Orientation = 'Axial'		#### SET TO AXIAL BY DEFAULT ###
		if Orientation == 'Radial': loc = electrodeloc[0]
		elif Orientation == 'Axial': loc = electrodeloc[1]
		Phaseaxis = GenerateAxis('Phase',Isym=Isymlist[l])

		#=============#

		SxLoc = list()
		#For all phases, process data and record for plotting.
		for k in range(0,len(Phaselist)):
			#Extract Ni and Ne variables for sheath processing.
			Ne = PhaseData[k][proclist[varlist.index('E')]]
			Ni = PhaseData[k][proclist[varlist.index('AR+')]]

			#Extract sheath width and record sheath width at electrodeloc
			Sx = SheathThickness(folder=l,Phase=Phaselist[k],Ne=Ne,Ni=Ni)
			for j in range(0,len(Sx)): 
				try: Sx[j] = ((SourceWidth[0]*dr[l])-Sx[j])*10	#Convert to mm
				except: Sx[j] = 0.0								#'NaN' = 0.0
			#endfor
			SxLoc.append(Sx[loc])
		#endfor
		SxDynRangeTrend.append(max(SxLoc)-min(SxLoc))		#Dynamic Range
		SxMeanExtTrend.append(sum(SxLoc)/len(SxLoc))		#Mean Extension
		SxMaxExtTrend.append(max(SxLoc))					#Max Extension

		#Calculate phase-averaged (mean) sheath velocity.
		#Assumes one sheath collapse and one sheath expansion per rf-cycle.
		RFPeriod = 1.0/FREQM[l]											#[s]
		SheathExtent = (sum(SxLoc)/len(SxLoc))/1E6						#[km]
		MeanSheathVelocity = (2*SheathExtent)/RFPeriod					#[km/s]
		SxMeanVelTrend.append( MeanSheathVelocity )						#[km/s]		

		#Calculate maximum instantaneous sheath velocity.
		#Assumes sheath collapse velocity > sheath expansion velocity
		Collapsed,CollapsedPhase = min(SxLoc),SxLoc.index(min(SxLoc))
		Extended,ExtendedPhase = max(SxLoc),SxLoc.index(max(SxLoc))

		SheathExtension = (Extended-Collapsed)/1000.0  					#[m]
		PhaseResolution = 1.0/(FREQM[l]*len(Phaseaxis))					#[s]
		SheathTime = (ExtendedPhase-CollapsedPhase)*PhaseResolution		#[s]
		try:
			MaxSheathVelocity = SheathExtension/SheathTime				#[m/s]
			SxMaxVelTrend.append(MaxSheathVelocity/1000.0)				#[km/s]
		except:
			SxMaxVelTrend.append(0.0)									#[km/s]
		#endtry

		#=============#

		#Increase sheath extension by required number of phasecycles.
		for m in range(1,phasecycles):
			for n in range(0,len(SxLoc)):
				SxLoc.append(SxLoc[n])
			#endfor
		#endfor

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
		fig,ax = figure(image_aspectratio,2)
		ax[0].plot(Phaseaxis,SxLoc, lw=2)
		Ylabel = 'Sheath Extension [mm]'
		ImageOptions(ax[0],Ylabel=Ylabel,Crop=False)

		#Plot Waveform onto Temporally collapsed PROES.
		ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
		ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

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
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Max Sheath Extension Trends'+ext)
	plt.close('all')

	#Plot mean sheath extension trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMeanExtTrend,VariedValuelist,NormFactor=0)
	Title='Mean Sheath Extension W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Mean Sheath Extension [mm]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Mean Sheath Extension Trends'+ext)
	plt.close('all')

	#Plot sheath dynamic range (max-min extension) trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxDynRangeTrend,VariedValuelist,NormFactor=0)
	Title='Sheath Dynamic Range W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Dynamic Range [mm]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Sheath Dynamic Range Trends'+ext)
	plt.close('all')

	#Plot maximum sheath velocity trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMaxVelTrend,VariedValuelist,NormFactor=0)
	Title='Maximum Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

	plt.savefig(DirSheath+'Max Sheath Velocity Trends'+ext)
	plt.close('all')

	#Plot mean sheath velocity trend for all folders
	fig,ax = figure(image_aspectratio)
	TrendPlotter(ax,SxMeanVelTrend,VariedValuelist,NormFactor=0)
	Title='Phase-averaged Sheath Velocity W.R.T '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Sheath Velocity [km/s]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Crop=False)

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

		#Create processlist for each folder as required. (Always get PPOT)
		PhaseData,Phaselist,proclist,varlist = ExtractPhaseData(folder=l,Variables=PhaseVariables+['E','AR+'])
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
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseData,PPOT)

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)

		ax.plot(Phaseaxis,ElectrodeWaveform, lw=2)
		for j in range(0,len(waveformlocs)): 
			ax.plot(Phaseaxis,VoltageWaveforms[j], lw=2)
			#ax.plot(Phaseaxis,WaveformBiases[j], 'k--', lw=2)
		#endfor
		#ax.plot(Phaseaxis,ElectrodeBias, 'k--', lw=2)

		Title = 'Phase-Resolved Voltage Waveforms for '+FolderNameTrimmer(Dirlist[l])
		Legend = ['rf self-bias: '+str(round(ElectrodeBias[0],2))+'V']
		Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)

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
				PROES = list()

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Phaselist)):
					#Refresh lists between each phasecycle.
					IntegratedDoFArray,DoFArrays = list(),list()

					#Collect each profile for stitching into a PROES image if required.
					if DoFWidth > 0:
						#Determine range of lineouts within the depth of field.
						DOFRegion = [(Lineouts[k]-DoFWidth),(Lineouts[k]+DoFWidth)]
						#Collect lineouts from DOF region and transpose to allow easy integration.
						for LineoutLoc in range(DOFRegion[0],DOFRegion[1]):
							if LineoutsOrientation[k] == 'Radial': DoFArrays.append(PlotRadialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc))
							elif LineoutsOrientation[k] == 'Axial': DoFArrays.append(PlotAxialProfile(PhaseData[j],proclist[i],varlist[i],LineoutLoc)[::-1])
							#endif
						#endfor
						DoFArrays = np.asarray(DoFArrays).transpose().tolist()

						#Integrate DoF lineouts to form a single phase point PROES lineout.
						for m in range(0,len(DoFArrays)):
							IntegratedDoFArray.append( sum(DoFArrays[m])/(DoFWidth*2+1) )
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
				#endfor

				#Increase PROES image by required number of phasecycles.
				for m in range(1,phasecycles):
					for n in range(0,len(PROES)):
						PROES.append(PROES[n])
					#endfor
				#endfor


				#Create figure and rotate PROES such that phaseaxis aligns with waveform.
				fig,ax = figure(image_aspectratio,2)
				PROES = ndimage.rotate(PROES, 90)
				#Choose correct axial or radial distance axis and create associated folder.
				x1,x2 = Phaseaxis[0],Phaseaxis[-1]
				if LineoutsOrientation[k] == 'Axial':
					lineoutstring = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm'
					NameString = varlist[i]+'_'+lineoutstring[2::]
					Ylabel = 'Axial Distance Z [cm]'
					Crop = [image_axialcrop[::-1],image_radialcrop] #Reversed accounting for rotation.
					y1,y2 = Zaxis[-1],Zaxis[0]						#Reversed accounting for top origin.
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 1:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					Crop = [image_radialcrop,image_axialcrop]
					y1,y2 = Raxis[-1],-Raxis[-1]
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 0:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = varlist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					Crop = [image_radialcrop,image_axialcrop]
					y1,y2 = Raxis[-1],0
				#endif
				DirPROESloc = CreateNewFolder(DirPROES,lineoutstring[3::])

				#Create PROES image along line of sight with phase-locked waveform.
				fig.suptitle( 'Simulated '+varlist[i]+' PROES for '+VariedValuelist[l]+lineoutstring+'\n DoF = '+str(round(DoFWidth*dz[l],2))+' cm', y=0.95, fontsize=18)
				im = ax[0].contour(PROES,extent=[x1,x2,y1,y2],origin='lower',aspect='auto')
				im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='bottom',aspect='auto')
				ImageOptions(ax[0],Xlabel='',Ylabel=Ylabel,Crop=Crop)
				ax[0].set_xticks([])
				ax[0].set_xlim(x1,x2)
				#Add Colourbar (Axis, Label, Bins)
				label = VariableLabelMaker(varlist)
				Colourbar(ax[0],label[i],5,Lim=CbarMinMax(PROES,LineoutsOrientation[k]))

				#Plot Waveform.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)
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
				fig,ax = figure(image_aspectratio,2)
				ax[0].plot(Phaseaxis,TemporalPROES, lw=2)
				Ylabel = 'Spatially Integrated '+varlist[i]
				ImageOptions(ax[0],Ylabel=Ylabel,Crop=False)

				#Plot Waveform onto Temporally collapsed PROES.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' TemporalPROES'+ext)
				plt.close('all')

				#Plot Temporally Collapsed PROES with required axis.
				fig,ax = figure(image_aspectratio,2)
#				try: ax[0].plot(Raxis,SpatialPROES, lw=2)		### HACKY ###
#				except: ax[0].plot(Zaxis,SpatialPROES, lw=2)	### HACKY ###
				Xlabel = 'Phase [$\omega$t/2$\pi$]'
				Ylabel = 'Temporally Integrated '+varlist[i]
				ImageOptions(ax[0],Xlabel,Ylabel,Crop=False)

				#Plot Waveform onto Spatially collapsed PROES.
#				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
#				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

				#plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' SpatialPROES'+ext)
				plt.close('all')

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

if any([savefig_trendphaseresolved, savefig_phaseresolve1D ,savefig_phaseresolve2D ,savefig_PROES]) == True:
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
#             Notes             #
#===============================#

# Disabled the following warning message regarding scalar assignment of 'arr' axis.
# /home/sjd549/.local/lib/python2.7/site-packages/numpy/ma/core.py:6385





	##### ANALYSIS NEEDS MADE INTO A GENERAL FUNCTION. #####
	##### ATTACH AS SEPERATE 1D/2D TREND (KNUDSEN LIKE) #####
	##### ALSO MAKE NORMALIZATION OF VZ/VR-NEUTRAL #####
#====================================================================#
				  	#SOUND SPEED - SHOCKWAVE ANALYSIS#
#====================================================================#

if True == False:
	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#For Argon, GasConstant in units of [J/kg K]
	AdiabaticIndex=1.66
	GasConstant=208

	# SoundVelocity=np.sqrt(k*Pressure/Density)

	SoundVelocity=np.sqrt(AdiabaticIndex*GasConstant*NeutralTemp)
#endif

#=====================================================================#
#=====================================================================#







#====================================================================#
				#MORE DETAILED GRAPHS/CONTOUR PLOTS#
#====================================================================#

if True == False:
	#SeaBorn Colour, as_cmap causes issue!
	#Set_Style should remove grids but doesn't.
	try:
		import seaborn as sns
		GlobalCmap = sns.color_palette("muted", as_cmap=True)
		sns.set_style("whitegrid", {'axes.grid' : False})
	except:
		GlobalCmap = matplotlib.cm.get_cmap('jet')
	#endtry

	#NB:
	#Need to replace plotting functions with:
	#im = ax.imshow(Image,extent=extent,cmap=GlobalCmap,origin="lower")
#endif

#=====================================================================#
#=====================================================================#









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








