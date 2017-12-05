#!/usr/bin/env python

#################################
#		Point of Contact		#
#								#
#	   Mr. Scott J. Doyle		#
#	   University of York		#
#	   York Plasma Institute	#
#	   1&2 Genesis Building		#
#	   North Yorkshire, UK		#
#	   sjd549@york.ac.uk		#
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
	os.system('sudo apt-get install python-pip')
	os.system('sudo apt-get install python-matplotlib')
	os.system('sudo apt-get install python-numpy')
	os.system('sudo apt-get install python-scipy')
	os.system('sudo apt-get install ffmpeg')
	os.system('pip install findtools')
	os.system('pip install tqdm')
	print ''
	print ''
#endif

import matplotlib.cm as cm
import numpy as np
import scipy as sp
import math as m
import os, sys
import os.path

from mpl_toolkits.axes_grid1 import make_axes_locatable
from findtools.find_files import (find_files, Match)
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *


#====================================================================#
				  		#DEFAULT PARAMETERS#
#====================================================================#

#Define Misc parameters, Do Not Change Unless Required.

#Enviroment variables.
Magmesh = 1				#initmesh.exe magnification factor. (almost obsolete)
ierr = 0				#OldDebugMode (almost obsolete)

#Create switchboard directory for GUI.
Switchboard = {}

#Activates various debug outputs.
DebugMode = False

#Tweaks and fixes for 'volitile' diagnostics.
AxialLine = 80 						#Z-axis line for thrust calculation.  YPR=80, ESCTest=42,
Manualbiasaxis = ''					#'Axial' or 'Radial'. (empty '' for auto)

#List of recognised atomic density sets, add new sets as required.
ArgonReduced = ['AR','AR+','AR*']
ArgonFull = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR+','AR2+','AR2*']
Oxygen = ['O','O+','O-','O*','O2','O2+','O2*']
AtomicSet = ['E']+ArgonReduced+ArgonFull+Oxygen

#List of recognized ground-state neutral species for fluid analysis.
NeutSpecies = ['AR','AR3S','O2']

#Commonly used variable sets.
Ar = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR','AR+','AR2+','AR2*','E','TE','P-POT','TG-AVE','RHO','PRESSURE','EF-TOT','POW-RF','POW-RF-E','S-AR+','S-AR4P','SEB-AR+','SEB-AR4P','EB-ESORC','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','FR-AR+','FZ-AR+']
O2 = ['O2','O2+','O','O+','O-','E','TE','P-POT','TG-AVE','PRESSURE','EF-TOT','POW-RF','POW-RF-E','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','VR-ION-','VZ-ION-','FR-O-','FZ-O-']
ArO2 = Ar+O2

MSHC_PCMC = ['AR^0.5S','EB-0.5S','ION-TOT0.5S','AR^1.1B','EB-1.1B','ION-TOT1.1B']
PR_PCMC = ['AR^0.25','EB-0.25','ION-TOT0.25']

MSHC_Phase = ['S-E','S-AR+','S-AR4P','TE','PPOT']
PR_Phase = ['S-E','S-AR+','S-AR4P','SEB-AR4P','TE','PPOT']
PR_Phase2 = ['S-E','S-AR+','SEB-AR+','SEB-AR4P','SRCE-2437','TE','PPOT']



#Paper Trend Locations
#SDoyle2017a: 
#dz(5.50/118), dr(2.55/102) height=[24,43], Trend=[19]

#SDoyle2017b: 
#PROES, (Z=1.52,2.25,3.23) radialineouts = [31,46,66]
#Dielectric locations: [16,31],[16,46],[16,66]

        
        
        








#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested IEDF/NEDF Variables.
IEDFVariables = PR_PCMC		#Requested iprofile_2d variables (no spaces)
NEDFVariables = []			#Requested nprofile_2d variables (no spaces)

#Requested movie1/movie_icp Variables.
IterVariables = ['E','S-E','PPOT','TE']		#Requested Movie_icp (iteration) Variables.
PhaseVariables = PR_Phase					#Requested Movie1 (phase) Variables.
electrodeloc = [30,46]						#Cell location of powered electrode [R,Z].
waveformlocs = [[16,31],[16,46],[16,66]]		#Cell locations of additional waveforms [R,Z].

phasecycles = 2								#Number of phase cycles to be plotted.
DoFWidth = 41								#PROES Depth of Field Cells
#electrodeloc	#YPR [30,46],[16,46] #SPR [0,107] #MSHC [0,12]
#waveformlocs 	#YPR [[16,31],[16,46],[16,66]]
#DOFWidth		#YPR 41=2cm   #MSHC 10=0.47cm

#Requested TECPLOT Variables
Variables = Ar
MultiVar = []						#Additional variables plotted ontop of [Variables]
radialineouts = [31,46,66] 					#Radial 1D-Profiles to be plotted (fixed Z-mesh) --
heightlineouts = []					#Axial 1D-Profiles to be plotted (fixed R-mesh) |
TrendLocation = [] 					#Cell location For Trend Analysis [R,Z], ([] = min/max)
#YPR H0;R[31,46,66] #MSHC H0,20;R20


#Requested diagnostics and plotting routines.
savefig_itermovie = False					#Requires movie_icp.pdt
savefig_plot2D = False						#Requires TECPLOT2D.PDT

savefig_monoprofiles = False				#Single Variables; fixed height/radius
savefig_multiprofiles = False				#Multi-Variables; same folder
savefig_comparelineouts = False				#Multi-Variables; all folders
savefig_trendcomparison = False				#Single Variables; fixed cell location (or max/min)

savefig_phaseresolve1D = False				#1D Phase Resolved Images
savefig_phaseresolve2D = False				#2D Phase Resolved Images
savefig_PROES = False						#Phase-Resolved 2D Images

savefig_IEDF = False						#IN DEVELOPMENT, WORKS BUT UNRELIABLE.
savefig_EEDF = False						#IN DEVELOPMENT, NO PLOTTING ROUTINE.


#Steady-State diagnostics terminal output toggles.
print_meshconvergence = False				#Make More General: <_numerictrendaxis>
print_generaltrends = False					#Verbose Min/Max Trend Output.
print_Knudsennumber = False
print_totalpower = False
print_DCbias = False
print_thrust = False


#Image plotting options.
image_extension = '.png'					#Extensions { '.png', '.jpg', '.eps' }
image_aspectratio = [10,10]					#[x,y] in cm [Doesn't rotate dynamically]
image_radialcrop = [0.6]					#[R,Z] in cm
image_axialcrop = [1,4]						#[R,Z] in cm
#YPR R[0.6];Z[1,4]   #MSHC R[0.0,1.0];Z[0.5,2.5]

image_plotsymmetry = True
image_contourplot = True
image_normalize = False
image_plotgrid = False
image_plotmesh = False						#### NOT IMPLIMENTED ####
image_logplot = False
image_rotate = True


#Write data to ASCII files.
write_trendcomparison = True
write_phaseresolve = True
write_lineouts = True
write_plot2D = True


#============================#

#Overrides the automatic image labelling.
titleoverride = []
legendoverride = []
xaxisoverride = []
xlabeloverride = []
ylabeloverride = []
cbaroverride = ['NotImplimented']

#Commonly Used:
#'0','30','60','90','120','150','180','210','240','270','300','330'
#'100','200','300','400','500','600','700','800','900','1000'
#'13.56MHz','27.12MHz','40.68MHz','54.24MHz','67.80MHz'
#'67.80MHz','54.24MHz','40.68MHz','27.12MHz','13.56MHz'
#'1.0mm','1.5mm','2.0mm','2.5mm','3.0mm'








#####TODO#####
#need to update the cbar function and introduce an override.
#need to get seaborn introduced into the program en-masse.
#need to introduce 'garbage collection' at the end of each diagnostic.
#meshconvergence --> numerictrendaxis needs made for comparisons.
#Generally rename the switchboard to increase clarity.

#need to clean up the phase data in general and functionalise as much as possible.
#need to clean up the IEDF/NEDF/EEDF sections and functionalise.
#need to add the ability to compare IEDF/NEDF profiles between folders.









#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()
Dirlist = list()
Variablelist = list()
Variablelists = list()
IEDFVariablelist = list()
MovieVariablelist = list()
MovieVariablelists = list()
MovieITERlist = list()
Moviephaselist = list()
Geometrylist = list()
Legendlist = list()

Globalvarlist = list()
Globalnumvars = list()

#Create mesh_size lists and SI conversion
R_mesh = list()
Z_mesh = list()
Raxis = list()
Zaxis = list()

Depth = list()
Radius = list()
Height = list()
dr = list()
dz = list()

Isymlist = list()

#Create lists to store data
rawdata_2D = list()
rawdata_kin = list()
rawdata_phasemovie = list()
rawdata_itermovie_icp = list()
rawdata_IEDF = list()
rawdata_mcs = list()

Data = list()					#Data[folder][Variable][Datapoint]
DataIEDF = list()				#Data[folder][Variable][Datapoint]
DataEEDF = list()				#Data[folder][Variable][Datapoint]
IterMovieData = list()			#ITERMovieData[folder][timestep][variable][datapoints]
PhaseMovieData = list()			#PhaseMovieData[folder][timestep][variable][datapoints]

Moviephaselist = list()			#'CYCL = n'
Movieiterlist = list()			#'ITER = n'
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
print '                                                            v0.10.4 '
print '--------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if savefig_plot2D == True:
	print'# 2D Steady-State Image Processing'
if savefig_itermovie == True:
	print'# 2D Convergence Movie Processing'
if True in [savefig_phaseresolve2D,savefig_PROES]:
	print'# 2D Phase Resolved Movie Processing'
if True in [savefig_phaseresolve1D]:
	print'# 1D Phase Resolved Profile Processing'
if True in [savefig_monoprofiles,savefig_multiprofiles,savefig_comparelineouts]:
	print'# 1D Steady-State Profile Processing'
if True in [print_generaltrends,print_Knudsennumber,print_totalpower,print_DCbias,print_thrust]:
	print'# 1D Specific Trend Analysis'
if savefig_trendcomparison == True:
	print'# 1D Steady-State Trend Processing'
if True in [savefig_EEDF,savefig_IEDF]:
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

#Find all files ending in dir recursively from current directory.
sh_files_pattern = Match(filetype='f', name='*.PDT')
found_files1 = find_files(path='./', match=sh_files_pattern)
sh_files_pattern = Match(filetype='f', name='*.pdt')
found_files2 = find_files(path='./', match=sh_files_pattern)
sh_files_pattern = Match(filetype='f', name='initmesh.out')
found_files3 = find_files(path='./', match=sh_files_pattern)
sh_files_pattern = Match(filetype='f', name='icp.nam')
found_files4 = find_files(path='./', match=sh_files_pattern)
sh_files_pattern = Match(filetype='f', name='icp.out')
found_files5 = find_files(path='./', match=sh_files_pattern)

#Organize the files into an array and sort them alphabetically.
for found_file in found_files1:
	Dir.append(found_file)
#endfor
for found_file in found_files2:
	Dir.append(found_file)
#endfor
for found_file in found_files3:
	Dir.append(found_file)
#endfor
for found_file in found_files4:
	Dir.append(found_file)
#endfor
for found_file in found_files5:
	Dir.append(found_file)
Dir.sort()

#Calculate the number of seperate simulations involved for plotting.
#Create preamble Dir list for saving plots back into relevant folders.
#Identifies the first '/' reading filename in reverse saving as 'n','m'
#These indices are then used to cut off file names and save the preamble.
try:
	Dirlist.append(Dir[0][:len(Dir[0])-Dir[0][::-1].index('/')])
except:
	print '#===============================#'
	print 'No data found, aborting analysis.'
	print '#===============================#'
	print ''
	exit()
#endtry
numfolders = 1
for i in range(0, len(Dir)-1):
	n = Dir[i][::-1].index('/')
	m = Dir[i+1][::-1].index('/')

	currentdir = Dir[i][:len(Dir[i])-n]
	nextdir = Dir[i+1][:len(Dir[i+1])-m]
	if currentdir != nextdir:
		Dirlist.append(nextdir)
		numfolders += 1
	#endif
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
		#endif

	except:
		#If data for current file exists, ask for manual input.
		if l <= len(TEC2D)-1:

			#Manual input of mesh sizes.
			ierr = 1

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
	try:
		SImeshdata = open(icpnam[l]).readlines()

		#Retrieve useful input variables from icp.nam.
#		NUMPHASE = int(filter(lambda x: x.isdigit(),filter(lambda x:'IMOVIE_FRAMES' in x,SImeshdata)[0]))
#		NUMMETALS = int(filter(lambda x: x.isdigit(),filter(lambda x:'IMETALS' in x,SImeshdata)[0]))+1
#		MATERIALS = filter(lambda x: 'CMETAL=' in x, SImeshdata)[0].split()[1:NUMMETALS]

		#Input frequencies/voltages/powers
#		VRFM = filter(lambda x: 'VRFM=' in x, SImeshdata)[0].split()[1:NUMMETALS]
#		VRFM2 = filter(lambda x: 'VRFM_2=' in x, SImeshdata)[0].split()[1:NUMMETALS]
#		FREQM = filter(lambda x: 'FREQM=' in x, SImeshdata)[0].split()[1:NUMMETALS]
#		FREQM2 = filter(lambda x: 'FREQM_2=' in x, SImeshdata)[0].split()[1:NUMMETALS]
#		FREQICP = float(filter(lambda x:'FREQ=' in x, SImeshdata)[0].strip(' \t\n\r,=FREQ'))
#		IRFPOW = float(filter(lambda x:'IRFPOW=' in x, SImeshdata)[0].strip(' \t\n\r,=IRFPOW'))

		#SI Conversion unit extraction.
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

		#clean up variables and assign required types.
#		for i in range(0,NUMMETALS-1):
#			MATERIALS[i] = MATERIALS[i].strip(',\'')
#			VRFM[i] = float(VRFM[i].strip(','))
#			VRFM2[i] = float(VRFM2[i].strip(','))
#			FREQM[i] = float(FREQM[i].strip(','))
#			FREQM2[i] = float(FREQM2[i].strip(','))
		#endfor

		#Obtain useful parameters for diagnostics
#		MinFreq = min(filter(lambda a: a != 0, FREQM+FREQM2+[FREQICP]))
#		MaxFreq = max(filter(lambda x: x != 0, FREQM+FREQM2+[FREQICP]))
		dr.append(Radius[-1]/(R_mesh[-1]-1))
		dz.append(Height[-1]/(Z_mesh[-1]-1))

	except:
		ierr = 2

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
		if "CYCL=" in Rawdata[i]: 
			break
		else: CurrentRow = Rawdata[i].split()

		#For all elements in the current row, convert to string and save in list.
		for j in range(0,len(CurrentRow)):
			try: DataArray1D.append(float(CurrentRow[j]))
			except: Avoids_String_Conversion_Error = 1
		#endfor
	#endfor

	#Seperate total 1D array into 2D array with data for each variable.
	#Offset data by a certain number of variable 'chunks' if requested.
	offset = (Rmesh*Zmesh)*offset	   #Convert from variables to rows.
	for i in range(0,numvariables):
		numstart = (Zmesh*Rmesh)*(i)+offset
		numend = (Zmesh*Rmesh)*(i+1)+offset
		CurrentFolderData.append(list(DataArray1D[numstart:numend]))
	#endfor

	return(CurrentFolderData)
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
		elif variablelist[i] == 'JZ-NET':
			Variable = 'Axial Net Current Density'
			VariableUnit = '[mA cm$^{-2}$]'
		elif variablelist[i] == 'JR-NET':
			Variable = 'Radial Net Current Density'
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
			Variable = Variablelist[i]
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
			profile[i] = profile[i]*(0.01)
		#endfor
	if IsStringInVariable(variable,['VR-ION+','VZ-ION+','VR-ION-','VZ-ION-']) == True:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*(0.01)*(0.001)
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

	#IEDF/NEDF file readin.
	if savefig_IEDF == True:

		#Define arguments and autorun conv_prof.exe if possible.
		IEDFVarArgs = ['1','1','1','1','1'] #### THIS IS HACKY, WON'T ALWAYS WORK ####
		args = ['pcmc.prof','title','1','1','1'] + IEDFVarArgs + ['0','0']
		DirAdditions = ['iprofile_tec2d.pdt','nprofile_tec2d.pdt','iprofile_tec1d.pdt', 'nprofile_tec1d.pdt']
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

	#Kinetics data readin - NOT CURRENTLY USED
	if True == False:

		#Load data from TECPLOT_KIN file and unpack into 1D array.
		rawdata, nn_kin = ExtractRawData(Dir,'TECPLOT_KIN.PDT',l)
		rawdata_kin.append(rawdata)
	#endif


#===================##===================#
#===================##===================#

	if savefig_itermovie == True:

		#Load data from movie_icp file and unpack into 1D array.
		try:
			itermovie_icp = filter(lambda x: 'movie_icp.pdt' in x, Dir)
			rawdata_itermovie_icp.append(open(itermovie_icp[l]).readlines())
			nn_itermovie = len(rawdata_itermovie_icp[l])
		except:
			ierr = 6
			print 'Unable to find movie_icp.pdt'
		#endtry

		#Identify length of variable section and save variables.
		for i in range(2,nn_itermovie):
			MovieVariablelist.append(str(rawdata_itermovie_icp[l][i][:-2].strip(' \t\n\r\"')))

			#Locate when variable section has stopped.
			if str(rawdata_itermovie_icp[l][i]).find('GEOMETRY') != -1:
				#Remove trailing value in variable list.
				MovieVariablelist = MovieVariablelist[:len(MovieVariablelist)-2]
				numvariables_movie = len(MovieVariablelist)
				break
			#endif
		#endfor

		#Identify length of header.
		for i in range(2,nn_itermovie):

			#Calculate headersize and identify beginning of data.
			if str(rawdata_itermovie_icp[l][i]).find('ITER') != -1:
				header_movie = i+1
				break
			#endif
		#endfor
		header_itermovie.append(header_movie)

		#Create Variablelists for each folder of data and refresh Variablelist
		MovieVariablelists.append(MovieVariablelist)
		MovieVariablelist = list()
		MovieITERlist_temp = list()
		data_array = list()

		#Unpack each row of 7 data points into single array of floats.
		#Removing 'spacing' between the floats and ignoring variables above data.
		for i in range(header_movie,nn_itermovie):
			numstart = 1
			for j in range(0,7):
				try:
					#Collect Iteration Details, then extract data as normal.
					if str(rawdata_itermovie_icp[l][i]).find('ITER') != -1:
						ITERstart = rawdata_itermovie_icp[l][i].find('ITER')
						MovieITERlist_temp.append(rawdata_itermovie_icp[l][i][ITERstart:ITERstart+9])
						break
					#endif
					data_array.append(float(rawdata_itermovie_icp[l][i][numstart:(numstart+10)]))

				except:
					This_means_there_was_a_space = 1
				#endtry
				numstart+= 11
			#endfor
		#endfor

		#Seperate total 1D array into sets of data for each variable.
		#Data is a 4D array of form (folder,timestep,variable,datapoints)
		tempdata,tempdata2 = list(),list()
		for j in range(1,len(MovieITERlist_temp)+1):

			#Collect data for each variable in turn, then reset per iteration.
			IterationStart = numvariables_movie*(j-1)
			IterationEnd = numvariables_movie*j
			for i in range(IterationStart,IterationEnd):
				#Offset of (Z_mesh[l]*R_mesh[l])*2 to avoid initial R and Z output.
				numstart = (Z_mesh[l]*R_mesh[l])*(i) + (Z_mesh[l]*R_mesh[l])*2
				numend = (Z_mesh[l]*R_mesh[l])*(i+1) + (Z_mesh[l]*R_mesh[l])*2
				tempdata.append(list(data_array[numstart:numend]))
			#endfor
			tempdata2.append(tempdata)
			tempdata = list()
		#endfor

		#Save all variables for folder[l] to Data and refresh lists.
		IterMovieData.append(tempdata2)
		MovieITERlist.append(MovieITERlist_temp)
		tempdata,tempdata2 = list(),list()
		data_array = list()
	#endif


#===================##===================#
#===================##===================#

	if True in [savefig_phaseresolve2D,savefig_phaseresolve1D,savefig_PROES]:

		#Load data from movie_icp file and unpack into 1D array.
		rawdata,nn_phasemovie = ExtractRawData(Dir,'movie1.pdt',l)
		rawdata_phasemovie.append(rawdata)

		#Read through all variables for each file and stop when list ends. 
		#Movie1 has geometry at top, therefore len(header) != len(variables).
		VariableEndMarker,HeaderEndMarker = 'GEOMETRY','ZONE'
#		Variablelist = ['Radius','Height']		#NewMethod
		Variablelist = list()					#OldMethod
		for i in range(2,nn_phasemovie):
			if HeaderEndMarker in str(rawdata_phasemovie[l][i]): 
				header_phase = i+2
				break
			if VariableEndMarker in str(rawdata_phasemovie[l][i]):
#				NumVariables = (i-1)	#Including R,Z		#NewMethod
				NumVariables = (i-3)	#Not Including R,Z	#OldMethod
			if len(rawdata_phasemovie[l][i]) > 1: 
				Variablelist.append(str(rawdata_phasemovie[l][i][:-2].strip(' \t\n\r\"')))
			#endif
		#endfor
		Variablelist = Variablelist[0:NumVariables]	#Keep Variables up till 'GEOMETRY'
		numvariables_phase = len(Variablelist)
		header_phasemovie.append(header_phase)









		#New functionalized method for extracting data, needs fixing.
		#The issue appears to be with incorrect saving of R/Z 'data'.
		#Old method skipped this, difficult to replicate in new method.
		#Would be nice to be able to extract it without the 'hacky' old method.
		OldMethod = True
		if OldMethod == False:

			#Fudged numbers until the icp.nam dictionary is fixed.
			NUMPHASE = 180

			#Rough method of obtaining the movie1.pdt cycle locations for data extraction.
			Cyclelocations = list()
			for j in range(0,len(rawdata_phasemovie[l])):
				if "CYCL=" in rawdata_phasemovie[l][j]:
					Cyclelocations.append(j+1)
				#endif
			#endfor

			#Cycle through all phases for current datafile, appending per cycle.
			CurrentFolderData,CurrentFolderPhaselist = list(),list()
			for i in range(0,NUMPHASE):
				CurrentPhaseData = SDFileFormatConvertorHPEM(rawdata_phasemovie[l],Cyclelocations[i],numvariables_phase,offset=2)

				CurrentFolderPhaselist.append('CYCL = '+str(i+1))
				CurrentFolderData.append(CurrentPhaseData)
			#endfor
			Moviephaselist.append(CurrentFolderPhaselist)
			PhaseMovieData.append(CurrentFolderData)
		#endif

			#========== NOTES ==========#
			#Inclusion of the optional 'offset' allows for saving of R,Z but offseting of data.
			#Maybe not offset the data and figure out how to change the diagnostics to reflect 
			#this change and allow the phasedata to work in the same way as the 2D data.

			Test = True
			if Test == True:
				#[Folder][Phase][Variable][Data]
				print len(PhaseMovieData), numfolders
				print len(PhaseMovieData[0]), NUMPHASE
				print len(PhaseMovieData[0][0]), len(Variablelist)
				print len(PhaseMovieData[0][0][0]), R_mesh[l]*Z_mesh[l]
	
				#Create empty 2D image of required size.
				for k in range(0,3):
					Data = PhaseMovieData[0][k][9-2]
					numrows = len(Data)/R_mesh[l]
					Image = np.zeros([numrows,R_mesh[l]])
	
					#Reshape data into 2D array for further processing.
					for j in range(0,numrows):
						for i in range(0,R_mesh[l]):
							Start = R_mesh[l]*j
							Row = Z_mesh[l]-1-j
							Image[Row,i] = Data[Start+i]
						#endfor
					#endfor
	
					plt.imshow(Image)
					plt.show()
					plt.close('all')
				#endfor
			#endif
		#endif










		if OldMethod == True:

			#Create Variablelists for each folder of data and refresh Variablelist
			MovieVariablelist,Moviephaselist_temp = list(),list()
			data_array = list()

			#Unpack each row of 7 data points into single array of floats.
			#Removing 'spacing' between the floats and ignoring variables above data.
			for i in range(header_phasemovie[l],nn_phasemovie):
				numstart = 1
				for j in range(0,7):
					try:
						#Collect Phase Details, then extract data as normal.
						if str(rawdata_phasemovie[l][i]).find('CYCL') != -1:
							phasestart = rawdata_phasemovie[l][i].find('CYCL')
							Moviephaselist_temp.append(rawdata_phasemovie[l][i][phasestart:phasestart+9])
							break
						#endif
						data_array.append(float(rawdata_phasemovie[l][i][numstart:(numstart+10)]))

					except:
						This_means_there_was_a_space = 1
					#endtry
					numstart+= 11
				#endfor
			#endfor

			#Seperate total 1D array into sets of data for each variable.
			#Data is a 4D array of form (folder,timestep,variable,datapoints)
			tempdata,tempdata2 = list(),list()
			for j in range(1,len(Moviephaselist_temp)+1):

				#Collect data for each variable in turn, then reset per iteration.
				IterationStart = numvariables_phase*(j-1)
				IterationEnd = numvariables_phase*j
				for i in range(IterationStart,IterationEnd):
					#Offset of (Z_mesh[l]*R_mesh[l])*2 to avoid initial R and Z output.
					numstart = (Z_mesh[l]*R_mesh[l])*(i) + (Z_mesh[l]*R_mesh[l])*2
					numend = (Z_mesh[l]*R_mesh[l])*(i+1) + (Z_mesh[l]*R_mesh[l])*2
					tempdata.append(list(data_array[numstart:numend]))
				#endfor
				tempdata2.append(tempdata)
				tempdata = list()
			#endfor

			#Save all variables for folder[l] to Data and refresh lists.
			PhaseMovieData.append(tempdata2)
			Moviephaselist.append(Moviephaselist_temp)
			tempdata,tempdata2 = list(),list()
			data_array = list()
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
del RADIUS,RADIUST,HEIGHT,HEIGHTT,DEPTH,SYM
del data_array,tempdata,tempdata2,templineout
del Energy,Fe,rawdata_mcs


#Alert user that readin process has ended and continue with selected diagnostics.
if any([savefig_plot2D, savefig_phaseresolve2D, savefig_itermovie, savefig_monoprofiles, savefig_multiprofiles, savefig_comparelineouts, savefig_phaseresolve1D, savefig_PROES, savefig_trendcomparison, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust, savefig_IEDF, savefig_EEDF]) == True:
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

#Returns a 2D array of inputted data with size [R_mesh] x [Z_mesh]
#Can optionally perform variable unit conversion if required.
#Image = ImageExtractor2D(Data,Variable=[]):
def ImageExtractor2D(Data,Variable=[],Rmesh=0,Zmesh=0):

	#WHY Rmesh=R_mesh[l],Zmesh=Z_mesh[l] NOT WORKING????
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
#fig,ax = figure(image_aspectratio,1)
def figure(aspectratio,subplots=1,shareX=False):
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
	else:
		fig, ax = plt.subplots(subplots, figsize=(9,9), sharex=shareX)
	#endif
	return(fig,ax)
#enddef



#=========================#
#=========================#



#Crops 2D image taking account of image rotation options.
#Takes image axis and returns nothing, use figure()
#Extent format: [ [Rmin,Rmax], [Zmin,Zmax] ] in cm.
#CropImage(ax[0],[[Rcrop],[Zcrop]]), assumes default axis.
def CropImage(ax=plt.gca(),Extent=[]):

		#Obtain default limits and rotate if needed.
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
		ax.set_xlim(R1,R2)
		ax.set_ylim(Z1,Z2)
	#endif
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
		ax.legend(Legend, loc=1, frameon=False)
	#endif

	#Set labels and ticksize.
	ax.set_xlabel(Xlabel, fontsize=24)
	ax.set_ylabel(Ylabel, fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#Force scientific notation for all axes, accounting for non-scalar x-axes.
	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')
	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
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
#ImagePlotter2D(Image,extent,image_aspectratio,fig,ax[0]):
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
		profile = Normalize(profile)
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
#ImagePlotter2D(Image,extent,image_aspectratio,fig,ax[0]):
def ImagePlotter2D(Image,extent,aspectratio,variable='N/A',fig=111,ax=111):

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

	#Plot image with or without contour plots, allowing for contour failures.
	if image_contourplot == True:
		try:
			im = ax.contour(Image,extent=extent,origin="lower")
			im = ax.imshow(Image,extent=extent,origin="lower")
		except:
			im = ax.imshow(Image,extent=extent,origin="lower")
		#endtry
	else:
		im = ax.imshow(Image,extent=extent,origin="lower")
	#endif
	return(fig,ax,im)
#enddef



#=========================#
#=========================#



#Creates a 1D image from an array of supplied points.
#Image plotted onto existing axes, figure() should be used.
#NormFactor = 0 will normalize to maximum of given profile.
#TrendPlotter(Image,Xaxis,0)
def TrendPlotter(TrendArray,Xaxis,NormFactor=0):

	#Normalize data to provided normalization factor if required.
	if image_normalize == True:
		TrendArray = Normalize(TrendArray,NormFactor)
	#endif

	#Choose how to plot the trends.
	if print_meshconvergence == True:
		#Plot results against number of cells in each mesh for convergence studies.
		numcells = list()
		for l in range(0,numfolders):
			numcells.append(Z_mesh[l]*R_mesh[l])
		#endfor
		plt.plot(numcells[::-1],TrendArray[::-1], 'o-', lw=2)
	else:
		#Plot results against strings pulled from folder names for batch studies.
		plt.plot(range(0,numfolders),TrendArray, 'o-', lw=2)
		if len(xaxisoverride) > 0:
			plt.xticks(np.arange(0,numfolders), xaxisoverride)
		else:
			plt.xticks(np.arange(0,numfolders), Xaxis)
		#endif
	#endif

	return()
#enddef



#=========================#
#=========================#



#Creates and plots a colourbar with given label and binsize.
#Takes colourbar axis, label string, number of colour bins
#Also takes normalization factors in form [min,max].
#Returns cbar axis if further changes are required.
def Colourbar(ax,Label,Bins,Norm=[]):

	#Colourbar plotting details
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)
	#Number of ticks, set scientific notation.
	tick_locator = ticker.MaxNLocator(nbins=Bins)
	cbar.locator = tick_locator
	cbar.set_label(Label, rotation=270,labelpad=30,fontsize=24)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.update_ticks()
	#Size of font
	cbar.ax.yaxis.offsetText.set(size=18)
	yticks(fontsize=18)

	#Apply normalized colourbar if limits specified.
	if len(Norm) == 2: im.set_clim(vmin=Norm[0], vmax=Norm[1])

	return(cbar)
#enddef



#=========================#
#=========================#


#Generates an SI axis for a 1D profile plot.
#Takes orientation, symmetry and phasecycle options.
#Returns 1D array in units of [cm] or [omega*t/2pi].
def GenerateAxis(Orientation,Isym=Isymlist[l],phasepoints=range(0,180)):
	
	#Extract number of phase datapoints and create axis list.
	phasepoints = len(phasepoints)
	axis = list()

	if Orientation == 'Radial':
		if Isym == 1:
			for i in range(0,2*R_mesh[l]):
				axis.append(i*dr[l])
		#endfor
		else:
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



#Obtains a radial 1D profile at a requested axial location.
#Returns a 1D array for plotting and performs unit conversion.
def PlotRadialProfile(Data,process,variable,lineout,Rmesh=0,Isym=0):

	#WHY Rmesh=R_mesh[l],Zmesh=Z_mesh[l] NOT WORKING????
	#If no mesh sizes supplied, collect sizes for current global folder.
	if Rmesh == 0 or Isym == 0:
		Rmesh,Isym = R_mesh[l],Isymlist[l]
	#endif

	#Obtain start location for requested data and perform SI conversion.
	ZStart = Rmesh*lineout
	ZEnd = Rmesh*(lineout+1)

	#Plot lines for each variable at each requested slice, ignoring ones that fail.
	try:
		#If mesh is symmetric, copy the data over and make full plot.
		if Isym == 1:
			Zend = len(Data[process])-ZStart
			Zstart = len(Data[process])-ZEnd
			Rlineout = Data[process][Zstart:Zend][::-1]
			for m in range(0,len(Rlineout)):
				Rlineout.append(Data[process][Zstart:Zend][m])
			#endfor

		#If the data isn't symmetric, just plot as is.
		elif Isym == 0:
			Zend = len(Data[process])-ZStart
			Zstart = len(Data[process])-ZEnd
			Rlineout = Data[process][Zstart:Zend]
		#endif
	except:
		This_Means_A_Variable_Was_Missing_From_The_Data = 1
	#endtry

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
		Trend,Min,Max = Normalize(Image,NormFactor=0)
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
		MaxValueTrend,MaxMin,MaxMax = Normalize(MaxValueTrend,NormFactor=0)
		MinValueTrend,MinMin,MinMax = Normalize(MinValueTrend,NormFactor=0)
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
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

		#Setting the radial ticks for beauty purposes ~HACKY~
		R = Radius[l]
		SymTicks = [round(-R, 1), round(-R/2, 1), 0, round(R/2, 1), round(R, 1)]
		NoSymTicks = [0, round(R/2, 1), round(R, 1)]

		#Reshape specific part of 1D Data array into 2D image for plotting.
		for k in tqdm(range(0,len(processlist))):

			#Extract full 2D image for further processing.
			Image = ImageExtractor2D(Data[l][processlist[k]],Variablelist[k])

			#Generate and rotate figure as requested.
			extent,aspectratio = DataExtent(l)
			fig,ax,im = ImagePlotter2D(Image,extent,aspectratio,Variablelist[k])

			#Define image beautification variables.
			if image_rotate == True:
				Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
			elif image_rotate == False:
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				plt.gca().invert_yaxis()
			#endif

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = '2D Steady State Plot of '+Variablelist[k]+' for \n'+Dirlist[l][2:-1]
			ImageOptions(ax,Xlabel,Ylabel,Title)

			#Add Colourbar (Axis, Label, Bins)
			label,bins = VariableLabelMaker(Variablelist),5
			cax = Colourbar(ax,label[k],bins)

			#Write data to ASCII files if requested.
			if write_plot2D == True:
				DirWrite = CreateNewFolder(Dir2Dplots, '2Dplots Data')
				WriteDataToFile(Image, DirWrite+Variablelist[k])
			#endif

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot '+Variablelist[k]+ext)
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
if savefig_itermovie == True:

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create new folder to keep convergence variable folders in.
		DirConvergence = CreateNewFolder(Dirlist[l],'Convergence/')

		#Create processlist for each folder as required.
		iterprocesslist,IterVariablelist = VariableEnumerator(IterVariables,rawdata_itermovie_icp[l],header_itermovie[l])
		#Skip over the R and Z processes as they are not saved properly in iterdata.
		for i in range(0,len(iterprocesslist)):
			iterprocesslist[i] = iterprocesslist[i]-2
		#endfor

		#Create list and x-axis for convergence trend plotting.
		ConvergenceTrends,Xaxis = list(),list()
		for i in range(0,len(MovieITERlist[l])):
			Xaxis.append(filter(lambda x: x.isdigit(), MovieITERlist[l][i]))
		#endfor

		#for all variables requested by the user.
		for i in tqdm(range(0,len(iterprocesslist))):

			#Create new folder to keep output plots.
			DirMovieplots = CreateNewFolder(DirConvergence,IterVariablelist[i]+'_2DConvergence/')

			#Create empty image array based on mesh size and symmetry options.
			try: 
				numrows = len(IterMovieData[l][0][0])/R_mesh[l]
				Image = np.zeros([numrows,R_mesh[l]])
			except: 
				print 'No Iteration Data Found For '+Dirlist[l]
				break
			#endtry
			
			#Append new list to convergenceTrends for each variable.
			ConvergenceTrends.append(list())

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for k in range(0,len(MovieITERlist[l])):

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(IterMovieData[l][k][iterprocesslist[i]],IterVariablelist[i])
				#Take Max value of image for general convergence trend.
				ConvergenceTrends[-1].append( sum(Image.flatten())/len(Image.flatten()) )

				#Generate and rotate figure as requested.
				extent,aspectratio = DataExtent(l)
				fig,ax,im = ImagePlotter2D(Image,extent,aspectratio,Variablelist[i])

				#Define image axis labels.
				if image_rotate == True:
					Xlabel,Ylabel = 'Axial Distance Z [cm]','Radial Distance R [cm]'
				elif image_rotate == False:
					Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
					plt.gca().invert_yaxis()
				#endif

				#Image plotting details.
				Title = str(MovieITERlist[l][k])
				ImageOptions(ax,Xlabel,Ylabel,Title)

				#Add Colourbar (Axis, Label, Bins)
				label,bins = VariableLabelMaker(IterVariablelist),5
				cax = Colourbar(ax,label[i],bins)

				#Save to seperate folders inside simulation folder.
				num1,num2,num3 = k % 10, k/10 % 10, k/100 % 10
				Number = str(num3)+str(num2)+str(num1)
				savefig(DirMovieplots+IterVariablelist[i]+'_'+Number+ext)
				plt.close('all')
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+IterVariablelist[i])
		#endfor


		#=================#


		#Plot a convergence check for all variables in each folder.
		Legend = VariableLabelMaker(IterVariablelist)
		fig, ax = plt.subplots(1, figsize=(10,10))

		#Normalize and plot each variable in ConvergenceTrends to single figure.
		for i in range(0,len(ConvergenceTrends)):
			ConvergenceTrends[i] = Normalize(ConvergenceTrends[i])[0]
			ax.plot(Xaxis,ConvergenceTrends[i], lw=2)
		#endfor

		#Image plotting details.
		Title = 'Convergence of '+str(IterVariablelist)+' for \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Simulation Iteration','Normalized Mesh-Average Value'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend,Crop=False)
		ax.set_ylim(0,1.02)

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
					if write_lineouts == True:
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
				plt.savefig(DirRlineouts+'1D_Radial_'+Variablelist[i]+' profiles'+ext)
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
					if write_lineouts == True:
						SaveString = '_Z='+str(round((heightlineouts[j])*dr[l], 2))+'cm'
						DirWrite = CreateNewFolder(DirZlineouts, 'Axial_Data')
						WriteDataToFile([Zaxis,Zlineout], DirWrite+Variablelist[i]+SaveString)
					#endif
				#endfor

				#Apply image options and axis labels.
				Title = 'Height Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]
				Xlabel,Ylabel = 'Axial Distance Z [cm]',VariableLabelMaker(Variablelist)[i]
				ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

				#Save profiles in previously created folder.
				plt.savefig(DirZlineouts+'1D_Height_'+Variablelist[i]+' profiles'+ext)
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
				if write_lineouts == True and l == 0:
					WriteFolder = 'Z='+str(round((radialineouts[j])*dz[l], 2))+'cm Data'
					DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
					try: os.remove(DirWrite+Variablelist[k])
					except: a=1
				if write_lineouts == True:
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
				if write_lineouts == True and l == 0:
					WriteFolder = 'R='+str(round((heightlineouts[j])*dr[l], 2))+'cm_Data'
					DirWrite = CreateNewFolder(DirComparisons, WriteFolder)
					try: os.remove(DirWrite+Variablelist[k])
					except: a=1
				if write_lineouts == True:
					WriteDataToFile(Zlineout+['\n'], DirWrite+Variablelist[k], 'a')
				#endif

				#Apply image options and axis labels.
				Title = 'Comparison of '+Variablelist[k]+' Profiles at Z='+str(round((heightlineouts[j])*dr[l], 2))+'cm for \n'+Dirlist[l][2:-1]
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


#=====================================================================#
#=====================================================================#

































#====================================================================#
				 #GENERAL ENERGY DISTRIBUTION ANALYSIS#
#====================================================================#

#====================================================================#
				#ION-NEUTRAL ANGULAR ENERGY DISTRIBUTIONS#
#====================================================================#

if savefig_IEDF == True:

	#For all simulation folders.
	for l in range(0,numfolders):

		#Create new folder for keeping EEDF/IEDF if required.
		DirEDF = CreateNewFolder(Dirlist[l],'EDFplots')

		#Create processlist for requested EDF species and extract images.
		processlist,Variablelist = VariableEnumerator(IEDFVariables,rawdata_IEDF[l],header_IEDFlist[l])

		#For all requested variables.
		for i in tqdm(range(0,len(processlist))):
			EDFprofile = list()

			#Extract image from required variable and create required profile lists.
			Image = ImageExtractor2D(DataIEDF[l][processlist[i]],Rmesh=EDFangle,Zmesh=EDFbins)

			#Flatten angular distribution across all angles to produce energy distribution.
			for j in range(0,len(Image)): EDFprofile.append(sum(Image[j]))

			#Transpose Image for plotting and reverse both lists due to reading error.
			Image, EDFprofile = Image[::-1].transpose(), EDFprofile[::-1]


			#Plot the angular distribution and EDF of the required species.
			fig,ax = figure([13,9],2)
			fig.suptitle(Dirlist[l][2::]+'\n'+Variablelist[i]+' Angular Energy Distribution Function', y=0.995, fontsize=16)
			Extent=[0,len(Image[0]), -len(Image)/2,len(Image)/2]

			im = ax[0].imshow(Image,extent=Extent)
			ImageOptions(ax[0],Ylabel='Angular Dispersion [deg]',Crop=False)				
			#Add Colourbar (Axis, Label, Bins)
			cax = Colourbar(ax[0],Variablelist[i]+' EDF($\\theta$)',5)

			ax[1].plot(EDFprofile, lw=2)
			Xlabel,Ylabel = 'Energy [eV]',Variablelist[i]+' EDF($\\theta$)'
			ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)
			ax[1].set_xlim(0,50)

			plt.tight_layout()
			plt.savefig(DirEDF+Variablelist[i]+'_EDF'+ext)
			plt.close('all')
		#endfor
	#endfor
#endif

#===============================#

if any([savefig_IEDF, savefig_EEDF]) == True:
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

if savefig_trendcomparison == True or print_generaltrends == True:

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
			TrendPlotter(Trend,Xaxis,NormFactor=0)
			Title='Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
			ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

			#Write data to ASCII format datafile if requested.
			if write_trendcomparison == True:
				DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
				DCASCII = [Xaxis,Trend]
				WriteDataToFile(Trend, DirASCII+Variablelist[k]+' Trends')
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
			TrendPlotter(MaxTrend,Xaxis,NormFactor=0)
			Title='Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Varied Property','Max '+YaxisLegend[k]
			ImageOptions(ax,Xlabel,Ylabel,Title,Legendlist,Crop=False)

			#Write data to ASCII format datafile if requested.
			if write_trendcomparison == True:
				DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
				DCASCII = [Xaxis,MaxTrend]
				WriteDataToFile(MaxTrend, DirASCII+Variablelist[k]+' Trends')
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



if savefig_trendcomparison == True or print_DCbias == True:

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
		if Manualbiasaxis == 'Radial':
			DCbias.append(RadialDCbias)
		elif Manualbiasaxis == 'Axial':
			DCbias.append(AxialDCbias)
		else:
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
	if write_trendcomparison == True:
		DirASCII = CreateNewFolder(DirTrends,'Trend_Data')
		DCASCII = [Xaxis,DCbias]
		WriteDataToFile(DCASCII, DirASCII+'DCbias trends')
	#endif

	#Plot and beautify the DCbias, applying normalization if requested.
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(DCbias,Xaxis,NormFactor=0)

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



if savefig_trendcomparison == True or print_totalpower == True:

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
		TrendPlotter(Powers[k],Xaxis,NormFactor=0)

		#Apply image options and axis labels.
		Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
		Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
		ImageOptions(ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

		plt.savefig(DirTrends+RequestedPowers[k]+' Deposition Trends'+ext)
		plt.close('all')
	#endfor

	#Plot a comparison of all power depositions requested.
	fig,ax = figure(image_aspectratio,1)
	for k in range(0,len(RequestedPowers)):TrendPlotter(Powers[k],Xaxis,NormFactor=0)

	#Apply image options and axis labels.
	Title = 'Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Power Deposited [W]'
	ImageOptions(ax,Xlabel,Ylabel,Title,Legend=RequestedPowers,Crop=False)

	#Write data to ASCII format datafile if requested.
	if write_trendcomparison == True:
		DirASCII, TotalPowerASCII = CreateNewFolder(DirTrends,'Trend_Data'), [Xaxis]
		for k in range(0,len(RequestedPowers)): TotalPowerASCII.append(Powers[k])
		WriteDataToFile(TotalPowerASCII, DirASCII+'TotalPower Trends')
	#endif

	plt.savefig(DirTrends+'Power Deposition Comparison'+ext)
	plt.close('all')
#endif



#====================================================================#
				  	#ION/NEUTRAL THRUST ANALYSIS#
#====================================================================#



if savefig_trendcomparison == True or print_thrust == True:

	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initiate lists required for storing data.
	NeutralThrustlist,IonThrustlist,Thrustlist = list(),list(),list()
	Xaxis = list()

	#For all folders.
	for l in range(0,numfolders):

		#Create extract data for the neutral flux and neutral velocity.
		processlist,Variablelist = VariableEnumerator(['AR3S','VZ-NEUTRAL','VZ-ION+','FZ-AR3S','FZ-AR+','PRESSURE'],rawdata_2D[l],header_2Dlist[l])

		#Update X-axis with folder information.
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Extract radial density, velocity and pressure profiles across the discharge plane.
		AbortDiagnostic = False
		try:
			Density = PlotRadialProfile(Data[l],processlist[0],Variablelist[0],AxialLine)
			NeutralVelocity = PlotRadialProfile(Data[l],processlist[1],Variablelist[1],AxialLine)
			IonVelocity = PlotRadialProfile(Data[l],processlist[2],Variablelist[2],AxialLine)
			NeutralAxialFlux = PlotRadialProfile(Data[l],processlist[3],Variablelist[3],AxialLine)
			IonAxialFlux = PlotRadialProfile(Data[l],processlist[4],Variablelist[4],AxialLine)
			Pressure = PlotRadialProfile(Data[l],processlist[5],Variablelist[5],AxialLine)
		except:
			NeutralAxialFlux = np.zeros(R_mesh[l]*2)
			IonAxialFlux = np.zeros(R_mesh[l]*2)
			Pressure = np.zeros(R_mesh[l]*2)
			AbortDiagnostic = True
		#endtry
		if len(NeutralAxialFlux) == 0 or AbortDiagnostic == True:
			Thrustlist = np.zeros([numfolders])
			break
		#endif

		#Define which gas is used and calculate neutral mass per atom.
		NeutralIsp,IonIsp = list(),list()
		Argon,Xenon = 39.948,131.29			 #amu
		NeutralMass = Argon*1.67E-27		 #Kg

		DefaultTechnique = True
		#Both techniques assume cylindrical geometry, cartesian geometry will be overestimated.
		if DefaultTechnique == True:

			#Thrust based on integration over concentric ion/neutral momentum loss rate.
			NeutralThrust,IonThrust = 0,0
			for i in range(0,R_mesh[l]):
				#Calculate radial plane area of a ring at radius [i], correcting for central r=0.
				Circumference = 2*np.pi*(i*(dr[l]/100))		#m
				CellArea = Circumference*(dr[l]/100)		#m^2
				if CellArea == 0:
					CellArea = np.pi*(dr[l]/100)**2			#m^2
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

			Thrust = NeutralThrust + IonThrust							#N
			IonIsp = (sum(IonIsp)/len(IonIsp))/9.81						#s
			NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81			#s

			NeutralThrustlist.append( round(NeutralThrust*1000,5) )		#mN
			IonThrustlist.append( round(IonThrust*1000,5) )				#mN
			Thrustlist.append( round(Thrust*1000,5) )					#mN

		else:
			#Integration over concentric neutral momentum and differential pressure.
			Thrust = 0
			for i in range(0,R_mesh[l]):
				#Calculate radial plane length at radius [i] to integrate over.
				Circumference = 2*np.pi*(i*(dr[l]/100))		#m
				MassDensity = NeutralMass*Density[i]		#Kg/m^3
				ExitVelocity = (NeutralVelocity[i])**2		#m^2/s^2
				Thrust += ( MassDensity*ExitVelocity + Pressure[i]*133 )*Circumference #N
			#endfor
			Thrustlist.append( round(Thrust*1000,5) )
		#endif

		#Display thrust to terminal if requested.
		if print_thrust == True:
			print Dirlist[l], '@ Z=',round(AxialLine*dz[l],2),'cm'
			print 'NeutralThrust', round(NeutralThrust*1000,2), 'mN @ ', round(NeutralIsp,2),'s'
			print 'IonThrust:', round(IonThrust*1000,4), 'mN @ ', round(IonIsp,2),'s'
			print 'Thrust:',round(Thrust*1000,4),'mN'
			print ''
		#endif
	#endfor

	#Plot requested thrusts to 1st or 2nd Yaxis as required.
	fig,ax1 = figure(image_aspectratio,1)
#	ax2 = ax1.twiny()
	TrendPlotter(Thrustlist,Xaxis,NormFactor=0)
#	TrendPlotter(IonThrustlist,Xaxis,NormFactor=0)

	#Apply image options and save figure.
	Title = 'Thrust with changing '+TrendVariable+' \n'+Dirlist[l][2:-1]
	Xlabel,Ylabel = 'Varied Property','Total Thrust [mN]'
	Legend = ['Total Thrust','Ion Thrust']
	ImageOptions(ax1,Xlabel,Ylabel,Title,Legend,Crop=False)

	plt.savefig(DirTrends+'Thrust Trends'+ext)
	plt.close('all')


	#=====#=====#
	#=====#=====#


	if True == False:
		try:
			ThrustEfficiency = list()
			for i in range(0,len(Thrustlist)):
				ThrustEfficiency.append(Thrustlist[i]/DepositedPowerList[i])
			#endfor

			fig,ax = figure(image_aspectratio,1)
			TrendPlottingOptions1D(ThrustEfficiency,Xaxis,NormFactor=0)

			plt.xlabel('Varied Property', fontsize=24)
			plt.ylabel('Thrust Efficiency [N/kW]', fontsize=24)
			ax.tick_params(axis='x', labelsize=18)
			ax.tick_params(axis='y', labelsize=18)
			plt.xticks(np.arange(0,numfolders), Xaxis)
			plt.show()
		except:
			a = 1
		#endtry
		plt.close('all')
	#endif
#endif




#====================================================================#
				  		#KNUDSEN NUMBER ANALYSIS#
#====================================================================#


if bool(set(NeutSpecies).intersection(Variables)) == True:
	if savefig_trendcomparison == True or print_Knudsennumber == True:

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
			fig,ax,im = ImagePlotter2D(Image,extent,aspectratio)

			#Image plotting details, invert Y-axis to fit 1D profiles.
			Title = 'Knudsen Number Image for \n'+Dirlist[l][2:-1]
			Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
			ImageOptions(ax,Xlabel,Ylabel,Title)
			plt.gca().invert_yaxis()

			#Add Colourbar (Axis, Label, Bins)
			label,bins = 'Knudsen Number',5
			cax = Colourbar(ax,label,bins)

			#Save Figure
			plt.savefig(Dir2Dplots+'KnudsenNumber'+ext)
			plt.close('all')
		#endfor

		#Plot a comparison of all average Knudsen numbers.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(KnudsenAverage,Xaxis,NormFactor=0)

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

if any([savefig_trendcomparison, print_generaltrends, print_Knudsennumber, print_totalpower, print_DCbias, print_thrust]) == True:
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

		#Create processlist for each folder as required. (Always get PPOT)
		Processlist,Variablelist = VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
		PPOT = VariableEnumerator(['PPOT'],rawdata_phasemovie[l],header_phasemovie[l])[0][0]

		#Subtract 2 from process as variables R&Z are not saved properly in phasedata.
		for i in range(0,len(Processlist)): Processlist[i] -= 2
		PPOT -= 2

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Moviephaselist[l])
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])


		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseMovieData[l],PPOT)

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
		if write_phaseresolve == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'2DPhase_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPhase+'VoltageWaveforms')
		#endif

		#===============#


		#for all variables requested by the user.
		for i in tqdm(range(0,len(Processlist))):

			#Create new folder to keep specific plots.
			DirMovieplots = CreateNewFolder(DirPhaseResolved,Variablelist[i]+'_2DPhaseResolved/')

			#Obtain maximum and minimum values of current variable over all phases.
			MaxNormalize,MinNormalize = list(),list()
			for j in range(0,phasecycles):
				Image = ImageExtractor2D(PhaseMovieData[l][j][Processlist[i]],Variablelist[i])
				if image_logplot == True: Image = np.log(Image)
				
				MinNormalize.append(Normalize(Image)[1])
				MaxNormalize.append(Normalize(Image)[2])
			#endfor
			MaxNormalize,MinNormalize = max(MaxNormalize),min(MinNormalize)

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for j in range(0,len(Moviephaselist[l])):

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(PhaseMovieData[l][j][Processlist[i]],Variablelist[i])

				#Obtain image extent and axis labels based on image symmetry and rotation.
				Xlabel,Ylabel = 'Radial Distance R [cm]','Axial Distance Z [cm]'
				if image_rotate == True: Xlabel,Ylabel = Ylabel,Xlabel
				extent,aspectratio = DataExtent(l)

				#Create figure and axes, plot image on top and waveform underneath.
				fig,ax = figure(image_aspectratio,2)
				Title = 'Phase-Resolved '+Variablelist[i]+'\n'+str(Moviephaselist[l][j])
				fig.suptitle(Title, y=0.97, fontsize=18)

				#Plot 2D image, applying image options and cropping as required.
				fig,ax[0],im = ImagePlotter2D(Image,extent,aspectratio,Variablelist[i],fig,ax[0])
				ImageOptions(ax[0],Xlabel,Ylabel,Crop=True)

				#Add Colourbar (Axis, Label, Bins)
				Ylabel = VariableLabelMaker(Variablelist)
				cax = Colourbar(ax[0],Ylabel[i],5,Norm=[MinNormalize,MaxNormalize])

				#Plot waveform and apply image options.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

				#Cleanup layout and save images.
				fig.tight_layout()
				plt.subplots_adjust(top=0.90)
				num1,num2,num3 = j % 10, j/10 % 10, j/100 % 10
				Number = str(num3)+str(num2)+str(num1)
				savefig(DirMovieplots+Variablelist[i]+'_'+Number+ext)
				plt.close('all')


				#Write Phase data in ASCII format if required.
				if write_phaseresolve == True:
					DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'2DPhase_Data')
					Cycle = str( Moviephaselist[l][j].replace(" ", "") )
					SaveString = DirASCIIPhase+Variablelist[i]+'_'+Cycle
					WriteDataToFile(Image, SaveString)
				#endif
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+Variablelist[i])
		#endfor
	#endfor

	print'---------------------------------------'
	print'# 2D Phase-Resolved Processing Complete'
	print'---------------------------------------'
#endif




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
		Processlist,Variablelist = VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
		PPOT = VariableEnumerator(['PPOT'],rawdata_phasemovie[l],header_phasemovie[l])[0][0]

		#Subtract 2 from process as variables R&Z are not saved properly in phasedata.
		for i in range(0,len(Processlist)): Processlist[i] -= 2
		PPOT -= 2

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Moviephaselist[l])
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])


		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseMovieData[l],PPOT)

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
		if write_phaseresolve == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'1DPhase_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPhase+'VoltageWaveforms')
		#endif

		#==============#


		#for all requested variables.
		for i in tqdm(range(0,len(Processlist))):

			#Create new folder to keep specific plots.
			DirMovieplots = CreateNewFolder(DirPhaseResolved,Variablelist[i]+'_1DPhaseResolved/')

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
					NameString= Variablelist[i]+'_'+str(round(Lineouts[k]*dr[l],2))+'cm[R]'
				if LineoutsOrientation[k] == 'Radial':
					NameString= Variablelist[i]+'_'+str(round(Lineouts[k]*dz[l],2))+'cm[Z]'
				if savefig_phaseresolve1D == True:
					Dir1DProfiles = CreateNewFolder(DirMovieplots,NameString)
				#endif

				#Collect Normalization data for plotting.
				for j in range(0,len(Moviephaselist[l])):
					#Record local maximum and minimum for each phase.
					if LineoutsOrientation[k] == 'Axial': Profile = PlotAxialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],Lineouts[k])
					elif LineoutsOrientation[k] == 'Radial': Profile = PlotRadialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],Lineouts[k])
					#endif
					Profile,Minimum,Maximum = Normalize(Profile)
					VariableMax.append(Maximum)
					VariableMin.append(Minimum)
				#endfor
				#Find global maximum and minimum for all phases.
				VariableMax = max(VariableMax)
				VariableMin = min(VariableMin)

				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Moviephaselist[l])):

					if LineoutsOrientation[k] == 'Axial':
						ZlineoutLoc,axis = Lineouts[k],Zaxis
						PhaseResolvedlineout = PlotAxialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])
						lineoutstring = ' @ R='+str(round(Lineouts[k]*dr[l],2))+'cm \n'
						Xlabel = 'Axial Distance Z [cm]'
					elif LineoutsOrientation[k] == 'Radial':
						RlineoutLoc,axis = Lineouts[k],Raxis
						PhaseResolvedlineout = PlotRadialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],RlineoutLoc,R_mesh[l],Isymlist[l])
						lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm \n'
						Xlabel = 'Radial Distance R [cm]'
					#endif

					#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
					fig,ax = figure(image_aspectratio,2)
					Ylabel = VariableLabelMaker(Variablelist)
					fig.suptitle('Phase-Resolved '+Variablelist[i]+' for '+VariedValuelist[l]+lineoutstring+str(Moviephaselist[l][j]), y=0.97, fontsize=16)

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
					if write_phaseresolve == True:
						DirASCIIPhase = CreateNewFolder(DirPhaseResolved,'1DPhase_Data')
						DirASCIIPhaseloc = CreateNewFolder(DirASCIIPhase,lineoutstring[3:-2])
						Cycle = str( Moviephaselist[l][j].replace(" ", "") )
						SaveString = DirASCIIPhaseloc+Variablelist[i]+'_'+Cycle
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
				#SHEATH DYNAMICS AND PROES IMAGES#
#====================================================================#

#Plot Phase-Resolved profiles with electrode voltage and requested variables.
if savefig_PROES == True:

	#Tnitiate any required lists.
	VoltageWaveforms,WaveformBiases,VariedValuelist = list(),list(),list()

	#for all folders.
	for l in range(0,numfolders):

		#Create global folders to keep output plots and collect graph title.
		DirPROES = CreateNewFolder(Dirlist[l],'PROES')

		#Update X-axis with folder information.
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get PPOT)
		Processlist,Variablelist = VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
		PPOT = VariableEnumerator(['PPOT'],rawdata_phasemovie[l],header_phasemovie[l])[0][0]

		#Subtract 2 from process as variables R&Z are not saved properly in phasedata.
		for i in range(0,len(Processlist)): Processlist[i] -= 2
		PPOT -= 2

		#Generate SI scale axes for lineout plots. ([omega*t/2pi] and [cm] respectively)
		Phaseaxis = GenerateAxis('Phase',Isymlist[l],Moviephaselist[l])
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])


		#=============#

		#Extract waveforms from desired electrode locations.
		for j in range(0,len(waveformlocs)):
			VoltageWaveforms.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[0])
			WaveformBiases.append(WaveformExtractor(PhaseMovieData[l],PPOT,waveformlocs[j])[1])
		#endfor
		ElectrodeWaveform,ElectrodeBias = WaveformExtractor(PhaseMovieData[l],PPOT)

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
		if write_phaseresolve == True:
			ASCIIWaveforms = [Phaseaxis,ElectrodeWaveform]
			for j in range(0,len(waveformlocs)):
				ASCIIWaveforms.append(VoltageWaveforms[j])
			#endfor
			DirASCIIPROES = CreateNewFolder(DirPROES,'PROES_Data')
			WriteDataToFile(ASCIIWaveforms, DirASCIIPROES+'VoltageWaveforms')
		#endif

		#==============#


		#for all requested variables.
		for i in tqdm(range(0,len(Processlist))):

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
				for j in range(0,len(Moviephaselist[l])):

					#Refresh lists between each phasecycle.
					IntegratedDoFArray,DoFArrays = list(),list()

					#Collect each profile for stitching into a PROES image if required.
					if DoFWidth > 0:

						#Determine range of lineouts within the depth of field.
						DOFRegion = [(Lineouts[k]-DoFWidth),(Lineouts[k]+DoFWidth)]
						#Collect lineouts from DOF region and transpose to allow easy integration.
						for RLineoutLoc in range(DOFRegion[0],DOFRegion[1]):
							DoFArrays.append(PlotRadialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],RLineoutLoc))
						#endfor
						DoFArrays = np.asarray(DoFArrays).transpose().tolist()

						#Integrate DoF lineouts to form a single phase point PROES lineout.
						for m in range(0,len(DoFArrays)):
							IntegratedDoFArray.append( sum(DoFArrays[m])/(DoFWidth*2+1) )
						#endif
						PROES.append(IntegratedDoFArray)

					#If no DoF then simply collect lineout from required location.
					elif DoFWidth == 0:
						PROES.append( PlotRadialProfile(PhaseMovieData[l][j],Processlist[i],Variablelist[i],Lineouts[k]) )
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
					NameString = Variablelist[i]+'_'+lineoutstring[2::]
					Ylabel = 'Axial Distance Z [cm]'
					y1,y2 = Zaxis[0],Zaxis[-1]
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 1:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = Variablelist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					y1,y2 = Raxis[-1],-Raxis[-1]
				elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 0:
					lineoutstring = ' @ Z='+str(round(Lineouts[k]*dz[l],2))+'cm'
					NameString = Variablelist[i]+lineoutstring[2::]
					Ylabel = 'Radial Distance R [cm]'
					y1,y2 = Raxis[-1],0
				#endif
				DirPROESloc = CreateNewFolder(DirPROES,lineoutstring[3::])

				#Create PROES image along line of sight with phase-locked waveform.
				fig.suptitle( 'Simulated '+Variablelist[i]+' PROES for '+VariedValuelist[l]+lineoutstring+'\n DoF = '+str(round(DoFWidth*dz[l],2))+' cm', y=0.95, fontsize=18)
				im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='bottom',aspect='auto')
				ImageOptions(ax[0],Xlabel='',Ylabel=Ylabel,Crop=True)
				ax[0].set_xticks([])
				ax[0].set_xlim(x1,x2)
				#Add Colourbar (Axis, Label, Bins)
				label = VariableLabelMaker(Variablelist)
				cax = Colourbar(ax[0],label[i],5)

				#Plot Waveform.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

				#Cleanup layout and save images.
				fig.tight_layout()
				plt.subplots_adjust(top=0.85)
				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' PROES'+ext)
				plt.close('all')

				#Write PROES data in ASCII format if required.
				if write_phaseresolve == True:
					DirASCIIPROES = CreateNewFolder(DirPROES,'PROES_Data')
					DirASCIIPROESloc = CreateNewFolder(DirASCIIPROES,lineoutstring[3::])
					WriteDataToFile(PROES, DirASCIIPROESloc+Variablelist[i]+'_PROES')
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

				#Plot Temporally collapsed PROES with phase axis.
				fig,ax = figure(image_aspectratio,2)
				ax[0].plot(Phaseaxis,TemporalPROES, lw=2)
				Ylabel = 'Spatially Integrated '+Variablelist[i]
				ImageOptions(ax[0],Ylabel=Ylabel,Crop=False)

				#Plot Waveform onto Temporally collapsed PROES.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
				Xlabel,Ylabel = 'Phase [$\omega$t/2$\pi$]','Electrode Potential [V]'
				ImageOptions(ax[1],Xlabel,Ylabel,Crop=False)

				plt.savefig(DirPROESloc+VariedValuelist[l]+' '+NameString+' TemporalPROES'+ext)
				plt.close('all')

				#Plot Spatially collapsed PROES with required axis.
				fig,ax = figure(image_aspectratio,2)
				ax[0].plot(Raxis,SpatialPROES, lw=2)
				Xlabel = 'Phase [$\omega$t/2$\pi$]'
				Ylabel = 'Temporally Integrated '+Variablelist[i]
				ImageOptions(ax[0],Xlabel,Ylabel,Crop=False)

				#Plot Waveform onto Spatially collapsed PROES.
				ax[1].plot(Phaseaxis, ElectrodeWaveform, lw=2)
				ax[1].plot(Phaseaxis, ElectrodeBias, 'k--', lw=2)
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

if any([savefig_phaseresolve1D ,savefig_phaseresolve2D ,savefig_PROES]) == True:
	print'----------------------------------'
	print'# Phase Resolved Profiles Complete'
	print'----------------------------------'
#endif

#=====================================================================#
#=====================================================================#




























































#====================================================================#
							  #CODE DUMP#
#====================================================================#



#ERROR LIST, OUTDATED -- HOWEVER STILL IN MAIN CODE.

# ierr 1  ==  Can't find mesh sizes.		Minor Error (manual input)
# ierr 2  ==  Can't find SI sizes.			Minor Error (manual input)
# ierr 3  ==  Can't find TECPLOT2D file.	CRITICAL ERROR (Exit)
# ierr 4  ==  Can't find TECPLOT_KIN file.	CRITICAL ERROR (Exit)
# ierr 5  ==  Can't find movie1 file  	    CRITICAL ERROR (Exit)
# ierr 6  ==  Can't find movie_icp file     CRITICAL ERROR (Exit)

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
	btnlist[1]["command"] = lambda: toggle('savefig_itermovie',1)
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



# WAVEFORM OFFSET, MAY BE REQUIRED FOR STUFF?
if True == False:
	 Offset = m.ceil((1-(MinFreq/MaxFreq))*len(VoltageWaveform)*phasecycles)
#endif


















