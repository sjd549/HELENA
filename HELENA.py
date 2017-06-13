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

import matplotlib.cm as cm
import numpy as np
import scipy as sp
import math as m
import os, sys
import os.path


from mpl_toolkits.axes_grid1 import make_axes_locatable
from findtools.find_files import (find_files, Match)
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
#from tqdm import tqdm	#Progress bar
from pylab import *





#====================================================================#
				  		#DEFAULT PARAMETERS#
#====================================================================#

#Define Misc parameters, Do Not Change Unless Required.

#Enviroment variables.
numfolders = 1			#Fudge
Magmesh = 1				#initmesh.exe Mag-Factor
ierr = 0				#OldDebugMode

#Activates various debug outputs.
DebugMode = False

#Create switchboard dictionary
Switchboard = {}

#Tweaks and fixes for 'volitile' diagnostics.
AxialLine = 80 						#Z-axis line for thrust calculation.  YPR=80
Manualbiasaxis = ''					#'Axial' or 'Radial'. (empty '' for auto)

#List of recognised atomic sets, add new sets as required.
ArgonReduced = ['AR','AR+','AR*']
ArgonFull = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR+','AR2+','AR2*']
Oxygen = ['O','O+','O-','O*','O2','O2+','O2*']
AtomicSet = ['E']+ArgonReduced+ArgonFull+Oxygen

#List of neutral species for fluid analysis.
NeutSpecies = ['AR','AR3S','O2']

#Commonly used variable sets.
ArReduced = ['AR','AR+','E','TE','P-POT','RHO','TG-AVE','PRESSURE','EF-TOT','POW-RF','POW-RF-E','S-AR+','SEB-AR+', 'VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','FR-AR+','FZ-AR+']
ArFull = ['AR3S','AR4SM','AR4SR','AR4SPM','AR4SPR','AR4P','AR4D','AR+','AR2+','AR2*','E','TE','P-POT','TG-AVE','RHO','PRESSURE','EF-TOT','POW-RF','POW-RF-E','S-AR+','SEB-AR+', 'VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','FR-AR+','FZ-AR+']
O2 = ['O2','O2+','O','O+','O-','E','TE','P-POT','TG-AVE','PRESSURE','EF-TOT','POW-RF','POW-RF-E','VR-NEUTRAL','VZ-NEUTRAL','VR-ION+','VZ-ION+','VR-ION-','VZ-ION-','FR-O-','FZ-O-']


#Paper Trend Locations
#MSHC dz(5.50/118), dr(2.55/102) height=[24,43], Trend=[19]











#====================================================================#
					#SWITCHBOARD AND DIAGNOSTICS#
#====================================================================#

#Requested movie1/movie_icp Variables.
IterVariables = ['S-E','E','PPOT','TE']		#Requested Movie_icp (iteration) Variables.
PhaseVariables = ['S-E','TE']				#Requested Movie1 (phase) Variables.
electrodeloc = [0,0]						#Centre Cell of powered electrode [R,Z]. (T,B,L,R)
phasecycles = 1								#Number of phase cycles to be plotted.
#YPR [30,47] #SPR [0,107] #MSHC [0,12]

#Requested TECPLOT Variables
Variables = ArFull
MultiVar = []						#Additional variables plotted ontop of [Variables]
radialineouts = [] 					#Radial 1D-Profiles to be plotted (fixed Z-mesh)
heightlineouts = [0]				#Axial 1D-Profiles to be plotted (fixed R-mesh)
TrendLocation = [] 					#Cell location For Trend Analysis [R,Z], ([] = min/max)
#YPR H0;R47 #MSHC H0,20;R20


#Requested plotting routines.
savefig_itermovie = False					#Requires movie_icp.pdt
savefig_plot2D = True

savefig_radialines = False
savefig_heightlines = False
savefig_multiprofiles = False
savefig_comparelineouts = True

savefig_phaseresolvelines = False			#1D Phase Resolved Images
savefig_phaseresolve2D = False				#2D Phase Resolved Images
savefig_sheathdynamics = False				#PROES style images

#Steady-State diagnostics and terminal outputs.
savefig_trendcomparison = True
print_meshconvergence = False
print_generaltrends = False
print_KnudsenNumber = False
print_totalpower = False
print_DCbias = False
print_thrust = False


#Image plotting options.
image_aspectratio = [10,10]					#[x,y] in inches
image_plotsymmetry = True
image_contourplot = True
image_normalize = False						#### NORMALIZES TO EACH PROFILE SEPERATELY ###
image_plotgrid = False
image_logplot = False
image_rotate = True

image_plotmesh = False						#### NOT IMPLIMENTED ####




#============================#

#Overrides for automatic image labelling. (NB - Currently only for ComparisionProfiles)
titleoverride = ['NotImplimented']
legendoverride = []
xlabeloverride = []						#Only for Trend Plotter
ylabeloverride = ['NotImplimented']
cbaroverride = ['NotImplimented']
gridoverride = ['NotImplimented']

#Commonly Used:
#['$0$','$\\frac{\pi}{6}$','$\\frac{\pi}{3}$','$\\frac{\pi}{2}$','$\\frac{2\pi}{3}$','$\\frac{3\pi}{4}$','$\pi$']
#'0','30','60','90','120','150','180','210','240','270','300','330'
#'100','200','300','400','500','600','700','800','900','1000'
#'13.56MHz','27.12MHz','40.68MHz','54.24MHz','67.80MHz'
#'67.80MHz','54.24MHz','40.68MHz','27.12MHz','13.56MHz'
#'1.0mm','1.5mm','2.0mm','2.5mm','3.0mm'










#====================================================================#
 						#INITIATE GLOBAL LISTS#
#====================================================================#

#Create lists for basic processing
Dir = list()
Dirlist = list()
Variablelist = list()
Variablelists = list()
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
rawdata_mcs = list()

Data = list()					#Data[folder][Variable][Datapoint]
IterMovieData = list()			#ITERMovieData[folder][timestep][variable][datapoints]
PhaseMovieData = list()			#PhaseMovieData[folder][timestep][variable][datapoints]

header_itermovie = list()
header_phasemovie = list()
header_2Dlist = list()

Lineoutlist = list()
Zlineout = list()
Rlineout = list()

EEDFeVarray = list()
EEDFarray = list()

ElectrodeBias = list()
ElectrodeVoltage = list()
PhaseResolvedData = list()
PhaseResolvedData_avg = list()





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
print '                                                             v0.9.2 '
print '--------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if savefig_plot2D == True:
	print'# 2D Steady-State Image Processing'
if savefig_itermovie == True:
	print'# 2D Convergence Movie Processing'
if savefig_phaseresolve2D == True:
	print'# 2D Phase Resolved Movie Processing'
if True in [savefig_phaseresolvelines,savefig_sheathdynamics]:
	print'# 1D Phase Resolved Profile Processing'
if True in [savefig_radialines,savefig_heightlines,savefig_multiprofiles]:
	print'# 1D Steady-State Profile Processing'
if True in [print_generaltrends,print_KnudsenNumber,print_totalpower,print_DCbias,print_thrust]:
	print'# 1D Specific Trend Analysis'
if savefig_trendcomparison == True:
	print'# 1D Steady-State Trend Processing'
if savefig_comparelineouts == True:
	print'# 1D Steady-State Profile Comparisons'
print '-----------------------------------------'
print ''





#====================================================================#
					#OBTAINING FILE DIRECTORIES#
#====================================================================#

#Obtain system RAM.
mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
mem_gib = mem_bytes/(1024.**3)

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
for i in range(0, len(Dir)-1):
	n = Dir[i][::-1].index('/')
	m = Dir[i+1][::-1].index('/')

	tempname1 = Dir[i][:len(Dir[i])-n]
	tempname2 = Dir[i+1][:len(Dir[i+1])-m]
	if tempname1 != tempname2:
		Dirlist.append(Dir[i+1][:len(Dir[i+1])-m])
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
		for i in range(0,len(meshdata)):
			if 'ZONE' in meshdata[i]:
				#Split zone size at comma, R&Z values are given by "I=,J=" respectively.
				R = filter(lambda x: x.isdigit(), meshdata[i].split(",")[0])
				Z = filter(lambda x: x.isdigit(), meshdata[i].split(",")[1])
				R_mesh.append( int(R) )
				Z_mesh.append( int(Z) )
			#endif
		#endfor

	except ValueError:
		#Identify mesh size from initmesh.out file. (Issues with Q-VT and Magmesh)
		meshdata = open(mesh[l]).readline()
		R_mesh.append([int(i) for i in meshdata.split()][1])
		if Magmesh == 1:
			Z_mesh.append([int(i)+1 for i in meshdata.split()][3])
		elif Magmesh == 2:
			Z_mesh.append([int(i)+3 for i in meshdata.split()][3])
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

		for i in range(0,len(SImeshdata)):
			#Identify geometry section for outread.
			if SImeshdata[i][0:17].strip(' ') == '$DATAIN_GEOMETRY':

				#Refresh geometry list between runs.
				Geometrylist = list()
				#Remove all wording, leaving only numbers in the geometry section.
				for j in range(1,8):
					Geometrylist.append(SImeshdata[i+j].strip(' \t\n\r,=AEIUDGHTRPSXZYM'))
				#endfor

				#Select if plasma size or total size has been used and output SI size.
				if float(Geometrylist[3]) > 0.0:
					Radius.append(float(Geometrylist[3]))
				elif float(Geometrylist[5]) > 0.0:
					Radius.append(float(Geometrylist[5]))
				else:
					ierr = 2
				#endif
				if float(Geometrylist[3]) > 0.0:
					Height.append(float(Geometrylist[4]))
				elif float(Geometrylist[5]) > 0.0:
					Height.append(float(Geometrylist[6]))
				else:
					ierr = 2
				#endif
				Depth.append(float(Geometrylist[2]))
				if image_plotsymmetry == True:
					Isymlist.append(int(Geometrylist[1]))
				else:
					Isymlist.append(0)
				#endif

				#Calculate SI units for plotting.
				dr.append(Radius[-1]/(R_mesh[-1]-1))
				dz.append(Height[-1]/(Z_mesh[-1]-1))
			#endif
		#endfor
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

#Takes a 1D or 2D array and writes to a datafile in ASCII format.
#Two imputs, first is data to be written, second is filename string.
#WriteDataToFile(Image, FolderNameTrimmer(Dirlist[l])+Variablelist[k])
def WriteDataToFile(Data,filename):

	#Determine dimensionality of profile.
	if isinstance(Data[0], (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, 'w')
		for m in range(0,len(Data)):
			for n in range(0,len(Data[m])):
				datafile.write(str(Data[m][n]))
				datafile.write(' ')
			#endfor
			datafile.write('\n')
		#endfor
		datafile.close()

	#Lowest dimention is still list.
	elif isinstance(Data, (list, np.ndarray) ) == True:
		#Open new textfile and output 2D image data.
		datafile = open(filename, 'w')
		for n in range(0,len(Data)):
			datafile.write(str(Data[n]))
			datafile.write(' ')
		#endfor
		datafile.close()

	return()
#enddef


#Reads 1D or 2D data from textfile in ASCII format.
#One input, filename string, returns data array.
def ReadDataFromFile(Filename):
	OutputData = list()

	#Determine dimensionality of profile.
	if isinstance(Data[0], (list, np.ndarray) ) == True:
		#Read in 2D data from ASCII formatted file.
		datafile = open(Filename)
		for m in range(0,len(datafile)):
			Row = datafile.readline().split()
			for n in range(0,len(Row)):
				Row[n] = float(Row[n])
			#endfor
			OutputData.append(Row)
		#endfor

	#Lowest dimention is still list.
	elif isinstance(Data, (list, np.ndarray) ) == True:
		#Read in 1D data from ASCII formatted file.
		datafile = open(Filename)
		Row = datafile.readline().split()
		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData)
#enddef


#VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
#Enumerates requested variables and produces a processlist for plotting.
def VariableEnumerator(Variables,Rawdata,Header):
	processlist = list()
	legendlist = list()

	#For all requested variables, in the requested data header, find which match.
	for j in range(0,len(Variables)):
		for i in range(0,Header):

			#Compare variables and if they match, add to the process list.
			#Default uses [1:-3] slice for the variable string.
			if Variables[j] == Rawdata[i].replace(" ", "")[1:-3]:
				processlist.append(i)
				legendlist.append(Rawdata[i].replace(" ", "")[1:-3])
				break
			#endif
		#endfor
	#endfor
	return processlist,legendlist
#enddef


#Takes array of strings and compares to variable string.
#Returns true if any element of stringarray is in variable.
def IsStringInVariable(variable,stringarray):

#	#Should perform the same task but returns True every time.
#	OutputType = ['EFLUX-Z','EFLUX-R','FZ-','FR-']
#	if (String in variable for String in OutputType):
#		print variable
#	#endif

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
	Ionizationlist = ['S-','SEB-']
	Velocitylist = ['VZ-','VR-']
	Fluxlist = ['FZ-','FR-','EFLUX-R','EFLUX-Z']

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
			Variable = 'Bulk Electron Excitation'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-E':
			Variable = 'Secondary Electron Excitation'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'S-AR+':
			Variable = 'Bulk Ar+ Ionization Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'
		elif variablelist[i] == 'SEB-AR+':
			Variable = 'Secondary Ar+ Ionization Rate'
			VariableUnit = '[m$^{-3}$s$^{-1}$]'

		#Explicit Species Temperatures.
		elif variablelist[i] == 'TE':
			Variable = 'Electron Temperature'
			VariableUnit = '[eV]'
		elif variablelist[i] == 'TG-AVE':
			Variable = 'Neutral Gas Temperature'
			VariableUnit = '[K]'
		elif variablelist[i] == 'T-AR':
			Variable = 'Ar Temperature'
			VariableUnit = '[K]'
		elif variablelist[i] == 'T-AR+':
			Variable = 'Ar+ Temperature'
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
		elif variablelist[i] == 'EAMB-R':
			Variable = 'Radial E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] == 'ER':
			Variable = 'Radial E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] == 'EAMB-Z':
			Variable = 'Axial E-Field Strength'
			VariableUnit = '[Vcm$^{-1}$]'
		elif variablelist[i] == 'EZ':
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
			VariableUnit = '[Wcm$^{-3}$]'
		elif variablelist[i] == 'POW-RF':
			Variable = 'RF-Power Deposited'
			VariableUnit = '[Wcm$^{-3}$]'
		elif variablelist[i] == 'POW-RF-E':
			Variable = 'RF-Power Deposited by e-'
			VariableUnit = '[Wcm$^{-3}$]'


		#Implicit Variables.
		elif IsStringInVariable(variablelist[i],Ionizationlist) == True:
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

	#For ionization rates, convert from [cm-3 s-1] to [m-3 s-1]
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
			profile[i] = profile[i]#*100	### STILL IN [V cm-1] ###
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

	#For densities, convert from [cm-3] to [m-3]. (AtomicSet is defined in default parameters)
	if variable in AtomicSet:
		for i in range(0,len(profile)):
			profile[i] = profile[i]*1E6
		#endfor
	#endif

	return(profile)
#enddef


#Takes folder names and returns anything after final underscore.
def FolderNameTrimmer(Dirlist):
	try:
		Underscoreloc = str(Dirlist[::-1]).index('_')
		cutoff = (len(Dirlist)-Underscoreloc)
		FinalItem = Dirlist[cutoff:-1]
	except:
		FinalItem = str(Dirlist[2:])
	#endtry
	return(FinalItem)
#enddef


#Creates a new folder if one does not already exist, returns new directory.
def CreateNewFolder(Dir,string):
	try:
		NewFolderDir = Dir+string+'/'
		os.mkdir(NewFolderDir, 0755);
	except:
		a = 1
	#endtry
	return(NewFolderDir)
#enddef


#===================##===================#
#===================##===================#


print'---------------------'
print'Beginning Data Readin'
print'---------------------'

#Extraction and organization of data from .PDT files.
for l in range(0,numfolders):

	#Load data from TECPLOT2D file and unpack into 1D array.
	try:
		PLOT2D = filter(lambda x: 'TECPLOT2D.PDT' in x, Dir)
		rawdata_2D.append(open(PLOT2D[l]).readlines())
		nn_2D = len(rawdata_2D[l])
	except:
		ierr = 3
		print 'Unable to find TECPLOT2D.PDT'
	#endtry


	#Read through all variables for each file and stop when list ends.
	Variablelist.append('Radius')
	Variablelist.append('Height')
	for i in range(2,nn_2D):
		Variablelist.append(str(rawdata_2D[l][i][:-2].strip(' \t\n\r\"')))
		#Locate when variable section has stopped.
		if str(rawdata_2D[l][i]).find('ZONE') != -1:
			#Calculate headersize and remove trailing value in variable list.
			Variablelist = Variablelist[:len(Variablelist)-1]
			numvariables_2D = len(Variablelist)
			header_2D = numvariables_2D + 2
			break
		#endif
	#endfor
	header_2Dlist.append(header_2D)

	#Create Variablelists for each folder of data and refresh Variablelist
	Variablelists.append(Variablelist)
	Variablelist = list()


	#Unpack each row of 7 data points into single array of floats.
	#Removing 'spacing' between the floats and ignoring variables above data.
	tempdata, data_array = list(),list()
	for i in range(header_2D,nn_2D):
		numstart = 1
		for j in range(0,7):
			try:
				data_array.append(float(rawdata_2D[l][i][numstart:(numstart+10)]))
			except:
				This_means_there_was_a_space = 1
			#endtry
			numstart+= 11
		#endfor
	#endfor

	#Seperate total 1D array into sets of data for each variable.
	#Data is a 3D array of form (folder,variable,datapoints)
	for i in range(0,numvariables_2D):
		numstart = (Z_mesh[l]*R_mesh[l])*(i)
		numend = (Z_mesh[l]*R_mesh[l])*(i+1)
		tempdata.append(list(data_array[numstart:numend]))
	#endfor

	#Save all variables for folder[l] to Data and refresh lists.
	Data.append(tempdata)
	tempdata = list()
	data_array = list()


#===================##===================#
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
#===================##===================#


	if True in [savefig_phaseresolve2D,savefig_phaseresolvelines,savefig_sheathdynamics]:

		#Load data from movie_icp file and unpack into 1D array.
		try:
			phasemovie_icp = filter(lambda x: 'movie1.pdt' in x, Dir)
			rawdata_phasemovie.append(open(phasemovie_icp[l]).readlines())
			nn_phasemovie = len(rawdata_phasemovie[l])
		except:
			ierr = 5
			print 'Unable to find movie1.pdt'
		#endtry


		#Identify length of variable section and save variables.
		for i in range(2,nn_phasemovie):
			MovieVariablelist.append(str(rawdata_phasemovie[l][i][:-2].strip(' \t\n\r\"')))

			if str(rawdata_phasemovie[l][i]).find('GEOMETRY') != -1:
				#Remove trailing value in variable list.
				MovieVariablelist = MovieVariablelist[:len(MovieVariablelist)-2]
				numvariables_movie = len(MovieVariablelist)
				break
			#endif
		#endfor

		#Identify length of header.
		for i in range(2,nn_phasemovie):

			#Calculate headersize and identify beginning of data.
			if str(rawdata_phasemovie[l][i]).find('CYCL') != -1:
				header_movie = i+1
				break
			#endif
		#endfor
		header_phasemovie.append(header_movie)

		#Create Variablelists for each folder of data and refresh Variablelist
		MovieVariablelists.append(MovieVariablelist)
		MovieVariablelist = list()
		Moviephaselist_temp = list()


		#Unpack each row of 7 data points into single array of floats.
		#Removing 'spacing' between the floats and ignoring variables above data.
		for i in range(header_movie,nn_phasemovie):
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
		PhaseMovieData.append(tempdata2)
		Moviephaselist.append(Moviephaselist_temp)
		tempdata,tempdata2 = list(),list()
		data_array = list()
	#endif


#===================##===================#
#===================##===================#
#===================##===================#


	#Species EDF data readin - NOT CURRENTLY USED
	if True == False:

		#Load data from MCS.PDT file and unpack into 1D array.
		try:
			mcs = filter(lambda x: 'MCS.PDT' in x, Dir)
			rawdata_mcs.append(open(mcs[l]).readlines())
			nn_mcs = len(rawdata_mcs[l])
		except:
			ierr = 7
			print 'Unable to find MCS.PDT'
		#endtry

		header_mcs = 2
		nn_mcs = 81
		#Unpack each row of data points into single array of floats.
		#Removing 'spacing' between the floats and ignoring variables above data.
		for i in range(header_mcs,nn_mcs):
			rowlength = len(rawdata_mcs[l][i].split())
			for j in range(0,rowlength):
				try:
					data_array.append(float( rawdata_mcs[l][i].split()[j] ))
				except:
					NaN_Value = 1
				#endtry
			#endfor
		#endfor
		plt.plot(data_array)
		plt.show()
	#endif


	#Kinetics data readin - NOT CURRENTLY USED
	if True == False:

		#Load data from TECPLOT_KIN file and unpack into 1D array.
		try:
			PLOTKIN = filter(lambda x: 'TECPLOT_KIN.PDT' in x, Dir)
			rawdata_kin.append(open(PLOTKIN[l]).readlines())
			nn_kin = len(rawdata_kin[l])
		except:
			ierr = 4
			print 'Unable to find TECPLOT_KIN.PDT'
		#endtry
	#endif


#===================##===================#
#===================##===================#
#===================##===================#


	#Percentage complete printout.
	oldpercentage = int(((float(l)/numfolders)*100.0))
	newpercentage = int(((float(l+1)/numfolders)*100.0))
	if oldpercentage != newpercentage:
		print newpercentage,'%'
	#endif
#endfor

#Empty and delete any non-global data lists.
tempdata,tempdata2 = list(),list()
data_array,templineout = list(),list()
del data_array,tempdata,tempdata2,templineout


#=====================================================================#
#=====================================================================#


#Create global list of all variable names and find shortest list.
for m in range(0,numfolders):
	#Alphabetize the Variablelist and keep global alphabetized list.
	tempvarlist = VariableEnumerator(Variables,rawdata_2D[m],header_2Dlist[m])[1]
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



#Alert user that readin process has ended and continue with selected diagnostics.
if any([savefig_plot2D, savefig_phaseresolve2D, savefig_itermovie, savefig_radialines, savefig_heightlines, savefig_comparelineouts, savefig_multiprofiles, savefig_phaseresolvelines, savefig_sheathdynamics, savefig_trendcomparison, print_generaltrends, print_KnudsenNumber, print_totalpower, print_DCbias, print_thrust]) == True:
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


#Identifies if variable exists in all simulations, rejects if not.
#Allows for the comparison of datasets with different icp.dat files.
def VariableInterpolator(processlist,Variablelist,Comparisonlist):

	#Return default if atomic physics is the same in all datasets.
	if len(list(set(Comparisonlist).symmetric_difference(Variablelist))) == 0 :
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
				if interpolation[i] in Variablelist[j]:
					del Variablelist[j]
					del processlist[j]
				else:
					j += 1
				#endif
			#endfor
		#endfor
	#endif

	return(processlist, Variablelist)
#enddef



#=========================#
#=========================#



#Returns a 2D array of inputted data with size [R_mesh] x [Z_mesh]
#Can optionally perform variable unit conversion if required.
def ImageExtractor2D(Data,R_mesh,Z_mesh,Variable=[]):

	#Create empty 2D image of required size.
	numrows = len(Data)/R_mesh
	Image = np.zeros([numrows,R_mesh])

	#Reshape data into 2D array for further processing.
	for j in range(0,numrows):
		for i in range(0,R_mesh):
			Start = R_mesh*j
			Row = Z_mesh-1-j
			Image[Row,i] = Data[Start+i]
		#endfor
	#endfor

	#Convert units if required.
	Image = VariableUnitConversion(Image,Variable)

	return(Image)
#enddef



#=========================#
#=========================#



#Create figure of desired size and with variable axes.
#Returns figure and axes seperately.
def figure(aspectratio,subplots=1):
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]))
	else:
		fig, ax = plt.subplots(subplots, figsize=(9,9))
	#endif
	return(fig,ax)
#enddef



#=========================#
#=========================#



#Applies plt.options to current figure based on user input.
#Returns nothing, current image is required, use figure().
def ImageOptions():

	#Log image if possible, else linear. Grid default is off.
	if image_logplot == True:
		try: ax.set_yscale('log')
		except: ax.set_yscale('linear')
	#endif
	if image_plotgrid == True: plt.grid(True)
	#endif

	return()
#enddef



#=========================#
#=========================#



#Takes 1D or 2D array and returns array normalized to maximum value.
def Normalize(profile):

	#determine dimensionality of profile and select normaliztion method.
	if isinstance(profile[0], (list, np.ndarray, np.generic) ) == True:

		#Normalize 2D array to local maximum.
		FlatImage = [item for sublist in profile for item in sublist]
		NormalizedImage,NormFactor = list(),max(FlatImage)
		for i in range(0,len(profile)):
			NormalizedImage.append( [x/NormFactor for x in profile[i]] )
		#endfor
		profile = NormalizedImage
		return(profile)

	#Lowest dimention is still list.
	elif isinstance(profile, (list, np.ndarray, np.generic) ) == True:

		#Fix for division by zero.
		if max(profile) != 0: normalize = max(profile)
		else: normalize = 1
		#endif

		#Normalize 1D array to local maximum.
		for i in range(0,len(profile)):
			profile[i] = profile[i]/normalize
		#endfor
	#endif

	return(profile)
#enddef



#=========================#
#=========================#



#Create figure and plot a 2D image with associated image plotting requirements.
#Returns plotted image, axes and figure.
def ImagePlotter(Image,extent,aspectratio):
	fig, ax = figure(aspectratio)

	#Apply any required numerical changes to the image.
	if image_logplot == True:
		Image = np.log(Image)
	elif image_normalize == True:
		Image = Normalize(Image)
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
def TrendPlotter(TrendArray,Xaxis,Normalize=1):

	#Normalize data to provided normalization factor if required.
	if image_normalize == True:
		for i in range(0,len(TrendArray)):
			TrendArray[i] = TrendArray[i]/Normalize
		#endfor
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
		if len(legendoverride) > 0:
			plt.xticks(np.arange(0,numfolders), legendoverride)
		else:
			plt.xticks(np.arange(0,numfolders), Xaxis)
		#endif
	#endif

	#Standard image_<string> user modifcations.
	ImageOptions()

	return()
#enddef



#=========================#
#=========================#



#Takes pre-produced 2D image and 'beautifies' it for plotting.
def SymmetryConverter2D(Image,Isym,R_mesh,Z_mesh,Radius,Height):

	#Obtain image standard (non-rotated) aspect ratio.
	if len(image_aspectratio) == 2:
		aspectratio = image_aspectratio
	else:
		aspectratio = [9,9]
	#endif

	#Rotate image by 90 degrees anticlockwise and plot.
	if image_rotate == True:
		#Flip image axes and aspect ratios.
		RotateImage,SymRotateImage = Image.swapaxes(0,1),list()
		aspectratio = aspectratio[::-1]

		#If mesh uses symmetry, modify Image to conform to this.
		if Isym == 1:
			#Create new image by reversing height profiles and adding to beginning of image.
			for i in range(R_mesh):
				SymRotateImage.append(RotateImage[(R_mesh-1)-i])
				if i == (R_mesh-1):
					for j in range(R_mesh):
						SymRotateImage.append(RotateImage[j])
					#endfor
				#endif
			#endfor

			extent=[0,Height, -Radius,Radius]
			fig,ax,im = ImagePlotter(SymRotateImage,extent,aspectratio)

		#If the mesh does not use symmetry, simply plot the data as is.
		elif Isym == 0:
			extent=[0,Height, 0,Radius]
			fig,ax,im = ImagePlotter(RotateImage,extent,aspectratio)
		#endif

#=========================#

	#plot image as default mesh orientation.
	elif image_rotate == False:
		numrows = len(Image)
		SymImage = np.zeros([numrows,2*R_mesh])

		#If mesh uses symmetry, modify Image to conform to this.
		if Isym == 1:
			#Create new image by reversing and adding itself on the LHS.
			for m in range(0,len(Image)):
				SymImage[m] = np.concatenate([Image[m][::-1],Image[m]])
			#endfor

			extent = [-Radius,Radius, 0,Height]
			fig,ax,im = ImagePlotter(SymImage,extent,aspectratio)

		#If the mesh does not use symmetry, simply plot the data as is.
		elif Isym == 0:
			extent=[0,Radius, 0,Height]
			fig,ax,im = ImagePlotter(Image,extent,aspectratio)
		#endif
	#endif

	return(fig,ax,im)
#enddef



#=========================#
#=========================#



#Creates and plots a colourbar with given label and binsize.
def Colourbar(ax,Label,Bins):
	#Colourbar plotting details
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)
	#Number of ticks
	tick_locator = ticker.MaxNLocator(nbins=Bins)
	cbar.locator = tick_locator
	cbar.update_ticks()
	cbar.set_label(Label, rotation=270,labelpad=30,fontsize=24)
	#Size of font
	cbar.ax.yaxis.offsetText.set(size=18)
	yticks(fontsize=18)
#enddef



#=========================#
#=========================#


#Generates an SI axis for a 1D profile plot.
#Takes orientation and symmetry options.
#Returns 1D array in units of [cm].
def GenerateAxis(Orientation,Isym):
	#Generate SI scale axes for lineout plots.
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
	#endif
	return(axis)
#enddef



#=========================#
#=========================#



#Obtains a radial 1D profile at a requested axial location.
#Returns a 1D array for plotting and performs unit conversion.
def PlotRadialProfile(Data,process,variable,lineout,R_mesh,Isym):

	#Obtain start location for requested data and perform SI conversion.
	ZStart = R_mesh*lineout
	ZEnd = R_mesh*(lineout+1)

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
def PlotAxialProfile(Data,process,variable,lineout,R_mesh,Z_mesh,Isym):

	#Refresh lineout data between lines.
	Zlineout = list()

	#Pull out Z-data point from each radial line of data and list them.
	for k in range(0,Z_mesh):
		datapoint = R_mesh*k + lineout
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



#Locate powered electrode for bias extraction.
#Returns mesh location based on input string or [R,Z] list.
def ElectrodeLoc(location,data):

	#Check for string reference to powered electrode location.
	if isinstance(location[0], basestring) == True:
		#Upper Central Electrode.
		if location[0] == 'T':
			RlineoutLoc = 0
			ZlineoutLoc = R_mesh[l]/2
		#RHSCenterEdge electrode.
		elif location[0] == 'R':
			RlineoutLoc = Z_mesh[l]/2
			ZlineoutLoc = R_mesh
		#LHSCenterEdge electrode.
		elif location[0] == 'L':
			RlineoutLoc = Z_mesh[l]/2
			ZlineoutLoc = 0
		#Bottom Central Electrode.
		if location[0] == 'B':
			RlineoutLoc = Z_mesh[l]-1
			ZlineoutLoc = R_mesh[l]/2

	#If integer reference found: convert, else default location (CentralTopCell)
	elif isinstance(location[0], basestring) == False:
		#if data is from movie1, Zaxis is back to front.
		if len(location) == 2 and data == 'Phase':
			RlineoutLoc = (Z_mesh[l]-location[1])-1
			ZlineoutLoc = location[0]
		#if data is from TECPLOT2D, Zaxis is correct orientation.
		elif len(location) == 2 and data == '2D':
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




#Takes folder directory and creates a movie from .png images contained within.
def Automovie(FolderDir,Output):

	#Correct file extention on name.
	CurrentDir = os.getcwd()
	Output = Output+'.mp4'
	Morph, FPS = 1, 24

	#Use ffmpeg to create the movies and save in relevent files.
	os.chdir(FolderDir)
	os.system("convert *.png -delay 1 -morph "+str(Morph)+" %05d.morph.jpg > /dev/null")
	os.system("ffmpeg -nostats -loglevel 0 -r "+str(FPS)+" -i %05d.morph.jpg "+Output)
	os.system("rm *.jpg")
	os.chdir(CurrentDir)
	return()
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
	Xaxis = list()
	Trend = list()

	#For all simulation folders.
	for l in range(0,numfolders):

		#Extract image with given process and variable name.
		Image = ImageExtractor2D(Data[l][process],R_mesh[l],Z_mesh[l],variable)

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
		Image = Image.flatten()
		Trend = Trend/max(Image)
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
	Xaxis = list()
	MaxValueTrend = list()
	MinValueTrend = list()
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
	#endfor

	#Display Min/Max value trends to terminal if requested.
	if print_generaltrends == True:
		print VariableLabelMaker(Variablelist)[p]
		print Orientation+'Maximum: ', round(max(MaxValueTrend), 5)
		print Orientation+'Minimum: ', round(min(MinValueTrend), 5)
		print ''
	#endif

	#Normalize to maximum value in each profile if required.
	if image_normalize == True:
		NormalizeMax = max(MaxValueTrend)
		NormalizeMin = min(MinValueTrend)
		for i in range(0,len(MaxValueTrend)):
			MaxValueTrend[i] = MaxValueTrend[i]/NormalizeMax
			try: MinValueTrend[i] = MinValueTrend[i]/NormalizeMin
			except: MinValueTrend[i] = 0
		#endfor
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
		electrodelocation = electrodeloc[1]
	elif len(PPOTlineout) in [R_mesh[l],R_mesh[l]*2]:
		electrodelocation = electrodeloc[0]
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

	#Grounded metal will have a PPOT of zero -- ##INCORRECT IF DC-BIAS == 0.0##
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

	#Refresh percentage counters.
	New, Old = 0.0, 1.0

	for l in range(0,numfolders):
		#Create new folder to keep output plots.
		Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')

		#Create processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

		#Reshape specific part of 1D Data array into 2D image for plotting.
		for k in range(0,len(processlist)):

			#Extract full 2D image for further processing.
			Image = ImageExtractor2D(Data[l][processlist[k]],R_mesh[l],Z_mesh[l],Variablelist[k])

			if image_rotate == True:
				#Label and save the horizontal 2D Plots.
				fig,ax,im = SymmetryConverter2D(Image,Isymlist[l],R_mesh[l],Z_mesh[l],Radius[l],Height[l])

				#Image plotting details
				plt.title('2D Steady State Plot of '+Variablelist[k]+' for \n'+Dirlist[l][2:-1],y=1.10)
				plt.xlabel('Axial Distance Z [cm]',fontsize=24)
				plt.ylabel('Radial Distance R [cm]',fontsize=24)
				if image_plotsymmetry == True:
					plt.yticks([round(-Radius[l], 1), round(-Radius[l]/2, 1), 0, round(Radius[l]/2, 1), round(Radius[l], 1)])
				elif image_plotsymmetry == False:
					plt.yticks([0, round(Radius[l]/2, 1), round(Radius[l], 1)])
				#endif
				xticks(fontsize=18)
				yticks(fontsize=18)

			else:
				#Label and save the vertical 2D Plots.
				fig,ax,im = SymmetryConverter2D(Image,Isymlist[l],R_mesh[l],Z_mesh[l],Radius[l],Height[l])
				#Image plotting details, invert Y-axis to fit 1D profiles.
				plt.title('2D Steady State Plot of '+Variablelist[k]+' for \n'+Dirlist[l][2:-1],y=1.03)
				plt.xlabel('Radial Distance R [cm]',fontsize=24)
				plt.ylabel('Axial Distance Z [cm]',fontsize=24)
				plt.gca().invert_yaxis()
				if image_plotsymmetry == True:
					plt.xticks([round(-Radius[l], 1), round(-Radius[l]/2, 1), 0, round(Radius[l]/2, 1), round(Radius[l], 1)])
				elif image_plotsymmetry == False:
					plt.xticks([0, round(Radius[l]/2, 1), round(Radius[l], 1)])
				#endif
				xticks(fontsize=18)
				yticks(fontsize=18)
			#endif

			#Add Colourbar (Axis, Label, Bins)
			Bins = 5
			label = VariableLabelMaker(Variablelist)
			cax = Colourbar(ax,label[k],Bins)

			#Save Figure
			plt.savefig(Dir2Dplots+'2DPlot '+Variablelist[k]+'.png')
#			plt.show()
			plt.clf()
			plt.close('all')

			#Percentage Complete Readout.
			New += 1.0
			Old += 1.0
			Total = len(processlist)*numfolders
			oldpercentage = int( (New/Total )*100.0 )
			newpercentage = int( (Old/Total )*100.0 )
			if round(oldpercentage,-1) != round(newpercentage,-1):
				print int(round(oldpercentage,-1)),'%'
				a=1
			#endif
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

	#Refresh percentage counters and calculate total loop cycles for percentage.
	New, Old, Total = 0.0, 1.0, 0
	for l in range(0,numfolders):
		Total += len(IterVariables)*len(MovieITERlist[l])
	#endfor

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
		ConvergenceTrend = list()
		Xaxis = list()
		for k in range(0,len(MovieITERlist[l])):
			Xaxis.append(filter(lambda x: x.isdigit(), MovieITERlist[l][k]))
		#endfor

		#for all variables requested by the user.
		for i in range(0,len(iterprocesslist)):

			#Create new folder to keep output plots.
			DirMovieplots = CreateNewFolder(DirConvergence,IterVariablelist[i]+'_2DConvergence/')

			#Create empty image array based on mesh size and symmetry options.
			numrows = len(IterMovieData[l][0][0])/R_mesh[l]
			Image = np.zeros([numrows,R_mesh[l]])

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for k in range(0,len(MovieITERlist[l])):

				#Extract full 2D image for further processing.
				Image = ImageExtractor2D(IterMovieData[l][k][iterprocesslist[i]],R_mesh[l],Z_mesh[l],IterVariablelist[i])
				#Take Max value of image for general convergence trend.
				ConvergenceTrend.append( sum(Image.flatten())/len(Image.flatten()) )

				#Label and save the 2D Plots.
				if image_rotate == True:
					#Label and save the horizontal 2D Plots.
					fig,ax,im = SymmetryConverter2D(Image,Isymlist[l],R_mesh[l],Z_mesh[l],Radius[l],Height[l])

					#Image plotting details
					plt.title(MovieITERlist[l][k], y=1.10, fontsize=26)
					plt.xlabel('Axial Distance Z [cm]',fontsize=24)
					plt.ylabel('Radial Distance R [cm]',fontsize=24)
					if image_plotsymmetry == True:
						plt.yticks([round(-Radius[l], 1), round(-Radius[l]/2, 1), 0, round(Radius[l]/2, 1), round(Radius[l], 1)])
					elif image_plotsymmetry == False:
						plt.yticks([0, round(Radius[l]/2, 1), round(Radius[l], 1)])
					#endif
					xticks(fontsize=18)
					yticks(fontsize=18)

				else:
					#Label and save the vertical 2D Plots.
					fig,ax,im = SymmetryConverter2D(Image,Isymlist[l],R_mesh[l],Z_mesh[l],Radius[l],Height[l])
					#Image plotting details, invert Y-axis to fit 1D profiles.
					plt.title(MovieITERlist[l][k], y=1.10, fontsize=26)
					plt.xlabel('Radial Distance R [cm]',fontsize=24)
					plt.ylabel('Axial Distance Z [cm]',fontsize=24)
					plt.gca().invert_yaxis()
					if image_plotsymmetry == True:
						plt.xticks([round(-Radius[l], 1), round(-Radius[l]/2, 1), 0, round(Radius[l]/2, 1), round(Radius[l], 1)])
					elif image_plotsymmetry == False:
						plt.xticks([0, round(Radius[l]/2, 1), round(Radius[l], 1)])
					#endif
					xticks(fontsize=18)
					yticks(fontsize=18)
				#endif

				#Add Colourbar (Axis, Label, Bins)
				Bins = 5
				label = VariableLabelMaker(IterVariablelist)
				cax = Colourbar(ax,label[i],Bins)

				#Save to seperate folders inside simulation folder.
				num1,num2,num3 = k % 10, k/10 % 10, k/100 % 10
				Number = str(num3)+str(num2)+str(num1)
				savefig(DirMovieplots+IterVariablelist[i]+'_'+Number+'.png')
#				plt.show()
				plt.clf()
				plt.close('all')

				#Percentage Complete Readout.
				New += 1.0
				Old += 1.0
				oldpercentage = int( (New/Total)*100.0 )
				newpercentage = int( (Old/Total)*100.0 )
				if round(oldpercentage,-1) != round(newpercentage,-1):
					print int(round(newpercentage,-1)),'%'
				#endif
			#endfor

			#Normalize current variable in convergence trend.
			normalize = max(ConvergenceTrend[i*len(Xaxis):(i+1)*len(Xaxis)])
			for m in range(0,len(Xaxis)):
				index = (i*len(Xaxis))+m
				ConvergenceTrend[index] = ConvergenceTrend[index]/normalize
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+IterVariablelist[i])
		#endfor

		#Plot a convergence check for all variables in each folder.
		legend = VariableLabelMaker(IterVariablelist)
		fig, ax = plt.subplots(1, figsize=(10,10))

		for j in range(0,len(iterprocesslist)):
			start = j*len(Xaxis)
			end = (j+1)*len(Xaxis)
			plt.plot(Xaxis,ConvergenceTrend[start:end], lw=2)
		#endfor
		plt.title('Convergence of '+str(IterVariablelist), fontsize=16, y=1.03)
		plt.xlabel('Simulation Iteration', fontsize=24)
		plt.ylabel('Normalized Mesh-Average Value', fontsize=24)
		plt.legend(legend, prop={'size':16}, loc=4)
		plt.ylim((0,1.02))
		xticks(fontsize=18)
		yticks(fontsize=18)
		savefig(DirConvergence+FolderNameTrimmer(Dirlist[l])+'_Convergence.png')
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
if savefig_radialines or savefig_heightlines == True:

	for l in range(0,numfolders):

		#Create processlist for each folder as required.
		processlist,Variablelist = VariableEnumerator(Variables,rawdata_2D[l],header_2Dlist[l])

		#Generate SI scale axes for lineout plots and refresh legend.
		Raxis = GenerateAxis('Radial',Isymlist[l])
		Zaxis = GenerateAxis('Axial',Isymlist[l])
		Legendlist = list()

		#Generate the radial (horizontal) lineouts for a specific height.
		if savefig_radialines == True and len(radialineouts) > 0:
			#Create folder to keep output plots.
			DirRlineouts = CreateNewFolder(Dirlist[l],'Radial_Profiles/')

			#Loop over all required variables and requested profile locations.
			for i in range(0,len(processlist)):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(radialineouts)):
					#Update legend with location of each lineout.
					if len(Legendlist) < len(radialineouts):
						Legendlist.append('Z='+str(round((radialineouts[j])*dz[l]*10, 2))+' mm')
					#endif
					ylabel = VariableLabelMaker(Variablelist)

					#Plot all requested radial lines on single image per variable.
					Rlineout = PlotRadialProfile(Data[l],processlist[i],Variablelist[i],radialineouts[j],R_mesh[l],Isymlist[l])

					#Plot lines for each variable at each requested slice.
					if image_logplot == True:
						plt.plot(Raxis,np.log(Rlineout), lw=2)
						#ax.set_yscale('log')
					elif image_normalize == True:
						Rlineout = Normalize(Rlineout)
						plt.plot(Raxis,Rlineout, lw=2)
					else:
						plt.plot(Raxis,Rlineout, lw=2)
					#endif
				#endfor

				#Save lines in previously created folder.
				plt.title('Radial Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]+' Simulation',position=(0.5,1.05))
				plt.legend(Legendlist,loc=1)
				plt.xlabel('Radial Distance R [cm]', fontsize=24)
				plt.ylabel(ylabel[i], fontsize=24)
				ax.tick_params(axis='x', labelsize=18)
				ax.tick_params(axis='y', labelsize=18)
				if image_plotgrid == True: plt.grid(True)

				plt.savefig(DirRlineouts+'1D_Radial_'+Variablelist[i]+' profiles.png')
#				plt.show()
				plt.close('fig')
			#endfor
			plt.close('all')
		#endif

#===================##===================#

		#Generate the vertical (height) lineouts for a given radius.
		if savefig_heightlines == True and len(heightlineouts) > 0:
			#Create folder to keep output plots.
			DirZlineouts = CreateNewFolder(Dirlist[l],'Axial_Profiles/')
			Legendlist = list()

			#Collect and plot required data.
			for i in range(0,len(processlist)):
				#Create fig of desired size.
				fig,ax = figure(image_aspectratio,1)

				for j in range(0,len(heightlineouts)):

					#Plot all requested radial lines on single image per variable.
					Zlineout = PlotAxialProfile(Data[l],processlist[i],Variablelist[i],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])

					#Plot lines for each variable at each requested slice.
					if image_logplot == True:
						plt.plot(Zaxis[0:len(Zlineout)],np.log(Zlineout[::-1]), lw=2)
						#ax.set_yscale('log')
					elif image_normalize == True:
						Zlineout = Normalize(Zlineout)
						plt.plot(Zaxis[0:len(Zlineout)],Zlineout[::-1], lw=2)
					else:
						plt.plot(Zaxis[0:len(Zlineout)],Zlineout[::-1], lw=2)
					#endif

					#Perform SI conversion and save to legend.
					if len(Legendlist) < len(heightlineouts):
						Rlegend = heightlineouts[j]*dr[l]*10
						Legendlist.append('R='+str(round(Rlegend, 2))+' mm')
					#endif
					ylabel = VariableLabelMaker(Variablelist)
				#endfor

				#Save lines in previously created folder.
				plt.title('Height Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1]+' Simulation',position=(0.5,1.05))
				plt.xlabel('Axial Distance Z [cm]', fontsize=24)
				plt.ylabel(ylabel[i], fontsize=24)
				ax.tick_params(axis='x', labelsize=18)
				ax.tick_params(axis='y', labelsize=18)
				plt.legend(Legendlist,loc=1)
				if image_plotgrid == True: plt.grid(True)

				plt.savefig(DirZlineouts+'1D_Height_'+Variablelist[i]+' profiles.png')
#				plt.show()
				plt.close('fig')
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
		for k in range(0,len(Variables)):

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
				ylabel = VariableLabelMaker(Variablelist)

				#Plot all radial profiles for all variables in one folder.
				Rlineout = PlotRadialProfile(Data[l],processlist[k],Variablelist[k],radialineouts[j],R_mesh[l],Isymlist[l])

				#Plot radial profile and allow for log y-axis if requested.
				if image_logplot == True:
					plt.plot(Raxis,np.log(Rlineout), lw=2)
					#ax.set_yscale('log')
				elif image_normalize == True:
					Rlineout = Normalize(Rlineout)
					plt.plot(Raxis,Rlineout, lw=2)
				else:
					plt.plot(Raxis,Rlineout, lw=2)
				#endif

				#Organize and beautify the plots.
				plt.title('Comparison of '+Variablelist[k]+' Profiles at Z='+str(round((radialineouts[j])*dz[l], 2))+'cm for \n'+Dirlist[l][2:-1]+' Simulation',position=(0.5,1.05))
				plt.xlabel('Radial Distance R [cm]', fontsize=24)
				plt.ylabel(ylabel[k], fontsize=24)
				ax.tick_params(axis='x', labelsize=18)
				ax.tick_params(axis='y', labelsize=18)
				if image_plotgrid == True: plt.grid(True)
				if len(legendoverride) > 0:
					plt.legend(legendoverride, prop={'size':16}, loc=1)
				else:
					plt.legend(Legendlist, prop={'size':16}, loc=1)
				#endif
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ Z='+str(round((radialineouts[j])*dz[l], 2))+'cm profiles.png')
	#		plt.show()
			plt.clf
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
		for k in range(0,len(Variables)):

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
				ylabel = VariableLabelMaker(Variablelist)

				#Obtain axial profile for each folder of the current variable.
				Zlineout = PlotAxialProfile(Data[l],processlist[k],Variablelist[k],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])

				#Plot axial profile and allow for log y-axis if requested.
				if image_logplot == True:
					plt.plot(Zaxis[0:len(Zlineout)],np.log(Zlineout[::-1]), lw=2)
					#ax.set_yscale('log')
				elif image_normalize == True:
					Zlineout = Normalize(Zlineout)
					plt.plot(Raxis[0:len(Zlineout)],Zlineout[::-1], lw=2)
				else:
					plt.plot(Zaxis[0:len(Zlineout)],Zlineout[::-1], lw=2)
				#endif

				#Beautify and label the plots, allowing for manual renaming if requested.
				plt.title('Comparison of '+Variablelist[k]+' Profiles at R='+str(round((heightlineouts[j])*dr[l], 2))+'cm for \n'+Dirlist[l][2:-1]+' Simulation',position=(0.5,1.05))
				plt.xlabel('Axial Distance Z [cm]', fontsize=22)
				plt.ylabel(ylabel[k], fontsize=22)
				ax.tick_params(axis='x', labelsize=18)
				ax.tick_params(axis='y', labelsize=18)
				if image_plotgrid == True: plt.grid(True)
				if len(legendoverride) > 0:
					plt.legend(legendoverride, prop={'size':16}, loc=1)
				else:
					plt.legend(Legendlist, prop={'size':16}, loc=1)
				#endif
			#endfor

			#Save one image per variable with data from all simulations.
			plt.savefig(DirProfile+Variablelist[k]+'@ R='+str(round((heightlineouts[j])*dr[l], 2))+'cm profiles.png')
	#		plt.show()
			plt.clf
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
		multiprocesslist, multiVariablelist = VariableEnumerator(MultiVar,rawdata_2D[l],header_2Dlist[l])

		#Create variable labels with SI unit conversions if required.
		ylabel = VariableLabelMaker(Variablelist)
		multiylabel = VariableLabelMaker(multiVariablelist)

		#Generate the vertical (height) lineouts for a given radius.
		if len(heightlineouts) > 0:

			#Generate SI scale axes for lineout plots.
			Zaxis = GenerateAxis('Axial',Isymlist[l])

			#Perform the plotting for all requested variables.
			for i in range(0,len(processlist)):

				#Extract the lineout data from the main data array.
				for j in range(0,len(heightlineouts)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((heightlineouts[j])*dr[l], 2))
					DirZlineouts = CreateNewFolder(Dirlineouts,'R='+Slice+'cm/')

					#Create legendlist
					legendlist = list()
					legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot the initial variable in processlist first.
					Zlineout = PlotAxialProfile(Data[l],processlist[i],Variablelist[i],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])
					if image_normalize == True: Zlineout = Normalize(Zlineout)
					plt.plot(Zaxis[0:len(Zlineout)],Zlineout[::-1], lw=2)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):

						#Plot profile for multiplot variables in compareprocesslist.
						Zlineout = PlotAxialProfile(Data[l],multiprocesslist[m],multiVariablelist[m],heightlineouts[j],R_mesh[l],Z_mesh[l],Isymlist[l])
						if image_logplot == True:
							plt.plot(Zaxis,np.log(Zlineout), lw=2)
							#ax.set_yscale('log')
						elif image_normalize == True:
							Zlineout = Normalize(Zlineout)
							plt.plot(Zaxis,Zlineout[::-1], lw=2)
						else:
							plt.plot(Zaxis,Zlineout[::-1], lw=2)
						#endif

						#Update legendlist with each variable compared.
						legendlist.append(VariableLabelMaker(multiVariablelist)[m])
					#endfor

					#Save figures in original folder.
					plt.title(str(round((heightlineouts[j])*dr[l], 2))+'cm Height profiles for '+Variablelist[i]+','' for \n'+Dirlist[l][2:-1],position=(0.5,1.05))
					plt.legend(legendlist, prop={'size':14}, loc=1)
					ax.tick_params(axis='x', labelsize=18)
					ax.tick_params(axis='y', labelsize=18)
					plt.xlabel('Axial Distance Z [cm]', fontsize=22)
					plt.ylabel(ylabel[i], fontsize=22)
					if image_plotgrid == True: plt.grid(True)

					R = 'R='+str(round((heightlineouts[j])*dr[l], 2))+'_'
					plt.savefig(DirZlineouts+R+Variablelist[i]+'_MultiProfiles.png')
					plt.close('fig')
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
			for i in range(0,len(processlist)):

				#Perform the plotting for all requested variables.
				for j in range(0,len(radialineouts)):
					#Create fig of desired size.
					fig,ax = figure(image_aspectratio,1)

					#Create folder to keep output plots.
					Slice = str(round((radialineouts[j])*dz[l], 2))
					DirRlineouts = CreateNewFolder(Dirlineouts,'Z='+Slice+'cm/')

					#Create legendlist
					legendlist = list()
					legendlist.append(VariableLabelMaker(Variablelist)[i])

					#Plot profile for initial variable in processlist.
					Rlineout = PlotRadialProfile(Data[l],processlist[i],Variablelist[i],radialineouts[j],R_mesh[l],Isymlist[l])
					if image_normalize == True: Rlineout = Normalize(Rlineout)
					plt.plot(Raxis,Rlineout,lw=2)

					#Plot all of the requested comparison variables for this plot.
					for m in range(0,len(multiprocesslist)):

						#Plot profile for multiplot variables in compareprocesslist.
						Rlineout = PlotRadialProfile(Data[l],multiprocesslist[m],multiVariablelist[m],radialineouts[j],R_mesh[l],Isymlist[l])
						if image_logplot == True:
							plt.plot(Raxis,np.log(Rlineout), lw=2)
							#ax.set_yscale('log')
						elif image_normalize == True:
							Rlineout = Normalize(Rlineout)
							plt.plot(Raxis,Rlineout, lw=2)
						else:
							plt.plot(Raxis,Rlineout, lw=2)
						#endif

						#Update legendlist with each variable compared.
						legendlist.append(VariableLabelMaker(multiVariablelist)[m])
					#endfor

					#Save lines in previously created folder.
					plt.title(str(round((radialineouts[j])*dz[l], 2))+'cm Radial Profiles for '+Variablelist[i]+' for \n'+Dirlist[l][2:-1],position=(0.5,1.05))
					plt.legend(legendlist, prop={'size':14}, loc=1)
					ax.tick_params(axis='x', labelsize=18)
					ax.tick_params(axis='y', labelsize=18)
					plt.xlabel('Radial Distance R [cm]', fontsize=22)
					plt.ylabel(ylabel[i], fontsize=22)
					if image_plotgrid == True: plt.grid(True)

					Z = 'Z='+str(round((radialineouts[j])*dz[l], 2))+'_'
					plt.savefig(DirRlineouts+Z+Variablelist[i]+'_MultiProfiles.png')
					plt.close('fig')
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
	for k in range(0,len(Variables)):

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
				Legendlist.append( 'R='+str(round((heightlineouts[j]*dr[l]*10), 2))+'mm' )
			#endif

			#Plot trends for each variable over all folders, applying image options.
			TrendPlotter(Trend,Xaxis,Normalize=1)

			#Beautify the plot before saving.
			plt.title('Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
			plt.ylabel('Max '+YaxisLegend[k], fontsize=24)
			ax.tick_params(axis='x', labelsize=18)
			ax.tick_params(axis='y', labelsize=18)
			plt.legend(Legendlist, prop={'size':16}, loc=1)
			if len(xlabeloverride) > 0:
				plt.xlabel(xlabeloverride[0], fontsize=24)
			else:
				plt.xlabel('Varied Property', fontsize=24)
			#endif

		#Save one image per variable with data from all simulations.
		if len(heightlineouts) > 0:
			plt.savefig(DirAxialTrends+'Axial Trends in '+Variablelist[k]+'.png')
#			plt.show()
			plt.clf
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
				Legendlist.append( 'Z='+str(round((radialineouts[j]*dz[l]*10), 2))+'mm' )
			#endif

			#Plot trends for each variable over all folders, applying image options.
			TrendPlotter(MaxTrend,Xaxis,Normalize=1)

			#Beautify the plot before saving.
			plt.title('Trend in max '+Variablelist[k]+' with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
			plt.ylabel('Max '+YaxisLegend[k], fontsize=24)
			ax.tick_params(axis='x', labelsize=18)
			ax.tick_params(axis='y', labelsize=18)
			plt.legend(Legendlist, prop={'size':16}, loc=1)
			if len(xlabeloverride) > 0:
				plt.xlabel(xlabeloverride[0], fontsize=24)
			else:
				plt.xlabel('Varied Property', fontsize=24)
			#endif

		#Save one image per variable with data from all simulations.
		if len(radialineouts) > 0:
			plt.savefig(DirRadialTrends+'Radial Trends in '+Variablelist[k]+'.png')
#			plt.show()
			plt.clf
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
		Rlineoutloc = ElectrodeLoc(electrodeloc,'2D')[0]
		Zlineoutloc = ElectrodeLoc(electrodeloc,'2D')[1]

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

	#Plot and beautify the DCbias, applying normalization if requested.
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(DCbias,Xaxis,Normalize=520)

	plt.title('Trend in DCbias with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
	plt.ylabel('DC bias [V]', fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	if len(xlabeloverride) > 0:
		plt.xlabel(xlabeloverride[0], fontsize=24)
	else:
		plt.xlabel('Varied Property', fontsize=24)
	#endif

	plt.savefig(DirTrends+'Powered Electrode DCbias.png')
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
	RequestedPowers = list()
	DepositedPowerList = list()
	Xaxis = list()

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

			#Update X-axis with folder information.
			Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

			#Extract full 2D power density image.
			PowerDensity = ImageExtractor2D(Data[l][processlist[k]],R_mesh[l],Z_mesh[l])

			#For power densities, convert from [Wcm-3] to [Wm-3].
			for i in range(0,len(PowerDensity)):
				PowerDensity[i] = PowerDensity[i]*1E6
			#endfor

			Power = 0
			#Integrates power per unit volume and produces total coupled power.
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
		Power = DepositedPowerList[k*numfolders:(k+1)*numfolders]
		TrendPlotter(Power,Xaxis,Normalize=1)

		plt.title('Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
		plt.ylabel('RF-Power Deposited [W]', fontsize=24)
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)
		if len(xlabeloverride) > 0:
			plt.xlabel(xlabeloverride[0], fontsize=24)
		else:
			plt.xlabel('Varied Property', fontsize=24)
		#endif

		plt.savefig(DirTrends+RequestedPowers[k]+' Deposition Trends.png')
		plt.close('all')
	#endfor

	#Plot a comparison of all power depositions requested.
	fig,ax = figure(image_aspectratio,1)
	for k in range(0,len(RequestedPowers)):
		Power = DepositedPowerList[k*numfolders:(k+1)*numfolders]
		TrendPlotter(Power,Xaxis,Normalize=1)
#		plt.show()
	#endfor

	plt.title('Power Deposition with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
	plt.ylabel('Power Deposited [W]', fontsize=24)
	plt.legend(RequestedPowers)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	if len(xlabeloverride) > 0:
		plt.xlabel(xlabeloverride[0], fontsize=24)
	else:
		plt.xlabel('Varied Property', fontsize=24)
	#endif


	plt.savefig(DirTrends+'Power Deposition Comparison.png')
	plt.close('all')
#endif



#====================================================================#
				  	#SIMPLE THRUST ANALYSIS#
#====================================================================#



if savefig_trendcomparison == True or print_thrust == True:

	#Create Trend folder to keep output plots.
	TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
	DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

	#Initiate lists required for storing data.
	Thrustlist = list()
	Xaxis = list()

	#For all folders.
	for l in range(0,numfolders):

		#Create extract data for the neutral flux and neutral velocity.
		processlist,Variablelist = VariableEnumerator(['AR','VZ-NEUTRAL','VZ-ION+','FZ-AR','FZ-AR+','PRESSURE'],rawdata_2D[l],header_2Dlist[l])

		#Update X-axis with folder information.
		Xaxis.append( FolderNameTrimmer(Dirlist[l]) )

		#Extract radial density, velocity and pressure profiles across the discharge plane.
		AbortDiagnostic = False
		try:
			Density = PlotRadialProfile(Data[l],processlist[0],Variablelist[0],AxialLine,R_mesh[l],Isymlist[l])
			NeutralVelocity = PlotRadialProfile(Data[l],processlist[1],Variablelist[1],AxialLine,R_mesh[l],Isymlist[l])
			IonVelocity = PlotRadialProfile(Data[l],processlist[2],Variablelist[2],AxialLine,R_mesh[l],Isymlist[l])
			NeutralAxialFlux = PlotRadialProfile(Data[l],processlist[3],Variablelist[3],AxialLine, R_mesh[l],Isymlist[l])
			IonAxialFlux = PlotRadialProfile(Data[l],processlist[4],Variablelist[4],AxialLine,R_mesh[l],Isymlist[l])
			Pressure = PlotRadialProfile(Data[l],processlist[5],Variablelist[5],AxialLine,R_mesh[l],Isymlist[l])
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
		Argon,Xenon = 39.948,131.29			 #amu
		NeutralMass = Argon*1.67E-27		 #Kg
		NeutralIsp = list()
		IonIsp = list()

		DefaultTechnique = True
		if DefaultTechnique == True:

			#Thrust based on integration over concentric ion/neutral momentum loss rate.
			NeutralThrust = 0
			IonThrust = 0
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
			Thrust = NeutralThrust + IonThrust
			Thrustlist.append( round(Thrust*1000,5) )
			try: IonIsp = (sum(IonIsp)/len(IonIsp))/9.81
			except: IonIsp = 0
			try: NeutralIsp = (sum(NeutralIsp)/len(NeutralIsp))/9.81
			except: NeutralIsp = 0

		else:
			#Thrust based on integration over concentric neutral momentum and pressure.
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
			print 'NeutralThrust', round(NeutralThrust*1000,2), 'mN @ ', round(NeutralIsp,2),'s-1'
			print 'IonThrust:', round(IonThrust*1000,4), 'mN @ ', round(IonIsp,2),'s-1'
			print 'Thrust:',round(Thrust*1000,4),'mN'
			print ''
		#endif
	#endfor

	#Plot and Beautify the thrust.
	fig,ax = figure(image_aspectratio,1)
	TrendPlotter(Thrustlist,Xaxis,Normalize=1)

	plt.title('Thrust with changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
	plt.ylabel('Thrust [mN]', fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	if len(xlabeloverride) > 0:
		plt.xlabel(xlabeloverride[0], fontsize=24)
	else:
		plt.xlabel('Varied Property', fontsize=24)
	#endif

	plt.savefig(DirTrends+'Thrust Trends.png')
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
			TrendPlottingOptions1D(ThrustEfficiency,Xaxis,Normalize=1)

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
	if savefig_trendcomparison == True or print_KnudsenNumber == True:

		#Create Trend folder to keep output plots.
		TrendVariable = filter(lambda x: x.isalpha(), FolderNameTrimmer(Dirlist[0]))
		DirTrends = CreateNewFolder(os.getcwd()+'/',TrendVariable+' Trends')

		#Initiate lists required for storing data.
		KnudsenAverage = list()
		Xaxis = list()

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
			Knudsen = np.zeros([Z_mesh[l],R_mesh[l]])

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
					Knudsen[Row,i] = KnudsenNumber
				#endfor
			#endfor

			#Create new folder to keep 2D output plots.
			Dir2Dplots = CreateNewFolder(Dirlist[l],'2Dplots')

			#Display average Knudsen number to terminal if requested.
			KnudsenAverage.append( sum(Knudsen)/(len(Knudsen[0])*len(Knudsen)) )
			if print_KnudsenNumber == True:
				print Dirlist[l]
				print 'Average Knudsen Number:', KnudsenAverage[l]
			#endif

			#Label and save the 2D Plots.
			fig,ax,im = SymmetryConverter2D(Knudsen,Isymlist[l],R_mesh[l],Z_mesh[l],Radius[l],Height[l])
			#Image plotting details
			plt.title('Knudsen Number Image for \n'+Dirlist[l][2:-1],y=1.03)
			plt.xlabel('Radial Distance R [cm]',fontsize=24)
			plt.ylabel('Axial Distance Z [cm]',fontsize=24)
			plt.gca().invert_yaxis()
			xticks(fontsize=18)
			yticks(fontsize=18)

			#Add Colourbar (Axis, Label, Bins)
			Bins = 5
			cax = Colourbar(ax,'Knudsen Number',Bins)

			#Save Figure
			plt.savefig(Dir2Dplots+'KnudsenNumber.png')
			plt.clf()
			plt.close('all')
		#endfor

		#Plot a comparison of all average Knudsen numbers.
		fig,ax = figure(image_aspectratio,1)
		TrendPlotter(KnudsenAverage,Xaxis,Normalize=1)

		plt.title('Average Knudsen Number with Changing '+TrendVariable+' \n'+Dirlist[l][2:-1] ,position=(0.5,1.05))
		plt.ylabel('Average Knudsen Number', fontsize=24)
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)
		if len(xlabeloverride) > 0:
			plt.xlabel(xlabeloverride[0], fontsize=24)
		else:
			plt.xlabel('Varied Property', fontsize=24)
		#endif

		plt.savefig(DirTrends+'KnudsenNumber Comparison.png')
		plt.close('all')

	#endif
#endif

#===============================#

if any([savefig_trendcomparison, print_generaltrends, print_KnudsenNumber, print_totalpower, print_DCbias, print_thrust]) == True:
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

	#Initialize required lists and refresh percentage counters.
	WaveformTitles = list()
	VoltageWaveform = list()
	New = 0.0
	Old = 1.0

	#for all folders being processed.
	for l in range(0,numfolders):

		#Create global folder to keep output plots.
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'2DPhase/')

		#Create processlist for each folder as required. (Always get PPOT)
		PhaseProcesslist,PhaseVariablelist = VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
		PPOT = VariableEnumerator(['PPOT'],rawdata_phasemovie[l],header_phasemovie[l])[0]

		#Subtract 2 from process as variables R&Z are not saved properly in phasedata.
		for i in range(0,len(PhaseProcesslist)): PhaseProcesslist[i] -= 2
		PPOT = PPOT[0]-2

		#Obtain applied voltage waveform, Refresh list between folders if needed.
		if len(VoltageWaveform) > 0: VoltageWaveform = list()
		for j in range(0,phasecycles):
			for i in range(0,len(Moviephaselist[l])):
				RElectrodeLoc = ElectrodeLoc(electrodeloc,'Phase')[0]
				ZElectrodeLoc = ElectrodeLoc(electrodeloc,'Phase')[1]
				VoltageWaveform.append( PlotAxialProfile(PhaseMovieData[l][i], PPOT,'PPOT',ZElectrodeLoc,R_mesh[l],Z_mesh[l],Isymlist[l])[RElectrodeLoc])
			#endfor
		#endfor

		#Generate a phase axis of units [omega*t/2pi] for plotting.
		Phaseaxis = list()
		for i in range(0,phasecycles*len(Moviephaselist[l])):
			Phaseaxis.append(  (np.pi*(i*2)/180)/(2*np.pi)  )
		#endfor

		#Generate SI scale axes.
		Raxis,Zaxis = list(),list()
		for i in range(0,Z_mesh[l]): Zaxis.append(i*dz[l])
 		for i in range(0,R_mesh[l]): Raxis.append(i*dr[l])
		#endfor

		#Calculate time averaged waveform bias, i.e. waveform symmetry.
		WaveformBias = list()
		for m in range(0,len(VoltageWaveform)):
			WaveformBias.append(sum(VoltageWaveform)/len(VoltageWaveform))
		#endfor

		#Extend the waveform to match requested number of phase cycles.
		WaveformTitles.append( FolderNameTrimmer(Dirlist[l]) )
		for m in range(0,(phasecycles-1)*len(VoltageWaveform)):
			VoltageWaveform.append(VoltageWaveform[m])
			WaveformBias.append(WaveformBias[m])
		#endfor

		#Plot the phase-resolved waveform.
		fig, ax = plt.subplots(1, figsize=(10,10))
		plt.plot(Phaseaxis,VoltageWaveform, lw=2)
		plt.plot(Phaseaxis,WaveformBias, 'k--', lw=2)
		plt.title('Phase-Resolved Voltage Waveform for '+WaveformTitles[l], y=1.03, fontsize=16)
		plt.xlabel('Phase [$\omega$t/2$\pi$]', fontsize=24)
		plt.ylabel('Potential [V]', fontsize=24)
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)

		plt.savefig(DirPhaseResolved+WaveformTitles[l]+' Waveform.png')
		plt.close('all')

		#===============#

		#for all variables requested by the user.
		for i in range(0,len(PhaseProcesslist)):

			#Create new folder to keep specific plots.
			DirMovieplots = CreateNewFolder(DirPhaseResolved,PhaseVariablelist[i]+'_2DPhaseResolved/')

			#Obtain maximum and minimum values of current variable over all phases.
			MaxNormalize,MinNormalize = list(),list()
			for j in range(0,phasecycles):
				MaxNormalize.append( max(ImageExtractor2D(PhaseMovieData[l][j][PhaseProcesslist[i]],R_mesh[l],Z_mesh[l],PhaseVariablelist[i]).flatten()) )
				MinNormalize.append( min(ImageExtractor2D(PhaseMovieData[l][j][PhaseProcesslist[i]],R_mesh[l],Z_mesh[l],PhaseVariablelist[i]).flatten()) )
			#endfor
			MaxNormalize = max(MaxNormalize)
			MinNormalize = min(MinNormalize)

			#Reshape specific part of 1D Data array into 2D image for plotting.
			for j in range(0,len(Moviephaselist[l])):

				#Extract full 2D image for further processing.  [ASSUMES SYMMETRIC IMAGE]
				Image = ImageExtractor2D(PhaseMovieData[l][j][PhaseProcesslist[i]],R_mesh[l],Z_mesh[l],PhaseVariablelist[i])
				#Create new image by reversing and adding itself on the LHS.
				SymImage = np.zeros([len(Image),2*R_mesh[l]])
				for m in range(0,len(Image)): SymImage[m] = np.concatenate([Image[m][::-1],Image[m]])
				Image = SymImage

				#Define image size and axis label names.
				extent = [-Raxis[-1],Raxis[-1],Zaxis[0],Zaxis[-1]]
				xlabel = 'Radial Distance R [cm]'
				ylabel = 'Axial Distance Z [cm]'

				#Ensure image is landscape for waveform to fit underneath.
				height,width = len(Image), len(Image[0])
				if height > width:
					extent=[Zaxis[0],Zaxis[-1],-Raxis[-1],Raxis[-1]]
					xlabel,ylabel = ylabel,xlabel
					Image = np.transpose(Image)
				#endif

				#Create figure and axes, plot image on top and waveform underneath.
				fig,ax = figure(image_aspectratio,2)
				fig.suptitle( 'Phase-Resolved '+PhaseVariablelist[i]+'\n'+str(Moviephaselist[l][j]), y=0.97, fontsize=18)

				if image_contourplot == True:
					im = ax[0].contour(Image,extent=extent,origin="lower", aspect='auto')
					im = ax[0].imshow(Image,extent=extent,origin="lower", aspect='auto')
				else:
					im = ax[0].imshow(Image,extent=extent,origin="lower", aspect='auto')
				#endif
				ax[0].set_xlabel('Axial Distance Z [cm]', fontsize=24)
				ax[0].set_ylabel('Radial Distance Z [cm]', fontsize=24)
				ax[0].tick_params(axis='x', labelsize=18)
				ax[0].tick_params(axis='y', labelsize=18)
				#Add Colourbar (Axis, Label, Bins)
				label = VariableLabelMaker(PhaseVariablelist)
				cax = Colourbar(ax[0],label[i],5)
				im.set_clim(vmin=MinNormalize, vmax=MaxNormalize)

				ax[1].plot(Phaseaxis, VoltageWaveform, lw=2)
				ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
				ax[1].set_xlabel('Phase [$\omega$t/2$\pi$]', fontsize=24)
				ax[1].set_ylabel('Potential [V]', fontsize=24)
				ax[1].tick_params(axis='x', labelsize=18)
				ax[1].tick_params(axis='y', labelsize=18)
				ax[1].set_xlim(0,phasecycles)


				#Calculate numbers used for filename ordering.
				num1,num2,num3 = j % 10, j/10 % 10, j/100 % 10
				Number = str(num3)+str(num2)+str(num1)

				#Cleanup layout and save images.
				fig.tight_layout()
				plt.subplots_adjust(top=0.90)
				savefig(DirMovieplots+PhaseVariablelist[i]+'_'+Number+'.png')
				plt.close('all')

				#Percentage Complete Readout.
				New += 1.0
				Old += 1.0
				Total = len(PhaseProcesslist)*numfolders*len(Moviephaselist[l])
				oldpercentage = int( (New/Total)*100.0 )
				newpercentage = int( (Old/Total)*100.0 )
				if round(oldpercentage,-1) != round(newpercentage,-1):
					print int(round(newpercentage,-1)),'%'
				#endif
			#endfor

			#Create .mp4 movie from completed images.
			Prefix = FolderNameTrimmer(Dirlist[l])
			Automovie(DirMovieplots,Prefix+'_'+PhaseVariablelist[i])
		#endfor
	#endfor

	print'---------------------------------------'
	print'# 2D Phase-Resolved Processing Complete'
	print'---------------------------------------'
#endif




#====================================================================#
			#PHASE RESOLVED PROFILES & SHEATH DYNAMICS (PROES)#
#====================================================================#

#Plot Phase-Resolved profiles with electrode voltage and requested variables.
if savefig_phaseresolvelines == True or savefig_sheathdynamics == True:

	#Refresh percentage counters.
	Total = numfolders*len(PhaseVariables)*(len(heightlineouts)+len(radialineouts))
	New = 0.0
	Old = 1.0
	#Initiate any required lists.
	VariedValuelist = list()

	#for all folders.
	for l in range(0,numfolders):

		#Create folders to keep output plots.
		DirPhaseResolved = CreateNewFolder(Dirlist[l],'PhaseResolved_profiles/')

		#Update X-axis with folder information.
		VariedValuelist.append( FolderNameTrimmer(Dirlist[l]) )

		#Create processlist for each folder as required. (Always get PPOT)
		PhaseProcesslist,PhaseVariablelist = VariableEnumerator(PhaseVariables,rawdata_phasemovie[l],header_phasemovie[l])
		PPOT = VariableEnumerator(['PPOT'],rawdata_phasemovie[l],header_phasemovie[l])[0]

		#Subtract 2 from process as variables R&Z are not saved properly in phasedata.
		for i in range(0,len(PhaseProcesslist)): PhaseProcesslist[i] -= 2
		PPOT = PPOT[0]-2

		#Generate a phase axis of units [omega*t/2pi] for plotting.
		Phaseaxis = list()
		for i in range(0,phasecycles*len(Moviephaselist[l])):
			Phaseaxis.append(  (np.pi*(i*2)/180)/(2*np.pi)  )
		#endfor

		#Generate SI scale axes for lineout plots.
		Zaxis = list()
		Raxis = list()
		for i in range(0,Z_mesh[l]): Zaxis.append(i*dz[l])
		if Isymlist[l] == 1:
			for i in range(-R_mesh[l],R_mesh[l]):
				Raxis.append(i*dr[l])
			#endfor
		elif Isymlist[l] == 0:
			for i in range(0,R_mesh[l]):
				Raxis.append(i*dr[l])
			#endfor
		#endif

		#==============#

		VoltageWaveform = list()
		#Obtain applied voltage waveform and normalization values seperately.
		for j in range(0,len(Moviephaselist[l])):

			#Obtain applied voltage waveform.
			RElectrodeLoc = ElectrodeLoc(electrodeloc,'Phase')[0]
			ZElectrodeLoc = ElectrodeLoc(electrodeloc,'Phase')[1]
			VoltageWaveform.append( PlotAxialProfile(PhaseMovieData[l][j], PPOT,'PPOT',ZElectrodeLoc,R_mesh[l],Z_mesh[l],Isymlist[l])[RElectrodeLoc])
		#endfor

		#Calculate time averaged waveform bias, i.e. waveform symmetry.
		WaveformBias = list()
		for m in range(0,len(VoltageWaveform)):
			WaveformBias.append(sum(VoltageWaveform)/len(VoltageWaveform))
		#endfor

		#Extend the waveform to match requested number of phase cycles.
		for m in range(0,(phasecycles-1)*len(VoltageWaveform)):
			VoltageWaveform.append(VoltageWaveform[m])
			WaveformBias.append(WaveformBias[m])
		#endfor

		#Plot the phase-resolved waveform.
		fig,ax = figure(image_aspectratio,1)
		plt.plot(Phaseaxis,VoltageWaveform, lw=2)
		plt.plot(Phaseaxis,WaveformBias, 'k--', lw=2)
		plt.title('Phase-Resolved Voltage Waveform for '+VariedValuelist[l], y=1.03, fontsize=16)
		plt.xlabel('Phase [$\omega$t/2$\pi$]', fontsize=24)
		plt.ylabel('Potential [V]', fontsize=24)
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)

		plt.savefig(DirPhaseResolved+VariedValuelist[l]+' Waveform.png')
		plt.close('all')

		#==============#


		#for all requested variables.
		for i in range(0,len(PhaseProcesslist)):

			#Refresh lineout lists between variables.
			LineoutsOrientation = list()
			Lineouts = list()

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
				VariableMax = list()
				VariableMin = list()
				PROES = list()

				#Create folders to keep output plots for each variable.
				if LineoutsOrientation[k] == 'Axial':
					NameString= PhaseVariablelist[i]+'_'+str(round(Lineouts[k]*dz[l],2))+'cm[Z]'
				elif LineoutsOrientation[k] == 'Radial':
					NameString= PhaseVariablelist[i]+'_'+str(round(Lineouts[k]*dr[l],2))+'cm[R]'
				Dir1DProfiles = CreateNewFolder(DirPhaseResolved,NameString+'_1Dprofiles/')

				#Collect Normalization data for plotting.
				for j in range(0,len(Moviephaselist[l])):
					#Record local maximum and minimum for each phase.
					if LineoutsOrientation[k] == 'Axial':
						ZlineoutLoc = Lineouts[k]
						VariableMax.append( max(PlotAxialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])) )
						VariableMin.append( min(PlotAxialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])) )
					elif LineoutsOrientation[k] == 'Radial':
						RlineoutLoc = Lineouts[k]
						VariableMax.append( max(PlotRadialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],RlineoutLoc,R_mesh[l],Isymlist[l])) )
						VariableMin.append( min(PlotRadialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],RlineoutLoc,R_mesh[l],Isymlist[l])) )
					#endif
				#endfor
				#Find global maximum and minimum for all phases.
				VariableMax = max(VariableMax)
				VariableMin = min(VariableMin)


				#for all recorded phases, plot spatially varying variable and waveform.
				for j in range(0,len(Moviephaselist[l])):

					if LineoutsOrientation[k] == 'Axial':
						#Axial 1D profiles, using first hightlineout as axial location.
						ZlineoutLoc = Lineouts[k]
						PhaseResolvedlineout = PlotAxialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],ZlineoutLoc,R_mesh[l],Z_mesh[l],Isymlist[l])
						lineoutstring = ' @ Z='+str(round(ZlineoutLoc*dz[l],2))+'cm \n'
						xlabel = 'Axial Distance Z [cm]'
						axis = Zaxis

					elif LineoutsOrientation[k] == 'Radial':
						#Radial 1D profiles, using first radialineout as axial location.
						RlineoutLoc = Lineouts[k]
						PhaseResolvedlineout = PlotRadialProfile(PhaseMovieData[l][j],PhaseProcesslist[i],PhaseVariablelist[i],RlineoutLoc,R_mesh[l],Isymlist[l])
						lineoutstring = ' @ R='+str(round(RlineoutLoc*dr[l],2))+'cm \n'
						xlabel = 'Radial Distance R [cm]'
						axis = Raxis
					#endif


					#Create figures and plot the 1D profiles. (ax[0]=variable, ax[1]=waveform)
					if savefig_phaseresolvelines == True:
						fig,ax = figure(image_aspectratio,2)
						ylabels = VariableLabelMaker(PhaseVariablelist)
						fig.suptitle( 'Phase-Resolved '+PhaseVariablelist[i]+' for '+VariedValuelist[l]+lineoutstring+str(Moviephaselist[l][j]), y=0.97, fontsize=16)

						ax[0].plot(axis, PhaseResolvedlineout[::-1], lw=2)
						ax[0].set_xlabel(xlabel, fontsize=24)
						ax[0].set_ylabel(ylabels[i], fontsize=24)
						ax[0].tick_params(axis='x', labelsize=18)
						ax[0].tick_params(axis='y', labelsize=18)
						ax[0].set_ylim(VariableMin,VariableMax*1.02)

						ax[1].plot(Phaseaxis, VoltageWaveform, lw=2)
						ax[1].axvline(Phaseaxis[j], color='k', linestyle='--', lw=2)
						ax[1].set_xlabel('Phase [$\omega$t/2$\pi$]', fontsize=24)
						ax[1].set_ylabel('Potential [V]', fontsize=24)
						ax[1].tick_params(axis='x', labelsize=18)
						ax[1].tick_params(axis='y', labelsize=18)
						ax[1].set_xlim(0,phasecycles)

						#Calculate numbers used for filename ordering.
						num1,num2,num3 = j % 10, j/10 % 10, j/100 % 10
						Number = str(num3)+str(num2)+str(num1)

						fig.tight_layout()
						plt.subplots_adjust(top=0.90)
						plt.savefig(Dir1DProfiles+NameString+'_'+Number+'.png')
						plt.close('all')
					#endif

					#Collect each profile for stitching into a PROES image if required.
					if savefig_sheathdynamics == True:
						PROES.append(PhaseResolvedlineout)
					#endif
				#endfor


				#==============#
				# 	  PROES    #
				#==============#

				#plot more detailed PROES sheath-dynamics if requested by user.
				if savefig_sheathdynamics == True:

					#Increase PROES image by required number of phasecycles.
					for m in range(1,phasecycles):
						for n in range(0,len(PROES)):
							PROES.append(PROES[n])
						#endfor
					#endfor

					#Create figure and rotate PROES such that phaseaxis aligns with waveform.
					fig,ax = figure(image_aspectratio,2)
					PROES = ndimage.rotate(PROES, 90)
					#Choose correct axial or radial distance axis.
					x1,x2 = Phaseaxis[0],Phaseaxis[-1]
					if LineoutsOrientation[k] == 'Axial':
						lineoutstring = ' @ R='+str(round(ZlineoutLoc*dr[l],2))+'cm'
						ylabel = 'Axial Distance Z [cm]'
						y1,y2 = Zaxis[0],Zaxis[-1]
					elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 1:
						lineoutstring = ' @ Z='+str(round(RlineoutLoc*dz[l],2))+'cm'
						ylabel = 'Radial Distance R [cm]'
						y1,y2 = -Raxis[-1],Raxis[-1]
					elif LineoutsOrientation[k] == 'Radial' and Isymlist[l] == 0:
						lineoutstring = ' @ Z='+str(round(RlineoutLoc*dz[l],2))+'cm'
						ylabel = 'Radial Distance R [cm]'
						y1,y2 = 0,Raxis[-1]
					#endif


					#Create PROES image along line of sight with phase-locked waveform.
					fig.suptitle( 'Simulated '+PhaseVariablelist[i]+' PROES for '+VariedValuelist[l]+lineoutstring, y=0.95, fontsize=18)
					im = ax[0].imshow(PROES,extent=[x1,x2,y1,y2],origin='bottom',aspect='auto')
					ax[0].set_ylabel(ylabel, fontsize=24)
					ax[0].tick_params(axis='x', labelsize=18)
					ax[0].tick_params(axis='y', labelsize=18)
					#Add Colourbar (Axis, Label, Bins)
#					label = VariableLabelMaker(PhaseVariablelist)
#					cax = Colourbar(ax[0],label[i],5)

					ax[1].plot(Phaseaxis, VoltageWaveform, lw=2)
					ax[1].set_xlabel('Phase [$\omega$t/2$\pi$]', fontsize=24)
					ax[1].set_ylabel('Potential [V]', fontsize=24)
					ax[1].tick_params(axis='x', labelsize=18)
					ax[1].tick_params(axis='y', labelsize=18)

					fig.tight_layout()
					plt.subplots_adjust(top=0.85)
					plt.savefig(DirPhaseResolved+VariedValuelist[l]+' '+NameString+' PROES.png')
					plt.close('all')


					#########################################
					#			UNDER CONSTRUCTION			#
					#########################################
					if True == False:
						#Collapse the 2D PROES image along the line of sight.
						fig,ax = figure(image_aspectratio,1)

						FlattenedPROES = list()
						for m in range(0,len(PROES)):
							FlattenedPROES.append( (sum(PROES[m][::])) )
						#endfor

						print y1,y2
						Filename = FolderNameTrimmer(Dirlist[l])+PhaseVariablelist[i]
						WriteDataToFile(FlattenedPROES,Filename)

						plt.plot(FlattenedPROES)
						plt.show()
					#endif
					#########################################
					#			UNDER CONSTRUCTION			#
					#########################################


					if l == numfolders and k == len(Lineouts):
						print'-------------------------------'
						print'# PROES Image analysis Complete'
						print'-------------------------------'
					#endif
				#endif

				#==============#
				# 	  PROES    #
				#==============#

				if savefig_phaseresolvelines == True:
					#Create .mp4 movie from completed images.
					Prefix = NameString+FolderNameTrimmer(Dirlist[l])
					Automovie(Dir1DProfiles,Prefix+'_'+PhaseVariablelist[i])
				#endif

				#Percentage Complete Readout.
				New += 1.0
				Old += 1.0
				oldpercentage = int( (New/Total)*100.0 )
				newpercentage = int( (Old/Total)*100.0 )
				if round(oldpercentage,-1) != round(newpercentage,-1):
					print int(round(oldpercentage, 0)),'%'
				#endif
			#endfor
		#endfor
	#endfor

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
