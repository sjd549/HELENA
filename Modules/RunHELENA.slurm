#!/bin/bash
# Set SLURM Variables
#SBATCH --job-name="HELENA3"					# Job Name
#SBATCH --mail-user=username@abdn.ac.uk			# User Mail
#SBATCH --mail-type="NONE"						# Mail if condition met
#SBATCH --output=""								# Output file location, default to dev/null
#SBATCH --nodes=1								# Default 1 compter/node
#SBATCH --ntasks=1								# Default 2 threads (Most nodes have a maximum of 40 cores)
#SBATCH --time=0-01:00:00						# Max Walltime, job is killed after this
#SBATCH --mem-per-cpu=6G						# Default 2G
#SBATCH --partition="spot-compute"				# Partition Name
 
date
hostname

./HELENA3.py									# Run simulation
 
date
exit 0







Add details of WinSCP
	- right click and select "open" to see images

Add details on how to identify and alter file permissions
	chmod +x "filename"
	Needed for ./comppar
	Needed for ./HELENA3.py

Add RunHELENA.slurm into training package
	TEST THIS FIRST...

Add details of how to build and use initmesh.exe
	Compilation gfortran initmesh.f -o initmesh.exe
	Best practices for use

Add details on how to restart a HPEM simulation
	Restarting_HPEM.txt

Add details of how to set modules and export automatically
	.bash_profile to be included:
	
# ########## ########## ##########
#
# Typically ~/.profile is executed by the command interpreter for login shells.
# This file is not read if ~/.bash_profile or ~/.bash_login exists.
# see /usr/share/doc/bash/examples/startup-files for examples.
# the files are located in the bash-doc package.
#
# the default umask is set in /etc/profile; for setting the umask
# for ssh logins, install and configure the libpam-umask package.
# umask 022
#
# ########## ########## ##########


#######################################
# User specific aliases and functions #
#######################################

# Load ifort libraries
module load compiler/latest > /dev/null
module load mpi/latest > /dev/null

# Set Environment Variables
export OMP_NUM_THREADS=2 > /dev/null
export MKL_NUM_THREADS=2 > /dev/null


#####################################
# Print Welcome Information To User #
#####################################

# Print welcome information
echo 'Batched and Running Jobs:'
squeue -u s09sd3

