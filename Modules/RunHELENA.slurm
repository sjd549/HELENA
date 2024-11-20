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
