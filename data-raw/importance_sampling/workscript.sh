#!/bin/bash
# This is a script for computing importance sampling estimates
# of the log partition function in parallel.

#SBATCH --ntasks=16
#SBATCH --account=uio --time=00:10:00 --mem-per-cpu=500M

source /cluster/bin/jobsetup

set -o errexit

module load R/3.5.0
module load openmpi.intel/1.10.2

mpirun -n 1 R --slave < estimate_z.R
