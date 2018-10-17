#!/bin/bash
# This is a script for computing importance sampling estimates
# of the log partition function in parallel.

#SBATCH --ntasks-per-node=16 --nodes=2
#SBATCH --account=uio --time=48:00:00 --mem-per-cpu=500M

source /cluster/bin/jobsetup

set -o errexit

module load R/3.5.0
module load openmpi.intel/1.10.2

mpirun -n 1 R --slave < estimate_z.R
