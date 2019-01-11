#!/bin/bash -l

#SBATCH -N 1         #Use 1 nodes
#SBATCH -t 02:00:00  #Set 10 minute time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH --output=submit.out 
#SBATCH --error=submit.err 

srun -n 1 -c 14 -A m1727 ./launch_pipe.sh

