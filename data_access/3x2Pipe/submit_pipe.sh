#!/bin/bash -l

#SBATCH -N 1         #Use 1 nodes
#SBATCH -n 16
#SBATCH -t 01:00:00  #Set 10 minute time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH --output=submit.out 
#SBATCH --error=submit.err 

srun -n 1 -c 16 ./launch_pipe.sh

