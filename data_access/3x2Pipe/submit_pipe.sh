#!/bin/bash -l

#SBATCH -N 1         #Use 1 nodes
#SBATCH -n 2
#SBATCH -t 00:40:00
#SBATCH --output=submit.out 
#SBATCH --error=submit.err 

srun -n 1 -c 2 ./launch_pipe.sh

