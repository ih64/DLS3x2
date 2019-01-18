#!/bin/bash -l

#SBATCH -N 1         #Use 1 nodes
#SBATCH -n 16
#SBATCH -t 00:40:00
#SBATCH --output=submit.out 
#SBATCH --error=submit.err 

srun -n 1 -c 16 ./launch_pipe.sh

