#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --time=00:10:00
#SBATCH --partition=blanca-dhl
#SBATCH --qos=blanca-dhl
#SBATCH --output=matlab_job_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=ezio.iacocca@colorado.edu

module purge
module load matlab

matlab -nodisplay -nodesktop -r "clear; YorkInitialState;clear all;exit();"
