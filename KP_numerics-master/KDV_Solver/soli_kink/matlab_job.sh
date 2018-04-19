#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --time=99:59:59
#SBATCH --partition=blanca-dhl
#SBATCH --qos=blanca-dhl
#SBATCH --output=/projects/mima6446/KP2/soli_kink/matlab_outs/matlab_job_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=michelle.maiden@colorado.edu

module purge
module load matlab

matlab -nodisplay -nodesktop -r "clear; hello; driver_KP_solver_param_space_segments_GPU; clear all;exit();"
