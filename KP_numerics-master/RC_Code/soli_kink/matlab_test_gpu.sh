#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --time=59:59:59
#SBATCH --partition=blanca-dhl
#SBATCH --qos=blanca-dhl
#SBATCH --output=/projects/mima6446/KP2/soli_kink/test_gpu.out
#SBATCH --mail-type=END
#SBATCH --mail-user=michelle.maiden@colorado.edu

module purge
module load matlab

matlab -nodisplay -nodesktop -r "clear; hello; test_gpu; clear all;exit();"
