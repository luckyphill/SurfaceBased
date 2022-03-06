#!/bin/bash 
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH --mem=16GB
#SBATCH --time=72:00:00
# NOTIFICATIONS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phillip.j.brown@adelaide.edu.au

module load arch/haswell
module load matlab

export FACEDIR='/hpcfs/users/a1738927/Research/SurfaceBased'

matlab -nodisplay -nodesktop -r "cd ../../; addpath(genpath([pwd,'/src']));addpath(genpath([pwd,'/HPC']));addpath(genpath([pwd,'/analysis'])); t = Tumour3DGrowth; t.LoadSimulationData; t.PlotData; quit();"
