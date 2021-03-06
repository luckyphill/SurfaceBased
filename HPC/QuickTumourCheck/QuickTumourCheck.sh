#!/bin/bash 
#SBATCH -p batch 
#SBATCH -N 1 
#SBATCH -n 4
#SBATCH --mem=8GB
#SBATCH --array=0-20
#SBATCH --time=72:00:00
# NOTIFICATIONS
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phillip.j.brown@adelaide.edu.au

module load arch/haswell
module load matlab

export FACEDIR='/hpcfs/users/a1738927/Research/SurfaceBased'

echo "array_job_index: $SLURM_ARRAY_TASK_ID"

i=1 
found=0 

while IFS=, read a b c d e f g h
do 
	if [ $i = $SLURM_ARRAY_TASK_ID ]; then
		echo "Running with [$a, $b, $c, $d, $e, $f, $g, $h]"
		found=1 

		break 
	fi 
	i=$((i + 1)) 
done < QuickTumourCheck.txt

if [ $found = 1 ]; then
	echo "matlab -nodisplay -nodesktop -r cd ../../; addpath(genpath([pwd,'/src']));addpath(genpath([pwd,'/HPC']));addpath(genpath([pwd,'/analysis'])); t = Tumour3D($a, $b, $c, $d, $e, $f, $g, $h); t.RunToConfluence(40, 200); v = Visualiser(t); v.ProduceMovie([],[],'Motion JPEG AVI'); quit()"
	matlab -nodisplay -nodesktop -r "cd ../../; addpath(genpath([pwd,'/src']));addpath(genpath([pwd,'/HPC']));addpath(genpath([pwd,'/analysis'])); t = Tumour3D($a, $b, $c, $d, $e, $f, $g, $h); t.RunToConfluence(40, 200); v = Visualiser(t); v.ProduceMovie([],[],'Motion JPEG AVI'); quit()" 
else 
  echo "SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is outside range of input file" 
fi
