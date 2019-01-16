#!/bin/bash
#SBATch --account=jlab
#SBATCH -N1
#SBATCH -n1
#SBATCH --time=00:25:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=2048       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=singularity-test
#SBATCH --output=out/slurm-%A_%a.out
#SBATCH --error=out/slurm-%A_%a.err

for arun in $(seq 6000 6100) ; do

  srun -l bash batch/slurm/run_singularity.sh bin/hc_coin -r $arun -n 100000 &

done

wait
