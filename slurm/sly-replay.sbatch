#!/bin/bash
#SBATCH --account=jlab
#SBATCH --time=12:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=2048       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=coin
#SBATCH --ntasks=2
#DISABLED_SBATCH --output=/home/sjia/out/coin-%A_%a.out
#DISABLED_SBATCH --error=/home/sjia/out/coin-%A_%a.err
#SBATCH --output=/lcrc/project/jlab/csv/log/coin-%A_%a.out
#SBATCH --error=/lcrc/project/jlab/csv/log/coin-%A_%a.err

# Load the default version of GNU parallel.
module purge
module load gcc/7.1.0-4bgguyp
module load parallel

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
ulimit -u 10000

# This specifies the options used to run srun. The "-N1 -n1" options are
# used to allocates a single core to each task.
srun="srun --exclusive -N1 -n1"
#JOBLOG=/home/sjia/out/coin-test.progress
JOBLOG=/home/sly2j/out/coin-test.progress
#JOBLOG=/lcrc/project/jlab/csv/log/coin-test.progress
LOGDIR=/lcrc/project/jlab/data/hallc/csv/replay/log

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog ${JOBLOG} --resume"

# Run the script, coin_replay.sh, using GNU parallel and srun. Parallel
# will run the runtask script for the numbers 1 through $SLURM_NTASKS 
$parallel "$srun ./coin_replay.sh {1} 2>&1 > $LOGDIR/coin_replay.sh.{1}" ::: {7593..7594}
