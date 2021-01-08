#!/bin/bash

module load singularity
module use /lcrc/project/jlab/local/etc/modulefiles
module load hallac_container/1.9.1

RUNGROUP=10*$1

export SINGULARITY_BINDPATH="/lcrc,/scratch" 
REPLAYDIR="/lcrc/project/jlab/csv/offline/online_csv"
ODIR="$REPLAYDIR/results"


echo "================================"
echo "do something for $RUNGROUP"
echo "================================"

cd ${REPLAYDIR}
root -q -b "shuo_analysis/check/skim_check.cxx+($RUNGROUP)" || exit $?
#root -q -b "shuo_analysis/simc_data/compare_sim_data_new.cxx+($RUNGROUP)" || exit$?
#root -q -b "shuo_analysis/Yield/statistic_runs_D2.cxx+($RUNGROUP)" || exit $?

echo "================================="
echo "end of do something $RUNGROUP"
echo "================================="
