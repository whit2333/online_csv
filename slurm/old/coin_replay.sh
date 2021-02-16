#!/bin/bash

module load singularity
module use /lcrc/project/jlab/local/etc/modulefiles
module load hallac_container/1.1.0

RUN=$1

export SINGULARITY_BINDPATH="/lcrc,/scratch" 
#export SINGULARITY_BINDPATH="/scratch" 
REPLAYDIR="/lcrc/project/jlab/csv/offline/online_csv"
TMPDIR="/scratch/replay-$RUN"
ODIR="$REPLAYDIR/ROOTfiles"


echo "================================"
echo "Running replay for COIN run $RUN"
echo "================================"

cd ${REPLAYDIR}
mkdir -p $TMPDIR/ROOTfiles
mkdir -p $TMPDIR/out
hcana -q -b "scripts/replay_production_coin_hElec_pProt.C($RUN,-1)"
rsync -va $TMPDIR/* $ODIR
rm -rf $TMPDIR

echo "================================="
echo "Finished processing COIN run $RUN"
echo "================================="
