#!/bin/bash

module load singularity
module use /lcrc/project/jlab/local/etc/modulefiles
module load hallac_container/1.1.0

RUN=$1

export SINGULARITY_BINDPATH="/lcrc,/scratch" 
REPLAYDIR="/lcrc/project/jlab/csv/offline/online_csv"
TMPDIR="/scratch/Calib-$RUN"
ODIR="$REPLAYDIR/Data"


echo "================================"
echo "Collect info for run $RUN"
echo "================================"

cd ${REPLAYDIR}
mkdir -p $TMPDIR/full
mkdir -p $TMPDIR/logs
hcana -q -b "jlabana/grepcurrent_to_json.c+($RUN)" || exit $?
rsync -va $TMPDIR/full/* $ODIR/full
mkdir -p $ODIR/log/log-$RUN
rsync -va $TMPDIR/logs/* $ODIR/log/log-$RUN
rm -rf $TMPDIR

echo "================================="
echo "Finished processing COIN run $RUN"
echo "================================="
