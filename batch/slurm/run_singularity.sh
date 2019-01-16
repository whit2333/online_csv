#!/bin/bash

module load singularity
singularity exec -B /lcrc:/lcrc /lcrc/project/jlab/images/cool_halls.simg bash "$1" "${@:2}"

