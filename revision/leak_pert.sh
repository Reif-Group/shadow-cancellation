#!/bin/bash
fracs="$1"
ROOT_FOLDER=$2
OUT_FOLDER=$3
CONC_FRAC=$4

for frac in $fracs; do
    `which python3` leak_pert.py $frac $ROOT_FOLDER $OUT_FOLDER $CONC_FRAC
done
