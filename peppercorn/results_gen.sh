#!/bin/bash

FOLDER=$1
TIME=$2
LABELS=$3

./pil.sh $1 $TIME orig_shadow_nocancel  "$LABELS"
./pil.sh $1 $TIME orig_shadow_cancel  "$LABELS"
./pil.sh $1 $TIME orig_shadow_pert_nocancel  "$LABELS"
./pil.sh $1 $TIME orig_shadow_pert_cancel  "$LABELS"
