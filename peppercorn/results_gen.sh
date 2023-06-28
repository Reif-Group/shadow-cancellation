#!/bin/bash

FOLDER=$1
TIME=$2
LABELS=$3

./pil.sh $1 $TIME orig_shadow_nocancel  "$LABELS" plot_o_s_nc
./pil.sh $1 $TIME orig_shadow_cancel  "$LABELS" plot_o_s_c
./pil.sh $1 $TIME orig_shadow_pert_nocancel  "$LABELS" plot_o_s_p_nc
./pil.sh $1 $TIME orig_shadow_pert_cancel  "$LABELS" plot_o_s_p_c
