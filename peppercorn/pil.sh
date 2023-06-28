#!/bin/bash
DIR=$1
TIME=$2
NAME=$3
LABELS=$4
PLOT_NAME=$5
mkdir -p $DIR/images
cat "$DIR/$NAME"_enum.pil | pilsimulator --t8 $TIME --pyplot "$DIR"/images/$NAME.png --labels $LABELS --labels-strict --nxy --header --labels-strict > $DIR/$5
