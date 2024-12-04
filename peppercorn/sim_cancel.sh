#!/bin/bash
DIR=$1
TIME=$2
NAME=$3
LABELS=$4
PREFIX=$5
mkdir -p $DIR/images
peppercorn -o "$DIR/$NAME"_enum.pil --reject-remote --k-slow 1e-5 --k-fast 1e-1 -c  --complex-prefix $PREFIX $DIR/$NAME.pil
cat "$DIR/$NAME"_enum.pil | pilsimulator --t8 $TIME --pyplot "$DIR"/images/$NAME.png --labels $LABELS --labels-strict --nxy --header --labels-strict > $DIR/plot.txt
