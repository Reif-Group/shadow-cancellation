#!/bin/bash
DIR=$1
TIME=$2
NAME=$3
LABELS=$4
PREFIX=$5
mkdir -p $DIR/images
peppercorn -o "$DIR/$NAME"_enum.pil  --max-complex-size 6 --complex-prefix $PREFIX  -L 7 -c  --ignore-branch-4way  $DIR/$NAME.pil
cat "$DIR/$NAME"_enum.pil | pilsimulator --t8 $TIME --pyplot "$DIR"/images/$NAME.png --labels $LABELS --labels-strict --nxy --header --labels-strict > $DIR/plot.txt
