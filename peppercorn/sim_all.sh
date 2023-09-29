#!/bin/bash

DIR=$1
NAME=$2

mkdir -p $DIR/original/images
peppercorn -o "$DIR/original/$NAME"_enum.pil  --max-complex-size 7 --complex-prefix M  --k-slow 1e-5 --k-slow 1e-2 -c --ignore-branch-4way  $DIR/original/$NAME.pil
peppercorn -o "$DIR/shadow/$NAME"_enum.pil  --max-complex-size 7 --complex-prefix MS  --k-slow 1e-5 --k-slow 1e-2 -c  --ignore-branch-4way  $DIR/shadow/$NAME.pil
peppercorn -o "$DIR/cancel/$NAME"_enum.pil  --max-complex-size 7 --complex-prefix EC  --k-slow 1e-5 --k-fast 1e-1 -c $DIR/cancel/$NAME.pil



