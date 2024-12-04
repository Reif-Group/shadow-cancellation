#!/bin/bash
ENUM_CODE=$3
TIME=$1
LABELS=$2
echo "$ENUM_CODE" | pilsimulator --t8 $TIME  --labels $LABELS --labels-strict --nxy --header
