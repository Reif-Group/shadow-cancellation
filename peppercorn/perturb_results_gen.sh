ROOT=$1 # rps
TIME=$2 # 10000
FILE=$3 #orig_shadow_pert_nocancel
LABELS=$4 # "A B C shA shB shC"

./pil.sh $ROOT/0_0 $TIME $FILE  "$LABELS" $FILE
./pil.sh $ROOT/0_4 $TIME $FILE  "$LABELS" $FILE
./pil.sh $ROOT/1 $TIME $FILE  "$LABELS" $FILE
./pil.sh $ROOT/4 $TIME $FILE  "$LABELS" $FILE
./pil.sh $ROOT/9 $TIME $FILE  "$LABELS" $FILE

