SUBFOLDER=$1
arr=$2 # "0 1 5 10 20 50 100"
TIME=$3
LABELS=$4
names="orig_shadow_cancel orig_shadow_nocancel"
echo $names

rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/$SUBFOLDER/$i $TIME $j "$LABELS" $j
    done
done

