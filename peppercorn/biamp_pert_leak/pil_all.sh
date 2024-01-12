names="orig_shadow_cancel orig_shadow_nocancel"
TIME=$1
arr=$2 # "__0_1 0_0 0_5 1 4 9"
LABELS=$3
mkdir -p images/
rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/pert_leak/$i $TIME $j "$LABELS" $j
    done
done

