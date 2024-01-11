SUBFOLDER=$1
arr=$2 # "__0_1 0_0 0_5 1 4 9"
TIME=$3
LABELS=$4
names="orig_shadow_cancel orig_shadow_nocancel orig_shadow_pert_cancel orig_shadow_pert_nocancel orig_shadow_pert_cancel_zeroconc orig_shadow_pert_nocancel_zeroconc"
echo $names

rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/$SUBFOLDER/$i $TIME $j "$LABELS" $j
    done
done

