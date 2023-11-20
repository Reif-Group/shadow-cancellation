names="orig_shadow_cancel orig_shadow_nocancel orig_shadow_pert_cancel orig_shadow_pert_nocancel orig_shadow_pert_nocancel_zeroconc orig_shadow_pert_cancel_zeroconc"
TIME=$1
arr=$2 # "__0_1 0_0 0_5 1 4 9"
LABELS=$3
mkdir -p images/
rm images/*
for i in $arr; do
    for j in $names; do
        nohup ../pil.sh `pwd`/$i $TIME $j "$LABELS" $j &
    done
done

