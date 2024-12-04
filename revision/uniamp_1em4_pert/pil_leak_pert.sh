#names="orig_shadow_cancel orig_shadow_nocancel orig_shadow_pert_cancel orig_shadow_pert_nocancel orig_shadow_pert_cancel_zeroconc orig_shadow_pert_nocancel_zeroconc"
#names="orig_shadow_pert_cancel orig_shadow_pert_nocancel"
names="orig_shadow_cancel orig_shadow_nocancel"
TIME=$1
arr=$2 # "__0_1 0_0 0_5 1 4 9"
LABELS="$3"
mkdir -p images/
rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/leak_pert/$i $TIME $j "$LABELS" $j
    done
done

