arr=$2 # "__0_1 0_0 0_5 1 4 9"
names="orig_shadow_pert_cancel orig_shadow_pert_nocancel"
TIME=$1

rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/pert/$i $TIME $j "Ap Br Cj shAp shBr shCj Aq Bs Ck shAq shBs shCk" $j
    done
done

