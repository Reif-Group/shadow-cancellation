arr="0_0 0_5 1 4 9"
names="orig_shadow_cancel orig_shadow_nocancel orig_shadow_pert_cancel orig_shadow_pert_nocancel orig_shadow_pert_nocancel_zeroconc orig_shadow_pert_cancel_zeroconc"

rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh ../rps/$i 30000 $j "Ap Br Cj shAp shBr shCj Aq Bs Ck shAq shBs shCk" $j
    done
done

