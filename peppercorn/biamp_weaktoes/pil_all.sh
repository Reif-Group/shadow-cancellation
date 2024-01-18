TIME=$1
arr=$2 # "__0_1 0_0 0_5 1 4 9"
LABELS=$3
mkdir -p images/
rm images/*
for i in $arr; do
    for j in $names; do
        ../pil.sh `pwd`/$i $TIME main "$LABELS" $i
    done
done

