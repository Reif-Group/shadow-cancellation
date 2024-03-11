ROOT=$1
OUT=$2
TIME=$3
FOLDERS="$4"
LABELS="$5"
CONC_FACT=$6

./leak_pert.sh "$FOLDERS" $ROOT $OUT $CONC_FACT
cd $ROOT
./leak_pert_pil_all.sh $TIME "$FOLDERS"  "$LABELS"
