FOLDER=$1
TIME=$2
LABELS="$3"

for j in original leaky occluded leaky_shadow; do
    ./pil.sh $FOLDER/$j $TIME main  "$LABELS" $j
done
