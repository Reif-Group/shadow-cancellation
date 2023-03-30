# ./sim.sh 1000 vanilla "Cj Br"
mkdir -p images
TIME=$1
peppercorn -o $2_enum.pil  $2.pil --condensed
cat $2_enum.pil | pilsimulator --t8 $TIME --pyplot images/$2.png --labels $3 --labels-strict --nxy --header --labels-strict > plot.txt
