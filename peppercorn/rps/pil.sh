# ./pil.sh 7200 rps_leak "Cj Br Ak"
 cat $2_enum.pil | pilsimulator --t8 $1 --pyplot images/$2.png --labels $3 --labels-strict --nxy --header --labels-strict > plot.txt
