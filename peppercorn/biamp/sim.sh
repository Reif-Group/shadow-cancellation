TIME=$1
peppercorn -o enum.pil main.pil
cat enum.pil | pilsimulator --t8 $TIME --pyplot images/biamp.png --labels Br  Cj ReactInt --labels-strict --nxy --header --labels-strict > plot.txt
