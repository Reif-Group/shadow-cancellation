Python3.10 might not work, due to the tkinter issue. Use 3.8. It works.

## PIL Simulator command

> cat enum.pil | pilsimulator --t8 2080 --pyplot images/biamp_leaky_shadow.png --labels Br sBr Cj sCj Ck sCk --labels-strict --nxy --header --labels-strict > plot.txt

> peppercorn -o enum.pil -c  -L 7 --ignore-branch-4way  main.pil

## To generate cancellation reactions

> peppercorn -o cancel_enum.pil --reject-remote --k-slow 1e-5 --k-fast 1e-1 -c  --complex-prefix EC cancel.pil

## Produce-Helper Leak mechanism

sLeakWaste = hcjr( fcr mcr scr + fcr( hckr( fcr( + ) ) ) ) sbr* @initial 0 nM
LeakWaste = sc mc fc hcj( + sb* ) fc*( hck*( fc*( + ) ) ) @initial 0 nM

This leak rate constant corresponds to 150 /M/s which is the highest leak that was observed.

reaction [condensed    = 2e-7 /nM/s ] ProduceBCjCk + HelperCCk -> LeakWaste + Ck
reaction [condensed    = 2e-7 /nM/s ] sProduceBCjCk + sHelperCCk -> sLeakWaste + sCk



##### Updated Instructions as of May 21, 2023

# Compile the main into main_enum

./sim.sh rps2/vanilla 10200 main "Ak Br Cj" M

# Compile the shadow circuit

./sim.sh rps2/shadow_vanilla 10200 main "sAk sBr sCj" MS

# Perturb the rate constants of the shadow circuit. 

1. Jupyter notebook. Run the perturb section
2. Check
 
# Check the performance
./results_gen.sh rps/0_7 102000

