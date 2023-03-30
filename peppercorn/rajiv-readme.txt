Python3.10 might not work, due to the tkinter issue. Use 3.8. It works.

## PIL Simulator command

> cat enum.pil | pilsimulator --t8 2080 --pyplot images/biamp_leaky_shadow.png --labels Br sBr Cj sCj Ck sCk --labels-strict --nxy --header --labels-strict > plot.txt

> peppercorn -o enum.pil -c  -L 7 --ignore-branch-4way  main.pil

## To generate cancellation reactions

> peppercorn -o cancel_enum.pil --reject-remote --k-slow 1e-5 --k-fast 1e-1 -c  --complex-prefix EC cancel.pil

## Produce-Helper Leak mechanism

sLeakWaste = hcjr( fcr mcr scr + fcr( hckr( fcr( + ) ) ) ) sbr* @initial 0 nM
LeakWaste = sc mc fc hcj( + sb* ) fc*( hck*( fc*( + ) ) ) @initial 0 nM

reaction [condensed    = 0.000001 /nM/s ] ProduceBCjCk + HelperCCk -> LeakWaste + Ck
reaction [condensed    = 0.000001 /nM/s ] sProduceBCjCk + sHelperCCk -> sLeakWaste + sCk
