#------
## R18: Ym -> W
# Resting complexes (8) 
R18_ReactYmW = hW sYm( mYm( + fYm* ) ) @initial 10000 nM
R18_6 = sYm( mYm( fYm( hYm + ) ) )
R18_7 = hW sYm mYm
R18_12 = hW( sYm( mYm + ) ) fW*
R18_ProduceYmW = sW mW fW( hW( + sYm* ) ) @initial 10000 nM

# Resting macrostates (8) 
macrostate R18_ReactYmW = [R18_ReactYmW]
macrostate R18_6 = [R18_6]
macrostate R18_7 = [R18_7]
macrostate R18_12 = [R18_12]
macrostate R18_ProduceYmW = [R18_ProduceYmW]

# Condensed reactions (3) 
reaction [condensed    =  0.0001 /nM/s ] R18_ReactYmW + Ym -> R18_6 + R18_7
reaction [condensed    = 0.001 /nM/s ] R18_7 + R18_ProduceYmW -> R18_12 + W