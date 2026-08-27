# The counterexample of least degree: the primitive affine group
# F_3^9 : (C_2^6 : D18) of degree 19683, with point stabiliser the
# monomial group of order 1152 generated below (entry (9,3,1,200) of the
# IRREDSOL library).  The regular vectors form a single orbit R of size
# 1152, exactly 96 vectors lie outside R + R, and the Saxl graph has
# diameter exactly 3.

Read("verify.g");

gens := [
  MonomialGen(9, [0,2,1,6,8,7,3,5,4], [1,1,1,1,1,1,1,1,1]),
  MonomialGen(9, [4,5,3,6,7,8,1,2,0], [1,1,1,1,1,1,1,1,1]),
  MonomialGen(9, [0,1,2,3,4,5,6,7,8], [1,1,1,1,1,1,2,1,2])
];;

VerifySaxlCounterexample(gens, 9, 1152, 96);
QUIT;
