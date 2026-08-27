# The primitive affine group F_3^10 : (C_2^5 : S_5) of degree 59049,
# with insoluble monomial point stabiliser of order 3840 generated
# below.  The regular vectors form a single orbit R of size 3840,
# exactly 64 vectors lie outside R + R, and the Saxl graph has diameter
# exactly 3.

Read("verify.g");

gens := [
  MonomialGen(10, [0,1,2,3,4,5,6,7,8,9], [2,1,1,1,2,1,1,1,2,2]),
  MonomialGen(10, [0,1,2,3,4,5,6,7,8,9], [1,2,1,1,2,1,1,2,1,2]),
  MonomialGen(10, [0,1,2,3,4,5,6,7,8,9], [1,1,2,1,1,1,2,2,2,1]),
  MonomialGen(10, [0,1,2,3,4,5,6,7,8,9], [1,1,1,2,2,1,2,2,1,1]),
  MonomialGen(10, [0,1,2,3,4,5,6,7,8,9], [1,1,1,1,1,2,2,1,2,2]),
  MonomialGen(10, [2,3,4,5,6,7,8,9,0,1], [1,1,1,1,1,1,1,1,1,1]),
  MonomialGen(10, [1,0,6,3,4,5,2,8,7,9], [1,1,1,1,1,1,1,1,1,1])
];;

VerifySaxlCounterexample(gens, 10, 3840, 64);
QUIT;
