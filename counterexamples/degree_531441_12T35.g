# The first of the two counterexamples of degree 531441: the primitive
# affine group F_3^12 : (C_2^8 : (12T35)), with monomial point
# stabiliser of order 18432 generated below.  The regular vectors form
# a single orbit R of size 18432, exactly 1600 vectors lie outside
# R + R, and the Saxl graph has diameter exactly 3.

Read("verify.g");

gens := [
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [2,1,1,1,1,1,1,1,2,1,1,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,2,1,1,1,1,1,1,1,2,1,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,2,1,1,1,1,1,1,1,2,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,1,2,1,1,1,1,1,1,1,2]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,1,1,2,1,1,1,2,1,1,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,1,1,1,2,1,1,1,2,1,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,1,1,1,1,2,1,1,1,2,1]),
  MonomialGen(12, [0,1,2,3,4,5,6,7,8,9,10,11], [1,1,1,1,1,1,1,2,1,1,1,2]),
  MonomialGen(12, [0,5,2,7,4,9,6,11,8,1,10,3], [1,1,1,1,1,1,1,1,1,1,1,1]),
  MonomialGen(12, [0,7,2,5,4,3,6,1,8,11,10,9], [1,1,1,1,1,1,1,1,1,1,1,1]),
  MonomialGen(12, [11,2,1,4,3,6,5,8,7,10,9,0], [1,1,1,1,1,1,1,1,1,1,1,1])
];;

VerifySaxlCounterexample(gens, 12, 18432, 1600);
QUIT;
