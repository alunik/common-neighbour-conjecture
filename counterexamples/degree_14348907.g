# The counterexample with perfect point stabiliser: the primitive
# affine group F_3^15 : AGL_4(2) of degree 14348907, where AGL_4(2)
# acts on the 15-dimensional deleted permutation module over GF(3) of
# its natural 2-transitive action on 16 points.  The regular vectors
# form a single orbit R of size 322560, exactly 32 vectors lie outside
# R + R, and the Saxl graph has diameter exactly 3.
#
# The deleted permutation module is the sum-zero subspace of GF(3)^16,
# with basis b_i = e_i - e_16 for i < 16; a vector in the subspace has
# coordinates c_i equal to its entries at the positions i < 16.

if LoadPackage("primgrp") = fail then Error("PrimGrp is required"); fi;
Read("verify.g");

P16 := AllPrimitiveGroups(NrMovedPoints, 16, Size, 322560);;
Require(Length(P16) = 1,
        "AGL(4,2) should be the unique primitive group of degree 16 and order 322560");
P16 := P16[1];;
Require(IsPerfectGroup(P16), "AGL(4,2) should be perfect");

F := GF(3);;
deleted := function(pi)
  local M, i, row;
  M := [];
  for i in [1..15] do
    row := ListWithIdenticalEntries(15, Zero(F));
    if i^pi < 16 then
      row[i^pi] := row[i^pi] + One(F);
    fi;
    if 16^pi < 16 then
      row[16^pi] := row[16^pi] - One(F);
    fi;
    ConvertToVectorRep(row, 3);
    M[i] := row;
  od;
  ConvertToMatrixRep(M, 3);
  return M;
end;;

gens := List(GeneratorsOfGroup(P16), deleted);;

VerifySaxlCounterexample(gens, 15, 322560, 32);
QUIT;
