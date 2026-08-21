# Construct exact stabilizers for all five Co3-orbits on the binary module.
# The stabilizers are exported as packed 22-bit matrix rows for the C++
# second-vector orbit census.
#
# The 11178-point dual orbit gives a faithful permutation copy of Co3.
# A primal vector v is embedded as its support
#   { f in Omega : f(v) = 1 }.
# Its vector stabilizer is therefore the set stabilizer of that support.
#
# Usage:
#   gap -q stabilisers.g > stabilisers.txt

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

F := GF(2);;
zero := Zero(F);;
d := 22;;
space_size := 2^d;;
expected_order := 495766656000;;

G := AtlasGroup(
  "Co3", IsMatrixGroup, true, Characteristic, 2, Dimension, d
);;
if G = fail then
  Error("ATLAS Co3 GF(2), dimension-22 representation is unavailable");
fi;
gens := GeneratorsOfGroup(G);;
if Size(G) <> expected_order then
  Error("unexpected Co3 order");
fi;

VectorFromMask := function(mask)
  return List([0 .. d - 1], i ->
    (QuoInt(mask, 2^i) mod 2) * One(F));
end;;

RowMask := function(row)
  local mask, i;
  mask := 0;
  for i in [1 .. d] do
    if not IsZero(row[i]) then
      mask := mask + 2^(i - 1);
    fi;
  od;
  return mask;
end;;

Dot := function(v, f)
  return Sum([1 .. d], i -> v[i] * f[i]);
end;;

SupportOf := function(v, orbit)
  local support, i;
  support := [];
  for i in [1 .. Length(orbit)] do
    if Dot(v, orbit[i]) <> zero then
      Add(support, i);
    fi;
  od;
  return support;
end;;

AllFixVector := function(matrices, v)
  local matrix;
  for matrix in matrices do
    if v * matrix <> v then
      return false;
    fi;
  od;
  return true;
end;;

# Maximal subgroup 2 is HS, of index 11178.  Its common fixed space on the
# dual module is one-dimensional and supplies a seed for the dual orbit.
HS := AtlasSubgroup(G, 2);;
if Size(HS) <> 44352000 then
  Error("unexpected order for ATLAS maximal subgroup 2");
fi;
hs_gens := GeneratorsOfGroup(HS);;
identity := IdentityMat(d, F);;
fixed_equations := List([1 .. d], i ->
  Concatenation(List(hs_gens, h -> TransposedMat(h)[i] - identity[i]))
);;
dual_fixed := NullspaceMat(fixed_equations);;
if Length(dual_fixed) <> 1 then
  Error("HS does not have a one-dimensional fixed space on the dual");
fi;
dual_seed := dual_fixed[1];;

dual_gens := List(gens, g -> TransposedMat(Inverse(g)));;
D := Group(dual_gens);;
dual_orbit := Orbit(D, dual_seed, OnRight);;
if Length(dual_orbit) <> 11178 then
  Error("unexpected dual orbit length");
fi;
action := ActionHomomorphism(D, dual_orbit, OnRight);;
P := Image(action);;
if Size(P) <> expected_order then
  Error("the dual-orbit action is not faithful Co3");
fi;

# The complete binary-vector orbit census is:
#   representative  orbit length  stabilizer order
#       0                 1          495766656000
#       1             37950             13063680
#       7           2608200               190080
#      13           1536975               322560
#    1023             11178             44352000
# Only representatives 7 and 13 can possibly have a regular stabilizer orbit
# on a set of 2^22 vectors.  All five are exported because their second-vector
# orbit representatives form a complete Co3-orbit transversal on GF(4)^22.
cases := [
  rec(rep := 0, orbit_length := 1,
      expected_stabilizer := 495766656000),
  rec(rep := 1, orbit_length := 37950,
      expected_stabilizer := 13063680),
  rec(rep := 7, orbit_length := 2608200, expected_stabilizer := 190080),
  rec(rep := 13, orbit_length := 1536975, expected_stabilizer := 322560),
  rec(rep := 1023, orbit_length := 11178,
      expected_stabilizer := 44352000)
];;

Print("CO3_F4_VECTOR_STABILIZERS_V1\n");
Print("ATLAS_REP Co3G1-f2r22B0\n");
Print("GROUP_ORDER ", Size(G), "\n");
Print("DUAL_ORBIT_LENGTH ", Length(dual_orbit), "\n");
Print("DUAL_ACTION_ORDER ", Size(P), "\n");
Print("CASES ", Length(cases), "\n");

for case in cases do
  v := VectorFromMask(case.rep);;
  if Length(Orbit(G, v, OnRight)) <> case.orbit_length then
    Error("binary-vector orbit length does not match C++ census");
  fi;
  support := SupportOf(v, dual_orbit);;
  HP := Stabilizer(P, support, OnSets);;
  if Size(HP) <> case.expected_stabilizer then
    Error("support set stabilizer has unexpected order");
  fi;

  # Pull a compact generating set back through the faithful dual action, then
  # undo the contragredient construction to recover primal matrices.
  hp_gens := SmallGeneratingSet(HP);;
  if hp_gens = fail then
    hp_gens := GeneratorsOfGroup(HP);
  fi;
  h_gens := List(hp_gens, p ->
    Inverse(TransposedMat(PreImagesRepresentative(action, p))));;
  H := Group(h_gens);;
  if Size(H) <> case.expected_stabilizer then
    Error("pulled-back primal generators have unexpected order");
  fi;
  if not AllFixVector(h_gens, v) then
    Error("a pulled-back generator does not fix the primal vector");
  fi;

  Print("CASE ", case.rep, " ", case.orbit_length, " ",
    Size(H), " ", Length(h_gens), " ", Length(support), "\n");
  for h in h_gens do
    Print("MATRIX");
    for row in h do
      Print(" ", RowMask(row));
    od;
    Print("\n");
  od;
  Print("ENDCASE\n");
od;
Print("END\n");

QUIT_GAP(0);
