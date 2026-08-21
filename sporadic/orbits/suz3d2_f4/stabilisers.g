# Construct exact stabilizers for all five 3.Suz.2-orbits on the binary
# 24-dimensional module and export their generators as packed row masks.
#
# The faithful 98280-point ATLAS action is aligned with an orbit of nonzero
# dual vectors.  Hence a primal vector v is encoded faithfully by
#   { i : dualVectors[i](v) = 1 },
# and its vector stabilizer is the corresponding set stabilizer.
#
# Usage:
#   gap -q -b stabilisers.g > stabilisers.txt

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;
Reset(GlobalMersenneTwister, 20260723);;

F := GF(2);;
zero := Zero(F);;
d := 24;;
expectedOrder := 2690072985600;;

matrixInfo := OneAtlasGeneratingSetInfo(
  "3.Suz.2", Characteristic, 2, Dimension, d
);;
matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
if matrixAtlas.repname <> "3Suzd2G1-f2r24B0" then
  Error("unexpected matrix representation");
fi;
matrixGenerators := matrixAtlas.generators;;
G := Group(matrixGenerators);;
SetSize(G, expectedOrder);;

permutationInfo := OneAtlasGeneratingSetInfo(
  "3.Suz.2", NrMovedPoints, 98280
);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
if permutationAtlas.repname <> "3Suzd2G1-p98280B0" then
  Error("unexpected permutation representation");
fi;
permutationGenerators := permutationAtlas.generators;;
P := Group(permutationGenerators);;
SetSize(P, expectedOrder);;

VectorFromMask := function(mask)
  return List([0 .. d - 1], position ->
    (QuoInt(mask, 2^position) mod 2) * One(F));
end;;

RowMask := function(row)
  local mask, position;
  mask := 0;
  for position in [1 .. d] do
    if not IsZero(row[position]) then
      mask := mask + 2^(position - 1);
    fi;
  od;
  return mask;
end;;

CanonicalSubspaceBasis := function(basis)
  return BasisVectors(CanonicalBasis(Subspace(F^d, basis)));
end;;

Dot := function(v, f)
  return Sum([1 .. d], position -> v[position] * f[position]);
end;;

SupportOf := function(v, orbit)
  local support, point;
  support := [];
  for point in [1 .. Length(orbit)] do
    if Dot(v, orbit[point]) <> zero then
      Add(support, point);
    fi;
  od;
  return support;
end;;

AllFixVector := function(matrices, v)
  return ForAll(matrices, matrix -> v * matrix = v);
end;;

# Maximal subgroup 4, (3 x U5(2)).2, fixes a two-dimensional subspace of the
# dual module.  Its three nonzero vectors correspond to the three points
# {1,35480,43655}; direct generator propagation fixes the alignment below.
maxProgram := AtlasProgram(
  "3.Suz.2", matrixAtlas.standardization, "maxes", 4
);;
maxMatrixGenerators := ResultOfStraightLineProgram(
  maxProgram.program, matrixGenerators
);;
dualGenerators := List(matrixGenerators,
  matrix -> TransposedMat(Inverse(matrix)));;
maxDualModule := GModuleByMats(
  List(maxMatrixGenerators,
    matrix -> TransposedMat(Inverse(matrix))), F
);;
series := MTX.BasesCompositionSeries(maxDualModule);;
if List(series, Length) <> [0, 2, 22, 24] then
  Error("unexpected composition series for maximal subgroup 4");
fi;
dualPlane := CanonicalSubspaceBasis(series[2]);;
dualSeed := dualPlane[1];;

dualVectors := [];;
dualVectors[1] := dualSeed;;
queue := [1];;
head := 1;;
while head <= Length(queue) do
  point := queue[head];;
  head := head + 1;;
  for generatorPosition in [1 .. Length(permutationGenerators)] do
    imagePoint := point ^ permutationGenerators[generatorPosition];;
    imageVector :=
      dualVectors[point] * dualGenerators[generatorPosition];;
    if not IsBound(dualVectors[imagePoint]) then
      dualVectors[imagePoint] := imageVector;
      Add(queue, imagePoint);
    elif dualVectors[imagePoint] <> imageVector then
      Error("matrix/permutation dual-orbit alignment conflict");
    fi;
  od;
od;
if Length(queue) <> 98280 or Length(dualVectors) <> 98280 then
  Error("unexpected aligned dual-orbit length");
fi;

# Both ATLAS representations use the same standard generators.  Recover a
# word for every permutation stabilizer generator and evaluate that word in
# the matrix generators.
epimorphism := EpimorphismFromFreeGroup(P);;
SetIsSurjective(epimorphism, true);;
freeGenerators := GeneratorsOfGroup(Source(epimorphism));;
if List(freeGenerators, generator -> Image(epimorphism, generator))
    <> permutationGenerators then
  Error("free presentation is not based on the ATLAS standard generators");
fi;

# Complete binary-vector orbit census, independently obtained by orbit_census.
cases := [
  rec(rep := 0, orbitLength := 1,
      expectedStabilizer := 2690072985600),
  rec(rep := 1, orbitLength := 1216215,
      expectedStabilizer := 2211840),
  rec(rep := 7, orbitLength := 7076160,
      expectedStabilizer := 380160),
  rec(rep := 96, orbitLength := 8386560,
      expectedStabilizer := 320760),
  rec(rep := 324, orbitLength := 98280,
      expectedStabilizer := 27371520)
];;

Print("SUZ3D2_F4_VECTOR_STABILIZERS_V1\n");
Print("MATRIX_ATLAS_REP ", matrixAtlas.repname, "\n");
Print("PERMUTATION_ATLAS_REP ", permutationAtlas.repname, "\n");
Print("GROUP_ORDER ", expectedOrder, "\n");
Print("DUAL_ORBIT_LENGTH ", Length(queue), "\n");
Print("CASES ", Length(cases), "\n");

for case in cases do
  v := VectorFromMask(case.rep);;
  support := SupportOf(v, dualVectors);;
  if case.rep = 0 then
    # Avoid rediscovering the order of the full matrix group: its order and
    # generators were already certified above.
    hGenerators := matrixGenerators;
  else
    HP := Stabilizer(P, support, OnSets);;
    if Size(HP) <> case.expectedStabilizer then
      Error("support stabilizer has unexpected order");
    fi;
    if expectedOrder / Size(HP) <> case.orbitLength then
      Error("support stabilizer contradicts binary orbit census");
    fi;
    hpGenerators := SmallGeneratingSet(HP);;
    if hpGenerators = fail then
      hpGenerators := GeneratorsOfGroup(HP);
    fi;
    hGenerators := List(hpGenerators, permutationElement ->
      MappedWord(
        PreImagesRepresentative(epimorphism, permutationElement),
        freeGenerators, matrixGenerators
      )
    );;
  fi;
  if not AllFixVector(hGenerators, v) then
    Error("a pulled-back matrix generator does not fix the primal vector");
  fi;

  Print("CASE ", case.rep, " ", case.orbitLength, " ",
    case.expectedStabilizer, " ", Length(hGenerators), " ",
    Length(support), "\n");
  for h in hGenerators do
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
