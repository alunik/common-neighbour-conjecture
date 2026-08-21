# Exact ATLAS export for the 3.Fi22 <= GL(27,4) result.
#
# The degree-61776 action of Fi22 is the action on conjugates of max4,
# O8+(2).3.2.  On the contragredient 27-dimensional module that subgroup
# fixes one projective line.  Starting from the common standard generators,
# this script transports the fixed line and the permutation point together;
# equality on every revisited edge checks the alignment.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

ROOT := "sporadic/large/fi22";;
OUT := Concatenation(ROOT, "/generated");;
Exec(Concatenation("mkdir -p ", OUT));

AssertFi22 := function(condition, message)
  if not condition then
    Error(Concatenation("FI22 EXPORT FAILURE: ", message));
  fi;
end;;

CanonicalLine := function(vector)
  local pivot;
  pivot := PositionProperty(vector, x -> not IsZero(x));
  AssertFi22(pivot <> fail, "attempt to normalize the zero vector");
  return vector / vector[pivot];
end;;

Code := function(element)
  if IsZero(element) then
    return 0;
  elif element = One(GF(4)) then
    return 1;
  elif element = Z(4) then
    return 2;
  elif element = Z(4)^2 then
    return 3;
  fi;
  Error("element is not in the reference GF(4) encoding");
end;;

Digits := function(vector)
  return Concatenation(List(vector, x -> String(Code(x))));
end;;

field := GF(4);;
matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
permutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22", NrMovedPoints, 61776);;
AssertFi22(matrixInfo <> fail and permutationInfo <> fail,
    "required ATLAS representations are unavailable");
matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
matrixGenerators := matrixAtlas.generators;;
permutationGenerators := permutationAtlas.generators;;
AssertFi22(matrixAtlas.repname = "3F22G1-f4r27aB0",
    "unexpected matrix representation");
AssertFi22(permutationAtlas.repname = "F22G1-p61776B0",
    "unexpected permutation representation");
AssertFi22(List(matrixGenerators, Order) = [2,13],
    "unexpected matrix generator orders");
AssertFi22(List(permutationGenerators, Order) = [2,13],
    "unexpected permutation generator orders");

max4 := AtlasProgram("Fi22", "maxes", 4);;
AssertFi22(max4 <> fail and max4.subgroupname = "O8+(2).3.2",
    "unexpected max4 program");
maxMatrixGenerators := ResultOfStraightLineProgram(
    max4.program, matrixGenerators);;
maxPermutationGenerators := ResultOfStraightLineProgram(
    max4.program, permutationGenerators);;
fixedPoints := Filtered([1..61776], point ->
    ForAll(maxPermutationGenerators, generator ->
      point ^ generator = point));;
AssertFi22(fixedPoints = [1], "max4 is not aligned at point 1");

# The max4 module has an invariant hyperplane.  Its annihilator is the
# required invariant line in the contragredient module.
compositionBases := MTX.BasesCompositionSeries(
    GModuleByMats(maxMatrixGenerators, field));;
AssertFi22(List(compositionBases, Length) = [0,26,27],
    "unexpected max4 composition dimensions");
dualSeed := CanonicalLine(
    NullspaceMat(TransposedMat(compositionBases[2]))[1]);;
dualGenerators := List(matrixGenerators,
    generator -> TransposedMat(generator^-1));;
dualMaxGenerators := List(maxMatrixGenerators,
    generator -> TransposedMat(generator^-1));;
AssertFi22(ForAll(dualMaxGenerators, generator ->
    RankMat([dualSeed, dualSeed * generator]) = 1),
    "annihilator line is not max4-invariant");

dualLines := [];;
dualLines[1] := dualSeed;;
queue := [1];;
head := 1;;
while head <= Length(queue) do
  point := queue[head];
  head := head + 1;
  for generatorPosition in [1..2] do
    imagePoint := point ^ permutationGenerators[generatorPosition];
    imageLine := CanonicalLine(
        dualLines[point] * dualGenerators[generatorPosition]);
    if not IsBound(dualLines[imagePoint]) then
      dualLines[imagePoint] := imageLine;
      Add(queue, imagePoint);
    else
      AssertFi22(dualLines[imagePoint] = imageLine,
          "matrix/permutation orbit alignment conflict");
    fi;
  od;
od;
AssertFi22(Length(queue) = 61776,
    "dual-line orbit does not have length 61776");
AssertFi22(Length(Set(dualLines)) = 61776,
    "dual-line orbit contains duplicate projective lines");
AssertFi22(RankMat(dualLines) = 27,
    "dual-line orbit does not span the dual module");

dualPath := Concatenation(OUT, "/dual_orbit.tsv");;
output := OutputTextFile(dualPath, false);;
SetPrintFormattingStatus(output, false);
for point in [1..61776] do
  AppendTo(output, point, "\t", Digits(dualLines[point]), "\n");
od;
CloseStream(output);

matrixPath := Concatenation(OUT, "/matrix_generators.txt");;
output := OutputTextFile(matrixPath, false);;
SetPrintFormattingStatus(output, false);
for generator in matrixGenerators do
  for row in generator do
    AppendTo(output, Digits(row), "\n");
  od;
od;
CloseStream(output);

form := MTX.InvariantSesquilinearForm(
    GModuleByMats(matrixGenerators, field));;
AssertFi22(form <> fail and RankMat(form) = 27,
    "nondegenerate invariant sesquilinear form not found");
formPath := Concatenation(OUT, "/hermitian_form.txt");;
output := OutputTextFile(formPath, false);;
SetPrintFormattingStatus(output, false);
for row in form do
  AppendTo(output, Digits(row), "\n");
od;
CloseStream(output);

# A reference regular seed from the earlier exact census.  The zero functionals
# span its full annihilator, so its zero-set stabilizer in the faithful
# 61776-point action is exactly its projective-line stabilizer.
decode := [Zero(field), One(field), Z(4), Z(4)^2];;
regularSeed := List("032133122122121333223313020",
    character -> decode[IntChar(character) - 47]);;
zeroSet := Set(Filtered([1..61776], point ->
    Sum([1..27], coordinate ->
      regularSeed[coordinate] * dualLines[point][coordinate]) =
        Zero(field)));;
AssertFi22(Length(zeroSet) = 15424,
    "reference regular seed has unexpected zero count");
AssertFi22(RankMat(dualLines{zeroSet}) = 26,
    "reference regular seed zeros do not span its annihilator");
permutationGroup := Group(permutationGenerators);;
SetSize(permutationGroup, permutationInfo.size);
regularStabilizer := Stabilizer(permutationGroup, zeroSet, OnSets);;
AssertFi22(Size(regularStabilizer) = 1,
    "reference regular seed is not projectively regular");

seedPath := Concatenation(OUT, "/regular_seed.txt");;
output := OutputTextFile(seedPath, false);;
SetPrintFormattingStatus(output, false);
AppendTo(output, Digits(regularSeed), "\n");
CloseStream(output);

Print("FI22_DUAL_ORBIT_EXPORT_V1\n");
Print("matrix_representation=", matrixAtlas.repname, "\n");
Print("permutation_representation=", permutationAtlas.repname, "\n");
Print("matrix_generator_orders=2,13\n");
Print("permutation_generator_orders=2,13\n");
Print("max4=", max4.subgroupname, " max4_order=", max4.size,
    " fixed_point=1\n");
Print("dual_orbit=61776 alignment_edges=", 2 * 61776,
    " distinct=61776 span=27\n");
Print("regular_seed=032133122122121333223313020",
    " zero_count=", Length(zeroSet),
    " zero_rank=26 projective_stabilizer=1\n");
Print("OUTPUT dual_orbit=", dualPath, "\n");
Print("OUTPUT matrix_generators=", matrixPath, "\n");
Print("OUTPUT hermitian_form=", formPath, "\n");
Print("OUTPUT regular_seed=", seedPath, "\n");
Print("Wrote the 61,776-point dual orbit and the module data.\n");
QUIT_GAP(0);
