# Export the two involution-centralizer actions that are not obtained from
# the general cyclic-class program in export_prime_order_spaces.g.
#
# For each class this script:
#   1. evaluates the ATLAS maximal-subgroup SLP on the common standard
#      generators of 3.Fi22 <= GL(27,4) and Fi22 <= Sym(3510);
#   2. computes the unique central involution of the evaluated maximal
#      subgroup and lifts the same word to the matrix representation;
#   3. proves that the evaluated subgroup is the full involution
#      centralizer by commutation and exact order;
#   4. recomputes the nonzero eigenspace basis and the restrictions of
#      every maximal-subgroup generator; and
#   5. writes the resulting eigenspace action for the C++ orbit program.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

ROOT := "sporadic/large/fi22";;
OUT := Concatenation(ROOT, "/generated");;
Exec(Concatenation("mkdir -p ", OUT));

Assert2AB := function(condition, message)
  if not condition then
    Error(Concatenation("FI22 2AB VERIFICATION FAILURE: ", message));
  fi;
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
  return Concatenation(List(vector, entry -> String(Code(entry))));
end;;

WriteLines := function(path, lines)
  local output, line;
  output := OutputTextFile(path, false);
  SetPrintFormattingStatus(output, false);
  for line in lines do
    AppendTo(output, line, "\n");
  od;
  CloseStream(output);
end;;

FixedPointCount := function(permutation, degree)
  return Number([1..degree], point -> point ^ permutation = point);
end;;

IsScalarMatrix := function(matrix, field)
  local dimension, diagonal, i, j;
  dimension := Length(matrix);
  diagonal := matrix[1][1];
  if IsZero(diagonal) then
    return false;
  fi;
  for i in [1..dimension] do
    if matrix[i][i] <> diagonal then
      return false;
    fi;
    for j in [1..dimension] do
      if i <> j and matrix[i][j] <> Zero(field) then
        return false;
      fi;
    od;
  od;
  return true;
end;;

field := GF(4);;
identity27 := IdentityMat(27, field);;
matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
permutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22", NrMovedPoints, 3510);;
Assert2AB(matrixInfo <> fail and permutationInfo <> fail,
    "required ATLAS representations are unavailable");
matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
matrixGenerators := matrixAtlas.generators;;
permutationGenerators := permutationAtlas.generators;;
permutationGroup := Group(permutationGenerators);;

Assert2AB(matrixAtlas.repname = "3F22G1-f4r27aB0",
    "unexpected matrix representation");
Assert2AB(permutationAtlas.repname = "F22G1-p3510B0",
    "unexpected degree-3510 permutation representation");
Assert2AB(matrixAtlas.standardization = permutationAtlas.standardization,
    "matrix and permutation representations use different standardizations");
Assert2AB(List(matrixGenerators, Order) = [2,13],
    "unexpected matrix standard-generator orders");
Assert2AB(List(permutationGenerators, Order) = [2,13],
    "unexpected permutation standard-generator orders");
Assert2AB(Size(permutationGroup) = 64561751654400,
    "degree-3510 generators do not generate Fi22");

table := CharacterTable("Fi22");;
classNames := AtlasClassNames(table);;
centralizerOrders := SizesCentralizers(table);;
classOrders := OrdersClassRepresentatives(table);;

# name, ATLAS maxes position, W1 cyclic position, W1 exponent,
# degree-3510 fixed points, expected nonzero eigenspace dimension.
specifications := [
  ["2A", 1, 1, 3, 694, 21],
  ["2B", 7, 6, 6, 182, 17]
];;

cyclic := AtlasProgram("Fi22", "cyclic");;
Assert2AB(cyclic <> fail, "ATLAS cyclic program is unavailable");
cyclicPermutations := ResultOfStraightLineProgram(
    cyclic.program, permutationGenerators);;
Assert2AB(Length(cyclicPermutations) = 27,
    "unexpected cyclic-program output length");

ExportSpecification := function(specification)
  local name, maximalPosition, cyclicPosition, cyclicExponent,
        expectedFixedPoints, expectedDimension, classPosition,
        expectedCentralizer, cyclicElement, maximal, maximalMatrixGenerators,
        maximalPermutationGenerators, maximalPermutationGroup, actualOrder,
        centre, centreGenerators, centreElement, freeEpimorphism,
        freeGenerators, centreWord, matrixElement, eigenvalues, bases,
        scalar, basisRows, dimensions, spaces, generatorPosition, lift,
        blockPermutation, blockMatrices, block, imageRows, target,
        coefficients, candidateTarget, lines, outputPath, commutator;

  name := specification[1];
  maximalPosition := specification[2];
  cyclicPosition := specification[3];
  cyclicExponent := specification[4];
  expectedFixedPoints := specification[5];
  expectedDimension := specification[6];
  classPosition := Position(classNames, name);
  Assert2AB(classPosition <> fail,
      Concatenation(name, ": class is absent from the Fi22 character table"));
  Assert2AB(classOrders[classPosition] = 2,
      Concatenation(name, ": character-table class is not involutory"));
  expectedCentralizer := centralizerOrders[classPosition];

  cyclicElement :=
      cyclicPermutations[cyclicPosition]^cyclicExponent;
  Assert2AB(Order(cyclicElement) = 2,
      Concatenation(name, ": W1 cyclic representative is not involutory"));
  Assert2AB(FixedPointCount(cyclicElement, 3510) = expectedFixedPoints,
      Concatenation(name,
          ": W1 cyclic representative has the wrong fixed-point count"));

  maximal := AtlasProgram("Fi22", "maxes", maximalPosition);
  Assert2AB(maximal <> fail,
      Concatenation(name, ": ATLAS maximal-subgroup program is unavailable"));
  Assert2AB(maximal.size = expectedCentralizer,
      Concatenation(name,
          ": ATLAS maximal-subgroup order is not the class centralizer order"));
  maximalMatrixGenerators := ResultOfStraightLineProgram(
      maximal.program, matrixGenerators);
  maximalPermutationGenerators := ResultOfStraightLineProgram(
      maximal.program, permutationGenerators);
  Assert2AB(Length(maximalMatrixGenerators) = 2 and
      Length(maximalPermutationGenerators) = 2,
      Concatenation(name,
          ": maximal-subgroup SLP did not return two generators"));

  maximalPermutationGroup := Group(maximalPermutationGenerators);
  actualOrder := Size(maximalPermutationGroup);
  Assert2AB(actualOrder = expectedCentralizer,
      Concatenation(name,
          ": evaluated maximal-subgroup generators have the wrong order"));
  centre := Centre(maximalPermutationGroup);
  Assert2AB(Size(centre) = 2,
      Concatenation(name,
          ": evaluated maximal subgroup does not have centre of order two"));
  centreGenerators := GeneratorsOfGroup(centre);
  Assert2AB(Length(centreGenerators) = 1,
      Concatenation(name, ": centre does not have one canonical generator"));
  centreElement := centreGenerators[1];
  Assert2AB(Order(centreElement) = 2,
      Concatenation(name, ": central element is not an involution"));
  Assert2AB(FixedPointCount(centreElement, 3510) = expectedFixedPoints,
      Concatenation(name,
          ": maximal-subgroup centre lies in the wrong Fi22 class"));
  Assert2AB(ForAll(maximalPermutationGenerators,
      generator -> centreElement * generator =
          generator * centreElement),
      Concatenation(name,
          ": maximal-subgroup generator does not centralize its centre"));

  # Since this subgroup centralizes the involution and its exact order is
  # the character-table centralizer order, it is the full centralizer.
  freeEpimorphism := EpimorphismFromFreeGroup(maximalPermutationGroup);
  SetIsSurjective(freeEpimorphism, true);
  freeGenerators := GeneratorsOfGroup(Source(freeEpimorphism));
  centreWord := PreImagesRepresentative(
      freeEpimorphism, centreElement);
  matrixElement := MappedWord(
      centreWord, freeGenerators, maximalMatrixGenerators);
  Assert2AB(Order(matrixElement) in [2,6],
      Concatenation(name,
          ": lifted projective involution has unexpected matrix order"));
  for lift in maximalMatrixGenerators do
    commutator := Comm(matrixElement, lift);
    Assert2AB(IsScalarMatrix(commutator, field),
        Concatenation(name,
            ": matrix lifts do not centralize the involution projectively"));
  od;

  eigenvalues := [One(field), Z(4), Z(4)^2];
  bases := [];
  for scalar in eigenvalues do
    basisRows := NullspaceMat(matrixElement - scalar * identity27);
    if Length(basisRows) > 0 then
      Add(bases, basisRows);
    fi;
  od;
  dimensions := List(bases, Length);
  Assert2AB(dimensions = [expectedDimension],
      Concatenation(name, ": unexpected nonzero eigenspace dimensions"));
  spaces := List(bases, rows -> VectorSpace(field, rows));

  blockMatrices := [];
  for generatorPosition in [1..Length(maximalMatrixGenerators)] do
    lift := maximalMatrixGenerators[generatorPosition];
    blockPermutation := [];
    Add(blockMatrices, []);
    for block in [1..Length(bases)] do
      imageRows := List(bases[block], row -> row * lift);
      target := fail;
      coefficients := fail;
      for candidateTarget in [1..Length(bases)] do
        if ForAll(imageRows, row -> row in spaces[candidateTarget]) then
          target := candidateTarget;
          coefficients := List(imageRows, row ->
              Coefficients(Basis(spaces[candidateTarget],
                  bases[candidateTarget]), row));
          break;
        fi;
      od;
      Assert2AB(target <> fail,
          Concatenation(name,
              ": centralizer generator does not preserve the eigenline set"));
      Assert2AB(dimensions[block] = dimensions[target] and
          RankMat(coefficients) = dimensions[block],
          Concatenation(name, ": singular restricted generator"));
      Add(blockPermutation, target);
      Add(blockMatrices[generatorPosition], coefficients);
    od;
    Assert2AB(blockPermutation = [1],
        Concatenation(name, ": unexpected eigenspace block permutation"));
  od;

  lines := [
    Concatenation("CLASS ", name),
    Concatenation("CENT ", String(expectedCentralizer)),
    "NBLOCKS 1",
    Concatenation("DIMS ", String(expectedDimension)),
    Concatenation("NGENS ",
        String(Length(maximalMatrixGenerators))),
    "BASIS 1"
  ];
  for basisRows in bases[1] do
    Add(lines, Digits(basisRows));
  od;
  for generatorPosition in [1..Length(maximalMatrixGenerators)] do
    Add(lines, Concatenation("GEN ", String(generatorPosition),
        " BLOCKPERM 1"));
    for basisRows in blockMatrices[generatorPosition][1] do
      Add(lines, Digits(basisRows));
    od;
  od;

  outputPath := Concatenation(OUT, "/cls_", name, ".txt");
  WriteLines(outputPath, lines);
  Print("CLASS ", name,
      " atlas_maxes=", maximalPosition,
      " subgroup=", maximal.subgroupname,
      " subgroup_order=", actualOrder,
      " centralizer_order=", expectedCentralizer,
      " fixed_points_3510=", expectedFixedPoints,
      " matrix_element_order=", Order(matrixElement),
      " eigendim=", expectedDimension,
      " generators=", Length(maximalMatrixGenerators),
      " output=", outputPath, "\n");
end;;

Print("Fi22 involution spaces\n");
Print("matrix_representation=", matrixAtlas.repname, "\n");
Print("permutation_representation=", permutationAtlas.repname, "\n");
for specification in specifications do
  ExportSpecification(specification);
od;
Print("Wrote the actions for 2A and 2B.\n");
QUIT_GAP(0);
