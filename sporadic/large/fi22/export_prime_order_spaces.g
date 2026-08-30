# Export centralizer actions on prime-order eigenline sets for
# 3.Fi22 <= GL(27,4). These files are consumed by
# enumerate_inner_orbits.cpp.
#
# The class representatives are obtained from ATLAS straight-line programs
# in common standard generators. A centralizer generator in the faithful
# degree-3510 Fi22 action is lifted by its word in the common standard
# generators.  The lift can commute with the chosen matrix representative
# only modulo the scalar centre, so the exporter explicitly records the
# induced permutation of the nonzero-eigenvalue blocks.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

ROOT := "sporadic/large/fi22";;
OUT := Concatenation(ROOT, "/generated");;
Exec(Concatenation("mkdir -p ", OUT));

AssertFi22 := function(condition, message)
  if not condition then
    Error(Concatenation("FI22 CLASS EXPORT FAILURE: ", message));
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
  return Concatenation(List(vector, x -> String(Code(x))));
end;;

field := GF(4);;
identity27 := IdentityMat(27, field);;
matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
permutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22", NrMovedPoints, 3510);;
AssertFi22(matrixInfo <> fail and permutationInfo <> fail,
    "required ATLAS representations are unavailable");
matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
matrixGenerators := matrixAtlas.generators;;
permutationGenerators := permutationAtlas.generators;;
permutationGroup := Group(permutationGenerators);;
SetSize(permutationGroup, permutationInfo.size);

cyclic := AtlasProgram("Fi22", "cyclic");;
AssertFi22(cyclic <> fail, "ATLAS cyclic program is unavailable");
cyclicMatrices := ResultOfStraightLineProgram(
    cyclic.program, matrixGenerators);;
cyclicPermutations := ResultOfStraightLineProgram(
    cyclic.program, permutationGenerators);;
AssertFi22(Length(cyclicMatrices) = 27 and
    Length(cyclicPermutations) = 27,
    "unexpected cyclic-program output length");

freeEpimorphism := EpimorphismFromFreeGroup(permutationGroup);;
SetIsSurjective(freeEpimorphism, true);
freeGenerators := GeneratorsOfGroup(Source(freeEpimorphism));;

# name, cyclic-output position, power, expected centralizer order,
# expected nonzero eigenspace dimensions, whether only bases are needed.
specifications := [
  ["2C",  2, 3,   1769472, [15],       false],
  ["3A",  2, 2,  19595520, [15,6,6],   false],
  ["3B",  6, 4,   2519424, [9,9,9],    false],
  ["3C",  1, 2,    139968, [9,9,9],    false],
  ["3D",  5, 3,     17496, [9,9,9],    false],
  ["5A", 22, 4,       600, [7],        false],
  ["7A", 17, 2,        42, [3],        true ],
  ["11A",24, 4,        22, [2],        true ],
  ["11B",24, 2,        22, [2],        true ],
  ["13A",16, 1,        13, [3],        true ],
  ["13B",16, 2,        13, [3],        true ]
];;

ExportSpecification := function(specification)
  local name, cyclicPosition, exponent, expectedCentralizer,
        expectedDimensions, tiny, matrixElement, permutationElement,
        centralizer, centralizerGenerators, matrixLifts, eigenvalues,
        bases, dimensions, scalar, basisRows, spaces, outputPath, output,
        generatorPosition, lift, blockPermutation, blockMatrices, block,
        target, imageRows, coefficients, candidateTarget, word, centralOrder;

  name := specification[1];
  cyclicPosition := specification[2];
  exponent := specification[3];
  expectedCentralizer := specification[4];
  expectedDimensions := specification[5];
  tiny := specification[6];
  matrixElement := cyclicMatrices[cyclicPosition]^exponent;
  permutationElement := cyclicPermutations[cyclicPosition]^exponent;
  AssertFi22(Order(permutationElement) =
      Int(Filtered(name, IsDigitChar)),
      Concatenation(name, ": wrong projective element order"));

  eigenvalues := [One(field), Z(4), Z(4)^2];
  bases := [];
  for scalar in eigenvalues do
    # GAP's NullspaceMat returns the left nullspace, exactly the row
    # eigenvectors used by the right matrix action.
    basisRows := NullspaceMat(
        matrixElement - scalar * identity27);
    if Length(basisRows) > 0 then
      Add(bases, basisRows);
    fi;
  od;
  dimensions := List(bases, Length);
  AssertFi22(dimensions = expectedDimensions,
      Concatenation(name, ": unexpected eigenspace dimensions"));
  spaces := List(bases, rows -> VectorSpace(field, rows));

  centralizer := Centralizer(permutationGroup, permutationElement);
  centralOrder := Size(centralizer);
  AssertFi22(centralOrder = expectedCentralizer,
      Concatenation(name, ": unexpected centralizer order"));
  centralizerGenerators := GeneratorsOfGroup(centralizer);
  matrixLifts := [];
  if not tiny then
    for generatorPosition in [1..Length(centralizerGenerators)] do
      word := PreImagesRepresentative(
          freeEpimorphism, centralizerGenerators[generatorPosition]);
      Add(matrixLifts, MappedWord(
          word, freeGenerators, matrixGenerators));
    od;
  fi;

  if tiny then
    outputPath := Concatenation(OUT, "/cls_", name, "_tiny.txt");
  else
    outputPath := Concatenation(OUT, "/cls_", name, ".txt");
  fi;
  output := OutputTextFile(outputPath, false);
  SetPrintFormattingStatus(output, false);
  AppendTo(output, "CLASS ", name, "\n");
  AppendTo(output, "CENT ", centralOrder, "\n");
  AppendTo(output, "NBLOCKS ", Length(bases), "\n");
  AppendTo(output, "DIMS ",
      JoinStringsWithSeparator(List(dimensions, String), " "), "\n");
  if not tiny then
    AppendTo(output, "NGENS ", Length(matrixLifts), "\n");
  fi;
  for block in [1..Length(bases)] do
    AppendTo(output, "BASIS ", block, "\n");
    for basisRows in bases[block] do
      AppendTo(output, Digits(basisRows), "\n");
    od;
  od;

  for generatorPosition in [1..Length(matrixLifts)] do
    lift := matrixLifts[generatorPosition];
    blockPermutation := [];
    blockMatrices := [];
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
      AssertFi22(target <> fail,
          Concatenation(name,
              ": a centralizer lift does not preserve the eigenline union"));
      AssertFi22(dimensions[block] = dimensions[target] and
          RankMat(coefficients) = dimensions[block],
          Concatenation(name, ": singular restricted generator"));
      Add(blockPermutation, target);
      Add(blockMatrices, coefficients);
    od;
    AssertFi22(Set(blockPermutation) = [1..Length(bases)],
        Concatenation(name, ": block action is not a permutation"));
    AppendTo(output, "GEN ", generatorPosition, " BLOCKPERM ",
        JoinStringsWithSeparator(List(blockPermutation, String), " "), "\n");
    for block in [1..Length(bases)] do
      for basisRows in blockMatrices[block] do
        AppendTo(output, Digits(basisRows), "\n");
      od;
    od;
  od;
  CloseStream(output);
  Print("CLASS ", name, " order=", Order(permutationElement),
      " centralizer=", centralOrder,
      " dims=", JoinStringsWithSeparator(List(dimensions, String), ","),
      " generators=", Length(matrixLifts),
      " output=", outputPath, "\n");
end;;

for specification in specifications do
  ExportSpecification(specification);
od;
Print("Wrote ", Length(specifications),
    " prime-order fixed-space actions.\n");
QUIT_GAP(0);
