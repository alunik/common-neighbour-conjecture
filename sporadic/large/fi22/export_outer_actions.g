# Exact fixed-space centralizer actions for the three outer involution
# classes of 3.Fi22.2 on its 54-dimensional GF(2) module.
#
# For an outer involution xbar in Fi22.2, an arbitrary matrix lift of a
# generator of C_Fi22.2(xbar) need only commute with a chosen lift x modulo
# the normal kernel C3.  Since x inverts C3, there is a unique kernel
# adjustment that commutes with x exactly.  These adjusted lifts generate
# C_H(x) for H = 3.Fi22.2 (the centralizer maps isomorphically onto the
# quotient centralizer), and hence their restrictions to Fix_V(x) give the exact
# action whose vector orbits meet every 3.Fi22.2-orbit stabilized by this
# outer class.
#
# The exporter also reconstructs the exact restriction-of-scalars bridge
# from the ATLAS 54-dimensional coordinates to the reference 27-dimensional
# GF(4) coordinates.  Each fixed-space basis row is written in both systems,
# so outer_orbit_bfs.cpp can emit directly consumable representatives.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

ROOT := "sporadic/large/fi22";;
OUT := Concatenation(ROOT, "/generated/outer_fixed");;
Exec(Concatenation("mkdir -p ", OUT));

AssertOuter := function(condition, message)
  if not condition then
    Error(Concatenation("FI22 OUTER FIXED-ACTION FAILURE: ", message));
  fi;
end;;

CanonicalRows := function(rows)
  return BasisVectors(CanonicalBasis(VectorSpace(GF(2), rows)));
end;;

BinaryCode := function(entry)
  if IsZero(entry) then
    return "0";
  fi;
  return "1";
end;;

BinaryDigits := function(vector)
  return Concatenation(List(vector, BinaryCode));
end;;

F4Code := function(entry)
  if IsZero(entry) then
    return 0;
  elif entry = One(GF(4)) then
    return 1;
  elif entry = Z(4) then
    return 2;
  elif entry = Z(4)^2 then
    return 3;
  fi;
  Error("entry is outside the reference GF(4) encoding");
end;;

F4Digits := function(vector)
  return Concatenation(List(vector, entry -> String(F4Code(entry))));
end;;

DecodeF4Digits := function(digits)
  local decode;
  decode := [Zero(GF(4)), One(GF(4)), Z(4), Z(4)^2];
  return List(digits, character -> decode[IntChar(character) - 47]);
end;;

field2 := GF(2);;
field4 := GF(4);;
identity54 := IdentityMat(54, field2);;

matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22.2", Characteristic, 2, Dimension, 54);;
permutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22.2", NrMovedPoints, 3510);;
coreMatrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
AssertOuter(ForAll([matrixInfo, permutationInfo, coreMatrixInfo],
    info -> info <> fail),
    "required ATLAS representations are unavailable");

matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
coreMatrixAtlas := AtlasGenerators(coreMatrixInfo.identifier);;
matrixGenerators := matrixAtlas.generators;;
permutationGenerators := permutationAtlas.generators;;
coreMatrixGenerators := coreMatrixAtlas.generators;;
AssertOuter(matrixAtlas.repname = "3F22d2G1-f2r54B0" and
    permutationAtlas.repname = "F22d2G1-p3510B0" and
    coreMatrixAtlas.repname = "3F22G1-f4r27aB0",
    "unexpected ATLAS representations");
AssertOuter(List(matrixGenerators, Order) = [2,18] and
    List(permutationGenerators, Order) = [2,18] and
    List(coreMatrixGenerators, Order) = [2,13] and
    Order(Product(matrixGenerators)) = 42 and
    Order(Product(permutationGenerators)) = 42,
    "unexpected standard-generator orders");

permutationGroup := Group(permutationGenerators);;
SetSize(permutationGroup, permutationInfo.size);;
AssertOuter(Size(permutationGroup) = 129123503308800,
    "unexpected Fi22.2 quotient order");
corePermutationGroup := DerivedSubgroup(permutationGroup);;
AssertOuter(Size(corePermutationGroup) = 64561751654400,
    "unexpected derived-subgroup order");

kernelProgram := AtlasProgram(
    "3.Fi22.2", matrixAtlas.standardization,
    "kernel", "Fi22.2");;
AssertOuter(kernelProgram <> fail,
    "normal-kernel program is unavailable");
kernelGenerator := ResultOfStraightLineProgram(
    kernelProgram.program, matrixGenerators)[1];;
AssertOuter(Order(kernelGenerator) = 3 and
    54 - RankMat(kernelGenerator - identity54) = 0 and
    54 - RankMat(kernelGenerator^2 - identity54) = 0,
    "normal C3 kernel is not fixed-point-free");

# Reconstruct the restriction-of-scalars bridge exactly, allowing the
# standard generators of the index-two core to differ by kernel powers.
maximalCore := AtlasProgram("3.Fi22.2", "maxes", 1);;
AssertOuter(maximalCore <> fail and
    maximalCore.subgroupname = "3.Fi22",
    "index-two core program is unavailable");
core54Generators := ResultOfStraightLineProgram(
    maximalCore.program, matrixGenerators);;
basis4over2 := CanonicalBasis(field4);;
blownCoreGenerators := List(coreMatrixGenerators,
    generator -> ImmutableMatrix(
      field2, BlownUpMat(basis4over2, generator)));;
sourceModule := GModuleByMats(blownCoreGenerators, field2);;
intertwiner := fail;;
coreAdjustment := fail;;
for firstPower in [0..2] do
  for secondPower in [0..2] do
    adjustedCore := [
      core54Generators[1] * kernelGenerator^firstPower,
      core54Generators[2] * kernelGenerator^secondPower
    ];
    if List(adjustedCore, Order) <> [2,13] then
      continue;
    fi;
    candidateIntertwiner := MTX.IsomorphismModules(
        sourceModule, GModuleByMats(adjustedCore, field2));
    if candidateIntertwiner <> fail then
      intertwiner := candidateIntertwiner;
      coreAdjustment := [firstPower, secondPower];
      core54Generators := adjustedCore;
      break;
    fi;
  od;
  if intertwiner <> fail then break; fi;
od;
AssertOuter(intertwiner <> fail and RankMat(intertwiner) = 54,
    "restriction-of-scalars modules are not isomorphic");
if ForAll([1..2], position ->
    blownCoreGenerators[position] * intertwiner =
      intertwiner * core54Generators[position]) then
  bridgeOrientation := "source_times_T";
elif ForAll([1..2], position ->
    core54Generators[position] * intertwiner =
      intertwiner * blownCoreGenerators[position]) then
  bridgeOrientation := "target_times_T";
else
  bridgeOrientation := fail;
fi;
AssertOuter(bridgeOrientation <> fail,
    "restriction bridge has no valid orientation");

# Build an explicit GF(2) coordinate inverse for BlownUpVector.  Row
# 2j-1 is the blown image of e_j and row 2j is that of Z(4)e_j.
blowRows := [];;
for coordinate in [1..27] do
  unit := List([1..27], position -> Zero(field4));
  unit[coordinate] := One(field4);
  Add(blowRows, BlownUpVector(basis4over2, unit));
  unit[coordinate] := Z(4);
  Add(blowRows, BlownUpVector(basis4over2, unit));
od;
blowMatrix := ImmutableMatrix(field2, blowRows);;
AssertOuter(RankMat(blowMatrix) = 54,
    "explicit GF(4)-to-GF(2) blow-up basis is singular");
blowInverse := blowMatrix^-1;;

TargetFromF4 := function(vector)
  local source;
  source := BlownUpVector(basis4over2, vector);
  if bridgeOrientation = "source_times_T" then
    return source * intertwiner;
  fi;
  return source * intertwiner^-1;
end;;

F4FromTarget := function(vector)
  local source, coefficients, result, coordinate, code;
  if bridgeOrientation = "source_times_T" then
    source := vector * intertwiner^-1;
  else
    source := vector * intertwiner;
  fi;
  coefficients := source * blowInverse;
  result := [];
  for coordinate in [1..27] do
    code := 0;
    if not IsZero(coefficients[2 * coordinate - 1]) then
      code := code + 1;
    fi;
    if not IsZero(coefficients[2 * coordinate]) then
      code := code + 2;
    fi;
    Add(result,
      [Zero(field4), One(field4), Z(4), Z(4)^2][code + 1]);
  od;
  return result;
end;;

testF4 := ListWithIdenticalEntries(27, Zero(field4));;
testF4[1] := One(field4);;
AssertOuter(F4FromTarget(TargetFromF4(testF4)) = testF4,
    "restriction bridge does not round trip");

freeEpimorphism := EpimorphismFromFreeGroup(permutationGroup);;
SetIsSurjective(freeEpimorphism, true);;
freeGenerators := GeneratorsOfGroup(Source(freeEpimorphism));;
AssertOuter(Length(freeGenerators) = 2 and
    List(freeGenerators, generator ->
      ImagesRepresentative(freeEpimorphism, generator)) =
        permutationGenerators,
    "free-group epimorphism is not aligned with the standard generators");

classes := ConjugacyClasses(permutationGroup);;
outerRecords := [];;
for classPosition in [1..Length(classes)] do
  quotientRepresentative := Representative(classes[classPosition]);
  if Order(quotientRepresentative) = 2 and
      not quotientRepresentative in corePermutationGroup then
    quotientCentralizer := Centralizer(
        permutationGroup, quotientRepresentative);
    Add(outerRecords, rec(
      position := classPosition,
      classSize := Size(classes[classPosition]),
      centralizerOrder := Size(quotientCentralizer),
      representative := quotientRepresentative,
      centralizer := quotientCentralizer
    ));
  fi;
od;
expectedPairs := [
  [61776,2090188800],
  [22239360,5806080],
  [19459440,6635520]
];;
AssertOuter(Length(outerRecords) = 3 and
    Set(List(outerRecords,
      record -> [record.classSize, record.centralizerOrder])) =
        Set(expectedPairs),
    Concatenation("unexpected outer involution inventory ",
      String(List(outerRecords, record ->
        [record.position, record.classSize,
          record.centralizerOrder]))));
orderedOuterRecords := [];;
for expectedPair in expectedPairs do
  matchingRecord := First(outerRecords, record ->
      record.classSize = expectedPair[1] and
      record.centralizerOrder = expectedPair[2]);
  AssertOuter(matchingRecord <> fail,
      "stable outer involution class pair was not found");
  Add(orderedOuterRecords, matchingRecord);
od;
outerPositions := List(orderedOuterRecords, record -> record.position);;

for outerId in [1..Length(orderedOuterRecords)] do
  outerRecord := orderedOuterRecords[outerId];
  classPosition := outerRecord.position;
  quotientRepresentative := outerRecord.representative;
  representativeWord := PreImagesRepresentative(
      freeEpimorphism, quotientRepresentative);
  rawMatrixRepresentative := MappedWord(
      representativeWord, freeGenerators, matrixGenerators);
  orderTwoLifts := Filtered([0..2], power ->
      Order(kernelGenerator^power * rawMatrixRepresentative) = 2);
  AssertOuter(Length(orderTwoLifts) >= 1,
      Concatenation("outer class ", String(outerId),
        " has no order-two matrix lift"));
  representativePower := orderTwoLifts[1];
  matrixRepresentative :=
    kernelGenerator^representativePower * rawMatrixRepresentative;
  AssertOuter(matrixRepresentative * kernelGenerator *
      matrixRepresentative^-1 = kernelGenerator^-1,
      Concatenation("outer class ", String(outerId),
        " does not invert the normal kernel"));

  fixedBasis := CanonicalRows(
      NullspaceMat(matrixRepresentative - identity54));
  AssertOuter(Length(fixedBasis) = 27 and
      ForAll(fixedBasis,
        row -> row * matrixRepresentative = row),
      Concatenation("outer class ", String(outerId),
        " has an invalid fixed-space basis"));
  fixedSpace := VectorSpace(field2, fixedBasis);
  fixedBasisObject := Basis(fixedSpace, fixedBasis);
  fixedF4Basis := List(fixedBasis, F4FromTarget);
  AssertOuter(ForAll([1..27], position ->
      TargetFromF4(fixedF4Basis[position]) = fixedBasis[position]),
      Concatenation("outer class ", String(outerId),
        " fixed basis fails the restriction-bridge round trip"));

  quotientCentralizer := outerRecord.centralizer;
  centralizerOrder := outerRecord.centralizerOrder;
  AssertOuter([outerRecord.classSize, centralizerOrder] =
      expectedPairs[outerId],
      Concatenation("outer class ", String(outerId),
        " stable class pair changed"));
  quotientCentralizerGenerators :=
    GeneratorsOfGroup(quotientCentralizer);
  adjustedMatrixGenerators := [];;
  adjustmentPowers := [];;
  restrictedMatrices := [];;
  for quotientGenerator in quotientCentralizerGenerators do
    generatorWord := PreImagesRepresentative(
        freeEpimorphism, quotientGenerator);
    rawMatrixGenerator := MappedWord(
        generatorWord, freeGenerators, matrixGenerators);
    commutingPowers := Filtered([0..2], power ->
        (kernelGenerator^power * rawMatrixGenerator) *
            matrixRepresentative =
          matrixRepresentative *
            (kernelGenerator^power * rawMatrixGenerator));
    AssertOuter(Length(commutingPowers) = 1,
        Concatenation("outer class ", String(outerId),
          " centralizer lift has nonunique kernel adjustment"));
    adjustmentPower := commutingPowers[1];
    adjustedMatrixGenerator :=
      kernelGenerator^adjustmentPower * rawMatrixGenerator;
    Add(adjustedMatrixGenerators, adjustedMatrixGenerator);
    Add(adjustmentPowers, adjustmentPower);
    images := List(fixedBasis,
        row -> row * adjustedMatrixGenerator);
    AssertOuter(ForAll(images, row -> row in fixedSpace),
        Concatenation("outer class ", String(outerId),
          " adjusted centralizer generator leaves the fixed space"));
    restrictedMatrix := List(images,
        row -> Coefficients(fixedBasisObject, row));
    AssertOuter(RankMat(restrictedMatrix) = 27 and
        ForAll([1..27], rowPosition ->
          Sum([1..27], coefficientPosition ->
              restrictedMatrix[rowPosition][coefficientPosition] *
                fixedBasis[coefficientPosition]) =
            images[rowPosition]),
        Concatenation("outer class ", String(outerId),
          " has an invalid restricted generator"));
    Add(restrictedMatrices,
        ImmutableMatrix(field2, restrictedMatrix));
  od;

  outputPath := Concatenation(
      OUT, "/outer_", String(outerId), ".txt");
  output := OutputTextFile(outputPath, false);
  SetPrintFormattingStatus(output, false);
  AppendTo(output, "SCHEMA FI22_OUTER_ACTION_V1\n");
  AppendTo(output, "CLASS OUT", outerId, "\n");
  AppendTo(output, "CLASS_POSITION ", classPosition, "\n");
  AppendTo(output, "CLASS_SIZE ", outerRecord.classSize, "\n");
  AppendTo(output, "CENT ", centralizerOrder, "\n");
  AppendTo(output, "DIM 27\n");
  AppendTo(output, "AMBIENT 54\n");
  AppendTo(output, "NGENS ", Length(restrictedMatrices), "\n");
  AppendTo(output, "REPRESENTATIVE_POWER ", representativePower, "\n");
  AppendTo(output, "REPRESENTATIVE_MATRIX\n");
  for row in matrixRepresentative do
    AppendTo(output, BinaryDigits(row), "\n");
  od;
  AppendTo(output, "BASIS_H\n");
  for row in fixedBasis do
    AppendTo(output, BinaryDigits(row), "\n");
  od;
  AppendTo(output, "BASIS_F4\n");
  for row in fixedF4Basis do
    AppendTo(output, F4Digits(row), "\n");
  od;
  for generatorPosition in [1..Length(restrictedMatrices)] do
    AppendTo(output, "GEN ", generatorPosition,
        " ADJUST ", adjustmentPowers[generatorPosition], "\n");
    for row in restrictedMatrices[generatorPosition] do
      AppendTo(output, BinaryDigits(row), "\n");
    od;
  od;
  CloseStream(output);

  Print("OUTER class_id=OUT", outerId,
      " class_position=", classPosition,
      " class_size=", outerRecord.classSize,
      " centralizer=", centralizerOrder,
      " representative_kernel_power=", representativePower,
      " fixed_dimension=", Length(fixedBasis),
      " centralizer_generators=", Length(restrictedMatrices),
      " adjustment_powers=", adjustmentPowers,
      " output=", outputPath, "\n");
od;

Print("Fi22 outer-involution actions\n");
Print("matrix_representation=", matrixAtlas.repname,
    " permutation_representation=", permutationAtlas.repname,
    " core_representation=", coreMatrixAtlas.repname, "\n");
Print("restriction_bridge_rank=", RankMat(intertwiner),
    " orientation=", bridgeOrientation,
    " core_adjustment=", coreAdjustment, "\n");
Print("normal_kernel_order=3 fixed_dimensions=0,0\n");
Print("outer_class_positions=", outerPositions,
    " outer_classes=", Length(outerPositions), "\n");
Print("Wrote all three outer-involution actions.\n");
QUIT_GAP(0);
