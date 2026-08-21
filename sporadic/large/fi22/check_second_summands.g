# Exact downstream checker for the H-native Fischer-color results, where
# H = 3.Fi22.2 <= GL(54,2).
#
# The C++ stage exports singleton cells of an H-equivariant refinement of the
# degree-3510 Fischer action.  Hence the quotient image of a vector stabilizer
# fixes every exported singleton.  This script independently:
#   * replays each first summand from the reference GF(4)^27 seed;
#   * checks target + first summand = second summand;
#   * applies the exact GF(4)^27 -> GF(2)^54 bridge;
#   * computes the pointwise singleton container in Fi22.2;
#   * lifts every element of that container to 3.Fi22.2 and tests all three
#     elements in its normal-C3 fibre on the original 54-bit vector.
#
# A row is accepted only when the exact vector stabilizer contains precisely
# the identity matrix.  The normal C3 is not assumed central in H.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

AssertWLExact := function(condition, message)
  if not condition then
    Error(Concatenation("FI22.2 WL EXACT FAILURE: ", message));
  fi;
end;;

root := "sporadic/large/fi22";;
if IsBound(GAPInfo.SystemEnvironment.FI22_REFINED_INPUT) then
  resultPath := GAPInfo.SystemEnvironment.FI22_REFINED_INPUT;;
else
  resultPath := Concatenation(
      root, "/generated/refined_candidates.tsv");;
fi;
if IsBound(GAPInfo.SystemEnvironment.FI22_STABILISER_OUTPUT) then
  outputPath := GAPInfo.SystemEnvironment.FI22_STABILISER_OUTPUT;;
else
  outputPath := Concatenation(
      root, "/generated/checked_sums.tsv");;
fi;
if IsBound(GAPInfo.SystemEnvironment.FI22_CONTAINER_LIMIT) then
  containerCap := Int(
      GAPInfo.SystemEnvironment.FI22_CONTAINER_LIMIT);;
else
  containerCap := 1000000;;
fi;
AssertWLExact(containerCap >= 1, "invalid pointwise-container cap");

field2 := GF(2);;
field4 := GF(4);;
identity54 := IdentityMat(54, field2);;

DecodeF4Digits := function(digits)
  local decode;
  AssertWLExact(Length(digits) = 27,
      "a GF(4) vector does not have length 27");
  decode := [Zero(field4), One(field4), Z(4), Z(4)^2];
  AssertWLExact(ForAll(digits, character ->
      IntChar(character) >= 48 and IntChar(character) <= 51),
      "a GF(4) digit is outside 0,1,2,3");
  return List(digits,
      character -> decode[IntChar(character) - 47]);
end;;

F4Code := function(entry)
  if IsZero(entry) then
    return 0;
  elif entry = One(field4) then
    return 1;
  elif entry = Z(4) then
    return 2;
  elif entry = Z(4)^2 then
    return 3;
  fi;
  Error("FI22.2 WL EXACT FAILURE: entry outside reference GF(4) encoding");
end;;

DecodeBits54 := function(digits)
  AssertWLExact(Length(digits) = 54 and
      ForAll(digits, character ->
        IntChar(character) = 48 or IntChar(character) = 49),
      "a binary bridge vector is malformed");
  return List(digits,
      character -> (IntChar(character) - 48) * One(field2));
end;;

matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22.2", Characteristic, 2, Dimension, 54);;
permutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22.2", NrMovedPoints, 3510);;
coreInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
AssertWLExact(ForAll([matrixInfo, permutationInfo, coreInfo],
    info -> info <> fail),
    "required ATLAS representations are unavailable");

matrixAtlas := AtlasGenerators(matrixInfo.identifier);;
permutationAtlas := AtlasGenerators(permutationInfo.identifier);;
coreAtlas := AtlasGenerators(coreInfo.identifier);;
h54 := matrixAtlas.generators;;
p3510Generators := permutationAtlas.generators;;
g27 := coreAtlas.generators;;
AssertWLExact(matrixAtlas.repname = "3F22d2G1-f2r54B0" and
    permutationAtlas.repname = "F22d2G1-p3510B0" and
    coreAtlas.repname = "3F22G1-f4r27aB0",
    "unexpected ATLAS representations");
AssertWLExact(List(h54, Order) = [2,18] and
    List(p3510Generators, Order) = [2,18] and
    List(g27, Order) = [2,13] and
    Order(Product(h54)) = 42 and
    Order(Product(p3510Generators)) = 42,
    "unexpected standard generators");

hGroup := Group(h54);;
SetSize(hGroup, matrixInfo.size);;
p3510 := Group(p3510Generators);;
SetSize(p3510, permutationInfo.size);;
AssertWLExact(Size(hGroup) = 387370509926400 and
    Size(p3510) = 129123503308800,
    "unexpected matrix or quotient group order");

kernelProgram := AtlasProgram(
    "3.Fi22.2", matrixAtlas.standardization,
    "kernel", "Fi22.2");;
AssertWLExact(kernelProgram <> fail,
    "normal-kernel program is unavailable");
kernelGenerator := ResultOfStraightLineProgram(
    kernelProgram.program, h54)[1];;
AssertWLExact(Order(kernelGenerator) = 3 and
    54 - RankMat(kernelGenerator - identity54) = 0 and
    54 - RankMat(kernelGenerator^2 - identity54) = 0,
    "normal C3 kernel is not fixed-point-free");

freeEpimorphism := EpimorphismFromFreeGroup(p3510);;
SetIsSurjective(freeEpimorphism, true);;
freeGenerators := GeneratorsOfGroup(Source(freeEpimorphism));;
AssertWLExact(Length(freeGenerators) = 2 and
    List(freeGenerators, generator ->
      ImagesRepresentative(freeEpimorphism, generator)) =
        p3510Generators,
    "free-group epimorphism is not aligned with standard generators");

# Load the bridge table produced by export_fischer_data.g.
bridgeRows := Filtered(SplitString(StringFile(Concatenation(
    root, "/generated/fischer/f4_to_h54.tsv")), "\n"),
    row -> Length(row) > 0);;
AssertWLExact(Length(bridgeRows) = 82 and
    bridgeRows[1] = "coordinate\tdigit\th54_bits",
    "restriction bridge table is missing or malformed");
bridgeContributions := [];;
for coordinate in [1..27] do
  bridgeContributions[coordinate] := [];;
od;
for row in bridgeRows{[2..Length(bridgeRows)]} do
  fields := SplitString(row, "\t");
  AssertWLExact(Length(fields) = 3,
      "malformed restriction bridge row");
  coordinate := Int(fields[1]);
  digit := Int(fields[2]);
  AssertWLExact(coordinate in [1..27] and digit in [1..3] and
      not IsBound(bridgeContributions[coordinate][digit]),
      "duplicate or out-of-range restriction bridge row");
  bridgeContributions[coordinate][digit] :=
      DecodeBits54(fields[3]);
od;
AssertWLExact(ForAll([1..27], coordinate ->
    ForAll([1..3], digit ->
      IsBound(bridgeContributions[coordinate][digit]))),
    "restriction bridge table is incomplete");
AssertWLExact(ForAll([1..27], coordinate ->
      bridgeContributions[coordinate][3] =
        bridgeContributions[coordinate][1] +
        bridgeContributions[coordinate][2]) and
    RankMat(Concatenation(List([1..27], coordinate ->
      [bridgeContributions[coordinate][1],
       bridgeContributions[coordinate][2]]))) = 54,
    "restriction bridge table is not additive and invertible");

BridgeF4 := function(vector)
  local result, coordinate, code;
  result := ListWithIdenticalEntries(54, Zero(field2));
  for coordinate in [1..27] do
    code := F4Code(vector[coordinate]);
    if code <> 0 then
      result := result + bridgeContributions[coordinate][code];
    fi;
  od;
  return result;
end;;

ReplayCoreWord := function(seed, word)
  local result, character, generatorPosition;
  result := seed;
  for character in word do
    generatorPosition := IntChar(character) - 48;
    AssertWLExact(generatorPosition in [1,2],
        "a witness word contains a symbol other than 1 or 2");
    result := result * g27[generatorPosition];
  od;
  return result;
end;;

seedLines := Filtered(SplitString(StringFile(Concatenation(
    root, "/generated/regular_seed.txt")), "\n", "\r"),
    line -> Length(line) > 0);;
AssertWLExact(Length(seedLines) = 1 and Length(seedLines[1]) = 27,
    "regular seed file is malformed");
referenceSeedDigits := seedLines[1];;
referenceSeed := DecodeF4Digits(referenceSeedDigits);;
coreMaxProgram := AtlasProgram("3.Fi22.2", "maxes", 1);;
AssertWLExact(coreMaxProgram <> fail,
    "index-two core maximal program is unavailable");
rawK54 := ResultOfStraightLineProgram(
    coreMaxProgram.program, h54);;
adjustedK54 := [];;
bridgeAdjustment := [];;
for generatorPosition in [1..2] do
  matchingPowers := [];;
  for kernelPower in [0..2] do
    candidateGenerator :=
        rawK54[generatorPosition] * kernelGenerator^kernelPower;
    bridgeMatches := true;
    for coordinate in [1..27] do
      for digit in [1..3] do
        basisVector := ListWithIdenticalEntries(27, Zero(field4));
        basisVector[coordinate] :=
            [One(field4), Z(4), Z(4)^2][digit];
        if BridgeF4(basisVector * g27[generatorPosition]) <>
            bridgeContributions[coordinate][digit] *
                candidateGenerator then
          bridgeMatches := false;
          break;
        fi;
      od;
      if not bridgeMatches then
        break;
      fi;
    od;
    if bridgeMatches then
      Add(matchingPowers, kernelPower);
    fi;
  od;
  AssertWLExact(Length(matchingPowers) = 1,
      "bridge/core-generator adjustment is not unique");
  Add(bridgeAdjustment, matchingPowers[1]);
  Add(adjustedK54, rawK54[generatorPosition] *
      kernelGenerator^matchingPowers[1]);
od;
AssertWLExact(List(adjustedK54, Order) = [2,13] and
    ForAll([1..2], position ->
      BridgeF4(referenceSeed * g27[position]) =
        BridgeF4(referenceSeed) * adjustedK54[position]),
    "bridge/core-generator alignment check failed");

resultRows := Filtered(
    SplitString(StringFile(resultPath), "\n"),
    row -> Length(row) > 0 and row[1] <> '#');;
AssertWLExact(Length(resultRows) >= 1,
    "the WL result has no data rows");

rowsByTarget := [];;
targetIds := [];;
for row in resultRows do
  fields := Filtered(SplitString(row, " "),
      field -> Length(field) > 0);
  AssertWLExact(Length(fields) = 19,
      "a WL result row does not have 19 fields");
  tid := Int(fields[1]);
  AssertWLExact(tid >= 0, "a target identifier is negative");
  if not IsBound(rowsByTarget[tid + 1]) then
    rowsByTarget[tid + 1] := [];;
    Add(targetIds, tid);
  fi;
  Add(rowsByTarget[tid + 1], fields);
od;
AssertWLExact(targetIds = Set(targetIds),
    "target rows are not grouped in increasing identifier order");
AssertWLExact(ForAll(targetIds, tid -> tid in [0..6105]),
    "a target identifier lies outside the expected range");
CandidateBetter := function(left, right)
  local leftSingletons, rightSingletons;
  leftSingletons := Int(left[17]);
  rightSingletons := Int(right[17]);
  if leftSingletons <> rightSingletons then
    return leftSingletons > rightSingletons;
  fi;
  if left[6] = "FAIL" then
    return false;
  elif right[6] = "FAIL" then
    return true;
  fi;
  return Int(left[6]) < Int(right[6]);
end;;
for tid in targetIds do
  Sort(rowsByTarget[tid + 1], CandidateBetter);
od;

output := OutputTextFile(outputPath, false);;
SetPrintFormattingStatus(output, false);
AppendTo(output,
    "tid\tclass\ttarget\tsurvivor\tword_length",
    "\tsingletons\tused_singletons",
    "\tquotient_container\tlifts_tested\texact_stabilizer\n");

accepted := 0;;
maxContainer := 0;;
totalLiftsTested := 0;;
for tid in targetIds do
  targetAccepted := false;
  referenceClass := rowsByTarget[tid + 1][1][2];
  referenceCoordinates := rowsByTarget[tid + 1][1][4];
  referenceTarget := rowsByTarget[tid + 1][1][5];
  for fields in rowsByTarget[tid + 1] do
    className := fields[2];
    targetDigits := fields[5];
    word := fields[10];
    rDigits := fields[11];
    uDigits := fields[12];
    singletonCount := Int(fields[17]);
    AssertWLExact(className = referenceClass and
        fields[4] = referenceCoordinates and
        targetDigits = referenceTarget,
        "candidate rows disagree on target metadata");
    if fields[6] = "FAIL" or singletonCount < 1 or
        fields[19] = "-" then
      Print("REJECT tid=", tid, " survivor=", fields[6],
          " reason=no_singleton\n");
      continue;
    fi;
    survivor := Int(fields[6]);
    singletonPoints := List(SplitString(fields[19], ","), Int);
    AssertWLExact(Length(singletonPoints) = singletonCount and
        Length(Set(singletonPoints)) = singletonCount and
        ForAll(singletonPoints, point -> point in [1..3510]),
        "the exported singleton list is malformed");

    targetF4 := DecodeF4Digits(targetDigits);
    if tid >= 5363 then
      AssertWLExact(className in ["OUT1", "OUT2", "OUT3"] and
          Length(fields[4]) = 54 and
          DecodeBits54(fields[4]) = BridgeF4(targetF4),
          Concatenation("outer target H54/F4 bridge mismatch at target ",
            String(tid)));
    else
      AssertWLExact(not className in ["OUT1", "OUT2", "OUT3"],
          "an outer class occurs in the inner-target range");
    fi;
    rF4 := DecodeF4Digits(rDigits);
    uF4 := DecodeF4Digits(uDigits);
    AssertWLExact(ReplayCoreWord(referenceSeed, word) = rF4,
        Concatenation("first-summand word replay failed at target ",
          String(tid)));
    AssertWLExact(targetF4 + rF4 = uF4,
        Concatenation("target decomposition failed at target ",
          String(tid)));
    uH := BridgeF4(uF4);
    AssertWLExact(ForAny(uH, entry -> not IsZero(entry)),
        "the second summand is zero");

    container := p3510;
    usedSingletons := 0;
    for point in singletonPoints do
      container := Stabilizer(container, point);
      usedSingletons := usedSingletons + 1;
      if Size(container) = 1 then
        break;
      fi;
    od;
    containerOrder := Size(container);
    if containerOrder > maxContainer then
      maxContainer := containerOrder;
    fi;
    if containerOrder > containerCap then
      Print("REJECT tid=", tid, " survivor=", survivor,
          " singleton_count=", singletonCount,
          " quotient_container=", containerOrder,
          " reason=container_cap\n");
      continue;
    fi;
    AssertWLExact(ForAll(GeneratorsOfGroup(container), generator ->
        ForAll(singletonPoints, point -> point ^ generator = point)),
        "computed container does not fix every singleton");

    exactStabilizer := 0;
    identityFixes := 0;
    liftsTested := 0;
    for quotientElement in Elements(container) do
      quotientWord := PreImagesRepresentative(
          freeEpimorphism, quotientElement);
      rawLift := MappedWord(
          quotientWord, freeGenerators, h54);
      for kernelPower in [0..2] do
        matrixLift := kernelGenerator^kernelPower * rawLift;
        liftsTested := liftsTested + 1;
        if uH * matrixLift = uH then
          exactStabilizer := exactStabilizer + 1;
          if matrixLift = identity54 then
            identityFixes := identityFixes + 1;
          fi;
        fi;
      od;
    od;
    AssertWLExact(liftsTested = 3 * containerOrder,
        "not every lift in the pointwise container was tested");
    totalLiftsTested := totalLiftsTested + liftsTested;
    if exactStabilizer <> 1 or identityFixes <> 1 then
      Print("REJECT tid=", tid, " survivor=", survivor,
          " singleton_count=", singletonCount,
          " used_singletons=", usedSingletons,
          " quotient_container=", containerOrder,
          " lifts_tested=", liftsTested,
          " exact_stabilizer=", exactStabilizer,
          " identity_fixes=", identityFixes,
          " reason=nonregular\n");
      continue;
    fi;

    targetAccepted := true;
    accepted := accepted + 1;
    AppendTo(output, tid, "\t", className, "\t", targetDigits, "\t",
        survivor, "\t", Length(word), "\t", singletonCount, "\t",
        usedSingletons, "\t", containerOrder, "\t", liftsTested,
        "\t", exactStabilizer, "\n");
    Print("ACCEPT tid=", tid, " class=", className,
        " survivor=", survivor,
        " singleton_count=", singletonCount,
        " used_singletons=", usedSingletons,
        " quotient_container=", containerOrder,
        " lifts_tested=", liftsTested,
        " exact_stabilizer=1\n");
    break;
  od;
  AssertWLExact(targetAccepted,
      Concatenation("no H-regular second summand at target ",
        String(tid)));
od;
CloseStream(output);

AssertWLExact(accepted = Length(targetIds),
    "not every target was accepted exactly once");
Print("SUMMARY accepted=", accepted,
    " max_quotient_container=", maxContainer,
    " total_lifts_tested=", totalLiftsTested, "\n");
Print("All supplied targets have a decomposition into two regular vectors.\n");
QUIT_GAP(0);
