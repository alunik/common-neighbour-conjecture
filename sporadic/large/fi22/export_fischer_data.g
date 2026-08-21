# Exact equivariant degree-3510 Fischer-action data for
# H = 3.Fi22.2 <= GL(54,2).
#
# The faithful degree-3510 action is an action of H/C3 = Fi22.2.  At each
# point its full preimage in H preserves a 12-dimensional dual subspace.
# We transport an ordered basis of these subspaces (never an unrelated RREF
# basis), and classify evaluation tuples under the full lifted stabilizer,
# including the normal C3 kernel.  Its two nonzero evaluation-tuple orbits
# have exact sizes 2016 and 2079.  Every standard-generator edge is then
# checked exactly: the induced inverse-transpose basis transition belongs to
# this full local value group and hence preserves the entire 0/2016/2079
# label table.  Thus the exported initial coloring, and every
# subsequent Fischer-graph refinement cell, is H-equivariant.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

ROOT := "sporadic/large/fi22";;
OUT := Concatenation(ROOT, "/generated/fischer");;
Exec(Concatenation("mkdir -p ", OUT));

AssertH := function(condition, message)
  if not condition then
    Error(Concatenation("FI22.2 FISCHER EXPORT FAILURE: ", message));
  fi;
end;;

CanonicalRows := function(rows)
  return BasisVectors(CanonicalBasis(VectorSpace(GF(2), rows)));
end;;

Bits := function(vector)
  return Concatenation(List(vector, entry ->
    String(IntFFE(entry))));
end;;

field2 := GF(2);;
field4 := GF(4);;
extensionMatrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22.2", Characteristic, 2, Dimension, 54);;
extensionPermutationInfo := OneAtlasGeneratingSetInfo(
    "Fi22.2", NrMovedPoints, 3510);;
coreMatrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
AssertH(ForAll([
    extensionMatrixInfo, extensionPermutationInfo, coreMatrixInfo],
    info -> info <> fail),
    "required ATLAS representations are unavailable");
extensionMatrixAtlas := AtlasGenerators(extensionMatrixInfo.identifier);;
extensionPermutationAtlas := AtlasGenerators(
    extensionPermutationInfo.identifier);;
coreMatrixAtlas := AtlasGenerators(coreMatrixInfo.identifier);;
h54 := extensionMatrixAtlas.generators;;
p3510Generators := extensionPermutationAtlas.generators;;
g27 := coreMatrixAtlas.generators;;
AssertH(List(h54, Order) = [2,18] and
    List(p3510Generators, Order) = [2,18] and
    Order(Product(h54)) = 42 and
    Order(Product(p3510Generators)) = 42,
    "unexpected extension standard generators");
AssertH(List(g27, Order) = [2,13],
    "unexpected core standard generators");

p3510 := Group(p3510Generators);;
SetSize(p3510, extensionPermutationInfo.size);;
expectedPointStabilizerOrder := extensionPermutationInfo.size / 3510;;

# The normal kernel is invisible in the 3510 action, but it changes every
# nonzero module vector and must be included in the local value-label group.
kernelProgram := AtlasProgram(
    "3.Fi22.2", extensionMatrixAtlas.standardization,
    "kernel", "Fi22.2");;
AssertH(kernelProgram <> fail,
    "normal-kernel program is unavailable");
normalKernelGenerator := ResultOfStraightLineProgram(
    kernelProgram.program, h54)[1];;
AssertH(Order(normalKernelGenerator) = 3 and
    54 - RankMat(normalKernelGenerator - IdentityMat(54, field2)) = 0 and
    54 - RankMat(normalKernelGenerator^2 -
        IdentityMat(54, field2)) = 0,
    "normal C3 kernel is not fixed-point-free");

# Locate the quotient maximal subgroup aligned with the degree-3510 action.
alignedMax := fail;;
for maximalPosition in [1..30] do
  maximal := AtlasProgram("Fi22.2", "maxes", maximalPosition);
  if maximal = fail then
    break;
  fi;
  if maximal.size = expectedPointStabilizerOrder then
    maxPermutationGenerators := ResultOfStraightLineProgram(
        maximal.program, p3510Generators);
    fixedPoints := Filtered([1..3510], point ->
        ForAll(maxPermutationGenerators,
            generator -> point ^ generator = point));
    if Length(fixedPoints) = 1 then
      alignedMax := [maximalPosition, maximal, fixedPoints[1],
          maxPermutationGenerators];
      break;
    fi;
  fi;
od;
AssertH(alignedMax <> fail,
    "aligned degree-3510 maximal subgroup was not found");
maximalPosition := alignedMax[1];;
maximal := alignedMax[2];;
basePoint := alignedMax[3];;
maxPermutationGenerators := alignedMax[4];;
maxMatrixGenerators := ResultOfStraightLineProgram(
    maximal.program, h54);;
fullMaxMatrixGenerators := Concatenation(
    maxMatrixGenerators, [normalKernelGenerator]);;
dualGenerators := List(h54,
    generator -> TransposedMat(generator^-1));;
dualFullMaxGenerators := List(fullMaxMatrixGenerators,
    generator -> TransposedMat(generator^-1));;

# Find the full-preimage-invariant dual 12-space, directly or as the
# annihilator of a primal 42-space.
primalCompositionBases := MTX.BasesCompositionSeries(
    GModuleByMats(fullMaxMatrixGenerators, field2));;
dualCompositionBases := MTX.BasesCompositionSeries(
    GModuleByMats(dualFullMaxGenerators, field2));;
seedCandidates := [];;
for rows in dualCompositionBases do
  if Length(rows) = 12 then
    AddSet(seedCandidates, CanonicalRows(rows));
  fi;
od;
for rows in primalCompositionBases do
  if Length(rows) = 42 then
    AddSet(seedCandidates,
        CanonicalRows(NullspaceMat(TransposedMat(rows))));
  fi;
od;
seedCandidates := Filtered(seedCandidates, seed ->
    Length(seed) = 12 and
    ForAll(dualFullMaxGenerators, generator ->
      ForAll(seed, row ->
        row * generator in VectorSpace(field2, seed))));;
AssertH(Length(seedCandidates) >= 1,
    "no invariant dual 12-space was found");

ValueDataForSeed := function(seed)
  local seedSpaceLocal, seedBasisLocal, restrictedDualLocal,
        imagesLocal, valueGeneratorsLocal, valueGroupLocal,
        nonzeroValuesLocal, valueOrbitsLocal;
  seedSpaceLocal := VectorSpace(field2, seed);
  seedBasisLocal := Basis(seedSpaceLocal, seed);
  restrictedDualLocal := [];
  for generator in dualFullMaxGenerators do
    imagesLocal := List(seed, row -> row * generator);
    AssertH(ForAll(imagesLocal, row -> row in seedSpaceLocal),
        "candidate is not full-preimage invariant");
    Add(restrictedDualLocal,
        List(imagesLocal,
          row -> Coefficients(seedBasisLocal, row)));
  od;
  valueGeneratorsLocal := List(restrictedDualLocal,
      matrix -> TransposedMat(matrix^-1));
  valueGroupLocal := Group(valueGeneratorsLocal);
  nonzeroValuesLocal := Filtered(
      Tuples([Zero(field2), One(field2)], 12),
      vector -> ForAny(vector, entry -> not IsZero(entry)));
  valueOrbitsLocal := ShallowCopy(
      Orbits(valueGroupLocal, nonzeroValuesLocal, OnRight));
  Sort(valueOrbitsLocal,
      function(left, right)
        return Length(left) < Length(right);
      end);
  return [restrictedDualLocal, valueGeneratorsLocal,
      nonzeroValuesLocal, valueOrbitsLocal];
end;;

matchingSeedData := [];;
for candidatePosition in [1..Length(seedCandidates)] do
  candidateSeed := seedCandidates[candidatePosition];
  candidateData := ValueDataForSeed(candidateSeed);
  candidateOrbitSizes := List(candidateData[4], Length);
  dualMatches := Number(dualCompositionBases, rows ->
      Length(rows) = 12 and
      CanonicalRows(rows) = candidateSeed);
  primalMatches := Number(primalCompositionBases, rows ->
      Length(rows) = 42 and
      CanonicalRows(NullspaceMat(TransposedMat(rows))) =
        candidateSeed);
  Print("STRUCTURAL seed_candidate=", candidatePosition,
      " of=", Length(seedCandidates),
      " dual_composition_matches=", dualMatches,
      " primal_annihilator_matches=", primalMatches,
      " value_orbit_sizes=", candidateOrbitSizes,
      " restricted_generator_orders=",
      List(candidateData[2], Order), "\n");
  if candidateOrbitSizes = [2016,2079] then
    Add(matchingSeedData,
        [candidatePosition, candidateSeed, candidateData]);
  fi;
od;
AssertH(Length(matchingSeedData) = 1,
    "the geometric dual 12-space was not selected uniquely");
selectedCandidatePosition := matchingSeedData[1][1];;
dualSeed := matchingSeedData[1][2];;
selectedValueData := matchingSeedData[1][3];;
restrictedDualGenerators := selectedValueData[1];;
valueGenerators := selectedValueData[2];;
valueGroup := Group(valueGenerators);;
nonzeroValues := selectedValueData[3];;
valueOrbits := selectedValueData[4];;
smallValueOrbit := Set(valueOrbits[1]);;
valueOrbitSets := List(valueOrbits, Set);;

# Transport an ordered basis along the first-discovery Schreier tree.
# Revisited edges are checked only for subspace equality here; their exact
# ordered-basis transition is checked against the label table below.
dualBases := [];;
dualBases[basePoint] := dualSeed;;
queue := [basePoint];;
head := 1;;
alignmentChecks := 0;;
while head <= Length(queue) do
  point := queue[head];
  head := head + 1;
  for generatorPosition in [1..2] do
    imagePoint := point ^ p3510Generators[generatorPosition];
    imageBasis := List(dualBases[point],
        row -> row * dualGenerators[generatorPosition]);
    if not IsBound(dualBases[imagePoint]) then
      dualBases[imagePoint] := imageBasis;
      Add(queue, imagePoint);
    else
      targetSpace := VectorSpace(field2, dualBases[imagePoint]);
      AssertH(RankMat(imageBasis) = 12 and
          ForAll(imageBasis, row -> row in targetSpace),
          "matrix/permutation dual-subspace alignment conflict");
      alignmentChecks := alignmentChecks + 1;
    fi;
  od;
od;
AssertH(Length(queue) = 3510,
    "dual 12-space orbit does not have length 3510");
canonicalDualSpaces := List(dualBases, CanonicalRows);;
AssertH(Length(Set(canonicalDualSpaces)) = 3510,
    "dual 12-space orbit contains duplicate subspaces");
AssertH(RankMat(Concatenation(dualBases)) = 54,
    "dual 12-spaces do not span the dual module");

# The uniquely selected seed has the full-preimage 2016/2079 value orbits.
# A basis transition A acts on evaluation tuples through A^(-T), not A.
AssertH(Length(nonzeroValues) = 4095 and
    List(valueOrbits, Length) = [2016,2079] and
    Sum(List(valueOrbits, Length)) = 4095,
    "selected value-label table is incomplete");

ValueLabel := function(value)
  if ForAll(value, IsZero) then
    return 0;
  elif value in smallValueOrbit then
    return 1;
  fi;
  return 2;
end;;

# Check the complete two-cell partition under every lifted
# point-stabilizer generator, including the normal-kernel generator.
fullGeneratorLabelChecks := 0;;
for valueGenerator in valueGenerators do
  for orbitPosition in [1..Length(valueOrbitSets)] do
    AssertH(ForAll(valueOrbitSets[orbitPosition], value ->
        value * valueGenerator in valueOrbitSets[orbitPosition]),
        "a lifted point-stabilizer generator changes a value label");
    fullGeneratorLabelChecks := fullGeneratorLabelChecks +
        Length(valueOrbitSets[orbitPosition]);
  od;
od;
AssertH(fullGeneratorLabelChecks =
    Length(valueGenerators) * 4095,
    "full-generator label-check count mismatch");

EvaluateOnBasis := function(vector, basis)
  return List(basis, functional ->
    Sum([1..54], coordinate ->
      vector[coordinate] * functional[coordinate]));
end;;

# Exact all-edge check.  If imageBasis = A * storedBasis, evaluation tuples
# change by A^(-T).  Exact membership in the full local value group proves
# preservation of every cell of the complete 2016/2079 orbit partition.
# We additionally check the actual vector/basis pairing on all 54 standard
# basis vectors on all 7020 directed generator edges.  Computing the complete
# 54-by-12 evaluation matrices at once is both exact and substantially faster
# than 379080 separate scalar-product loops.
edgeTransitionChecks := 0;;
edgeValueGroupChecks := 0;;
edgeExhaustiveLabelChecks := 0;;
checkedValueTransitions := [];;
directValueChecks := 0;;
for point in [1..3510] do
  for generatorPosition in [1..2] do
    imagePoint := point ^ p3510Generators[generatorPosition];
    imageBasis := List(dualBases[point],
        row -> row * dualGenerators[generatorPosition]);
    targetSpace := VectorSpace(field2, dualBases[imagePoint]);
    targetBasisObject := Basis(targetSpace, dualBases[imagePoint]);
    transition := List(imageBasis,
        row -> Coefficients(targetBasisObject, row));
    AssertH(RankMat(transition) = 12,
        "singular ordered-basis edge transition");
    valueTransition := TransposedMat(transition^-1);
    AssertH(valueTransition in valueGroup,
        "an edge transition is outside the full local value group");
    checkedPosition := Position(
        checkedValueTransitions, valueTransition);
    if checkedPosition = fail then
      for orbitPosition in [1..Length(valueOrbitSets)] do
        AssertH(ForAll(valueOrbitSets[orbitPosition], value ->
            value * valueTransition in
                valueOrbitSets[orbitPosition]),
            "an edge transition changes a value label");
      od;
      Add(checkedValueTransitions, valueTransition);
    fi;
    edgeTransitionChecks := edgeTransitionChecks + 1;
    edgeValueGroupChecks := edgeValueGroupChecks + 1;
    edgeExhaustiveLabelChecks :=
        edgeExhaustiveLabelChecks + 4095;
    sourceEvaluations := TransposedMat(dualBases[point]);
    targetEvaluations := h54[generatorPosition] *
        TransposedMat(dualBases[imagePoint]);
    AssertH(targetEvaluations =
        sourceEvaluations * valueTransition,
        "direct standard-vector tuple transport failed");
    AssertH(ForAll([1..54], coordinate ->
        ValueLabel(sourceEvaluations[coordinate]) =
          ValueLabel(targetEvaluations[coordinate])),
        "direct standard-vector label check failed");
    directValueChecks := directValueChecks + 54;
  od;
od;
AssertH(edgeTransitionChecks = 2 * 3510 and
    edgeValueGroupChecks = 2 * 3510 and
    edgeExhaustiveLabelChecks = 2 * 3510 * 4095 and
    directValueChecks = 2 * 3510 * 54,
    "edge-check count mismatch");

# The 693-suborbit orbital is the undirected Fischer graph.  Its rows are
# transported and checked along the same standard-generator edges.
pointStabilizer := Stabilizer(p3510, basePoint);;
SetSize(pointStabilizer, expectedPointStabilizerOrder);;
suborbits := Orbits(pointStabilizer, [1..3510]);;
AssertH(Set(List(suborbits, Length)) = [1,693,2816],
    "unexpected degree-3510 subdegrees");
baseNeighbors := First(suborbits, orbit -> Length(orbit) = 693);;
neighbors := [];;
neighbors[basePoint] := Set(baseNeighbors);;
queue := [basePoint];;
head := 1;;
graphChecks := 0;;
while head <= Length(queue) do
  point := queue[head];
  head := head + 1;
  for generatorPosition in [1..2] do
    imagePoint := point ^ p3510Generators[generatorPosition];
    imageNeighbors := Set(List(neighbors[point],
        neighbor -> neighbor ^ p3510Generators[generatorPosition]));
    if not IsBound(neighbors[imagePoint]) then
      neighbors[imagePoint] := imageNeighbors;
      Add(queue, imagePoint);
    else
      AssertH(neighbors[imagePoint] = imageNeighbors,
          "Fischer-graph alignment conflict");
      graphChecks := graphChecks + 1;
    fi;
  od;
od;
AssertH(Length(queue) = 3510 and
    ForAll(neighbors, row -> Length(row) = 693),
    "Fischer graph is incomplete");
AssertH(ForAll([1..3510], point ->
    ForAll(neighbors[point],
      neighbor -> point in neighbors[neighbor])),
    "Fischer graph is not undirected");

# Deterministic point base for strict singleton results.
base := [];;
container := p3510;;
while Size(container) > 1 do
  point := SmallestMovedPoint(container);;
  Add(base, point);
  container := Stabilizer(container, point);;
od;
AssertH(Size(container) = 1,
    "exported quotient point base is not a base");

# Exact restriction-of-scalars bridge.  Rather than asking C++ to reproduce
# GAP's basis convention, export the 54-bit H-coordinate contribution of
# each of the three nonzero F4 symbols in every one of the 27 coordinates.
basis4over2 := CanonicalBasis(field4);;
blownGenerators := List(g27,
    generator -> ImmutableMatrix(
        field2, BlownUpMat(basis4over2, generator)));;
sourceModule := GModuleByMats(blownGenerators, field2);;
coreMax := AtlasProgram("3.Fi22.2", "maxes", 1);;
AssertH(coreMax <> fail,
    "index-two core maximal program is unavailable");
k54 := ResultOfStraightLineProgram(coreMax.program, h54);;
intertwiner := fail;;
adjustment := fail;;
for firstKernelPower in [0..2] do
  for secondKernelPower in [0..2] do
    adjustedK54 := [
      k54[1] * normalKernelGenerator^firstKernelPower,
      k54[2] * normalKernelGenerator^secondKernelPower
    ];
    if List(adjustedK54, Order) <> [2,13] then
      continue;
    fi;
    targetModule := GModuleByMats(adjustedK54, field2);
    candidateIntertwiner := MTX.IsomorphismModules(
        sourceModule, targetModule);
    if candidateIntertwiner <> fail then
      intertwiner := candidateIntertwiner;
      adjustment := [firstKernelPower, secondKernelPower];
      k54 := adjustedK54;
      break;
    fi;
  od;
  if intertwiner <> fail then break; fi;
od;
AssertH(intertwiner <> fail and RankMat(intertwiner) = 54,
    "restriction-of-scalars module bridge was not found");
if ForAll([1..2], position ->
    blownGenerators[position] * intertwiner =
      intertwiner * k54[position]) then
  bridgeMatrix := intertwiner;
  bridgeOrientation := "source_times_T";
elif ForAll([1..2], position ->
    k54[position] * intertwiner =
      intertwiner * blownGenerators[position]) then
  bridgeMatrix := intertwiner^-1;
  bridgeOrientation := "source_times_T_inverse";
else
  Error("FI22.2 FISCHER EXPORT FAILURE: bridge orientation check failed");
fi;
AssertH(ForAll([1..2], position ->
    blownGenerators[position] * bridgeMatrix =
      bridgeMatrix * k54[position]),
    "normalized bridge conjugation check failed");

decode4 := [Zero(field4), One(field4), Z(4), Z(4)^2];;
bridgeContributions := [];;
for coordinate in [1..27] do
  bridgeContributions[coordinate] := [];;
  for digit in [1..3] do
    coreVector := ListWithIdenticalEntries(27, Zero(field4));;
    coreVector[coordinate] := decode4[digit + 1];;
    bridgeContributions[coordinate][digit] :=
        BlownUpVector(basis4over2, coreVector) * bridgeMatrix;
  od;
od;

path := Concatenation(OUT, "/bases.txt");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
for point in [1..3510] do
  AppendTo(output,
      Concatenation(List(dualBases[point], Bits)), "\n");
od;
CloseStream(output);

path := Concatenation(OUT, "/functionals.txt");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
for value in nonzeroValues do
  AppendTo(output, Bits(value), " ", ValueLabel(value), "\n");
od;
CloseStream(output);

path := Concatenation(OUT, "/neighbors.txt");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
for point in [1..3510] do
  AppendTo(output,
      JoinStringsWithSeparator(
          List(neighbors[point], neighbor -> String(neighbor - 1)), " "),
      "\n");
od;
CloseStream(output);

path := Concatenation(OUT, "/baseB.txt");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
AppendTo(output, JoinStringsWithSeparator(List(base, String), " "), "\n");
CloseStream(output);

path := Concatenation(OUT, "/p3510_generators.tsv");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
AppendTo(output, "point\tgenerator_1\tgenerator_2\n");
for point in [1..3510] do
  AppendTo(output, point, "\t",
      point ^ p3510Generators[1], "\t",
      point ^ p3510Generators[2], "\n");
od;
CloseStream(output);

path := Concatenation(OUT, "/f4_to_h54.tsv");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
AppendTo(output, "coordinate\tdigit\th54_bits\n");
for coordinate in [1..27] do
  for digit in [1..3] do
    AppendTo(output, coordinate, "\t", digit, "\t",
        Bits(bridgeContributions[coordinate][digit]), "\n");
  od;
od;
CloseStream(output);

path := Concatenation(OUT, "/h54_generators.txt");;
output := OutputTextFile(path, false);;
SetPrintFormattingStatus(output, false);
for generator in h54 do
  for row in generator do
    AppendTo(output, Bits(row), "\n");
  od;
od;
CloseStream(output);

Print("Fi22 Fischer-action data\n");
Print("matrix_representation=", extensionMatrixAtlas.repname,
    " permutation_representation=", extensionPermutationAtlas.repname,
    " core_representation=", coreMatrixAtlas.repname, "\n");
Print("extension_order=", 3 * extensionPermutationInfo.size,
    " quotient_order=", extensionPermutationInfo.size,
    " kernel_order=3 kernel_fixed_dimension=0\n");
Print("maximal_position=", maximalPosition,
    " maximal_name=", maximal.subgroupname,
    " quotient_maximal_order=", maximal.size,
    " base_point=", basePoint, "\n");
Print("composition_dimensions_primal=",
    List(primalCompositionBases, Length),
    " dual=", List(dualCompositionBases, Length), "\n");
Print("dual_spaces=3510 dimension=12 distinct=3510 span=54",
    " alignment_checks=", alignmentChecks, "\n");
Print("value_orbits=0,2016,2079 full_preimage_includes_kernel=1\n");
Print("edge_transition_checks=", edgeTransitionChecks,
    " edge_value_group_membership_checks=", edgeValueGroupChecks,
    " edge_exhaustive_label_checks=", edgeExhaustiveLabelChecks,
    " distinct_edge_transitions=", Length(checkedValueTransitions),
    " direct_standard_vector_checks=", directValueChecks, "\n");
Print("fischer_graph=3510 valency=693 graph_checks=", graphChecks, "\n");
Print("point_base=", base, " length=", Length(base), "\n");
Print("restriction_bridge orientation=", bridgeOrientation,
    " kernel_adjustment=", adjustment, " rank=54\n");
Print("Wrote the Fischer graph and local actions.\n");
QUIT_GAP(0);
