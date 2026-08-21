# Check that the vector used as the first summand is regular for
# 3.Fi22.2 on its 54-dimensional module over GF(2).

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

AssertExt := function(condition, message)
  if not condition then
    Error(Concatenation("FI22 EXTENSION PROBE FAILURE: ", message));
  fi;
end;;

field2 := GF(2);;
field4 := GF(4);;
root := "sporadic/large/fi22";;

extensionMatrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22.2", Characteristic, 2, Dimension, 54);;
extensionPermutation61776Info := OneAtlasGeneratingSetInfo(
    "Fi22.2", NrMovedPoints, 61776);;
coreMatrixInfo := OneAtlasGeneratingSetInfo(
    "3.Fi22", Characteristic, 2, Dimension, 27);;
AssertExt(ForAll([
    extensionMatrixInfo, extensionPermutation61776Info,
    coreMatrixInfo], info -> info <> fail),
    "required ATLAS representations are unavailable");

extensionMatrixAtlas := AtlasGenerators(
    extensionMatrixInfo.identifier);;
extensionPermutation61776Atlas := AtlasGenerators(
    extensionPermutation61776Info.identifier);;
coreMatrixAtlas := AtlasGenerators(coreMatrixInfo.identifier);;
h54 := extensionMatrixAtlas.generators;;
p61776gens := extensionPermutation61776Atlas.generators;;
g27 := coreMatrixAtlas.generators;;

Print("REP matrix=", extensionMatrixAtlas.repname,
    " p61776=", extensionPermutation61776Atlas.repname, "\n");
Print("ORDERS matrix=", List(h54, Order),
    " p61776=", List(p61776gens, Order), "\n");
Print("PRODUCT_ORDERS matrix=", Order(Product(h54)),
    " p61776=", Order(Product(p61776gens)), "\n");

# The first maximal subgroup is the index-two core 3.Fi22.
max1 := AtlasProgram("3.Fi22.2", "maxes", 1);;
AssertExt(max1 <> fail, "max1 program unavailable");
k54 := ResultOfStraightLineProgram(max1.program, h54);;
Print("MAX1 name=", max1.subgroupname,
    " size=", max1.size,
    " generator_orders=", List(k54, Order), "\n");

kernelProgram := AtlasProgram(
    "3.Fi22.2", extensionMatrixAtlas.standardization,
    "kernel", "Fi22.2");;
AssertExt(kernelProgram <> fail, "normal-kernel program unavailable");
centralGenerator := ResultOfStraightLineProgram(
    kernelProgram.program, h54)[1];;
AssertExt(Order(centralGenerator) = 3 and
    54 - RankMat(centralGenerator - IdentityMat(54, field2)) = 0 and
    54 - RankMat(centralGenerator^2 - IdentityMat(54, field2)) = 0,
    "normal C3 kernel is not fixed-point-free");
Print("NORMAL_KERNEL order=3 fixed_dimension=",
    54 - RankMat(centralGenerator - IdentityMat(54, field2)), "\n");

# Blow up the F4 module to F2 and construct an exact module isomorphism.
basis4over2 := CanonicalBasis(field4);;
blownGenerators := List(g27,
    generator -> ImmutableMatrix(
        field2, BlownUpMat(basis4over2, generator)));;
sourceModule := GModuleByMats(blownGenerators, field2);;
intertwiner := fail;;
adjustment := fail;;
for firstCentralPower in [0..2] do
  for secondCentralPower in [0..2] do
    adjustedK54 := [
      k54[1] * centralGenerator^firstCentralPower,
      k54[2] * centralGenerator^secondCentralPower
    ];
    if List(adjustedK54, Order) <> [2,13] then
      continue;
    fi;
    targetModule := GModuleByMats(adjustedK54, field2);
    candidateIntertwiner := MTX.IsomorphismModules(
        sourceModule, targetModule);
    if candidateIntertwiner <> fail then
      intertwiner := candidateIntertwiner;
      adjustment := [firstCentralPower, secondCentralPower];
      k54 := adjustedK54;
      break;
    fi;
  od;
  if intertwiner <> fail then break; fi;
od;
AssertExt(intertwiner <> fail and RankMat(intertwiner) = 54,
    "restriction-of-scalars modules are not isomorphic");
orientation := fail;;
if ForAll([1..2], position ->
    blownGenerators[position] * intertwiner =
      intertwiner * k54[position]) then
  orientation := "source_times_T";
elif ForAll([1..2], position ->
    k54[position] * intertwiner =
      intertwiner * blownGenerators[position]) then
  orientation := "target_times_T";
fi;
AssertExt(orientation <> fail, "intertwiner conjugation check failed");
Print("RESTRICTION_BRIDGE rank=54 orientation=", orientation,
    " central_adjustment=", adjustment,
    " adjusted_orders=", List(k54, Order), "\n");

# Find the maximal subgroup aligned with the degree-61776 action.
alignedMax := fail;;
for maximalPosition in [1..20] do
  maximal := AtlasProgram("Fi22.2", "maxes", maximalPosition);
  if maximal = fail then break; fi;
  if maximal.size = extensionPermutation61776Info.size / 61776 then
    alignedMax := [maximalPosition, maximal];
    break;
  fi;
od;
AssertExt(alignedMax <> fail, "degree-61776 maximal not identified");
maximalPosition := alignedMax[1];;
maximal := alignedMax[2];;
maxP61776 := ResultOfStraightLineProgram(
    maximal.program, p61776gens);;
fixedPoints := Filtered([1..61776], point ->
    ForAll(maxP61776, generator -> point ^ generator = point));;
Print("ALIGNED_MAX position=", maximalPosition,
    " name=", maximal.subgroupname,
    " size=", maximal.size,
    " fixed_points=", fixedPoints, "\n");

# The corresponding 54-dimensional maximal subgroup must preserve a
# codimension-two subspace; its annihilator is the dual 2-space attached
# to point 1.
maxH54 := ResultOfStraightLineProgram(maximal.program, h54);;
compositionBases := MTX.BasesCompositionSeries(
    GModuleByMats(maxH54, field2));;
dualCompositionBases := MTX.BasesCompositionSeries(
    GModuleByMats(List(maxH54,
        generator -> TransposedMat(generator^-1)), field2));;
Print("MAX_COMPOSITION primal=", List(compositionBases, Length),
    " dual=", List(dualCompositionBases, Length), "\n");

# Exact H-regularity of the reference core seed, after transport to the ATLAS
# 54-dimensional coordinates, using the aligned dual-plane zero set.
codimensionTwoPosition := PositionProperty(
    compositionBases, basis -> Length(basis) = 52);;
AssertExt(codimensionTwoPosition <> fail,
    "maximal subgroup has no invariant codimension-two subspace");
dualPlaneSeed := NullspaceMat(
    TransposedMat(compositionBases[codimensionTwoPosition]));;
AssertExt(Length(dualPlaneSeed) = 2,
    "annihilator is not two-dimensional");
dualH54 := List(h54,
    generator -> TransposedMat(generator^-1));;
dualPlanes := [];;
dualPlaneBasePoint := fixedPoints[1];;
dualPlanes[dualPlaneBasePoint] := CanonicalBasis(
    VectorSpace(field2, dualPlaneSeed));;
dualPlanes[dualPlaneBasePoint] :=
    BasisVectors(dualPlanes[dualPlaneBasePoint]);;
queue := [dualPlaneBasePoint];;
head := 1;;
while head <= Length(queue) do
  point := queue[head];
  head := head + 1;
  for generatorPosition in [1..2] do
    imagePoint := point ^ p61776gens[generatorPosition];
    imagePlane := CanonicalBasis(VectorSpace(field2,
        List(dualPlanes[point],
            vector -> vector * dualH54[generatorPosition])));
    imagePlane := BasisVectors(imagePlane);
    if not IsBound(dualPlanes[imagePoint]) then
      dualPlanes[imagePoint] := imagePlane;
      Add(queue, imagePoint);
    else
      AssertExt(dualPlanes[imagePoint] = imagePlane,
          "dual-plane alignment conflict");
    fi;
  od;
od;
AssertExt(Length(queue) = 61776,
    "dual-plane orbit does not have length 61776");
Print("DUAL_PLANES count=61776 alignment_edges=123552\n");

seedLines := Filtered(SplitString(StringFile(Concatenation(
    root, "/generated/regular_seed.txt")), "\n", "\r"),
    line -> Length(line) > 0);;
AssertExt(Length(seedLines) = 1 and Length(seedLines[1]) = 27,
    "regular seed file is malformed");
referenceF4 := List(seedLines[1], character ->
    [Zero(field4), One(field4), Z(4), Z(4)^2][
      IntChar(character) - 47]);;
referenceBlown := BlownUpVector(basis4over2, referenceF4);;
if orientation = "source_times_T" then
  reference54 := referenceBlown * intertwiner;
else
  reference54 := referenceBlown * intertwiner^-1;
fi;
referenceZeros := Set(Filtered([1..61776], point ->
    ForAll(dualPlanes[point], functional ->
      Sum([1..54], coordinate ->
        reference54[coordinate] * functional[coordinate]) =
          Zero(field2))));;
planeRows := Concatenation(dualPlanes{referenceZeros});;
Print("REFERENCE_EXT zero_count=", Length(referenceZeros),
    " zero_rank=", RankMat(planeRows), "\n");
p61776 := Group(p61776gens);;
SetSize(p61776, extensionPermutation61776Info.size);
referenceSetStabilizer := Stabilizer(p61776, referenceZeros, OnSets);;
AssertExt(Size(referenceSetStabilizer) = 1,
    "reference seed does not have trivial quotient zero-set stabilizer");
Print("REFERENCE_EXT set_stabilizer=", Size(referenceSetStabilizer), "\n");
Print("The chosen first summand is regular for 3.Fi22.2.\n");
QUIT_GAP(0);
