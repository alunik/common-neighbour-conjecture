# 3.Suz on a 12-dimensional module over GF(16).
#
# The AtlasRep module is defined over GF(4) and then extended to GF(16).  The
# calculation below exhibits three inequivalent projective regular orbits.
# Their combined size is greater than 16^11 - 1, so the projective-density
# lemma gives R + R = V.

if LoadPackage("atlasrep") <> true then
  Error("This program requires the AtlasRep package");
fi;

Require := function(condition, message)
  if not condition then
    Error(message);
  fi;
end;

RunSuzCalculation := function()
  local matrixInfo, matrixAtlas, matrixGenerators,
        permutationInfo, permutationAtlas, permutationGenerators,
        field4, field16, eta, nonzeroScalars, matrixModule,
        dualModule, frobeniusGenerators, frobeniusModule,
        moduleIsomorphism, kernelProgram, centralGenerators,
        centralGenerator, maxProgram, maxMatrixGenerators,
        maxPermutationGenerators, maxDualModule, series, canonical,
        dualSeed, dualGenerators, dualPoints, queue, head, point,
        generatorPosition, imagePoint, imageVector, permutationGroup,
        epimorphism, freeGenerators, decode, candidateCodes, candidates,
        zeroSets, candidatePosition, candidate, firstVector, secondVector,
        zeroSet, setStabilizer, vector16, projectiveStabilizerOrder,
        stabilizerElement, word, matrixLift, matrixLift16, scalar,
        otherPosition, suzOrder, maximalGroupOrder, regularLowerBound,
        densityThreshold;

  # Standard matrix representation of 3.Suz and the compatible action of
  # Suz on 32760 points (the cosets of U5(2)).
  matrixInfo := OneAtlasGeneratingSetInfo(
    "3.Suz", Characteristic, 2, Dimension, 12
  );
  Require(matrixInfo <> fail,
    "The 12-dimensional representation of 3.Suz is missing");
  matrixAtlas := AtlasGenerators(matrixInfo.identifier);
  Require(matrixAtlas <> fail, "Could not load the 3.Suz matrices");
  matrixGenerators := matrixAtlas.generators;

  permutationInfo := OneAtlasGeneratingSetInfo(
    "Suz", NrMovedPoints, 32760
  );
  Require(permutationInfo <> fail,
    "The 32760-point action of Suz is missing");
  permutationAtlas := AtlasGenerators(permutationInfo.identifier);
  Require(permutationAtlas <> fail,
    "Could not load the Suz permutations");
  permutationGenerators := permutationAtlas.generators;

  Require(List(matrixGenerators, Order) = [2, 3],
    "Unexpected matrix generators");
  Require(List(permutationGenerators, Order) = [2, 3],
    "Unexpected permutation generators");
  Require(Order(matrixGenerators[1] * matrixGenerators[2]) = 13,
    "Unexpected product of the matrix generators");
  Require(Order(permutationGenerators[1] * permutationGenerators[2]) = 13,
    "Unexpected product of the permutation generators");

  field4 := GF(4);
  field16 := GF(16);
  eta := Z(16);
  nonzeroScalars := Filtered(Elements(field16), x -> not IsZero(x));

  matrixModule := GModuleByMats(matrixGenerators, field4);
  Require(MTX.IsIrreducible(matrixModule), "The GF(4) module is reducible");

  # The dual is the Frobenius conjugate.  Thus the same calculation covers
  # both choices of the 12-dimensional module.
  dualModule := MTX.DualModule(matrixModule);
  frobeniusGenerators := List(matrixGenerators, generator ->
    ImmutableMatrix(field4,
      List(generator, row -> List(row, entry -> entry^2)))
  );
  frobeniusModule := GModuleByMats(frobeniusGenerators, field4);
  moduleIsomorphism := MTX.IsomorphismModules(dualModule, frobeniusModule);
  Require(moduleIsomorphism <> fail and RankMat(moduleIsomorphism) = 12,
    "The dual and Frobenius-conjugate modules are not isomorphic");

  # The kernel of 3.Suz -> Suz is scalar, so projective stabilisers can be
  # computed in the smaller permutation action.
  kernelProgram := AtlasProgram(
    "3.Suz", matrixAtlas.standardization, "kernel", "Suz"
  );
  Require(kernelProgram <> fail, "The quotient program is missing");
  centralGenerators := ResultOfStraightLineProgram(
    kernelProgram.program, matrixGenerators
  );
  Require(Length(centralGenerators) = 1,
    "Unexpected kernel generators");
  centralGenerator := centralGenerators[1];
  Require(Order(centralGenerator) = 3 and IsScalar(centralGenerator),
    "The central subgroup is not represented by scalars of order three");

  # The third maximal subgroup is 3 x U5(2).  Its fixed line on the dual
  # module has a projective orbit of length 32760.
  maxProgram := AtlasProgram(
    "3.Suz", matrixAtlas.standardization, "maxes", 3
  );
  Require(maxProgram <> fail, "The U5(2) subgroup program is missing");
  maxMatrixGenerators := ResultOfStraightLineProgram(
    maxProgram.program, matrixGenerators
  );
  maxPermutationGenerators := ResultOfStraightLineProgram(
    maxProgram.program, permutationGenerators
  );
  Require(Filtered([1 .. 32760], point ->
    ForAll(maxPermutationGenerators, generator -> point ^ generator = point)
  ) = [1], "The subgroup does not fix the expected point");

  maxDualModule := GModuleByMats(
    List(maxMatrixGenerators,
      generator -> TransposedMat(generator^-1)),
    field4
  );
  series := MTX.BasesCompositionSeries(maxDualModule);
  Require(List(series, Length) = [0, 1, 11, 12],
    "Unexpected composition series for the restricted dual module");

  canonical := function(vector)
    local firstNonzero;
    firstNonzero := PositionProperty(vector, entry -> not IsZero(entry));
    Require(firstNonzero <> fail, "Cannot projectivise the zero vector");
    return vector / vector[firstNonzero];
  end;

  dualSeed := canonical(series[2][1]);
  dualGenerators := List(matrixGenerators,
    generator -> TransposedMat(generator^-1));

  # Align the matrix orbit with the permutation action.  This makes the
  # later set-stabiliser calculations much faster than using matrices alone.
  dualPoints := [];
  dualPoints[1] := dualSeed;
  queue := [1];
  head := 1;
  while head <= Length(queue) do
    point := queue[head];
    head := head + 1;
    for generatorPosition in [1 .. Length(permutationGenerators)] do
      imagePoint := point ^ permutationGenerators[generatorPosition];
      imageVector := canonical(
        dualPoints[point] * dualGenerators[generatorPosition]
      );
      if not IsBound(dualPoints[imagePoint]) then
        dualPoints[imagePoint] := imageVector;
        Add(queue, imagePoint);
      else
        Require(dualPoints[imagePoint] = imageVector,
          "The matrix and permutation actions do not align");
      fi;
    od;
  od;
  Require(Length(queue) = 32760,
    "The dual projective orbit has the wrong length");

  permutationGroup := Group(permutationGenerators);
  SetSize(permutationGroup, permutationInfo.size);
  epimorphism := EpimorphismFromFreeGroup(permutationGroup);
  SetIsSurjective(epimorphism, true);
  freeGenerators := GeneratorsOfGroup(Source(epimorphism));

  # Coordinates in GF(4) are written as 0, 1, z, z^2.  Each pair (a,b)
  # represents the GF(16) vector a + eta*b.
  decode := function(code)
    if code = 0 then
      return Zero(field4);
    elif code = 1 then
      return One(field4);
    elif code = 2 then
      return Z(4);
    elif code = 3 then
      return Z(4)^2;
    fi;
    Error("Invalid GF(4) coordinate");
  end;

  candidateCodes := [
    [
      [3,3,0,0,0,2,1,2,2,2,0,3],
      [0,2,1,3,0,0,3,1,0,2,1,1]
    ],
    [
      [0,3,1,1,1,1,1,1,0,2,3,1],
      [0,2,0,1,2,2,1,2,1,3,2,0]
    ],
    [
      [1,3,1,0,1,1,2,0,2,0,1,1],
      [2,2,2,2,1,1,2,3,1,2,0,0]
    ]
  ];
  candidates := List(candidateCodes, pair ->
    List(pair, vector -> List(vector, decode))
  );
  Require(ForAll(candidates, pair -> RankMat(pair) = 2),
    "A candidate pair is linearly dependent over GF(4)");

  zeroSets := [];
  Print("3.Suz on GF(16)^12\n");
  Print("candidate  zero-set size  set stabiliser  projective stabiliser\n");

  for candidatePosition in [1 .. Length(candidates)] do
    candidate := candidates[candidatePosition];
    firstVector := candidate[1];
    secondVector := candidate[2];

    # A projective stabiliser of a + eta*b must preserve the dual points
    # which vanish on both a and b.
    zeroSet := Set(Filtered([1 .. 32760], point ->
      Sum([1 .. 12], position ->
        firstVector[position] * dualPoints[point][position]) = Zero(field4)
      and
      Sum([1 .. 12], position ->
        secondVector[position] * dualPoints[point][position]) = Zero(field4)
    ));
    Add(zeroSets, zeroSet);

    setStabilizer := Stabilizer(permutationGroup, zeroSet, OnSets);
    vector16 := firstVector + eta * secondVector;

    # Lift each element of the small set stabiliser back to the matrix
    # representation and test whether it fixes the projective line.
    projectiveStabilizerOrder := 0;
    for stabilizerElement in Elements(setStabilizer) do
      word := PreImagesRepresentative(epimorphism, stabilizerElement);
      matrixLift := MappedWord(
        word, freeGenerators, matrixGenerators
      );
      matrixLift16 := ImmutableMatrix(field16, matrixLift);
      if ForAny(nonzeroScalars, scalar ->
        vector16 * matrixLift16 = scalar * vector16
      ) then
        projectiveStabilizerOrder := projectiveStabilizerOrder + 1;
      fi;
    od;

    Require(projectiveStabilizerOrder = 1,
      "A candidate line has a nontrivial projective stabiliser");
    Print(candidatePosition, String(Length(zeroSet), 17),
      String(Size(setStabilizer), 16),
      String(projectiveStabilizerOrder, 23), "\n");
  od;

  # Distinct zero-set orbits imply distinct projective vector orbits.
  for candidatePosition in [1 .. Length(zeroSets)] do
    for otherPosition in [candidatePosition + 1 .. Length(zeroSets)] do
      Require(
        RepresentativeAction(
          permutationGroup,
          zeroSets[candidatePosition],
          zeroSets[otherPosition],
          OnSets
        ) = fail,
        "Two candidate lines lie in the same orbit"
      );
    od;
  od;

  suzOrder := permutationInfo.size;
  maximalGroupOrder := 15 * suzOrder;
  regularLowerBound := Length(candidates) * maximalGroupOrder;
  densityThreshold := 16^11 - 1;

  Print("\ndistinct regular projective orbits: ",
    Length(candidates), "\n");
  Print("lower bound for regular vectors:    ",
    regularLowerBound, "\n");
  Print("density threshold:                  ",
    densityThreshold, "\n");

  Require(regularLowerBound > densityThreshold,
    "The density inequality was not obtained");
  Print("Conclusion: every vector is a sum of two regular vectors.\n");
end;

RunSuzCalculation();
QUIT_GAP(0);
