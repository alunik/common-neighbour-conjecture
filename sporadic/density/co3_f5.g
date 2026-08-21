# Co3 on its 23-dimensional module over GF(5).
#
# Every nonregular nonzero vector is an eigenvector of some prime-order
# element.  We sum the sizes of all relevant eigenspaces class by class.  If
# the resulting lower bound for regular vectors is greater than 5^22 - 1,
# the projective-density lemma gives R + R = V.

if LoadPackage("atlasrep") <> true then
  Error("This program requires the AtlasRep package");
fi;
if LoadPackage("ctbllib") <> true then
  Error("This program requires the CTblLib package");
fi;

Require := function(condition, message)
  if not condition then
    Error(message);
  fi;
end;

q := 5;
d := 23;
field := GF(q);
identity := IdentityMat(d, field);

# Load the standard matrix representation and the ordinary class table.
table := CharacterTable("Co3");
info := OneAtlasGeneratingSetInfo(
  "Co3", Characteristic, q, Dimension, d
);
Require(info <> fail, "The 23-dimensional GF(5) representation is missing");

atlas := AtlasGenerators(info.identifier);
Require(atlas <> fail, "Could not load the Co3 matrices");
generators := atlas.generators;
group := Group(generators);

Require(Size(group) = 495766656000, "The matrices do not generate Co3");
Require(
  MTX.IsIrreducible(GModuleByMats(generators, field)),
  "The module is reducible"
);

# AtlasRep supplies one word in the standard generators for each conjugacy
# class.  Only the classes of prime order are needed.
classProgram := AtlasProgram("Co3", atlas.standardization, "classes");
Require(classProgram <> fail, "The class-representative program is missing");

representatives := ResultOfStraightLineProgram(
  classProgram.program, generators
);
classNames := AtlasClassNames(table);
classOrders := OrdersClassRepresentatives(table);
classSizes := SizesConjugacyClasses(table);
Require(classProgram.outputs = classNames,
  "The class names in AtlasRep and the character table do not align");
Require(
  List(representatives, Order) = classOrders,
  "The matrix representatives and class table do not align"
);

primeClasses := Filtered(
  [1 .. Length(classOrders)],
  position -> IsPrimeInt(classOrders[position])
);

Print("Co3 on GF(5)^23\n");
Print("class  size          eigenspace dimensions over GF(5)\n");

nonregularUpperBound := 0;
for position in primeClasses do
  representative := representatives[position];
  dimensions := List(
    [1 .. q - 1],
    scalar -> d - RankMat(
      representative - (scalar * One(field)) * identity
    )
  );

  # There are 5^e - 1 nonzero vectors in an e-dimensional eigenspace.
  classContribution := classSizes[position] * Sum(
    dimensions,
    dimension -> q^dimension - 1
  );
  nonregularUpperBound := nonregularUpperBound + classContribution;

  Print(
    classNames[position], "  ",
    String(classSizes[position], 12), "  ",
    JoinStringsWithSeparator(List(dimensions, String), ","), "\n"
  );
od;

nonzeroVectors := q^d - 1;
regularLowerBound := nonzeroVectors - nonregularUpperBound;
densityThreshold := q^(d - 1) - 1;

Print("\nupper bound for nonregular vectors: ",
  nonregularUpperBound, "\n");
Print("lower bound for regular vectors:    ",
  regularLowerBound, "\n");
Print("density threshold:                  ",
  densityThreshold, "\n");

Require(
  regularLowerBound > densityThreshold,
  "The density inequality was not obtained"
);
Print("Conclusion: every vector is a sum of two regular vectors.\n");

QUIT_GAP(0);
