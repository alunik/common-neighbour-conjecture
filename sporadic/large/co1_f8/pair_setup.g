# Compute the stabilizer of one representative from each Co1-vector orbit.
# The resulting matrices are consumed by pair_orbits.cpp.

LoadPackage("atlasrep");;

WORK := GetEnv("CO1_WORK_DIR");;
if WORK = fail then Error("set CO1_WORK_DIR"); fi;

AssertTrue := function(condition, message)
  if not condition then
    Error(Concatenation("PAIR-LEVEL PREPARE FAILURE: ", message));
  fi;
end;;
MaskOfRow := function(row)
  return Sum([1..Length(row)],
      i -> Int(row[i]) * 2^(i - 1));
end;;

M := AtlasGroup("Co1", Dimension, 24, Characteristic, 2);;
P := AtlasGroup("Co1", NrMovedPoints, 98280);;
mgens := GeneratorsOfGroup(M);;
pgens := GeneratorsOfGroup(P);;
AssertTrue(Size(M) = 4157776806543360000 and
    Size(P) = 4157776806543360000, "Co1 orders");
hom := GroupHomomorphismByImagesNC(P, M, pgens, mgens);;

orbit_lines := Filtered(
    SplitString(StringFile(Concatenation(WORK,
        "/vectors_orbit.tsv")), "\n"),
    line -> Length(line) > 0);;
functional_masks := List(orbit_lines{[2..98281]},
    line -> Int(SplitString(line, "\t")[2]));;
F := GF(2);;
one := One(F);;
Bits := mask -> List([0..23],
    bit -> (QuoInt(mask, 2^bit) mod 2) * one);;
evaluation_matrix := TransposedMat(
    ImmutableMatrix(F, List(functional_masks, Bits)));;
AssertTrue(RankMat(evaluation_matrix) = 24, "dual orbit rank");
SupportMask := function(mask)
  return Set(Positions(Bits(mask) * evaluation_matrix, one));
end;;

vector_lines := Filtered(
    SplitString(StringFile(Concatenation(WORK,
        "/vectors_vector_orbits.tsv")), "\n"),
    line -> Length(line) > 0);;
AssertTrue(vector_lines[1] = "orbit\trepresentative_mask\tsize",
    "vector-orbit header");
vector_rows := List(vector_lines{[2..Length(vector_lines)]},
    line -> SplitString(line, "\t"));;
AssertTrue(Length(vector_rows) = 4, "expected four vector orbits");
first_masks := List(vector_rows, row -> Int(row[2]));;
first_orbit_sizes := List(vector_rows, row -> Int(row[3]));;
AssertTrue(Sum(first_orbit_sizes) = 2^24,
    "vector-orbit sizes do not sum to 2^24");
expected_stabilizer_orders := List(first_orbit_sizes,
    size -> QuoInt(Size(P), size));;

output := OutputTextFile(
    Concatenation(WORK, "/pair_input.txt"), false);;
SetPrintFormattingStatus(output, false);;
PrintTo(output, "CO1_F8_PAIR_LEVEL_V1\n");
for index in [1..Length(first_masks)] do
  a := first_masks[index];;
  if a = 0 then
    K := P;;
  else
    K := Stabilizer(P, SupportMask(a), OnSets);;
  fi;
  AssertTrue(Size(K) = expected_stabilizer_orders[index],
      Concatenation("first stabilizer order at mask ", String(a)));
  kgens := GeneratorsOfGroup(K);;
  kmats := List(kgens, generator -> Image(hom, generator));;
  AssertTrue(ForAll(kmats, matrix -> Bits(a) * matrix = Bits(a)),
      Concatenation("matrix image does not fix mask ", String(a)));

  AppendTo(output, "first ", a,
      " orbit_size ", first_orbit_sizes[index],
      " stabilizer_order ", Size(K),
      " generator_count ", Length(kmats), "\n");
  for matrix in kmats do
    AppendTo(output, "matrix");
    for row in matrix do
      AppendTo(output, " ", MaskOfRow(row));
    od;
    AppendTo(output, "\n");
  od;
  Print("First-vector representative ", a,
      ": orbit size ", first_orbit_sizes[index],
      ", stabilizer order ", Size(K), "\n");
od;
CloseStream(output);;
QUIT;
