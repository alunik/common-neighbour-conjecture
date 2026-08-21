# Compute the exact Co1 stabilizer of one pair-orbit representative.
# The aligned permutation model makes set stabilizers quick, and the resulting
# elements are mapped back to 24 by 24 matrices for the C++ orbit enumerator.

LoadPackage("atlasrep");;

WORK := GetEnv("CO1_WORK_DIR");;
if WORK = fail then Error("set CO1_WORK_DIR"); fi;
PAIR_STABILIZERS := Concatenation(WORK, "/pair_stabilizers");;

AssertTrue := function(condition, message)
  if not condition then
    Error(Concatenation("PAIR-STABILIZER EXPORT FAILURE: ", message));
  fi;
end;;

pair_index_text := GetEnv("CO1_PAIR_INDEX");;
AssertTrue(pair_index_text <> fail, "set CO1_PAIR_INDEX");
pair_index := Int(pair_index_text);;

MaskOfRow := function(row)
  return Sum([1..Length(row)], i -> Int(row[i]) * 2^(i - 1));
end;;

pair_lines := Filtered(
    SplitString(StringFile(Concatenation(
        WORK, "/pair_orbits.tsv")), "\n"),
    line -> Length(line) > 0);;
AssertTrue(Length(pair_lines) = 47, "expected 46 pair-orbit rows");
AssertTrue(pair_lines[1] =
    "pair_index\ta\tb\tfirst_orbit_size\tb_orbit_size\tpair_stabilizer_order",
    "pair-orbit header");
AssertTrue(pair_index >= 1 and pair_index < Length(pair_lines),
    "pair index out of range");
fields := SplitString(pair_lines[pair_index + 1], "\t");;
AssertTrue(Length(fields) = 6 and Int(fields[1]) = pair_index,
    "malformed or misaligned pair row");
a := Int(fields[2]);;
b := Int(fields[3]);;
first_orbit_size := Int(fields[4]);;
b_orbit_size := Int(fields[5]);;
expected_pair_stabilizer_order := Int(fields[6]);;

P := AtlasGroup("Co1", NrMovedPoints, 98280);;
M := AtlasGroup("Co1", Dimension, 24, Characteristic, 2);;
AssertTrue(P <> fail and M <> fail, "ATLAS representations unavailable");
AssertTrue(Size(P) = 4157776806543360000 and
    Size(M) = 4157776806543360000, "Co1 order");
pgens := GeneratorsOfGroup(P);;
mgens := GeneratorsOfGroup(M);;
AssertTrue(List(pgens, Order) = [2,3] and List(mgens, Order) = [2,3],
    "standard-generator orders");
hom := GroupHomomorphismByImagesNC(P, M, pgens, mgens);;

orbit_lines := Filtered(
    SplitString(StringFile(Concatenation(
        WORK, "/vectors_orbit.tsv")), "\n"),
    line -> Length(line) > 0);;
AssertTrue(Length(orbit_lines) = 98281, "dual-orbit row count");
AssertTrue(orbit_lines[1] = "point\tfunctional_mask",
    "dual-orbit header");
functional_masks := List(orbit_lines{[2..98281]},
    line -> Int(SplitString(line, "\t")[2]));;
AssertTrue(Length(Set(functional_masks)) = 98280,
    "dual-orbit masks are not distinct");

F := GF(2);;
one := One(F);;
Bits := mask -> List([0..23],
    bit -> (QuoInt(mask, 2^bit) mod 2) * one);;
evaluation_matrix := TransposedMat(
    ImmutableMatrix(F, List(functional_masks, Bits)));;
AssertTrue(RankMat(evaluation_matrix) = 24, "dual-orbit rank");
SupportMask := function(mask)
  return Set(Positions(Bits(mask) * evaluation_matrix, one));
end;;

support_a := SupportMask(a);;
support_b := SupportMask(b);;
K := Stabilizer(P, support_a, OnSets);;
AssertTrue(Size(K) * first_orbit_size = Size(P),
    "first stabilizer/order mismatch");
L := Stabilizer(K, support_b, OnSets);;
AssertTrue(Size(L) = expected_pair_stabilizer_order,
    "pair stabilizer order");
AssertTrue(first_orbit_size * b_orbit_size * Size(L) = Size(P),
    "pair orbit-stabilizer identity");

lgens := SmallGeneratingSet(L);;
if lgens = fail then lgens := GeneratorsOfGroup(L); fi;
lgens := Filtered(lgens, generator -> not IsOne(generator));;
if Length(lgens) = 0 then
  AssertTrue(Size(L) = 1, "empty generator list for nontrivial group");
else
  AssertTrue(Size(Group(lgens)) = Size(L),
      "small generators do not generate pair stabilizer");
fi;
lmats := List(lgens, generator -> Image(hom, generator));;
AssertTrue(ForAll(lmats,
    matrix -> Bits(a) * matrix = Bits(a) and
              Bits(b) * matrix = Bits(b)),
    "exported matrix does not fix the pair");

Exec(Concatenation("mkdir -p ", PAIR_STABILIZERS));;
output_path := Concatenation(PAIR_STABILIZERS,
    "/pair-", String(pair_index), ".txt");;
output := OutputTextFile(output_path, false);;
SetPrintFormattingStatus(output, false);;
PrintTo(output, "CO1_F8_PAIR_STABILIZER_V1\n");
AppendTo(output, "pair_index ", pair_index, "\n");
AppendTo(output, "a ", a, "\n");
AppendTo(output, "b ", b, "\n");
AppendTo(output, "first_orbit_size ", first_orbit_size, "\n");
AppendTo(output, "b_orbit_size ", b_orbit_size, "\n");
AppendTo(output, "pair_stabilizer_order ", Size(L), "\n");
AppendTo(output, "generator_count ", Length(lmats), "\n");
for matrix in lmats do
  AppendTo(output, "matrix");
  for row in matrix do
    AppendTo(output, " ", MaskOfRow(row));
  od;
  AppendTo(output, "\n");
od;
CloseStream(output);;

Print("PAIR_STABILIZER_EXPORT",
    " pair=", pair_index,
    " order=", Size(L),
    " generators=", Length(lmats), "\n");
QUIT;
