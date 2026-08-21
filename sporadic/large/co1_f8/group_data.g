# Export aligned matrix and permutation generators for Co1.
#
# The C++ programs use 24-bit integers for vectors in the natural
# F_2 Co1-module.  GAP's ATLASRep package supplies both representations;
# the fixed Co2 point aligns the two without any choices in later stages.

LoadPackage("atlasrep");;

WORK := GetEnv("CO1_WORK_DIR");;
if WORK = fail then Error("set CO1_WORK_DIR"); fi;

AssertTrue := function(condition, message)
  if not condition then Error(Concatenation("PREPARE FAILURE: ", message)); fi;
end;;

MaskOfRow := function(row)
  return Sum([1..Length(row)],
      i -> Int(row[i]) * 2^(i - 1));
end;;

matrix_group := AtlasGroup("Co1", Dimension, 24, Characteristic, 2);;
permutation_group := AtlasGroup("Co1", NrMovedPoints, 98280);;
matrix_generators := GeneratorsOfGroup(matrix_group);;
permutation_generators := GeneratorsOfGroup(permutation_group);;
AssertTrue(List(matrix_generators, Order) = [2,3],
    "matrix standard-generator orders");
AssertTrue(List(permutation_generators, Order) = [2,3],
    "permutation standard-generator orders");

max1 := AtlasProgram("Co1", "maxes", 1);;
AssertTrue(max1 <> fail, "Atlas max1 program unavailable");
AssertTrue(max1.subgroupname = "Co2", "Atlas max1 is not Co2");
AssertTrue(max1.size = 42305421312000, "Co2 order mismatch");
max1_matrix_generators := ResultOfStraightLineProgram(
    max1.program, matrix_generators);;
max1_permutation_generators := ResultOfStraightLineProgram(
    max1.program, permutation_generators);;
fixed_points := Filtered([1..98280],
    point -> ForAll(max1_permutation_generators,
        generator -> point ^ generator = point));;
AssertTrue(fixed_points = [86703], "Co2 fixed-point alignment mismatch");

dual_max1_generators := List(max1_matrix_generators,
    generator -> TransposedMat(generator^-1));;
dual_module := GModuleByMats(dual_max1_generators, GF(2));;
composition_bases := MTX.BasesCompositionSeries(dual_module);;
AssertTrue(List(composition_bases, Length) = [0,1,23,24],
    "dual Co2 composition dimensions");
dual_seed := composition_bases[2][1];;
AssertTrue(ForAll(dual_max1_generators,
    generator -> dual_seed * generator = dual_seed),
    "dual seed is not Co2-fixed");

dual_generators := List(matrix_generators,
    generator -> TransposedMat(generator^-1));;
input_path := Concatenation(WORK, "/group_data.txt");;
output := OutputTextFile(input_path, false);;
SetPrintFormattingStatus(output, false);
PrintTo(output, "CO1_F8_ATLAS_INPUT_V1\n");
AppendTo(output, "seed_point ", fixed_points[1], "\n");
AppendTo(output, "seed_mask ", MaskOfRow(dual_seed), "\n");
for generator in dual_generators do
  AppendTo(output, "dual_matrix");
  for row in generator do
    AppendTo(output, " ", MaskOfRow(row));
  od;
  AppendTo(output, "\n");
od;
for generator in matrix_generators do
  AppendTo(output, "primal_matrix");
  for row in generator do
    AppendTo(output, " ", MaskOfRow(row));
  od;
  AppendTo(output, "\n");
od;
for generator in permutation_generators do
  AppendTo(output, "permutation");
  for point in [1..98280] do
    AppendTo(output, " ", point ^ generator);
  od;
  AppendTo(output, "\n");
od;
CloseStream(output);

Print("Co1 generators: 24-dimensional matrices and 98280-point permutations\n");
Print("Aligned through the Co2 fixed point 86703\n");
QUIT;
