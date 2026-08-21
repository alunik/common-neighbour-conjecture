LoadPackage("atlasrep");;

SOURCE_DIR := GetEnv("SUZ6_F19_SOURCE");;
WORK_DIR := GetEnv("SUZ6_F19_WORK");;
if SOURCE_DIR = fail or WORK_DIR = fail then
  Error("run this program through run.sh");
fi;

Read(Concatenation(SOURCE_DIR, "/matrices.g"));;
Read(Concatenation(SOURCE_DIR, "/vectors.g"));;

Check := function(condition, message)
  if not condition then
    Error(Concatenation("CHECK FAILED: ", message));
  fi;
end;;

NormalizeRow := function(vector)
  local first, scale;
  first := PositionProperty(vector, x -> x mod 19 <> 0);
  Check(first <> fail, "attempt to normalize the zero vector");
  scale := PowerModInt(vector[first] mod 19, 17, 19);
  return List(vector, x -> (x * scale) mod 19);
end;;

ImageRow := function(vector, matrix)
  return List([1..12],
      j -> Sum([1..12], i -> vector[i] * matrix[i][j]) mod 19);
end;;

ReadOrbit := function(path)
  local lines;
  lines := SplitString(StringFile(path), "\n");
  Check(lines[1] = "index\tfunctional", "unexpected orbit-file header");
  lines := Filtered(lines{[2..Length(lines)]}, line -> Length(line) > 0);
  return List(lines,
      line -> List(SplitString(SplitString(line, "\t")[2], ","), Int));
end;;

IsProjectivelyFixed := function(vector, matrices)
  return ForAll(matrices, matrix -> RankMat([vector, vector * matrix]) = 1);
end;;

CheckModule := function(label, functional_matrices, invariant_seed,
                        representatives)
  local orbit, integer_matrices, generators, projective_group,
        point, k, image, expected, vector, zero_set, stabilizer,
        zero_sizes;

  orbit := ReadOrbit(Concatenation(WORK_DIR, "/", label, "_orbit.tsv"));
  Check(Length(orbit) = 32760,
      Concatenation("module ", label, ": wrong dual-orbit length"));
  Check(Length(Set(orbit)) = 32760,
      Concatenation("module ", label, ": repeated dual-orbit point"));
  Check(ForAll(orbit, vector -> NormalizeRow(vector) = vector),
      Concatenation("module ", label, ": non-normalized orbit point"));
  Check(orbit[1] = invariant_seed,
      Concatenation("module ", label, ": wrong invariant line"));

  SUZ6_F19_PERM_A := ();;
  SUZ6_F19_PERM_B := ();;
  Read(Concatenation(WORK_DIR, "/", label, "_perms.g"));
  generators := [SUZ6_F19_PERM_A, SUZ6_F19_PERM_B];
  projective_group := Group(generators);
  Check(IsTransitive(projective_group, [1..32760]),
      Concatenation("module ", label, ": projective action is not transitive"));
  Check(Size(projective_group) = 448345497600,
      Concatenation("module ", label, ": wrong projective-group order"));

  # Recompute every edge of the two generated permutations from the matrices.
  integer_matrices := List(functional_matrices,
      matrix -> List(matrix, row -> List(row, Int)));
  for point in [1..32760] do
    for k in [1..2] do
      image := NormalizeRow(ImageRow(orbit[point], integer_matrices[k]));
      expected := orbit[point ^ generators[k]];
      Check(image = expected,
          Concatenation("module ", label, ": incorrect action at point ",
              String(point)));
    od;
  od;

  Print("\nModule ", label, ": projective group of order ",
      Size(projective_group), " on 32760 points.\n");
  zero_sizes := [];
  for point in [1..Length(representatives)] do
    vector := representatives[point];
    Check(NormalizeRow(vector) = vector,
        Concatenation("module ", label, ": non-normalized representative"));
    zero_set := Set(Filtered([1..32760],
        i -> Sum([1..12], j -> vector[j] * orbit[i][j]) mod 19 = 0));
    stabilizer := Stabilizer(projective_group, zero_set, OnSets);
    Check(Size(stabilizer) = 1,
        Concatenation("module ", label,
            ": representative has a nontrivial projective stabilizer"));
    Add(zero_sizes, Length(zero_set));
    Print("  representative ", point,
        ": trivial projective stabilizer; zero-set size ",
        Length(zero_set), "\n");
  od;

  Check(Length(representatives) = 15,
      Concatenation("module ", label, ": expected 15 representatives"));
  Check(Length(Set(zero_sizes)) = 15,
      Concatenation("module ", label,
          ": zero-set sizes do not distinguish the 15 orbits"));
  Print("Module ", label,
      ": 15 distinct regular projective orbits.\n");
end;;

Check(List(SUZ6_F19_GENERATORS, Order) = [4,3],
    "standard-generator orders are not 4 and 3");

maximal_subgroup_program := AtlasProgram("6.Suz", "maxes", 3);;
Check(maximal_subgroup_program <> fail,
    "AtlasRep does not provide the required maximal-subgroup program");
Check(maximal_subgroup_program.subgroupname = "6xU5(2)",
    "unexpected third maximal subgroup");
Check(maximal_subgroup_program.size = 82114560,
    "unexpected maximal-subgroup order");
maximal_subgroup_generators := ResultOfStraightLineProgram(
    maximal_subgroup_program.program, SUZ6_F19_GENERATORS);;

central_kernel_program := AtlasProgram("6.Suz", "kernel", "Suz");;
Check(central_kernel_program <> fail,
    "AtlasRep does not provide the central-kernel program");
central_generator := ResultOfStraightLineProgram(
    central_kernel_program.program, SUZ6_F19_GENERATORS)[1];;
Check(Order(central_generator) = 6, "the scalar centre does not have order 6");
Check(ForAll([1..12],
    i -> ForAll([1..12],
        j -> (i = j) or central_generator[i][j] = Zero(GF(19)))),
    "the central generator is not scalar");
Check(Length(Set(List([1..12], i -> central_generator[i][i]))) = 1,
    "the central generator has unequal diagonal entries");

seed_A := [1,18,0,1,8,0,7,11,18,12,0,0];;
seed_B := [1,4,14,5,8,12,0,16,2,5,14,14];;
field_seed_A := List(seed_A, x -> (x mod 19) * One(GF(19)));;
field_seed_B := List(seed_B, x -> (x mod 19) * One(GF(19)));;
functional_A := List(SUZ6_F19_GENERATORS,
    generator -> TransposedMat(generator^-1));;
functional_B := SUZ6_F19_GENERATORS;;
Check(IsProjectivelyFixed(field_seed_A,
    List(maximal_subgroup_generators,
        generator -> TransposedMat(generator^-1))),
    "the module-A seed is not fixed by 6xU5(2)");
Check(IsProjectivelyFixed(field_seed_B, maximal_subgroup_generators),
    "the module-B seed is not fixed by 6xU5(2)");

CheckModule("A", functional_A, seed_A, SUZ6_F19_VECTORS_A);
CheckModule("B", functional_B, seed_B, SUZ6_F19_VECTORS_B);

order_6Suz := 2690072985600;;
order_maximal_scalar_group := 3 * order_6Suz;;
regular_lower_bound := 15 * order_maximal_scalar_group;;
projective_threshold := 19^11 - 1;;
Check(regular_lower_bound > projective_threshold,
    "regular-vector lower bound does not exceed the projective threshold");

Print("\nEach module has at least 15 * ", order_maximal_scalar_group,
    " = ", regular_lower_bound, " regular nonzero vectors.\n");
Print(regular_lower_bound, " > 19^11 - 1 = ", projective_threshold, ".\n");
Print("Therefore R + R = V for both modules and every allowed scalar subgroup.\n");
QUIT;
