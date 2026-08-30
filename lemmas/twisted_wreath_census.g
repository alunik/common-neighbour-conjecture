# Independent exact census for the A5^6 twisted-wreath input.
#
# If L <= P0, an L-fixed point in the induced-function model is determined
# independently on every double coset LxQ0.  Its value on LxQ0 must be fixed
# by phi0(Q0 cap x^-1 L x), so
#
# |Fix_{A5^6}(L)| =
#   Product_{x in L\P0/Q0} |C_A5(phi0(Q0 cap x^-1 L x))|.
#
# Inverting the table of marks decomposes the A5^6 action into transitive
# P0-sets.  Entry 1 is the number of regular P0-orbits.  This computation is
# independent of the coordinate-cycle calculation used for the witness.

SizeScreen([1000, 24]);

FixCount := function(P, Q, T, L)
  local double_cosets, double_coset, x, intersection, centraliser, answer;
  double_cosets := DoubleCosets(P, L, Q);
  answer := 1;
  for double_coset in double_cosets do
    x := Representative(double_coset);
    intersection := Intersection(Q, L^x);
    centraliser := Filtered(
      Elements(T),
      t -> ForAll(
        GeneratorsOfGroup(intersection),
        h -> t*h = h*t
      )
    );
    answer := answer * Length(centraliser);
  od;
  return answer;
end;

Census := function(P, Q, T, label, expected_regular, expected_nonregular)
  local tom, marks, fixed, i, L, orbit_coefficients, degree;
  tom := TableOfMarks(P);
  marks := MatTom(tom);
  fixed := [];
  for i in [1..Length(marks)] do
    L := RepresentativeTom(tom, i);
    fixed[i] := FixCount(P, Q, T, L);
  od;
  orbit_coefficients := fixed * Inverse(marks);
  if orbit_coefficients * marks <> fixed then
    Error("table-of-marks inversion failed");
  fi;
  if ForAny(orbit_coefficients, x -> not IsInt(x) or x < 0) then
    Error("invalid orbit decomposition");
  fi;
  degree := Size(T)^6;
  if orbit_coefficients[1] <> expected_regular then
    Error("regular-orbit count mismatch");
  fi;
  if degree - expected_regular * Size(P) <> expected_nonregular then
    Error("nonregular-point count mismatch");
  fi;
  Print(
    "CENSUS case=", label,
    " top_order=", Size(P),
    " degree=", degree,
    " subgroup_classes=", Length(marks),
    " regular_orbits=", orbit_coefficients[1],
    " regular_points=", orbit_coefficients[1] * Size(P),
    " nonregular_points=", expected_nonregular,
    " nonregular_orbits=", Sum(orbit_coefficients{[2..Length(marks)]}),
    "\n"
  );
end;

P := SymmetricGroup(6);
Q := Stabilizer(P, 1);
T := DerivedSubgroup(Q);
Census(P, Q, T, "s6", 64790243, 7025040);

QUIT;
