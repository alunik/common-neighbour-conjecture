///////////////////////////////////////////////////////////////////////////
// Independent exact quotient inventory for simple-diagonal actions in
// 10^7 < degree <= 10^8.
//
// For socle T^k, the full primitive diagonal normalizer has quotient
// Out(T) x S_k.  Its preimage is primitive precisely when the projection
// to S_k is transitive.  Huang's thesis closes the rows whose projected
// top is neither A_k nor S_k; the latter rows remain computational targets.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);

OuterGroup := function(tag)
    if tag eq "1" then
        return sub<Sym(1) | >;
    elif tag eq "C2" then
        return CyclicGroup(2);
    elif tag eq "C4" then
        return CyclicGroup(4);
    elif tag eq "C6" then
        return CyclicGroup(6);
    elif tag eq "V4" then
        ambient := Sym(4);
        return sub<ambient |
            ambient!(1,2)(3,4), ambient!(1,3)(2,4)>;
    end if;
    error "unknown outer group";
    return sub<Sym(1) | >;
end function;

// name, |T|, k, Out(T), expected selected, expected A_k/S_k residual,
// expected raw top-transitive classes.  Huang's theorem removes the
// non-A_k/S_k classes except in its exceptional A5,k=5 row.
rows := [
    <"PSL(2,19)", 3420, 3, "C2", 5, 5, 5>,
    <"PSL(2,16)", 4080, 3, "C4", 8, 8, 8>,
    <"PSL(3,3)", 5616, 3, "C2", 5, 5, 5>,
    <"PSU(3,3)", 6048, 3, "C2", 5, 5, 5>,
    <"PSL(2,23)", 6072, 3, "C2", 5, 5, 5>,
    <"PSL(2,25)", 7800, 3, "V4", 16, 16, 16>,
    <"M11", 7920, 3, "1", 2, 2, 2>,
    <"PSL(2,27)", 9828, 3, "C6", 12, 12, 12>,
    <"A6", 360, 4, "V4", 16, 16, 68>,
    <"A5", 60, 5, "C2", 13, 5, 13>
];

total := 0;
residual := 0;
for row in rows do
    name := row[1];
    t_order := row[2];
    k := row[3];
    outer := OuterGroup(row[4]);
    top := Sym(k);
    Q, inclusions, projections := DirectProduct(outer, top);
    assert Order(Q) eq Order(outer) * Factorial(k);

    degree := t_order^(k - 1);
    assert 10000000 lt degree and degree le 100000000;
    count := 0;
    residual_count := 0;
    raw_count := 0;
    for subgroup_record in Subgroups(Q) do
        B := subgroup_record`subgroup;
        projected_top := sub<top |
            [g @ projections[2] : g in Generators(B)]>;
        if not IsTransitive(projected_top) then continue; end if;

        // Necessary order gate for a base of size two.
        stabilizer_order := t_order * Order(B);
        if stabilizer_order gt degree - 1 then continue; end if;

        raw_count +:= 1;
        is_residual := Order(projected_top) in
                       {Factorial(k) div 2, Factorial(k)};
        selected := k eq 3 or name eq "A5" or is_residual;
        if not selected then continue; end if;
        count +:= 1;
        if is_residual then residual_count +:= 1; end if;
        printf "SD_QUOTIENT|%o|k=%o|degree=%o|outer=%o|B=%o|top=%o|residual=%o\n",
               name, k, degree, row[4], Order(B), Order(projected_top),
               is_residual;
    end for;
    assert count eq row[5] and residual_count eq row[6] and
           raw_count eq row[7];
    total +:= count;
    residual +:= residual_count;
    printf "SD_PAIR_COMPLETE|%o|k=%o|degree=%o|selected=%o|residual=%o|raw=%o\n",
           name, k, degree, count, residual_count, raw_count;
end for;

assert total eq 87 and residual eq 79;
print "SD8_QUOTIENT_AUDIT_PASS|pairs=10|classes=87|residual=79";
