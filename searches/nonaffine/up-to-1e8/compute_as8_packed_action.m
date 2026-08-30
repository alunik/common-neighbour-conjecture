///////////////////////////////////////////////////////////////////////////
// Exact consumer for one AS8_INTRINSIC_ACTION_PACK_V1 payload.
//
// Load the immutable AS8_INTRINSIC_ACTION_PACK_V1 payload immediately before
// this file on the Magma command line.  The payload reconstructs the exact
// permutation groups G and H.  This
// consumer verifies the intrinsic binding, computes the prime-order FPR
// bound, and only when that bound is inconclusive computes H\G/H exactly.
///////////////////////////////////////////////////////////////////////////

require assigned AS8PackedSimpleId and assigned AS8PackedOrdinal and
        assigned AS8PackedMultiplicity and assigned AS8PackedEnumerationLabel and
        assigned AS8PackedDegree and assigned AS8PackedGroupOrder and
        assigned AS8PackedStabilizerOrder and assigned AS8PackedAmbientDegree and
        assigned AS8PackedG and assigned AS8PackedH:
    "incomplete AS8 intrinsic action pack";

SetColumns(0);

G := AS8PackedG;
H := AS8PackedH;
label := Sprintf("AS8_pack_sid%o_d%o_G%o_H%o_i%oof%o",
                 AS8PackedSimpleId, AS8PackedDegree, AS8PackedGroupOrder,
                 AS8PackedStabilizerOrder, AS8PackedOrdinal,
                 AS8PackedMultiplicity);

require Type(G) eq GrpPerm and Type(H) eq GrpPerm:
    "packed action is not a permutation-group pair";
require Degree(G) eq AS8PackedAmbientDegree and H subset G:
    "packed ambient/subgroup binding mismatch";
require Order(G) eq AS8PackedGroupOrder and
        Order(H) eq AS8PackedStabilizerOrder and
        Index(G, H) eq AS8PackedDegree:
    "packed intrinsic invariants mismatch";

AS8EmitFPR := function(G, H, label)
    degree := Index(G, H);
    printf "AS8_PACKED_ACTION_V1|%o|sid=%o|ordinal=%o|actions=%o|degree=%o|G=%o|H=%o|ambient_degree=%o|enumeration_label=%o\n",
           label, AS8PackedSimpleId, AS8PackedOrdinal,
           AS8PackedMultiplicity, degree, Order(G), Order(H), Degree(G),
           AS8PackedEnumerationLabel;

    started := Cputime();
    hclasses := Classes(H);
    prime_classes := [c : c in hclasses | IsPrime(c[1])];
    g_representatives := [G | ];
    intersections := [];
    g_class_sizes := [];
    for c in prime_classes do
        x := G!c[3];
        position := 0;
        for i in [1 .. #g_representatives] do
            if IsConjugate(G, x, g_representatives[i]) then
                position := i;
                break;
            end if;
        end for;
        if position eq 0 then
            Append(~g_representatives, x);
            Append(~intersections, c[2]);
            Append(~g_class_sizes, Order(G) div Order(Centralizer(G, x)));
        else
            intersections[position] +:= c[2];
        end if;
    end for;

    bound := &+[Rationals() | intersections[i]^2 / g_class_sizes[i]
                : i in [1 .. #intersections]];
    for i in [1 .. #intersections] do
        term := Rationals()!(intersections[i]^2) / g_class_sizes[i];
        printf "FPR_CLASS|sequence=%o|intersection=%o|G_class=%o|term_num=%o|term_den=%o\n",
               i, intersections[i], g_class_sizes[i], Numerator(term),
               Denominator(term);
    end for;
    printf "AS8_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o|cpu_ms=%o\n",
           label, #prime_classes, #intersections, Numerator(bound),
           Denominator(bound), bound lt 1/2,
           Round(1000 * Cputime(started));
    return bound;
end function;

AS8EmitDoubleCosets := function(G, H, label)
    degree := Index(G, H);
    h_order := Order(H);
    started := Cputime();
    representatives, double_coset_sizes :=
        DoubleCosetRepresentatives(G, H, H);
    elapsed := Cputime(started);
    assert #representatives eq #double_coset_sizes;
    assert #representatives ge 1 and representatives[1] eq Identity(G);

    subdegree_sum := 0;
    regular_orbits := 0;
    printf "AS8_DOUBLE_COSET_V1|%o|degree=%o|G=%o|H=%o|rank=%o|cpu_ms=%o\n",
           label, degree, Order(G), h_order, #representatives,
           Round(1000 * elapsed);
    for i in [1 .. #representatives] do
        double_size := double_coset_sizes[i];
        assert double_size mod h_order eq 0;
        subdegree := double_size div h_order;
        subdegree_sum +:= subdegree;
        if double_size eq h_order^2 then regular_orbits +:= 1; end if;
        printf "DOUBLE_COSET|sequence=%o|double_size=%o|subdegree=%o|regular=%o\n",
               i, double_size, subdegree, double_size eq h_order^2;
    end for;
    assert subdegree_sum eq degree;
    printf "AS8_DOUBLE_COSET_COMPLETE|%o|rank=%o|regular_orbits=%o|subdegree_sum=%o|degree=%o|cpu_ms=%o\n",
           label, #representatives, regular_orbits, subdegree_sum, degree,
           Round(1000 * elapsed);
    return regular_orbits, #representatives;
end function;

bound := AS8EmitFPR(G, H, label);
if bound lt 1/2 then
    disposition := "fpr_density";
    regular_orbits := 0;
    regular_points := 0;
    rank := 0;
else
    regular_orbits, rank := AS8EmitDoubleCosets(G, H, label);
    regular_points := regular_orbits * Order(H);
    if regular_orbits eq 0 then
        disposition := "base_gt_2";
    elif 2 * regular_points gt AS8PackedDegree then
        disposition := "exact_density";
    else
        disposition := "graph_residual";
    end if;
end if;

printf "AS8_PACKED_ACTION_COMPLETE_V1|%o|sid=%o|ordinal=%o|actions=%o|degree=%o|G=%o|H=%o|disposition=%o|rank=%o|regular_orbits=%o|regular_points=%o\n",
       label, AS8PackedSimpleId, AS8PackedOrdinal, AS8PackedMultiplicity,
       AS8PackedDegree, AS8PackedGroupOrder, AS8PackedStabilizerOrder,
       disposition, rank, regular_orbits, regular_points;
