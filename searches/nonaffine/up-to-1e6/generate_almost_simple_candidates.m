///////////////////////////////////////////////////////////////////////////
// Database-independent almost-simple candidate generator for the Saxl
// search through degree ASXMaxDegree (at most 10^6).
//
// The key reduction is rigorous.  If an almost-simple primitive action of
// degree d has base size two, then
//
//        |Soc(G)| <= |G| = d |G_1| <= d(d-1).
//
// Thus for d <= 10^6 only nonabelian simple socles of order < 10^12 need
// be constructed.  Magma's simple-group list, automorphism constructors,
// and maximal-subgroup algorithms are used; PrimitiveGroup(s) is never
// called unless ASXValidateKnown is explicitly set for a validation run.
//
// Required globals:
//   ASXFirstSimpleId, ASXLastSimpleId, ASXMaxDegree
// Optional globals:
//   ASXValidateKnown (default false)
//   ASXSelectedLabel (default unset; emit only the matching action)
//   ASXExportOrbitIndexMap (default false; append the exact H-orbit index
//                           of every point for compressed orbital-transition
//                           calculations)
//
// Standard output is PRIMITIVE_SAXL_V1 for the C++ engine.  Progress and
// optional known-catalogue identifications are written to standard error.
///////////////////////////////////////////////////////////////////////////

require assigned ASXFirstSimpleId and assigned ASXLastSimpleId and
        assigned ASXMaxDegree:
    "set ASXFirstSimpleId, ASXLastSimpleId, ASXMaxDegree";
if Type(ASXFirstSimpleId) eq MonStgElt then
    ASXFirstSimpleId := StringToInteger(ASXFirstSimpleId);
end if;
if Type(ASXLastSimpleId) eq MonStgElt then
    ASXLastSimpleId := StringToInteger(ASXLastSimpleId);
end if;
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
if not assigned ASXValidateKnown then ASXValidateKnown := false; end if;
if Type(ASXValidateKnown) eq MonStgElt then
    ASXValidateKnown := ASXValidateKnown in {"1", "true", "True"};
end if;
if not assigned ASXBaseTwoOnly then
    ASXBaseTwoOnly := not ASXValidateKnown;
end if;
if not assigned ASXExportOrbitIndexMap then
    ASXExportOrbitIndexMap := false;
end if;
if Type(ASXExportOrbitIndexMap) eq MonStgElt then
    ASXExportOrbitIndexMap :=
        ASXExportOrbitIndexMap in {"1", "true", "True"};
end if;
if Type(ASXBaseTwoOnly) eq MonStgElt then
    ASXBaseTwoOnly := ASXBaseTwoOnly in {"1", "true", "True"};
end if;

require 1 le ASXFirstSimpleId and ASXFirstSimpleId le ASXLastSimpleId and
        ASXLastSimpleId le NumberOfSimpleGroups(): "bad simple-group interval";
require 5 le ASXMaxDegree and ASXMaxDegree le 1000000: "bad degree cap";
SetColumns(0);

ASXPrintPermutation := procedure(g, d)
    // A degree-d image needs only ceil(log_94(d)) printable bytes.  This is
    // substantially smaller than space-separated decimal images while
    // remaining streamable and backward-compatible in the C++ parser.
    width := 1;
    capacity := 94;
    while capacity lt d do
        width +:= 1;
        capacity *:= 94;
    end while;
    printf "packed_gen %o ", width;
    for point in [1 .. d] do
        value := point ^ g - 1;
        encoded := "";
        for digit_index in [1 .. width] do
            encoded cat:= CodeToString(33 + (value mod 94));
            value := value div 94;
        end for;
        assert value eq 0;
        printf "%o", encoded;
    end for;
    print "";
end procedure;

// Neighbourhood transport only requires graph automorphisms acting
// transitively; the exported permutations need not generate all of G.  Try
// short random sequences, certify transitivity exactly in the coset action,
// and fall back to a guaranteed prefix of the original generators.
ASXTransitiveActionGenerators := function(P, G, action_map)
    d := Degree(P);
    for length in [2, 3] do
        for attempt in [1 .. 3] do
            source_trial := [Random(G) : i in [1 .. length]];
            trial := [g @ action_map : g in source_trial];
            if #Orbit(sub<P | trial>, 1) eq d then
                return trial;
            end if;
        end for;
    end for;
    chosen := [P | ];
    for g in Setseq(Generators(G)) do
        Append(~chosen, g @ action_map);
        if #Orbit(sub<P | chosen>, 1) eq d then
            return chosen;
        end if;
    end for;
    assert false;
    return chosen;
end function;

ASXPSL2CohortOrderObstructed := function(q, simple_order, degree_cap)
    if q lt 1415 then return false; end if;
    factorisation := Factorization(q);
    assert #factorisation eq 1;
    prime := factorisation[1][1];
    field_exponent := factorisation[1][2];

    // Besides Borel, torus-normaliser, and exceptional subgroups, Dickson's
    // list has subfield groups for prime divisors r of the field exponent.
    // Test both PSL and PGL orders conservatively; if either gives an action
    // inside the cap which is not order-obstructed, retain the whole cohort.
    for r in PrimeDivisors(field_exponent) do
        q0 := prime^(field_exponent div r);
        psl_order := q0 * (q0^2 - 1) div GCD(2, q0 - 1);
        pgl_order := q0 * (q0^2 - 1);
        // Dickson's square-field case has PGL(2,q0) maximal; PSL(2,q0)
        // lies inside it.  For odd prime field index the maximal subfield
        // subgroup is PSL(2,q0).
        possible_orders := r eq 2 select [pgl_order] else [psl_order];
        for subgroup_order in possible_orders do
            if simple_order mod subgroup_order ne 0 then continue; end if;
            action_degree := simple_order div subgroup_order;
            if action_degree le degree_cap and
               subgroup_order le action_degree - 1 then
                return false;
            end if;
        end for;
    end for;
    return true;
end function;

ASXEmit := procedure(G, H, simple_id, extension_number, maximal_number)
    d := Index(G, H);
    h_order := Order(H);
    label := Sprintf("AS_sid%o_ext%o_max%o_d%o",
                     simple_id, extension_number, maximal_number, d);
    if assigned ASXSelectedLabel and label ne ASXSelectedLabel then
        return;
    end if;
    classification := h_order gt d - 1 select "order_obstruction" else "compute";
    regular_orbits := 0;
    regular_points := [Integers() | ];
    orbit_representatives := [Integers() | ];
    orbit_sizes := [Integers() | ];
    orbit_index_map := [Integers() | 0 : point in [1 .. d]];
    if classification eq "compute" then
        action_map, P := CosetAction(G, H);
        assert Degree(P) eq d;
        assert IsPrimitive(P);
        assert Order(P) eq Order(G);
        assert Order(Kernel(action_map)) eq 1;
        permutation_stabilizer := Stabilizer(P, 1);
        assert Order(permutation_stabilizer) eq h_order;
        for orbit in Orbits(permutation_stabilizer) do
            Append(~orbit_representatives, Min(orbit));
            Append(~orbit_sizes, #orbit);
            if ASXExportOrbitIndexMap then
                orbit_index := #orbit_representatives;
                for point in orbit do
                    orbit_index_map[point] := orbit_index;
                end for;
            end if;
            if #orbit eq h_order then
                regular_orbits +:= 1;
                regular_points cat:= Sort(Setseq(orbit));
            end if;
        end for;
        if #regular_points eq 0 then
            classification := "no_regular_orbit";
        elif #regular_points eq d - 1 then
            classification := "complete";
        elif 2 * #regular_points gt d then
            classification := "density";
        else
            classification := "graph";
        end if;
    end if;

    // The order obstruction needs no degree-d permutation construction.
    // A validation run still constructs the action for independent catalogue
    // identification.
    if classification eq "order_obstruction" and ASXValidateKnown and
       d le PrimitiveGroupDatabaseLimit() then
        action_map, P := CosetAction(G, H);
        assert Degree(P) eq d;
        assert IsPrimitive(P);
        assert Order(P) eq Order(G);
        assert Order(Kernel(action_map)) eq 1;
    end if;

    print "action";
    printf "label %o\n", label;
    printf "degree %o\n", d;
    printf "stabilizer_order %o\n", h_order;
    printf "classification %o\n", classification;
    printf "regular_orbits %o\n", regular_orbits;
    printf "regular_count %o\n", #regular_points;
    if classification eq "graph" then
        printf "regular ";
        for point in regular_points do printf "%o ", point; end for;
        print "";
        printf "orbit_representatives_count %o\n", #orbit_representatives;
        printf "orbit_representatives ";
        for index in [1 .. #orbit_representatives] do
            printf "%o %o ", orbit_representatives[index], orbit_sizes[index];
        end for;
        print "";
        if ASXExportOrbitIndexMap then
            assert not 0 in orbit_index_map;
            printf "orbit_index_map ";
            for orbit_index in orbit_index_map do
                printf "%o ", orbit_index;
            end for;
            print "";
        end if;
        gens := ASXTransitiveActionGenerators(P, G, action_map);
        assert #Orbit(sub<P | gens>, 1) eq d;
        printf "gens %o\n", #gens;
        for g in gens do ASXPrintPermutation(g, d); end for;
    end if;
    print "end";

    if ASXValidateKnown and d le PrimitiveGroupDatabaseLimit() then
        known_id, known_degree := PrimitiveGroupIdentification(P);
        assert known_degree eq d;
        fprintf "/dev/stderr",
                "ASX_KNOWN label=%o degree=%o id=%o order=%o\n",
                label, d, known_id, Order(P);
    end if;
end procedure;

print "PRIMITIVE_SAXL_V1";
ASXActions := 0;
ASXEligible := 0;
for simple_id in [ASXFirstSimpleId .. ASXLastSimpleId] do
    simple_tuple, simple_order := SimpleGroupId(simple_id);
    if simple_order gt ASXMaxDegree * (ASXMaxDegree - 1) then
        break;
    end if;

    // Dickson's maximal-subgroup list gives a complete production sieve for
    // PSL(2,q).  For q >= 1415, the split/nonsplit torus
    // normalisers have index at least q(q-1)/2 > 10^6, and the exceptional
    // A4/S4/A5 indices are larger still.  The only remaining action within
    // the degree cap is the projective-line action of degree q+1, whose
    // Borel stabiliser has order q(q-1)/gcd(2,q-1) > q.  The same conclusion
    // holds for the PGL extension.  Subfield subgroups are checked separately
    // and conservatively by ASXPSL2CohortOrderObstructed.
    if ASXBaseTwoOnly and simple_tuple[1] eq 1 and
       simple_tuple[2] eq 1 and
       ASXPSL2CohortOrderObstructed(simple_tuple[3], simple_order,
                                    ASXMaxDegree) then
        fprintf "/dev/stderr",
                "ASX_SIMPLE_COHORT_ORDER_OBSTRUCTION id=%o name=%o order=%o reason=PSL2_Dickson_q_ge_1415\n",
                simple_id, SimpleGroupName(simple_id), simple_order;
        continue;
    end if;

    // Magma's AutomorphismGroupSimpleGroup constructor for O'Nan depends on
    // an optional local ATLAS representation file.  This production-only
    // exclusion is independent of that file.  The complete maximal-subgroup
    // lists in GAP's Character Table Library give, for O'N, precisely two
    // classes below the degree cap, both of index 122760 and stabilizer order
    // 3753792 > 122759.  For O'N.2 every core-free maximal subgroup has index
    // greater than 10^6.  The companion GAP certificate reproduces every
    // maximal-subgroup order and index used here.
    if ASXBaseTwoOnly and simple_tuple eq <18, 21, 0> then
        fprintf "/dev/stderr",
                "ASX_SIMPLE_COHORT_ORDER_OBSTRUCTION id=%o name=%o order=%o reason=ON_GAP_CTBLLIB_MAXES\n",
                simple_id, SimpleGroupName(simple_id), simple_order;
        continue;
    end if;

    A := AutomorphismGroupSimpleGroup(simple_tuple);
    T := Socle(A);
    assert Order(T) eq simple_order;
    Q, quotient_map := quo<A | T>;
    extension_records := Subgroups(Q);
    extension_number := 0;

    for extension_record in extension_records do
        extension_number +:= 1;
        G := extension_record`subgroup @@ quotient_map;
        assert T subset G;
        normalizer_G := Normalizer(A, G);
        maximals := MaximalSubgroups(G : IndexLimit := ASXMaxDegree);
        representatives := [];
        maximal_number := 0;

        for maximal_record in maximals do
            H := maximal_record`subgroup;
            if T subset H then continue; end if;
            d := Index(G, H);
            if d lt 5 or d gt ASXMaxDegree then continue; end if;

            // Fuse G-conjugacy classes which are interchanged by an outer
            // automorphism normalising G.  This is exactly the equivalence
            // needed for permutation isomorphism of the coset actions.
            duplicate := false;
            for old_H in representatives do
                if Index(G, old_H) eq d and
                   IsConjugate(normalizer_G, H, old_H) then
                    duplicate := true;
                    break;
                end if;
            end for;
            if duplicate then continue; end if;
            Append(~representatives, H);

            maximal_number +:= 1;
            ASXEmit(G, H, simple_id, extension_number, maximal_number);
            ASXActions +:= 1;
            if Order(H) le d - 1 then ASXEligible +:= 1; end if;
        end for;
    end for;
    fprintf "/dev/stderr",
            "ASX_SIMPLE_COMPLETE id=%o name=%o order=%o total_actions=%o eligible_actions=%o\n",
            simple_id, SimpleGroupName(simple_id), simple_order,
            ASXActions, ASXEligible;
end for;
fprintf "/dev/stderr",
        "ASX_GENERATOR_COMPLETE first=%o last=%o max_degree=%o actions=%o eligible=%o\n",
        ASXFirstSimpleId, ASXLastSimpleId, ASXMaxDegree,
        ASXActions, ASXEligible;
quit;
