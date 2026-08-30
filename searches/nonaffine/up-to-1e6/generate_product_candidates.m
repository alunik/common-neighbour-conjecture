///////////////////////////////////////////////////////////////////////////
// Complete base-two-candidate generator for primitive product actions of
// degree at most 10^6.
//
// A product action has degree n^k and socle order s^k, where s is the
// socle order of its primitive almost-simple or diagonal component.  Base
// size two forces s^k <= n^k(n^k-1), hence s < n^2 and, since n divides s,
// s <= n(n-1).  Components failing this bound are therefore excluded as
// whole cohorts.  The component degree is at most 1000, so the established
// primitive catalogue is used only to obtain the unique largest member of
// each component cohort, exactly as in Stratford's construction.
///////////////////////////////////////////////////////////////////////////

require assigned ASXMaxDegree: "set ASXMaxDegree";
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
if not assigned ASXValidateKnown then ASXValidateKnown := false; end if;
if Type(ASXValidateKnown) eq MonStgElt then
    ASXValidateKnown := ASXValidateKnown in {"1", "true", "True"};
end if;
if not assigned ASXBaseTwoOnly then ASXBaseTwoOnly := false; end if;
if Type(ASXBaseTwoOnly) eq MonStgElt then
    ASXBaseTwoOnly := ASXBaseTwoOnly in {"1", "true", "True"};
end if;
if not assigned ASXValidateComponentNormalizers then
    ASXValidateComponentNormalizers := false;
end if;
if Type(ASXValidateComponentNormalizers) eq MonStgElt then
    ASXValidateComponentNormalizers :=
        ASXValidateComponentNormalizers in {"1", "true", "True"};
end if;
require 25 le ASXMaxDegree and ASXMaxDegree le 1000000: "bad degree cap";
SetColumns(0);

ASXPrintPermutation := procedure(g, d)
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

ASXTransitiveGenerators := function(G)
    d := Degree(G);
    for length in [2, 3] do
        for attempt in [1 .. 3] do
            trial := [Random(G) : i in [1 .. length]];
            if #Orbit(sub<G | trial>, 1) eq d then return trial; end if;
        end for;
    end for;
    chosen := [G | ];
    for g in Setseq(Generators(G)) do
        Append(~chosen, g);
        if #Orbit(sub<G | chosen>, 1) eq d then return chosen; end if;
    end for;
    assert false;
    return chosen;
end function;

print "PRIMITIVE_SAXL_V1";
action_count := 0;
eligible_count := 0;
component_count := 0;
order_filtered_layers := 0;
maximum_component_degree := Floor(Sqrt(ASXMaxDegree));

for n in [5 .. maximum_component_degree] do
    candidates := [* *];
    if NumberOfPrimitiveAlmostSimpleGroups(n) gt 0 then
        for P in PrimitiveGroups(n : Filter := "AlmostSimple") do
            Append(~candidates, P);
        end for;
    end if;
    if NumberOfPrimitiveDiagonalGroups(n) gt 0 then
        for P in PrimitiveGroups(n : Filter := "Diagonal") do
            Append(~candidates, P);
        end for;
    end if;
    if #candidates eq 0 then continue; end if;
    normalizers := [* *];
    eligible_candidates := [* P : P in candidates |
                                  Order(Socle(P)) le n * (n - 1) *];

    // The complete primitive catalogue already contains the full symmetric
    // normalizer N_Sym(n)(T): it is primitive and of the same O'Nan--Scott
    // type as T.  Therefore, within each Sym(n)-conjugacy class of component
    // socles, the catalogue member of greatest order is exactly that
    // normalizer (up to permutation conjugacy).  Selecting it from the
    // catalogue avoids a very expensive generic normalizer computation in
    // degrees near 1000 while preserving the constructed wreath products up
    // to permutation isomorphism.
    while #eligible_candidates gt 0 do
        component_socle := Socle(eligible_candidates[1]);
        cohort_positions := [i : i in [1 .. #eligible_candidates] |
                                  Socle(eligible_candidates[i]) eq component_socle or
                                  IsConjugate(Sym(n),
                                              Socle(eligible_candidates[i]),
                                              component_socle)];
        largest_position := cohort_positions[1];
        for position in cohort_positions do
            if Order(eligible_candidates[position]) gt
               Order(eligible_candidates[largest_position]) then
                largest_position := position;
            end if;
        end for;
        component_normalizer := eligible_candidates[largest_position];
        assert IsPrimitive(component_normalizer);

        if ASXValidateComponentNormalizers then
            computed_normalizer := Normalizer(Sym(n), component_socle);
            assert Order(computed_normalizer) eq Order(component_normalizer);
            assert IsConjugate(Sym(n), computed_normalizer,
                               component_normalizer);
        end if;

        Append(~normalizers, component_normalizer);
        for position in Reverse(Sort(cohort_positions)) do
            Remove(~eligible_candidates, position);
        end for;
    end while;

    for component_normalizer in normalizers do
        component_count +:= 1;
        component_id, component_degree :=
            PrimitiveGroupIdentification(component_normalizer);
        assert component_degree eq n;
        component_socle_order := Order(Socle(component_normalizer));
        degree := n * n;
        exponent := 2;
        while degree le ASXMaxDegree do
            top := Sym(exponent);
            W := PrimitiveWreathProduct(component_normalizer, top);
            assert Degree(W) eq degree;
            S := Socle(W);
            assert Order(S) eq component_socle_order^exponent;
            Q, quotient_map := quo<W | S>;
            subgroup_number := 0;

            if ASXBaseTwoOnly then
                maximum_quotient_order := degree * (degree - 1) div Order(S);
                allowed_orders := [order : order in Divisors(Order(Q)) |
                                           order le maximum_quotient_order and
                                           order mod exponent eq 0];
                subgroup_records := [];
                for allowed_order in allowed_orders do
                    subgroup_records cat:=
                        Subgroups(Q : OrderEqual := allowed_order);
                end for;
                order_filtered_layers +:= 1;
            else
                subgroup_records := Subgroups(Q);
            end if;

            for subgroup_record in subgroup_records do
                quotient_subgroup := subgroup_record`subgroup;
                stabilizer_order := Order(S) * Order(quotient_subgroup) div degree;
                assert stabilizer_order * degree eq
                       Order(S) * Order(quotient_subgroup);
                // In production search mode, primitivity of an already
                // order-obstructed preimage is immaterial: no such preimage
                // can have base size two.  Avoid constructing it in degree
                // up to 10^6.  Validation mode retains the full enumeration.
                if ASXBaseTwoOnly and stabilizer_order gt degree - 1 then
                    error "order-restricted subgroup enumeration failed";
                end if;
                G := quotient_subgroup @@ quotient_map;
                if not IsPrimitive(G) then continue; end if;
                subgroup_number +:= 1;
                assert Order(G) div degree eq stabilizer_order;
                classification := stabilizer_order gt degree - 1
                                  select "order_obstruction" else "compute";
                regular_orbits := 0;
                regular_points := [Integers() | ];
                orbit_representatives := [Integers() | ];
                orbit_sizes := [Integers() | ];
                if classification eq "compute" then
                    H := Stabilizer(G, 1);
                    assert Order(H) eq stabilizer_order;
                    for orbit in Orbits(H) do
                        Append(~orbit_representatives, Min(orbit));
                        Append(~orbit_sizes, #orbit);
                        if #orbit eq stabilizer_order then
                            regular_orbits +:= 1;
                            regular_points cat:= Sort(Setseq(orbit));
                        end if;
                    end for;
                    if #regular_points eq 0 then
                        classification := "no_regular_orbit";
                    elif #regular_points eq degree - 1 then
                        classification := "complete";
                    elif 2 * #regular_points gt degree then
                        classification := "density";
                    else
                        classification := "graph";
                    end if;
                end if;
                label := Sprintf("PA_comp%o_%o_k%o_sub%o_d%o",
                                 n, component_id, exponent,
                                 subgroup_number, degree);
                print "action";
                printf "label %o\n", label;
                printf "degree %o\n", degree;
                printf "stabilizer_order %o\n", stabilizer_order;
                printf "classification %o\n", classification;
                printf "regular_orbits %o\n", regular_orbits;
                printf "regular_count %o\n", #regular_points;
                if classification eq "graph" then
                    printf "regular ";
                    for point in regular_points do printf "%o ", point; end for;
                    print "";
                    printf "orbit_representatives_count %o\n",
                           #orbit_representatives;
                    printf "orbit_representatives ";
                    for index in [1 .. #orbit_representatives] do
                        printf "%o %o ", orbit_representatives[index],
                                          orbit_sizes[index];
                    end for;
                    print "";
                    gens := ASXTransitiveGenerators(G);
                    assert #Orbit(sub<G | gens>, 1) eq degree;
                    printf "gens %o\n", #gens;
                    for g in gens do ASXPrintPermutation(g, degree); end for;
                end if;
                if classification ne "order_obstruction" then eligible_count +:= 1; end if;
                print "end";
                action_count +:= 1;

                if ASXValidateKnown and degree le PrimitiveGroupDatabaseLimit() then
                    known_id, known_degree := PrimitiveGroupIdentification(G);
                    assert known_degree eq degree;
                    fprintf "/dev/stderr",
                            "ASX_PA_KNOWN label=%o degree=%o id=%o order=%o\n",
                            label, degree, known_id, Order(G);
                end if;
            end for;
            fprintf "/dev/stderr",
                    "ASX_PRODUCT_COMPLETE component_degree=%o component_id=%o exponent=%o degree=%o total_actions=%o eligible=%o\n",
                    n, component_id, exponent, degree,
                    action_count, eligible_count;
            if degree gt ASXMaxDegree div n then break; end if;
            degree *:= n;
            exponent +:= 1;
        end while;
    end for;
end for;
fprintf "/dev/stderr",
        "ASX_PRODUCT_GENERATOR_COMPLETE max_degree=%o components=%o actions=%o eligible=%o order_filtered_layers=%o base_two_only=%o\n",
        ASXMaxDegree, component_count, action_count, eligible_count,
        order_filtered_layers, ASXBaseTwoOnly;
quit;
