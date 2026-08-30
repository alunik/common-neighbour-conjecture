///////////////////////////////////////////////////////////////////////////
// Complete diagonal-type generator through degree 10^6.
//
// For T nonabelian simple and m >= 2 the degree is |T|^(m-1), and the
// normalizer quotient is Out(T) x S_m.  For m=2 every preimage is primitive
// but its point stabilizer contains a diagonal T of order d, so base size
// two is impossible; those actions are emitted symbolically.  The very few
// m>=3 cases are constructed explicitly and sent to the Saxl engine.
///////////////////////////////////////////////////////////////////////////

require assigned ASXMaxDegree: "set ASXMaxDegree";
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
require 5 le ASXMaxDegree and ASXMaxDegree le 1000000: "bad degree cap";
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

ASXBuildDiagonalData := function(A, T, m)
    r := Degree(A);
    ambient := Sym(r * m);
    EmbedBlock := function(g, block)
        images := [1 .. r * m];
        offset := (block - 1) * r;
        for point in [1 .. r] do
            images[offset + point] := offset + point ^ g;
        end for;
        return ambient ! images;
    end function;
    EmbedDiagonal := function(g)
        images := [Integers() | ];
        for block in [1 .. m] do
            offset := (block - 1) * r;
            images cat:= [offset + point ^ g : point in [1 .. r]];
        end for;
        return ambient ! images;
    end function;
    EmbedTop := function(g)
        images := [Integers() | ];
        for block in [1 .. m] do
            images cat:= [(block ^ g - 1) * r + point : point in [1 .. r]];
        end for;
        return ambient ! images;
    end function;

    top := Sym(m);
    n_generators := [EmbedBlock(g, block) : g in Setseq(Generators(T)),
                                              block in [1 .. m]];
    diagonal_generators := [EmbedDiagonal(g) : g in Setseq(Generators(A))];
    top_generators := [EmbedTop(g) : g in Setseq(Generators(top))];
    N := sub<ambient | n_generators>;
    D := sub<ambient | diagonal_generators cat top_generators>;
    K := sub<ambient | n_generators cat diagonal_generators cat top_generators>;
    assert Order(N) eq Order(T)^m;
    assert Order(D) eq Order(A) * Factorial(m);
    assert Order(K) eq Order(T)^m * (Order(A) div Order(T)) * Factorial(m);
    return K, N, D;
end function;

print "PRIMITIVE_SAXL_V1";
action_count := 0;
eligible_count := 0;
for simple_id in [1 .. NumberOfSimpleGroups()] do
    simple_tuple, simple_order := SimpleGroupId(simple_id);
    if simple_order gt ASXMaxDegree then break; end if;
    A := AutomorphismGroupSimpleGroup(simple_tuple);
    T := Socle(A);
    assert Order(T) eq simple_order;

    degree := simple_order;
    m := 2;
    while degree le ASXMaxDegree do
        K, N, D := ASXBuildDiagonalData(A, T, m);
        Q, quotient_map := quo<K | N>;

        if m eq 2 then
            // Every subgroup of K/N gives a primitive diagonal group.  No
            // degree-d permutation construction is needed for the order
            // obstruction |G_1| >= |T| = d.
            subgroup_number := 0;
            for subgroup_record in Subgroups(Q) do
                subgroup_number +:= 1;
                quotient_order := Order(subgroup_record`subgroup);
                stabilizer_order := simple_order * quotient_order;
                assert stabilizer_order gt degree - 1;
                print "action";
                printf "label SD_sid%o_m%o_sub%o_d%o\n",
                       simple_id, m, subgroup_number, degree;
                printf "degree %o\n", degree;
                printf "stabilizer_order %o\n", stabilizer_order;
                print "classification order_obstruction";
                print "regular_orbits 0";
                print "regular_count 0";
                print "end";
                action_count +:= 1;
            end for;
        else
            action_map, P := CosetAction(K, D);
            assert Degree(P) eq degree;
            assert Order(Kernel(action_map)) eq 1;
            S := Socle(P);
            assert Order(S) eq simple_order^m;
            QP, map_P := quo<P | S>;
            subgroup_number := 0;
            for subgroup_record in Subgroups(QP) do
                G := subgroup_record`subgroup @@ map_P;
                if not IsPrimitive(G) then continue; end if;
                subgroup_number +:= 1;
                stabilizer_order := Order(G) div degree;
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
                print "action";
                printf "label SD_sid%o_m%o_sub%o_d%o\n",
                       simple_id, m, subgroup_number, degree;
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
            end for;
        end if;

        fprintf "/dev/stderr",
                "ASX_DIAGONAL_COMPLETE id=%o name=%o m=%o degree=%o total_actions=%o eligible=%o\n",
                simple_id, SimpleGroupName(simple_id), m, degree,
                action_count, eligible_count;
        if degree gt ASXMaxDegree div simple_order then break; end if;
        degree *:= simple_order;
        m +:= 1;
    end while;
end for;
fprintf "/dev/stderr",
        "ASX_DIAGONAL_GENERATOR_COMPLETE max_degree=%o actions=%o eligible=%o\n",
        ASXMaxDegree, action_count, eligible_count;
quit;
