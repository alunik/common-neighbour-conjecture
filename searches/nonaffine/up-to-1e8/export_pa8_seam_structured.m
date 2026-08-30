///////////////////////////////////////////////////////////////////////////
// Exact compact structured-pack export for PA component actions beyond the
// installed primitive-group catalogue: 8192 <= n <= 10000.
//
// The base-two order sieve gives |T| <= n(n-1).  We enumerate every simple
// T under that bound, every almost-simple extension T <= G <= Aut(T), and
// every maximal subgroup of index in the seam.  Actions are grouped by the
// Aut(T)-class of H meet T and the largest preserving extension is retained.
// The final quotient census is exact in (N/T) wr S2.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);

PAXPrintPermutation := procedure(prefix, g, d)
    width := 1;
    capacity := 94;
    while capacity lt d do width +:= 1; capacity *:= 94; end while;
    printf "%o %o ", prefix, width;
    block_size := 4096;
    for first in [1 .. d by block_size] do
        pieces := [""];
        for point in [first .. Min(d, first + block_size - 1)] do
            value := point ^ g - 1;
            encoded := "";
            for digit_index in [1 .. width] do
                encoded cat:= CodeToString(33 + (value mod 94));
                value div:= 94;
            end for;
            assert value eq 0;
            Append(~pieces, encoded);
        end for;
        printf "%o", &cat pieces;
    end for;
    print "";
end procedure;

PAXPrintTopPermutation := procedure(g, d)
    printf "top_gen";
    for point in [1 .. d] do printf " %o", point ^ g; end for;
    print "";
end procedure;

PAXDecompose := function(element, small_quotient, top_projection,
                         top_inclusion, base_inclusions, coordinate_tuples)
    top_element := element @ top_projection;
    top_lift := top_element @ top_inclusion;
    base_element := element * top_lift^-1;
    assert base_element in Kernel(top_projection);
    for coordinate_tuple in coordinate_tuples do
        candidate := Identity(small_quotient);
        for coordinate in [1 .. #base_inclusions] do
            candidate *:= coordinate_tuple[coordinate] @
                          base_inclusions[coordinate];
        end for;
        if candidate eq base_element then
            return coordinate_tuple, top_element;
        end if;
    end for;
    error "could not decompose seam quotient generator";
end function;

minimum_degree := 8192;
maximum_degree := 10000;
maximum_simple_order := maximum_degree * (maximum_degree - 1);

// sid, component degree, |T|, |N|, |N/T|, quotient candidates.
expected := [
    <32, 8555, 102660, 205320, 2, 3>,
    <33, 9455, 113460, 226920, 2, 3>,
    <48, 8585, 515100, 515100, 1, 1>,
    <58, 8646, 1123980, 2247960, 2, 4>,
    <58, 8515, 1123980, 2247960, 2, 4>,
    <59, 9453, 1285608, 2571216, 2, 4>,
    <59, 9316, 1285608, 2571216, 2, 4>,
    <60, 9730, 1342740, 2685480, 2, 4>,
    <60, 9591, 1342740, 2685480, 2, 4>,
    <67, 8256, 2097024, 14679168, 7, 4>,
    <74, 8505, 3265920, 26127360, 8, 38>,
    <80, 9750, 4680000, 9360000, 2, 4>,
    <87, 8424, 6065280, 24261120, 4, 24>,
    <97, 9920, 9999360, 19998720, 2, 3>
];

// Each entry stores sid, n, T, H meet T, current largest G, and its H.
cohorts := [* *];
raw_action_count := 0;

for simple_id in [1 .. NumberOfSimpleGroups()] do
    simple_tuple, simple_order := SimpleGroupId(simple_id);
    if simple_order gt maximum_simple_order then break; end if;

    A := AutomorphismGroupSimpleGroup(simple_tuple);
    T := Socle(A);
    assert Order(T) eq simple_order;
    outer, quotient_map := quo<A | T>;

    for extension_record in Subgroups(outer) do
        G := extension_record`subgroup @@ quotient_map;
        normalizer_G := Normalizer(A, G);
        representatives := [];

        for maximal_record in MaximalSubgroups(G :
                                               IndexLimit := maximum_degree) do
            H := maximal_record`subgroup;
            if T subset H then continue; end if;
            n := Index(G, H);
            if n lt minimum_degree or n gt maximum_degree then continue; end if;
            if Order(H) gt n - 1 then continue; end if;

            duplicate := false;
            for old_H in representatives do
                if Index(G, old_H) eq n and
                   IsConjugate(normalizer_G, H, old_H) then
                    duplicate := true;
                    break;
                end if;
            end for;
            if duplicate then continue; end if;
            Append(~representatives, H);
            raw_action_count +:= 1;

            T_stabilizer := H meet T;
            match := 0;
            for index in [1 .. #cohorts] do
                old := cohorts[index];
                if old[1] eq simple_id and old[2] eq n and
                   IsConjugate(A, T_stabilizer, old[4]) then
                    match := index;
                    break;
                end if;
            end for;
            if match eq 0 then
                Append(~cohorts,
                       <simple_id, n, T, T_stabilizer, G, H>);
            elif Order(G) gt Order(cohorts[match][5]) then
                cohorts[match] :=
                    <simple_id, n, T, T_stabilizer, G, H>;
            end if;
        end for;
    end for;
end for;

assert raw_action_count eq 29 and #cohorts eq 14;

global_histogram := AssociativeArray(Integers());
total_candidates := 0;
for cohort in cohorts do
    simple_id := cohort[1];
    n := cohort[2];
    T := cohort[3];
    G := cohort[5];
    H := cohort[6];
    action_map, N := CosetAction(G, H);
    assert Order(Kernel(action_map)) eq 1 and Degree(N) eq n and
           IsPrimitive(N);
    S := Socle(N);
    assert Order(S) eq Order(T);
    socle_primitive := IsPrimitive(S);

    expected_matches := [row : row in expected |
        row[1] eq simple_id and row[2] eq n];
    assert #expected_matches eq 1;
    expected_row := expected_matches[1];
    assert Order(S) eq expected_row[3] and
           Order(N) eq expected_row[4];

    component_outer, component_map := quo<N | S>;
    assert Order(component_outer) eq expected_row[5];
    small_quotient, base_inclusions, top_inclusion, top_projection :=
        WreathProduct(component_outer, Sym(2));
    socle_order := Order(S)^2;
    degree := n^2;
    maximum_quotient_order := degree * (degree - 1) div socle_order;
    allowed_orders := [order : order in Divisors(Order(small_quotient)) |
        order le maximum_quotient_order and order mod 2 eq 0];

    candidate_count := 0;
    subgroup_classes_checked := 0;
    top_transitive_count := 0;
    selected_subgroups := [sub<small_quotient | > : index in [1 .. expected_row[6]]];
    selected_orders := [0 : index in [1 .. expected_row[6]]];
    selected_horders := [0 : index in [1 .. expected_row[6]]];
    selected_checks := [0 : index in [1 .. expected_row[6]]];
    local_histogram := AssociativeArray(Integers());
    for allowed_order in allowed_orders do
        for subgroup_record in Subgroups(small_quotient :
                                          OrderEqual := allowed_order) do
            subgroup_classes_checked +:= 1;
            B := subgroup_record`subgroup;
            if not IsTransitive(B @ top_projection) then continue; end if;
            top_transitive_count +:= 1;
            if not socle_primitive then
                first_block := {1 .. Degree(component_outer)};
                block_stabilizer := Stabilizer(B, first_block);
                component_generators := [component_outer |
                    component_outer![j^g : j in [1 .. Degree(component_outer)]]
                    : g in Generators(block_stabilizer)];
                coordinate_outer := sub<component_outer | component_generators>;
                coordinate_component := coordinate_outer @@ component_map;
                if not IsPrimitive(coordinate_component) then continue; end if;
            end if;
            candidate_count +:= 1;
            assert candidate_count le expected_row[6];
            selected_subgroups[candidate_count] := B;
            selected_orders[candidate_count] := allowed_order;
            selected_horders[candidate_count] :=
                Order(S)^2 * allowed_order div degree;
            selected_checks[candidate_count] := subgroup_classes_checked;
            if not IsDefined(local_histogram, allowed_order) then
                local_histogram[allowed_order] := 0;
            end if;
            local_histogram[allowed_order] +:= 1;
            if not IsDefined(global_histogram, allowed_order) then
                global_histogram[allowed_order] := 0;
            end if;
            global_histogram[allowed_order] +:= 1;
        end for;
    end for;
    printf "PA_SEAM_COHORT|sid=%o|n=%o|T=%o|N=%o|outer=%o|degree=%o|socle_primitive=%o|candidates=%o|provisional=%o\n",
           simple_id, n, Order(S), Order(N), Order(component_outer), degree,
           socle_primitive, candidate_count, expected_row[6];
    printf "PA_SEAM_LOCAL_HIST|sid=%o|n=%o|outer_id=%o|",
           simple_id, n, IdentifyGroup(component_outer);
    for key in Sort(Setseq(Keys(local_histogram))) do
        printf "%o:%o,", key, local_histogram[key];
    end for;
    print "";
    assert candidate_count eq expected_row[6];

    layer_label := Sprintf("PA8_seam_sid%o_n%o_k2_d%o",
                           simple_id, n, degree);
    outer_elements := [element : element in component_outer];
    coordinate_tuples := [[component_outer | ]];
    for coordinate in [1 .. 2] do
        next_tuples := [];
        for prefix in coordinate_tuples do
            for outer_element in outer_elements do
                Append(~next_tuples, prefix cat [outer_element]);
            end for;
        end for;
        coordinate_tuples := next_tuples;
    end for;
    assert #coordinate_tuples eq Order(component_outer)^2;
    socle_generators := Setseq(Generators(S));
    socle_point_stabilizer := Stabilizer(S, 1);
    socle_stabilizer_generators := Setseq(Generators(socle_point_stabilizer));
    assert #socle_generators gt 0 and sub<S | socle_generators> eq S and
           #Orbit(sub<S | socle_generators>, 1) eq n and
           sub<S | socle_stabilizer_generators> eq socle_point_stabilizer;

    printf "PA8_SEAM_PACK_BEGIN|%o|sid=%o|checked=%o|top_transitive=%o\n",
           layer_label, simple_id, subgroup_classes_checked,
           top_transitive_count;
    print "PRIMITIVE_SAXL_PA_STRUCTURED_V1";
    print "layer";
    printf "component_degree %o\n", n;
    print "exponent 2";
    printf "socle_order %o\n", Order(S);
    printf "outer_order %o\n", Order(component_outer);
    printf "socle_gens %o\n", #socle_generators;
    for generator in socle_generators do
        PAXPrintPermutation("component_gen", generator, n);
    end for;
    printf "socle_stabilizer_gens %o\n", #socle_stabilizer_generators;
    for generator in socle_stabilizer_generators do
        PAXPrintPermutation("component_gen", generator, n);
    end for;
    printf "actions %o\n", candidate_count;
    for index in [1 .. candidate_count] do
        B := selected_subgroups[index];
        subgroup_generators := Setseq(Generators(B));
        assert #subgroup_generators gt 0 and
               sub<small_quotient | subgroup_generators> eq B;
        print "action";
        printf "label %o_sub%o\n", layer_label, index;
        printf "stabilizer_order %o\n", selected_horders[index];
        printf "quotient_order %o\n", selected_orders[index];
        printf "smallq_check_index %o\n", selected_checks[index];
        printf "quotient_gens %o\n", #subgroup_generators;
        for generator in subgroup_generators do
            coordinate_tuple, top_element := PAXDecompose(
                generator, small_quotient, top_projection, top_inclusion,
                base_inclusions, coordinate_tuples);
            print "qgen";
            PAXPrintTopPermutation(top_element, 2);
            print "component_gens 2";
            for coordinate in [1 .. 2] do
                outer_element := coordinate_tuple[coordinate];
                component_lift := outer_element eq Identity(component_outer)
                    select Identity(N) else outer_element @@ component_map;
                assert component_lift @ component_map eq outer_element;
                PAXPrintPermutation("component_gen", component_lift, n);
            end for;
            print "end_qgen";
        end for;
        print "end";
    end for;
    print "end_layer";
    printf "PA8_SEAM_PACK_END|%o\n", layer_label;
    total_candidates +:= candidate_count;
end for;

expected_histogram := AssociativeArray(Integers());
for pair in [<2,10>, <4,28>, <8,30>, <14,2>, <16,16>, <32,10>,
             <64,6>, <98,1>, <128,1>] do
    expected_histogram[pair[1]] := pair[2];
end for;

printf "PA_SEAM_GLOBAL_HIST|";
for key in Sort(Setseq(Keys(global_histogram))) do
    printf "%o:%o,", key, global_histogram[key];
end for;
print "";
printf "PA8_SEAM_CORRECTED_CENSUS|raw=29|cohorts=14|candidates=%o\n",
       total_candidates;
assert total_candidates eq 104;
assert Sort(Setseq(Keys(global_histogram))) eq
       Sort(Setseq(Keys(expected_histogram)));
assert &and[global_histogram[key] eq expected_histogram[key]
            : key in Keys(expected_histogram)];
print "PA8_SEAM_AUDIT_PASS|raw=29|cohorts=14|candidates=104";
