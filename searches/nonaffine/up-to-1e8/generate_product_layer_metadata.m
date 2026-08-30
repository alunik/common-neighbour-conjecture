///////////////////////////////////////////////////////////////////////////
// One exact PA component/layer metadata task.
//
// The base-two order sieve is applied in the small quotient W/S before a
// candidate preimage is tested for primitivity.  No point-stabiliser orbits
// are enumerated here.  This isolates expensive degree-d graph work behind
// a complete, auditable survivor ledger.
//
// Optional globals:
//   PAXMode               ("metadata" or "engine", default "metadata")
//   PAXEngineMode         (boolean/numeric CLI-safe alias for engine mode)
//   PAXSelectedCandidate  (required in engine mode)
//   PAXBatchEngineMode    (emit every candidate in one engine stream)
//   PAXExpectedCandidates (required exact count in batch engine mode)
///////////////////////////////////////////////////////////////////////////

require assigned PAXComponentDegree and assigned PAXComponentId and
        assigned PAXExponent and assigned PAXExpectedSocleOrder and
        assigned PAXExpectedOuterOrder and assigned PAXExpectedMetadataMethod:
    "set component degree, ID, exponent, and certified inventory facts";
if Type(PAXComponentDegree) eq MonStgElt then
    PAXComponentDegree := StringToInteger(PAXComponentDegree);
end if;
if Type(PAXComponentId) eq MonStgElt then
    PAXComponentId := StringToInteger(PAXComponentId);
end if;
if Type(PAXExponent) eq MonStgElt then
    PAXExponent := StringToInteger(PAXExponent);
end if;
if Type(PAXExpectedSocleOrder) eq MonStgElt then
    PAXExpectedSocleOrder := StringToInteger(PAXExpectedSocleOrder);
end if;
if Type(PAXExpectedOuterOrder) eq MonStgElt then
    PAXExpectedOuterOrder := StringToInteger(PAXExpectedOuterOrder);
end if;
if not assigned PAXMinDegree then PAXMinDegree := 10000001; end if;
if not assigned PAXMaxDegree then PAXMaxDegree := 100000000; end if;
if Type(PAXMinDegree) eq MonStgElt then
    PAXMinDegree := StringToInteger(PAXMinDegree);
end if;
if Type(PAXMaxDegree) eq MonStgElt then
    PAXMaxDegree := StringToInteger(PAXMaxDegree);
end if;
if not assigned PAXMode then PAXMode := "metadata"; end if;
if assigned PAXEngineMode then
    if Type(PAXEngineMode) eq MonStgElt then
        PAXEngineMode := PAXEngineMode in {"1", "true", "True"};
    end if;
    if Type(PAXEngineMode) eq RngIntElt then
        PAXEngineMode := PAXEngineMode ne 0;
    end if;
    if PAXEngineMode then PAXMode := "engine"; end if;
end if;
if not assigned PAXBatchEngineMode then PAXBatchEngineMode := false; end if;
if Type(PAXBatchEngineMode) eq MonStgElt then
    PAXBatchEngineMode := PAXBatchEngineMode in {"1", "true", "True"};
end if;
if Type(PAXBatchEngineMode) eq RngIntElt then
    PAXBatchEngineMode := PAXBatchEngineMode ne 0;
end if;
if PAXBatchEngineMode then PAXMode := "engine"; end if;
if not assigned PAXSelectedCandidate then
    PAXSelectedCandidate := 0;
elif Type(PAXSelectedCandidate) eq MonStgElt then
    PAXSelectedCandidate := StringToInteger(PAXSelectedCandidate);
end if;
if not assigned PAXExpectedCandidates then
    PAXExpectedCandidates := 0;
elif Type(PAXExpectedCandidates) eq MonStgElt then
    PAXExpectedCandidates := StringToInteger(PAXExpectedCandidates);
end if;
require 5 le PAXComponentDegree and
        PAXComponentDegree le PrimitiveGroupDatabaseLimit():
    "component degree outside primitive catalogue";
require PAXComponentId ge 1 and PAXExponent ge 2: "bad component/layer";
require PAXExpectedSocleOrder gt 0 and PAXExpectedOuterOrder gt 0 and
        PAXExpectedMetadataMethod in {
            "quotient_first_primitive_socle",
            "quotient_first_component_and_top"
        }: "bad certified inventory facts";
require PAXMode eq "metadata":
    "PA8 action export uses the separate compact structured-pack source";
SetColumns(0);

PAXPrintPermutation := procedure(g, d)
    width := 1;
    capacity := 94;
    while capacity lt d do
        width +:= 1;
        capacity *:= 94;
    end while;
    printf "packed_gen %o ", width;
    block_size := 4096;
    for first in [1 .. d by block_size] do
        pieces := [""];
        for point in [first .. Min(d, first + block_size - 1)] do
            value := point ^ g - 1;
            encoded := "";
            for digit_index in [1 .. width] do
                encoded cat:= CodeToString(33 + (value mod 94));
                value := value div 94;
            end for;
            assert value eq 0;
            Append(~pieces, encoded);
        end for;
        printf "%o", &cat pieces;
    end for;
    print "";
end procedure;

PAXTransitiveGenerators := function(P)
    d := Degree(P);
    for length in [2, 3] do
        for attempt in [1 .. 3] do
            trial := [Random(P) : index in [1 .. length]];
            if #Orbit(sub<P | trial>, 1) eq d then return trial; end if;
        end for;
    end for;
    chosen := [P | ];
    for g in Setseq(Generators(P)) do
        Append(~chosen, g);
        if #Orbit(sub<P | chosen>, 1) eq d then return chosen; end if;
    end for;
    assert false;
    return chosen;
end function;

PAXEmitEngine := procedure(G, stabilizer_order, label)
    d := Degree(G);
    H := Stabilizer(G, 1);
    assert Order(H) eq stabilizer_order;
    hgens := Setseq(Generators(H));
    assert #hgens gt 0;
    gens := PAXTransitiveGenerators(G);
    assert #Orbit(sub<G | gens>, 1) eq d;

    print "action";
    printf "label %o\n", label;
    printf "degree %o\n", d;
    printf "stabilizer_order %o\n", stabilizer_order;
    print "classification compute";
    print "regular_orbits 0";
    print "regular_count 0";
    printf "hgens %o\n", #hgens;
    for g in hgens do PAXPrintPermutation(g, d); end for;
    printf "gens %o\n", #gens;
    for g in gens do PAXPrintPermutation(g, d); end for;
    print "end";
end procedure;

n := PAXComponentDegree;
k := PAXExponent;
degree := n^k;
require PAXMinDegree le degree and degree le PAXMaxDegree:
    "layer degree outside requested interval";
component_normalizer := PrimitiveGroup(n, PAXComponentId);
require IsPrimitive(component_normalizer): "component is not primitive";
component_socle := Socle(component_normalizer);
socle_order := Order(component_socle);
require socle_order le n * (n - 1):
    "component cohort fails the base-two order sieve";
outer_order := Order(component_normalizer) div socle_order;
require socle_order eq PAXExpectedSocleOrder and
        outer_order eq PAXExpectedOuterOrder:
    "component disagrees with certified inventory";

// Recheck that the selected catalogue member is the largest member in its
// Sym(n)-socle cohort; this binds the inventory ID to N_Sym(n)(T).
component_outer, component_to_outer :=
    quo<component_normalizer | component_socle>;
require Order(component_outer) eq outer_order:
    "component quotient order disagrees with inventory";
small_quotient, base_inclusions, top_inclusion, top_projection :=
    WreathProduct(component_outer, Sym(k));
quotient_order := outer_order^k * Factorial(k);
require Order(small_quotient) eq quotient_order:
    "small PA quotient has the wrong order";
socle_primitive := IsPrimitive(component_socle);
metadata_method := socle_primitive select
    "quotient_first_primitive_socle" else
    "quotient_first_component_and_top";
require metadata_method eq PAXExpectedMetadataMethod:
    "component socle primitivity disagrees with inventory";
maximum_quotient_order := degree * (degree - 1) div (socle_order^k);
allowed_orders := [order : order in Divisors(quotient_order) |
                           order le maximum_quotient_order and
                           order mod k eq 0];

task_label := Sprintf("PA8_comp%o_%o_k%o_d%o",
                      n, PAXComponentId, k, degree);
printf "PA_LAYER|%o|%o|%o|%o|%o|%o|%o|%o|%o\n",
       task_label, n, PAXComponentId, k, degree, socle_order,
       outer_order, quotient_order, maximum_quotient_order;

candidate_count := 0;
subgroup_classes_checked := 0;
top_transitive_count := 0;
for allowed_order in allowed_orders do
    for subgroup_record in Subgroups(small_quotient :
                                      OrderEqual := allowed_order) do
        subgroup_classes_checked +:= 1;
        quotient_subgroup := subgroup_record`subgroup;
        if not IsTransitive(quotient_subgroup @ top_projection) then
            continue;
        end if;
        top_transitive_count +:= 1;
        if not socle_primitive then
            first_block := {1 .. Degree(component_outer)};
            block_stabilizer := Stabilizer(quotient_subgroup, first_block);
            component_generators := [component_outer |
                component_outer![j^g : j in [1 .. Degree(component_outer)]]
                : g in Generators(block_stabilizer)];
            coordinate_outer := sub<component_outer | component_generators>;
            coordinate_component := coordinate_outer @@ component_to_outer;
            if not IsPrimitive(coordinate_component) then continue; end if;
        end if;
        candidate_count +:= 1;
        stabilizer_order := socle_order^k * allowed_order div degree;
        assert stabilizer_order le degree - 1;
        label := Sprintf("%o_sub%o", task_label, candidate_count);
        full_quotient := allowed_order eq quotient_order;
        printf "PA_CANDIDATE|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o\n",
               label, task_label, n, PAXComponentId, k, degree,
               socle_order, allowed_order, stabilizer_order,
               full_quotient, subgroup_classes_checked;
    end for;
end for;

printf "PA_LAYER_COMPLETE|%o|allowed_orders=%o|subgroup_classes_checked=%o|candidates=%o\n",
       task_label, #allowed_orders, subgroup_classes_checked,
       candidate_count;
fprintf "/dev/stderr", "PAX8_METADATA_METHOD|%o|%o|top_transitive=%o\n",
        task_label, metadata_method, top_transitive_count;
quit;
