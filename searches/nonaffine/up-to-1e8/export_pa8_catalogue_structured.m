///////////////////////////////////////////////////////////////////////////
// Emit one compact exact PA8 structured pack from a certified catalogue
// layer.  Only the component action and Q=(N/T) wr S_k are constructed;
// the degree-n^k wreath product is never materialised in Magma.
///////////////////////////////////////////////////////////////////////////

require assigned PAXComponentDegree and assigned PAXComponentId and
        assigned PAXExponent and assigned PAXExpectedSocleOrder and
        assigned PAXExpectedOuterOrder and assigned PAXExpectedMetadataMethod and
        assigned PAXExpectedLabels and assigned PAXExpectedAllowedOrders and
        assigned PAXExpectedHOrders and assigned PAXExpectedCheckIndices and
        assigned PAXExpectedSubgroupClassesChecked and
        assigned PAXExpectedTopTransitiveCount:
    "missing certified PA8 structured-pack globals";

PAXInteger := function(value)
    return Type(value) eq MonStgElt select StringToInteger(value) else value;
end function;

n := PAXInteger(PAXComponentDegree);
component_id := PAXInteger(PAXComponentId);
k := PAXInteger(PAXExponent);
expected_socle_order := PAXInteger(PAXExpectedSocleOrder);
expected_outer_order := PAXInteger(PAXExpectedOuterOrder);
expected_checked := PAXInteger(PAXExpectedSubgroupClassesChecked);
expected_top_transitive := PAXInteger(PAXExpectedTopTransitiveCount);
if Type(PAXExpectedLabels) eq MonStgElt then
    PAXExpectedLabels := Split(PAXExpectedLabels, ",");
end if;
if Type(PAXExpectedAllowedOrders) eq MonStgElt then
    PAXExpectedAllowedOrders :=
        [StringToInteger(x) : x in Split(PAXExpectedAllowedOrders, ",")];
end if;
if Type(PAXExpectedHOrders) eq MonStgElt then
    PAXExpectedHOrders :=
        [StringToInteger(x) : x in Split(PAXExpectedHOrders, ",")];
end if;
if Type(PAXExpectedCheckIndices) eq MonStgElt then
    PAXExpectedCheckIndices :=
        [StringToInteger(x) : x in Split(PAXExpectedCheckIndices, ",")];
end if;

expected_count := #PAXExpectedLabels;
require 5 le n and n le PrimitiveGroupDatabaseLimit() and component_id ge 1 and
        k ge 2 and expected_socle_order gt 0 and expected_outer_order gt 0 and
        expected_count gt 0 and #PAXExpectedAllowedOrders eq expected_count and
        #PAXExpectedHOrders eq expected_count and
        #PAXExpectedCheckIndices eq expected_count and
        #{x : x in PAXExpectedLabels} eq expected_count and
        expected_checked gt 0 and expected_top_transitive ge expected_count:
    "invalid certified PA8 structured-pack globals";
require PAXExpectedMetadataMethod in {
            "quotient_first_primitive_socle",
            "quotient_first_component_and_top"
        }: "bad metadata method";
SetColumns(0);

PAXJoin := function(pieces, separator)
    if #pieces eq 0 then return ""; end if;
    answer := pieces[1];
    for index in [2 .. #pieces] do
        answer cat:= separator cat pieces[index];
    end for;
    return answer;
end function;

PAXPrintPermutation := procedure(prefix, g, d)
    width := 1;
    capacity := 94;
    while capacity lt d do
        width +:= 1;
        capacity *:= 94;
    end while;
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

component_normalizer := PrimitiveGroup(n, component_id);
require IsPrimitive(component_normalizer): "component is not primitive";
component_socle := Socle(component_normalizer);
socle_order := Order(component_socle);
outer_order := Order(component_normalizer) div socle_order;
require socle_order eq expected_socle_order and outer_order eq expected_outer_order:
    "component disagrees with certified inventory";
socle_primitive := IsPrimitive(component_socle);
metadata_method := socle_primitive select
    "quotient_first_primitive_socle" else
    "quotient_first_component_and_top";
require metadata_method eq PAXExpectedMetadataMethod:
    "component metadata method drift";
component_outer, component_to_outer :=
    quo<component_normalizer | component_socle>;
require Order(component_outer) eq outer_order:
    "component quotient order drift";
small_quotient, base_inclusions, top_inclusion, top_projection :=
    WreathProduct(component_outer, Sym(k));
degree := n^k;
quotient_order := outer_order^k * Factorial(k);
require Order(small_quotient) eq quotient_order:
    "small quotient order drift";
maximum_quotient_order := degree * (degree - 1) div (socle_order^k);
allowed_orders := [order : order in Divisors(quotient_order) |
    order le maximum_quotient_order and order mod k eq 0];
task_label := Sprintf("PA8_comp%o_%o_k%o_d%o", n, component_id, k, degree);

selected_subgroups := [sub<small_quotient | > : index in [1 .. expected_count]];
actual_check_indices := [];
candidate_count := 0;
subgroup_classes_checked := 0;
top_transitive_count := 0;
for allowed_order in allowed_orders do
    for subgroup_record in Subgroups(small_quotient : OrderEqual := allowed_order) do
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
            coordinate_component := coordinate_outer @@ component_to_outer;
            if not IsPrimitive(coordinate_component) then continue; end if;
        end if;
        candidate_count +:= 1;
        require candidate_count le expected_count:
            "pack census exceeds metadata ledger";
        stabilizer_order := socle_order^k * allowed_order div degree;
        require PAXExpectedLabels[candidate_count] eq
                    Sprintf("%o_sub%o", task_label, candidate_count) and
                PAXExpectedAllowedOrders[candidate_count] eq allowed_order and
                PAXExpectedHOrders[candidate_count] eq stabilizer_order:
            "pack candidate disagrees with metadata ledger";
        Append(~actual_check_indices, subgroup_classes_checked);
        selected_subgroups[candidate_count] := B;
    end for;
end for;
require candidate_count eq expected_count and
        subgroup_classes_checked eq expected_checked and
        top_transitive_count eq expected_top_transitive:
    "pack terminal census disagrees with metadata ledger";

// Subgroups() may permute conjugacy classes within a fixed order cohort when
// the same abstract quotient is reconstructed in a fresh Magma process.  The
// labels in the upstream ledger are sequential identifiers, not canonical
// subgroup fingerprints.  Preserve the exact fresh-to-ledger bijection while
// retaining the upstream labels and check indices in the structured pack.
mismatch_count := #[index : index in [1 .. expected_count] |
    actual_check_indices[index] ne PAXExpectedCheckIndices[index]];
printf "PA8_RENUMBERING_V2|layer=%o|mismatches=%o|actual=%o|expected=%o\n",
       task_label, mismatch_count,
       PAXJoin([IntegerToString(x) : x in actual_check_indices], ","),
       PAXJoin([IntegerToString(x) : x in PAXExpectedCheckIndices], ",");

outer_elements := [element : element in component_outer];
coordinate_tuples := [[component_outer | ]];
for coordinate in [1 .. k] do
    next_tuples := [];
    for prefix in coordinate_tuples do
        for outer_element in outer_elements do
            Append(~next_tuples, prefix cat [outer_element]);
        end for;
    end for;
    coordinate_tuples := next_tuples;
end for;
require #coordinate_tuples eq outer_order^k:
    "small quotient base enumeration failed";

PAXDecompose := function(element)
    top_element := element @ top_projection;
    top_lift := top_element @ top_inclusion;
    base_element := element * top_lift^-1;
    assert base_element in Kernel(top_projection);
    for coordinate_tuple in coordinate_tuples do
        candidate := Identity(small_quotient);
        for coordinate in [1 .. k] do
            candidate *:= coordinate_tuple[coordinate] @ base_inclusions[coordinate];
        end for;
        if candidate eq base_element then
            return coordinate_tuple, top_element;
        end if;
    end for;
    error "could not decompose small quotient generator";
end function;

socle_generators := Setseq(Generators(component_socle));
socle_point_stabilizer := Stabilizer(component_socle, 1);
socle_stabilizer_generators := Setseq(Generators(socle_point_stabilizer));
require #socle_generators gt 0 and
        sub<component_socle | socle_generators> eq component_socle and
        #Orbit(sub<component_socle | socle_generators>, 1) eq n and
        sub<component_socle | socle_stabilizer_generators> eq
            socle_point_stabilizer:
    "component generator certification failed";

print "PRIMITIVE_SAXL_PA_STRUCTURED_V1";
print "layer";
printf "component_degree %o\n", n;
printf "exponent %o\n", k;
printf "socle_order %o\n", socle_order;
printf "outer_order %o\n", outer_order;
printf "socle_gens %o\n", #socle_generators;
for generator in socle_generators do
    PAXPrintPermutation("component_gen", generator, n);
end for;
printf "socle_stabilizer_gens %o\n", #socle_stabilizer_generators;
for generator in socle_stabilizer_generators do
    PAXPrintPermutation("component_gen", generator, n);
end for;
printf "actions %o\n", expected_count;
for index in [1 .. expected_count] do
    B := selected_subgroups[index];
    subgroup_generators := Setseq(Generators(B));
    require #subgroup_generators gt 0 and
            sub<small_quotient | subgroup_generators> eq B:
        "small quotient generator certification failed";
    print "action";
    printf "label %o\n", PAXExpectedLabels[index];
    printf "stabilizer_order %o\n", PAXExpectedHOrders[index];
    printf "quotient_order %o\n", PAXExpectedAllowedOrders[index];
    printf "smallq_check_index %o\n", PAXExpectedCheckIndices[index];
    printf "quotient_gens %o\n", #subgroup_generators;
    for generator in subgroup_generators do
        coordinate_tuple, top_element := PAXDecompose(generator);
        print "qgen";
        PAXPrintTopPermutation(top_element, k);
        printf "component_gens %o\n", k;
        for coordinate in [1 .. k] do
            outer_element := coordinate_tuple[coordinate];
            component_lift := outer_element eq Identity(component_outer)
                select Identity(component_normalizer)
                else outer_element @@ component_to_outer;
            require component_lift @ component_to_outer eq outer_element:
                "component outer lift certification failed";
            PAXPrintPermutation("component_gen", component_lift, n);
        end for;
        print "end_qgen";
    end for;
    print "end";
end for;
print "end_layer";
fprintf "/dev/stderr",
    "PAX8_STRUCTURED_COMPLETE layer=%o actions=%o checked=%o top_transitive=%o\n",
    task_label, expected_count, subgroup_classes_checked, top_transitive_count;
quit;
