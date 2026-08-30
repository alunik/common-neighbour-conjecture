///////////////////////////////////////////////////////////////////////////
// Export an exact, compact product-action description for selected PA V3
// quotient descriptors.  This script deliberately never constructs the
// degree-n^k primitive wreath product.  Its output is expanded by the C++
// pa_structured_expander into the existing PRIMITIVE_SAXL_V1 contract.
//
// Required globals:
//   PAXComponentDegree, PAXComponentId, PAXExponent
//   PAXStructuredLabels
//   PAXStructuredAllowedOrders
//   PAXStructuredHOrders
//   PAXStructuredCheckIndices
//   PAXExpectedQDegree, PAXExpectedQOrder, PAXExpectedQGeneratorCount
//   PAXExpectedQGeneratorImages
//   PAXStructuredDescriptorGeneratorImages
//   PAXExpectedSocleOrder, PAXExpectedOuterOrder
///////////////////////////////////////////////////////////////////////////

require assigned PAXComponentDegree and assigned PAXComponentId and
        assigned PAXExponent and assigned PAXStructuredLabels and
        assigned PAXStructuredAllowedOrders and assigned PAXStructuredHOrders and
        assigned PAXStructuredCheckIndices and assigned PAXExpectedQDegree and
        assigned PAXExpectedQOrder and assigned PAXExpectedQGeneratorCount and
        assigned PAXExpectedQGeneratorImages and
        assigned PAXStructuredDescriptorGeneratorImages and
        assigned PAXExpectedSocleOrder and assigned PAXExpectedOuterOrder:
    "missing structured PA export globals";

PAXInteger := function(value)
    return Type(value) eq MonStgElt select StringToInteger(value) else value;
end function;

n := PAXInteger(PAXComponentDegree);
component_id := PAXInteger(PAXComponentId);
k := PAXInteger(PAXExponent);
expected_socle_order := PAXInteger(PAXExpectedSocleOrder);
expected_outer_order := PAXInteger(PAXExpectedOuterOrder);
expected_qdegree := PAXInteger(PAXExpectedQDegree);
expected_qorder := PAXInteger(PAXExpectedQOrder);
expected_qgenerator_count := PAXInteger(PAXExpectedQGeneratorCount);

if Type(PAXStructuredLabels) eq MonStgElt then
    PAXStructuredLabels := Split(PAXStructuredLabels, ",");
end if;
if Type(PAXStructuredAllowedOrders) eq MonStgElt then
    PAXStructuredAllowedOrders :=
        [StringToInteger(x) : x in Split(PAXStructuredAllowedOrders, ",")];
end if;
if Type(PAXStructuredHOrders) eq MonStgElt then
    PAXStructuredHOrders :=
        [StringToInteger(x) : x in Split(PAXStructuredHOrders, ",")];
end if;
if Type(PAXStructuredCheckIndices) eq MonStgElt then
    PAXStructuredCheckIndices :=
        [StringToInteger(x) : x in Split(PAXStructuredCheckIndices, ",")];
end if;

selected_count := #PAXStructuredLabels;
require selected_count gt 0 and
        #PAXStructuredAllowedOrders eq selected_count and
        #PAXStructuredHOrders eq selected_count and
        #PAXStructuredCheckIndices eq selected_count and
        #{x : x in PAXStructuredLabels} eq selected_count and
        #{x : x in PAXStructuredCheckIndices} eq selected_count and
        &and[x gt 0 : x in PAXStructuredAllowedOrders] and
        &and[x gt 1 : x in PAXStructuredHOrders] and
        &and[x gt 0 : x in PAXStructuredCheckIndices]:
    "invalid structured PA selection arrays";
require Type(PAXExpectedQGeneratorImages) eq MonStgElt and
        Type(PAXStructuredDescriptorGeneratorImages) eq MonStgElt and
        expected_qdegree gt 0 and expected_qorder gt 0 and
        expected_qgenerator_count gt 0:
    "invalid pinned quotient descriptor";

SetColumns(0);

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

PAXPermutationFromCSV := function(csv, ambient, expected_degree)
    images := [StringToInteger(x) : x in Split(csv, ",")];
    assert #images eq expected_degree;
    return ambient ! images;
end function;

PAXTupleToPoint := function(tuple, alphabet_size)
    return 1 + &+[(tuple[index] - 1) * alphabet_size^(index - 1)
                  : index in [1 .. #tuple]];
end function;

PAXPointToTuple := function(point, alphabet_size, tuple_length)
    quotient := point - 1;
    tuple := [];
    for index in [1 .. tuple_length] do
        Append(~tuple, (quotient mod alphabet_size) + 1);
        quotient div:= alphabet_size;
    end for;
    assert quotient eq 0;
    return tuple;
end function;

PAXCoordinateLift := function(component_element, alphabet_size,
                              tuple_length, coordinate)
    product_degree := alphabet_size^tuple_length;
    images := [];
    for point in [1 .. product_degree] do
        tuple := PAXPointToTuple(point, alphabet_size, tuple_length);
        tuple[coordinate] := tuple[coordinate]^component_element;
        Append(~images, PAXTupleToPoint(tuple, alphabet_size));
    end for;
    return Sym(product_degree) ! images;
end function;

PAXTopLift := function(top_element, alphabet_size, tuple_length)
    product_degree := alphabet_size^tuple_length;
    images := [];
    for point in [1 .. product_degree] do
        tuple := PAXPointToTuple(point, alphabet_size, tuple_length);
        image_tuple := [tuple[index^(top_element^-1)]
                        : index in [1 .. tuple_length]];
        Append(~images, PAXTupleToPoint(image_tuple, alphabet_size));
    end for;
    return Sym(product_degree) ! images;
end function;

component_normalizer := PrimitiveGroup(n, component_id);
require IsPrimitive(component_normalizer): "component is not primitive";
component_socle := Socle(component_normalizer);
socle_order := Order(component_socle);
outer_order := Order(component_normalizer) div socle_order;
require socle_order eq expected_socle_order and
        outer_order eq expected_outer_order:
    "component disagrees with certified inventory";
require IsTransitive(component_socle): "component socle is not transitive";

component_outer, component_to_outer :=
    quo<component_normalizer | component_socle>;
require Order(component_outer) eq outer_order:
    "component quotient has the wrong order";
small_quotient, base_inclusions, top_inclusion, top_projection :=
    WreathProduct(component_outer, Sym(k));
degree := n^k;
quotient_order := outer_order^k * Factorial(k);
require Order(small_quotient) eq quotient_order:
    "small PA quotient has the wrong order";

// The base has order at most 256 in the certified degree-10^7 inventory.
// Enumerating its coordinate tuples gives a representation-independent exact
// inverse to the canonical base inclusions.
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
    "failed to enumerate the small quotient base";

PAXDecompose := function(element)
    top_element := element @ top_projection;
    top_lift := top_element @ top_inclusion;
    base_element := element * top_lift^-1;
    assert base_element in Kernel(top_projection);
    for coordinate_tuple in coordinate_tuples do
        candidate := Identity(small_quotient);
        for coordinate in [1 .. k] do
            candidate *:= coordinate_tuple[coordinate] @
                          base_inclusions[coordinate];
        end for;
        if candidate eq base_element then
            assert candidate * top_lift eq element;
            return coordinate_tuple, top_element;
        end if;
    end for;
    error "could not decompose a small quotient element";
end function;

// Recover the exact pinned descriptor subgroups, rather than attempting to
// replay a subgroup-enumeration index.  Magma's class enumeration order is not
// a cross-session identity key.  The quotient representation and every
// descriptor generator image are already sealed.  Rebuild the original
// structure-preserving base/top map once for the whole layer and pull the
// pinned generators back through it.  This avoids both subgroup-enumeration
// drift and the normalizer ambiguity of a bare permutation conjugacy.
quotient_ambient := Sym(expected_qdegree);
expected_qgenerator_strings := Split(PAXExpectedQGeneratorImages, ";");
require #expected_qgenerator_strings eq expected_qgenerator_count:
    "pinned quotient generator count mismatch";
pinned_quotient_generators := [quotient_ambient |
    PAXPermutationFromCSV(csv, quotient_ambient, expected_qdegree)
    : csv in expected_qgenerator_strings];
pinned_quotient := sub<quotient_ambient | pinned_quotient_generators>;
require Order(pinned_quotient) eq expected_qorder and
        expected_qorder eq quotient_order:
    "pinned quotient has the wrong order";
// Rebuild the exact structure-preserving map once for this entire layer.  The
// expensive degree-n^k construction is therefore paid once per layer pack,
// never once per action.  This is the same map certified by the V3 descriptor
// prepass, including its base-coordinate and top-generator conventions.
top := Sym(k);
large_wreath := PrimitiveWreathProduct(component_normalizer, top);
large_socle := Socle(large_wreath);
large_quotient, large_quotient_map := quo<large_wreath | large_socle>;
require Order(large_socle) eq socle_order^k and
        Degree(large_quotient) eq expected_qdegree and
        Order(large_quotient) eq expected_qorder:
    "rebuilt large quotient has the wrong certified facts";
large_qgenerators := Setseq(Generators(large_quotient));
require #large_qgenerators eq #pinned_quotient_generators and
        &and[quotient_ambient ! large_qgenerators[index] eq
             pinned_quotient_generators[index]
             : index in [1 .. #large_qgenerators]]:
    "rebuilt large quotient disagrees with the pinned fingerprint";

structure_pairs := [];
large_base_generators := [large_quotient | ];
for coordinate in [1 .. k] do
    for outer_generator in Generators(component_outer) do
        component_lift := outer_generator @@ component_to_outer;
        ambient_lift := PAXCoordinateLift(component_lift, n, k,
                                          coordinate);
        quotient_lift := (large_wreath ! ambient_lift) @ large_quotient_map;
        Append(~structure_pairs,
               <outer_generator @ base_inclusions[coordinate],
                quotient_lift>);
        Append(~large_base_generators, quotient_lift);
    end for;
end for;
large_top_generators := [large_quotient | ];
for top_generator in Generators(top) do
    ambient_lift := PAXTopLift(top_generator, n, k);
    quotient_lift := (large_wreath ! ambient_lift) @ large_quotient_map;
    Append(~structure_pairs,
           <top_generator @ top_inclusion, quotient_lift>);
    Append(~large_top_generators, quotient_lift);
end for;
small_to_pinned := hom<small_quotient -> large_quotient | structure_pairs>;
small_base := Kernel(top_projection);
mapped_base := sub<large_quotient |
    [small_to_pinned(generator) : generator in Generators(small_base)]>;
explicit_large_base := sub<large_quotient | large_base_generators>;
explicit_large_top := sub<large_quotient | large_top_generators>;
require #Kernel(small_to_pinned) eq 1 and
        Image(small_to_pinned) eq large_quotient and
        mapped_base eq explicit_large_base and
        Order(explicit_large_base) eq outer_order^k and
        Order(explicit_large_top) eq Factorial(k) and
        Order(explicit_large_base meet explicit_large_top) eq 1:
    "exact structure-preserving small-to-pinned map failed";

descriptor_strings := Split(PAXStructuredDescriptorGeneratorImages, "#");
require #descriptor_strings eq selected_count:
    "pinned descriptor action count mismatch";
selected_subgroups := [];
for selected_index in [1 .. selected_count] do
    generator_strings := Split(descriptor_strings[selected_index], ";");
    pinned_generators := [quotient_ambient |
        PAXPermutationFromCSV(csv, quotient_ambient, expected_qdegree)
        : csv in generator_strings];
    require &and[generator in pinned_quotient
                 : generator in pinned_generators]:
        "descriptor generator is outside the pinned quotient";
    small_generators := [small_quotient |
        (large_quotient ! generator) @@ small_to_pinned
        : generator in pinned_generators];
    subgroup := sub<small_quotient | small_generators>;
    require Order(subgroup) eq PAXStructuredAllowedOrders[selected_index]:
        "pulled-back descriptor subgroup has the wrong order";
    top_subgroup := sub<Sym(k) |
        [generator @ top_projection : generator in Generators(subgroup)]>;
    require #Orbit(top_subgroup, 1) eq k:
        "pulled-back descriptor subgroup has intransitive top";
    expected_horder := socle_order^k * Order(subgroup) div degree;
    require expected_horder eq PAXStructuredHOrders[selected_index] and
            expected_horder le degree - 1:
        "pulled-back descriptor has the wrong stabilizer order";
    Append(~selected_subgroups, subgroup);
end for;

socle_generators := Setseq(Generators(component_socle));
socle_point_stabilizer := Stabilizer(component_socle, 1);
socle_stabilizer_generators :=
    Setseq(Generators(socle_point_stabilizer));
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
printf "actions %o\n", selected_count;

for selected_index in [1 .. selected_count] do
    subgroup := selected_subgroups[selected_index];
    subgroup_generators := Setseq(Generators(subgroup));
    require #subgroup_generators gt 0 and
            sub<small_quotient | subgroup_generators> eq subgroup:
        "small quotient generator certification failed";
    print "action";
    printf "label %o\n", PAXStructuredLabels[selected_index];
    printf "stabilizer_order %o\n", PAXStructuredHOrders[selected_index];
    printf "quotient_order %o\n", PAXStructuredAllowedOrders[selected_index];
    printf "smallq_check_index %o\n", PAXStructuredCheckIndices[selected_index];
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
    "PAX_STRUCTURED_COMPLETE n=%o id=%o k=%o actions=%o pinned_qdegree=%o pinned_qorder=%o\n",
    n, component_id, k, selected_count, expected_qdegree, expected_qorder;
