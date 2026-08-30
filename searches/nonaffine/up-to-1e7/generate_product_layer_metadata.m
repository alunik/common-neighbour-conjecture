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
//   PAXExpectedLabels, PAXExpectedAllowedOrders, PAXExpectedHOrders
//     (required exact manifest arrays in batch engine mode)
//   PAXChunkFirst, PAXChunkLast
//     (required canonical inclusive interval in batch engine mode)
//   PAXDescriptorMode (emit exact quotient-subgroup generator descriptors)
//   PAXDescriptorEngineMode (reconstruct and emit only pinned descriptors)
//   PAXDescriptorSequences, PAXDescriptorAllowedOrders, PAXDescriptorHOrders,
//   PAXDescriptorGeneratorCounts, PAXDescriptorGeneratorImages
//     (selected canonical rows; generator images use # / ; / , delimiters)
//   PAXExpectedQDegree, PAXExpectedQOrder, PAXExpectedQGeneratorCount,
//   PAXExpectedQGeneratorImages (exact embedded-quotient fingerprint)
//   PAXCompressEngineGenerators, PAXCompressionSeed,
//   PAXCompressionAttempts, PAXCompressionMaxSize
//     (descriptor-engine-only exact bounded generator compression)
//   PAXCertifiedInventoryMode (opt-in boolean; default false)
//   PAXExpectedSocleOrder, PAXExpectedOuterOrder
//     (required positive exact inventory facts in certified mode)
//   PAXExpectedMetadataMethod (required exact adaptive method in certified
//     metadata and descriptor modes)
//   PAXExpectedSubgroupClassesChecked, PAXExpectedTopTransitiveCount
//     (required exact terminal metadata summary in descriptor mode)
///////////////////////////////////////////////////////////////////////////

require assigned PAXComponentDegree and assigned PAXComponentId and
        assigned PAXExponent: "set component degree, ID, and exponent";
if Type(PAXComponentDegree) eq MonStgElt then
    PAXComponentDegree := StringToInteger(PAXComponentDegree);
end if;
if Type(PAXComponentId) eq MonStgElt then
    PAXComponentId := StringToInteger(PAXComponentId);
end if;
if Type(PAXExponent) eq MonStgElt then
    PAXExponent := StringToInteger(PAXExponent);
end if;
if not assigned PAXMinDegree then PAXMinDegree := 1000001; end if;
if not assigned PAXMaxDegree then PAXMaxDegree := 10000000; end if;
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
if not assigned PAXDescriptorMode then PAXDescriptorMode := false; end if;
if not assigned PAXDescriptorEngineMode then PAXDescriptorEngineMode := false; end if;
if Type(PAXDescriptorMode) eq MonStgElt then
    PAXDescriptorMode := PAXDescriptorMode in {"1", "true", "True"};
end if;
if Type(PAXDescriptorEngineMode) eq MonStgElt then
    PAXDescriptorEngineMode := PAXDescriptorEngineMode in {"1", "true", "True"};
end if;
if Type(PAXDescriptorMode) eq RngIntElt then PAXDescriptorMode := PAXDescriptorMode ne 0; end if;
if Type(PAXDescriptorEngineMode) eq RngIntElt then PAXDescriptorEngineMode := PAXDescriptorEngineMode ne 0; end if;
require not (PAXDescriptorMode and PAXDescriptorEngineMode):
    "descriptor producer and consumer modes are exclusive";
if PAXDescriptorMode then PAXMode := "descriptor"; end if;
if PAXDescriptorEngineMode then PAXMode := "engine"; end if;
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
if not assigned PAXChunkFirst then PAXChunkFirst := 0; end if;
if not assigned PAXChunkLast then PAXChunkLast := 0; end if;
if Type(PAXChunkFirst) eq MonStgElt then
    PAXChunkFirst := StringToInteger(PAXChunkFirst);
end if;
if Type(PAXChunkLast) eq MonStgElt then
    PAXChunkLast := StringToInteger(PAXChunkLast);
end if;
if not assigned PAXCertifiedInventoryMode then
    PAXCertifiedInventoryMode := false;
end if;
if Type(PAXCertifiedInventoryMode) eq MonStgElt then
    PAXCertifiedInventoryMode :=
        PAXCertifiedInventoryMode in {"1", "true", "True"};
end if;
if Type(PAXCertifiedInventoryMode) eq RngIntElt then
    PAXCertifiedInventoryMode := PAXCertifiedInventoryMode ne 0;
end if;
if not assigned PAXExpectedSocleOrder then PAXExpectedSocleOrder := 0; end if;
if not assigned PAXExpectedOuterOrder then PAXExpectedOuterOrder := 0; end if;
if not assigned PAXExpectedMetadataMethod then PAXExpectedMetadataMethod := ""; end if;
if not assigned PAXExpectedSubgroupClassesChecked then PAXExpectedSubgroupClassesChecked := 0; end if;
if not assigned PAXExpectedTopTransitiveCount then PAXExpectedTopTransitiveCount := 0; end if;
if Type(PAXExpectedSocleOrder) eq MonStgElt then
    PAXExpectedSocleOrder := StringToInteger(PAXExpectedSocleOrder);
end if;
if Type(PAXExpectedOuterOrder) eq MonStgElt then
    PAXExpectedOuterOrder := StringToInteger(PAXExpectedOuterOrder);
end if;
if Type(PAXExpectedSubgroupClassesChecked) eq MonStgElt then
    PAXExpectedSubgroupClassesChecked := StringToInteger(PAXExpectedSubgroupClassesChecked);
end if;
if Type(PAXExpectedTopTransitiveCount) eq MonStgElt then
    PAXExpectedTopTransitiveCount := StringToInteger(PAXExpectedTopTransitiveCount);
end if;
if PAXCertifiedInventoryMode then
    require PAXExpectedSocleOrder gt 0 and PAXExpectedOuterOrder gt 0:
        "certified inventory mode needs positive socle and outer orders";
    if PAXMode in {"metadata", "descriptor"} then
        require PAXExpectedMetadataMethod in
                {"quotient_first_primitive_socle", "quotient_first_component_and_top"}:
            "certified metadata needs an expected method";
    end if;
else
    require PAXExpectedSocleOrder eq 0 and PAXExpectedOuterOrder eq 0:
        "inventory facts require certified inventory mode";
    require PAXExpectedMetadataMethod eq "":
        "metadata method requires certified inventory mode";
end if;
if not assigned PAXExpectedLabels then PAXExpectedLabels := []; end if;
if not assigned PAXExpectedAllowedOrders then PAXExpectedAllowedOrders := []; end if;
if not assigned PAXExpectedHOrders then PAXExpectedHOrders := []; end if;
if not assigned PAXDescriptorSequences then PAXDescriptorSequences := []; end if;
if not assigned PAXDescriptorAllowedOrders then PAXDescriptorAllowedOrders := []; end if;
if not assigned PAXDescriptorHOrders then PAXDescriptorHOrders := []; end if;
if not assigned PAXDescriptorGeneratorCounts then PAXDescriptorGeneratorCounts := []; end if;
if not assigned PAXDescriptorGeneratorImages then PAXDescriptorGeneratorImages := ""; end if;
if not assigned PAXExpectedQDegree then PAXExpectedQDegree := 0; end if;
if not assigned PAXExpectedQOrder then PAXExpectedQOrder := 0; end if;
if not assigned PAXExpectedQGeneratorCount then PAXExpectedQGeneratorCount := 0; end if;
if not assigned PAXExpectedQGeneratorImages then PAXExpectedQGeneratorImages := ""; end if;
if not assigned PAXCompressEngineGenerators then PAXCompressEngineGenerators := false; end if;
if not assigned PAXCompressionSeed then PAXCompressionSeed := 0; end if;
if not assigned PAXCompressionAttempts then PAXCompressionAttempts := 0; end if;
if not assigned PAXCompressionMaxSize then PAXCompressionMaxSize := 0; end if;
if Type(PAXCompressEngineGenerators) eq MonStgElt then
    PAXCompressEngineGenerators := PAXCompressEngineGenerators in {"1", "true", "True"};
end if;
if Type(PAXCompressEngineGenerators) eq RngIntElt then
    PAXCompressEngineGenerators := PAXCompressEngineGenerators ne 0;
end if;
if Type(PAXCompressionSeed) eq MonStgElt then PAXCompressionSeed := StringToInteger(PAXCompressionSeed); end if;
if Type(PAXCompressionAttempts) eq MonStgElt then PAXCompressionAttempts := StringToInteger(PAXCompressionAttempts); end if;
if Type(PAXCompressionMaxSize) eq MonStgElt then PAXCompressionMaxSize := StringToInteger(PAXCompressionMaxSize); end if;
if Type(PAXExpectedQDegree) eq MonStgElt then PAXExpectedQDegree := StringToInteger(PAXExpectedQDegree); end if;
if Type(PAXExpectedQOrder) eq MonStgElt then PAXExpectedQOrder := StringToInteger(PAXExpectedQOrder); end if;
if PAXBatchEngineMode or PAXDescriptorMode then
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
end if;
if PAXDescriptorEngineMode then
    if Type(PAXExpectedLabels) eq MonStgElt then
        PAXExpectedLabels := Split(PAXExpectedLabels, ",");
    end if;
    if Type(PAXDescriptorSequences) eq MonStgElt then
        PAXDescriptorSequences := [StringToInteger(x) : x in Split(PAXDescriptorSequences, ",")];
    end if;
    if Type(PAXDescriptorAllowedOrders) eq MonStgElt then
        PAXDescriptorAllowedOrders := [StringToInteger(x) : x in Split(PAXDescriptorAllowedOrders, ",")];
    end if;
    if Type(PAXDescriptorHOrders) eq MonStgElt then
        PAXDescriptorHOrders := [StringToInteger(x) : x in Split(PAXDescriptorHOrders, ",")];
    end if;
    if Type(PAXDescriptorGeneratorCounts) eq MonStgElt then
        PAXDescriptorGeneratorCounts := [StringToInteger(x) : x in Split(PAXDescriptorGeneratorCounts, ",")];
    end if;
    require #PAXDescriptorSequences gt 0 and
            #PAXDescriptorSequences eq #PAXExpectedLabels and
            #PAXDescriptorSequences eq #PAXDescriptorAllowedOrders and
            #PAXDescriptorSequences eq #PAXDescriptorHOrders and
            #PAXDescriptorSequences eq #PAXDescriptorGeneratorCounts:
        "descriptor engine arrays have unequal lengths";
    if Type(PAXExpectedQGeneratorCount) eq MonStgElt then
        PAXExpectedQGeneratorCount := StringToInteger(PAXExpectedQGeneratorCount);
    end if;
    require PAXExpectedQDegree gt 0 and PAXExpectedQOrder gt 0 and
            PAXExpectedQGeneratorCount gt 0 and
            Type(PAXExpectedQGeneratorImages) eq MonStgElt and
            #PAXExpectedQGeneratorImages gt 0:
        "descriptor engine needs a positive quotient fingerprint";
    require PAXCompressEngineGenerators and PAXCompressionSeed gt 0 and
            PAXCompressionAttempts gt 0 and
            2 le PAXCompressionMaxSize and PAXCompressionMaxSize le 6:
        "descriptor engine needs the pinned exact generator-compression policy";
else
    require PAXExpectedQDegree eq 0 and PAXExpectedQOrder eq 0 and
            PAXExpectedQGeneratorCount eq 0 and
            PAXExpectedQGeneratorImages eq "":
        "quotient fingerprint requires descriptor engine mode";
    require not PAXCompressEngineGenerators and PAXCompressionSeed eq 0 and
            PAXCompressionAttempts eq 0 and PAXCompressionMaxSize eq 0:
        "generator compression is only accepted in descriptor engine mode";
end if;
require 5 le PAXComponentDegree and
        PAXComponentDegree le PrimitiveGroupDatabaseLimit():
    "component degree outside primitive catalogue";
require PAXComponentId ge 1 and PAXExponent ge 2: "bad component/layer";
require PAXMode in {"metadata", "engine", "descriptor"}: "bad PAXMode";
if PAXMode eq "engine" then
    if PAXDescriptorEngineMode then
        require PAXBatchEngineMode eq false and PAXSelectedCandidate eq 0 and
                PAXExpectedCandidates eq 0 and
                PAXChunkFirst eq 0 and PAXChunkLast eq 0 and
                Type(PAXExpectedLabels) eq SeqEnum:
            "descriptor engine mode conflicts with traversal selectors";
    elif PAXBatchEngineMode then
        require PAXSelectedCandidate eq 0 and PAXExpectedCandidates ge 1 and
                Type(PAXExpectedLabels) eq SeqEnum and
                Type(PAXExpectedAllowedOrders) eq SeqEnum and
                Type(PAXExpectedHOrders) eq SeqEnum and
                #PAXExpectedLabels eq PAXExpectedCandidates and
                #PAXExpectedAllowedOrders eq PAXExpectedCandidates and
                #PAXExpectedHOrders eq PAXExpectedCandidates and
                #{x : x in PAXExpectedLabels} eq PAXExpectedCandidates and
                &and[Type(x) eq MonStgElt and #x gt 0 : x in PAXExpectedLabels] and
                &and[Type(x) eq RngIntElt and x gt 0 : x in PAXExpectedAllowedOrders] and
                &and[Type(x) eq RngIntElt and x gt 1 : x in PAXExpectedHOrders]:
            "batch engine mode needs valid exact manifest arrays";
        require 1 le PAXChunkFirst and PAXChunkFirst le PAXChunkLast and
                PAXChunkLast le PAXExpectedCandidates:
            "batch engine mode needs a valid canonical chunk interval";
    else
        require PAXSelectedCandidate ge 1 and PAXExpectedCandidates eq 0:
            "single engine mode needs only positive PAXSelectedCandidate";
        require PAXChunkFirst eq 0 and PAXChunkLast eq 0:
            "chunk selectors require batch engine mode";
    end if;
elif PAXDescriptorMode then
    require PAXExpectedCandidates ge 1 and
            Type(PAXExpectedAllowedOrders) eq SeqEnum and
            Type(PAXExpectedHOrders) eq SeqEnum and
            #PAXExpectedAllowedOrders eq PAXExpectedCandidates and
            #PAXExpectedHOrders eq PAXExpectedCandidates and
            PAXExpectedSubgroupClassesChecked gt 0 and
            PAXExpectedTopTransitiveCount ge PAXExpectedCandidates and
            &and[Type(x) eq RngIntElt and x gt 0 : x in PAXExpectedAllowedOrders] and
            &and[Type(x) eq RngIntElt and x gt 1 : x in PAXExpectedHOrders]:
        "descriptor producer needs the exact candidate signature multiset";
else
    require PAXExpectedSubgroupClassesChecked eq 0 and
            PAXExpectedTopTransitiveCount eq 0:
        "terminal metadata summary is only accepted in descriptor mode";
end if;
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

PAXJoin := function(pieces, separator)
    if #pieces eq 0 then return ""; end if;
    answer := pieces[1];
    for index in [2 .. #pieces] do
        answer cat:= separator cat pieces[index];
    end for;
    return answer;
end function;

PAXPermutationCSV := function(g, d)
    return PAXJoin([Sprint(i^g) : i in [1 .. d]], ",");
end function;

// This is deliberately representation-specific.  Descriptor production and
// consumption must reconstruct byte-identical ordered quotient generators;
// degree and order alone do not identify the embedded quotient action.
PAXGroupGeneratorFingerprint := function(P)
    generator_strings :=
        [PAXPermutationCSV(g, Degree(P)) : g in Setseq(Generators(P))];
    assert #generator_strings gt 0;
    return #generator_strings, PAXJoin(generator_strings, ";");
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

PAXTransitiveGenerators := function(P)
    chosen := Setseq(Generators(P));
    assert #chosen gt 0 and sub<P | chosen> eq P and
           #Orbit(sub<P | chosen>, 1) eq Degree(P);
    return chosen;
end function;

PAXExactEngineGenerators := function(H)
    all_generators := Setseq(Generators(H));
    assert #all_generators gt 0;
    if not PAXCompressEngineGenerators then return all_generators; end if;
    SetSeed(PAXCompressionSeed);
    for generator_count in [2 .. PAXCompressionMaxSize] do
        for attempt in [1 .. PAXCompressionAttempts] do
            trial := [H | Random(H) : index in [1 .. generator_count]];
            trial_subgroup := sub<H | trial>;
            if Order(trial_subgroup) eq Order(H) then
                // Because trial_subgroup <= H and the groups are finite, equal
                // orders prove exact equality; this is not a heuristic gate.
                assert trial_subgroup eq H;
                return trial;
            end if;
        end for;
    end for;
    return all_generators;
end function;

PAXEmitEngine := procedure(G, stabilizer_order, label)
    d := Degree(G);
    H := Stabilizer(G, 1);
    assert Order(H) eq stabilizer_order;
    hgens := PAXExactEngineGenerators(H);
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
if PAXCertifiedInventoryMode then
    require socle_order eq PAXExpectedSocleOrder and
            outer_order eq PAXExpectedOuterOrder:
        "selected component disagrees with certified inventory facts";
end if;

// Recheck that the selected catalogue member is the largest member in its
// Sym(n)-socle cohort; this binds the inventory ID to N_Sym(n)(T).
// In opt-in certified mode the multi-root composer has independently replayed
// the exact promoted inventory payload and cohort proof, pinned its audit hash,
// and passed the exact catalogue ID/order facts.  Repeating the cohort scan is
// redundant and can be prohibitively slow; default mode retains the old check.
if not PAXCertifiedInventoryMode then
    for candidate in PrimitiveGroups(n : Filter := "AlmostSimple") do
        candidate_socle := Socle(candidate);
        if Order(candidate_socle) eq socle_order and
           (candidate_socle eq component_socle or
            IsConjugate(Sym(n), candidate_socle, component_socle)) then
            assert Order(candidate) le Order(component_normalizer);
        end if;
    end for;
end if;

task_label := Sprintf("PA7_comp%o_%o_k%o_d%o",
                      n, PAXComponentId, k, degree);

// Certified metadata is computed entirely in Q=(N/T) wr S_k.  For B<=Q, the
// preimage in N wr S_k is primitive in product action iff (i) its top image is
// transitive and (ii) its coordinate component is primitive.  When T itself
// is primitive, (ii) is automatic.  For novelty actions, recover the exact
// coordinate component A<=N/T from the stabilizer of the first imprimitive
// block and test the preimage A@@component_to_outer in the original degree n.
if PAXCertifiedInventoryMode and PAXMode eq "metadata" then
    socle_primitive := IsPrimitive(component_socle);
    expected_method := socle_primitive select
        "quotient_first_primitive_socle" else
        "quotient_first_component_and_top";
    require PAXExpectedMetadataMethod eq expected_method:
        "component socle primitivity disagrees with method manifest";
    component_outer, component_to_outer :=
        quo<component_normalizer | component_socle>;
    require Order(component_outer) eq outer_order:
        "component quotient order disagrees with certified outer order";
    small_quotient, base_inclusions, top_inclusion, top_projection :=
        WreathProduct(component_outer, Sym(k));
    quotient_order := outer_order^k * Factorial(k);
    require Order(small_quotient) eq quotient_order:
        "small PA quotient has the wrong order";
    maximum_quotient_order :=
        degree * (degree - 1) div (socle_order^k);
    allowed_orders := [order : order in Divisors(quotient_order) |
                               order le maximum_quotient_order and
                               order mod k eq 0];
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
                block_stabilizer := Stabilizer(quotient_subgroup,
                                               first_block);
                component_generators := [component_outer |
                    component_outer![j^g :
                        j in [1 .. Degree(component_outer)]]
                    : g in Generators(block_stabilizer)];
                coordinate_outer := sub<component_outer |
                                        component_generators>;
                coordinate_component :=
                    coordinate_outer @@ component_to_outer;
                if not IsPrimitive(coordinate_component) then continue; end if;
            end if;
            candidate_count +:= 1;
            stabilizer_order :=
                socle_order^k * allowed_order div degree;
            require stabilizer_order le degree - 1:
                "small-quotient candidate fails stabilizer order sieve";
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
    fprintf "/dev/stderr", "PAX_METADATA_METHOD|%o|%o|top_transitive=%o\n",
            task_label, expected_method, top_transitive_count;
    quit;
end if;

top := Sym(k);
W := PrimitiveWreathProduct(component_normalizer, top);
assert Degree(W) eq degree and IsPrimitive(W);
S := Socle(W);
assert Order(S) eq socle_order^k;
Q, quotient_map := quo<W | S>;
maximum_quotient_order := degree * (degree - 1) div Order(S);
allowed_orders := [order : order in Divisors(Order(Q)) |
                           order le maximum_quotient_order and
                           order mod k eq 0];
q_generator_count, q_generator_images := PAXGroupGeneratorFingerprint(Q);

if PAXDescriptorEngineMode then
    require Degree(Q) eq PAXExpectedQDegree and Order(Q) eq PAXExpectedQOrder and
            q_generator_count eq PAXExpectedQGeneratorCount and
            q_generator_images eq PAXExpectedQGeneratorImages:
        "reconstructed quotient disagrees with descriptor fingerprint";
    subgroup_strings := Split(PAXDescriptorGeneratorImages, "#");
    require #subgroup_strings eq #PAXDescriptorSequences:
        "descriptor subgroup string count mismatch";
    print "PRIMITIVE_SAXL_V1";
    for descriptor_index in [1 .. #PAXDescriptorSequences] do
        generator_strings := Split(subgroup_strings[descriptor_index], ";");
        require #generator_strings eq PAXDescriptorGeneratorCounts[descriptor_index] and
                #generator_strings gt 0:
            "descriptor generator count mismatch";
        quotient_generators := [Q | ];
        for generator_string in generator_strings do
            images := [StringToInteger(x) : x in Split(generator_string, ",")];
            require #images eq Degree(Q):
                "descriptor permutation has the wrong degree";
            ambient_generator := Sym(Degree(Q)) ! images;
            require ambient_generator in Q:
                "descriptor permutation is not in the reconstructed quotient";
            Append(~quotient_generators, Q ! ambient_generator);
        end for;
        quotient_subgroup := sub<Q | quotient_generators>;
        allowed_order := PAXDescriptorAllowedOrders[descriptor_index];
        stabilizer_order := PAXDescriptorHOrders[descriptor_index];
        require Order(quotient_subgroup) eq allowed_order and
                allowed_order in allowed_orders:
            "descriptor subgroup order is wrong";
        G := quotient_subgroup @@ quotient_map;
        require IsPrimitive(G) and Order(G) eq Order(S) * allowed_order and
                stabilizer_order eq Order(S) * allowed_order div degree and
                stabilizer_order le degree - 1:
            "descriptor preimage fails the exact action checks";
        sequence := PAXDescriptorSequences[descriptor_index];
        require PAXExpectedLabels[descriptor_index] eq
                Sprintf("%o_sub%o", task_label, sequence):
            "descriptor canonical label disagrees with sequence";
        PAXEmitEngine(G, stabilizer_order,
                      PAXExpectedLabels[descriptor_index]);
    end for;
    fprintf "/dev/stderr",
            "PAX_DESCRIPTOR_CHUNK_COMPLETE label=%o qdegree=%o qorder=%o first=%o last=%o emitted=%o generator_compression=exact_seeded_v1 seed=%o attempts=%o max_size=%o\n",
            task_label, Degree(Q), Order(Q), PAXDescriptorSequences[1],
            PAXDescriptorSequences[#PAXDescriptorSequences],
            #PAXDescriptorSequences, PAXCompressionSeed,
            PAXCompressionAttempts, PAXCompressionMaxSize;
    quit;
end if;

// The descriptor prepass is the only action-export stage that traverses
// subgroup classes.  It repeats the already-certified small-quotient census,
// constructs one exact isomorphism small_Q -> Q, and serializes the images in
// the fixed quotient permutation representation Q.  In particular, it never
// constructs or tests all large-degree preimages.  The terminal metadata
// labels are arbitrary slots within each exact (|B|,|H|) signature; assigning
// the smallest unused slot gives a proved global bijection while preserving
// every quotient action exactly once.
if PAXDescriptorMode then
    socle_primitive := IsPrimitive(component_socle);
    expected_method := socle_primitive select
        "quotient_first_primitive_socle" else
        "quotient_first_component_and_top";
    require PAXExpectedMetadataMethod eq expected_method:
        "component socle primitivity disagrees with descriptor method";
    component_outer, component_to_outer :=
        quo<component_normalizer | component_socle>;
    require Order(component_outer) eq outer_order:
        "component quotient order disagrees with certified outer order";
    small_quotient, base_inclusions, top_inclusion, top_projection :=
        WreathProduct(component_outer, Sym(k));
    require Order(small_quotient) eq Order(Q):
        "small and embedded quotient orders differ";

    // Build the map on the canonical base-coordinate and top generators.
    // A bare IsIsomorphic(small_quotient,Q) is deliberately not used: an
    // arbitrary abstract automorphism need not preserve the wreath base or
    // the coordinate-component primitivity criterion.
    structure_pairs := [];
    large_base_generators := [Q | ];
    for coordinate in [1 .. k] do
        for outer_generator in Generators(component_outer) do
            component_lift := outer_generator @@ component_to_outer;
            ambient_lift := PAXCoordinateLift(component_lift, n, k,
                                              coordinate);
            require ambient_lift in W:
                "coordinate lift is not in the primitive wreath product";
            quotient_lift := (W ! ambient_lift) @ quotient_map;
            Append(~structure_pairs,
                   <outer_generator @ base_inclusions[coordinate],
                    quotient_lift>);
            Append(~large_base_generators, quotient_lift);
        end for;
    end for;
    large_top_generators := [Q | ];
    for top_generator in Generators(top) do
        ambient_lift := PAXTopLift(top_generator, n, k);
        require ambient_lift in W:
            "top lift is not in the primitive wreath product";
        quotient_lift := (W ! ambient_lift) @ quotient_map;
        Append(~structure_pairs,
               <top_generator @ top_inclusion, quotient_lift>);
        Append(~large_top_generators, quotient_lift);
    end for;
    small_to_large := hom<small_quotient -> Q | structure_pairs>;
    small_base := Kernel(top_projection);
    mapped_base := sub<Q |
        [small_to_large(g) : g in Generators(small_base)]>;
    explicit_large_base := sub<Q | large_base_generators>;
    explicit_large_top := sub<Q | large_top_generators>;
    require #Kernel(small_to_large) eq 1 and
            Image(small_to_large) eq Q and
            mapped_base eq explicit_large_base and
            Order(explicit_large_base) eq outer_order^k and
            Order(explicit_large_top) eq Factorial(k) and
            Order(explicit_large_base meet explicit_large_top) eq 1 and
            sub<Q | Setseq(Generators(explicit_large_base)) cat
                    Setseq(Generators(explicit_large_top))> eq Q and
            &and[(pair[1] @ small_to_large) eq pair[2]
                 : pair in structure_pairs]:
        "structure-preserving small-to-embedded quotient map failed";

    descriptor_allowed := [0 : index in [1 .. PAXExpectedCandidates]];
    descriptor_horders := [0 : index in [1 .. PAXExpectedCandidates]];
    descriptor_generator_counts := [0 : index in [1 .. PAXExpectedCandidates]];
    descriptor_generator_images := ["" : index in [1 .. PAXExpectedCandidates]];
    descriptor_check_indices := [0 : index in [1 .. PAXExpectedCandidates]];
    subgroup_classes_checked := 0;
    top_transitive_count := 0;
    candidate_count := 0;
    for allowed_order in allowed_orders do
        for subgroup_record in Subgroups(small_quotient :
                                          OrderEqual := allowed_order) do
            subgroup_classes_checked +:= 1;
            small_subgroup := subgroup_record`subgroup;
            if not IsTransitive(small_subgroup @ top_projection) then
                continue;
            end if;
            top_transitive_count +:= 1;
            if not socle_primitive then
                first_block := {1 .. Degree(component_outer)};
                block_stabilizer := Stabilizer(small_subgroup, first_block);
                component_generators := [component_outer |
                    component_outer![j^g :
                        j in [1 .. Degree(component_outer)]]
                    : g in Generators(block_stabilizer)];
                coordinate_outer := sub<component_outer | component_generators>;
                coordinate_component := coordinate_outer @@ component_to_outer;
                if not IsPrimitive(coordinate_component) then continue; end if;
            end if;

            stabilizer_order := Order(S) * allowed_order div degree;
            require stabilizer_order le degree - 1:
                "descriptor candidate fails stabilizer order sieve";
            eligible := [index : index in [1 .. PAXExpectedCandidates] |
                descriptor_generator_counts[index] eq 0 and
                PAXExpectedAllowedOrders[index] eq allowed_order and
                PAXExpectedHOrders[index] eq stabilizer_order];
            require #eligible gt 0:
                "descriptor candidate signature multiset differs from metadata";
            canonical_index := eligible[1];
            small_generators := Setseq(Generators(small_subgroup));
            require #small_generators gt 0:
                "candidate small quotient descriptor has no generators";
            large_generators := [Q | small_to_large(g) : g in small_generators];
            large_subgroup := sub<Q | large_generators>;
            require Order(large_subgroup) eq allowed_order:
                "mapped quotient descriptor has the wrong order";
            descriptor_allowed[canonical_index] := allowed_order;
            descriptor_horders[canonical_index] := stabilizer_order;
            descriptor_generator_counts[canonical_index] := #large_generators;
            descriptor_generator_images[canonical_index] := PAXJoin(
                [PAXPermutationCSV(g, Degree(Q)) : g in large_generators], ";");
            descriptor_check_indices[canonical_index] := subgroup_classes_checked;
            candidate_count +:= 1;
        end for;
    end for;
    require candidate_count eq PAXExpectedCandidates and
            subgroup_classes_checked eq PAXExpectedSubgroupClassesChecked and
            top_transitive_count eq PAXExpectedTopTransitiveCount and
            &and[count gt 0 : count in descriptor_generator_counts]:
        "descriptor census differs from terminal metadata summary";
    printf "PAX_DESCRIPTOR_V2|%o|qdegree=%o|qorder=%o|qgenerator_count=%o|qgenerator_images=%o\n",
           task_label, Degree(Q), Order(Q), q_generator_count,
           q_generator_images;
    for sequence in [1 .. PAXExpectedCandidates] do
        printf "PAX_DESCRIPTOR|%o|%o|%o|%o|%o|smallq_check_index=%o\n",
               sequence, descriptor_allowed[sequence],
               descriptor_horders[sequence],
               descriptor_generator_counts[sequence],
               descriptor_generator_images[sequence],
               descriptor_check_indices[sequence];
    end for;
    printf "PAX_DESCRIPTOR_COMPLETE|%o|candidates=%o|subgroup_classes_checked=%o|top_transitive=%o\n",
           task_label, candidate_count, subgroup_classes_checked,
           top_transitive_count;
    quit;
end if;

if PAXMode eq "metadata" then
    printf "PA_LAYER|%o|%o|%o|%o|%o|%o|%o|%o|%o\n",
           task_label, n, PAXComponentId, k, degree, socle_order,
           outer_order, Order(Q),
           maximum_quotient_order;
else
    print "PRIMITIVE_SAXL_V1";
end if;
candidate_count := 0;
subgroup_classes_checked := 0;
selected_count := 0;
batch_assignment := [0 : i in [1 .. PAXExpectedCandidates]];
for allowed_order in allowed_orders do
    for subgroup_record in Subgroups(Q : OrderEqual := allowed_order) do
        subgroup_classes_checked +:= 1;
        quotient_subgroup := subgroup_record`subgroup;
        G := quotient_subgroup @@ quotient_map;
        if not IsPrimitive(G) then continue; end if;
        candidate_count +:= 1;
        stabilizer_order := Order(S) * allowed_order div degree;
        assert stabilizer_order le degree - 1;
        assert Order(G) eq Order(S) * allowed_order;
        label := Sprintf("%o_sub%o", task_label, candidate_count);
        full_quotient := allowed_order eq Order(Q);
        if PAXMode eq "metadata" then
            printf "PA_CANDIDATE|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o|%o\n",
                   label, task_label, n, PAXComponentId, k, degree,
                   socle_order, allowed_order, stabilizer_order,
                   full_quotient, subgroup_classes_checked;
        elif PAXBatchEngineMode then
            eligible := [expected_index : expected_index in [1 .. PAXExpectedCandidates] |
                batch_assignment[expected_index] eq 0 and
                PAXExpectedAllowedOrders[expected_index] eq allowed_order and
                PAXExpectedHOrders[expected_index] eq stabilizer_order];
            require #eligible gt 0:
                "batch PA candidate signature multiset differs from manifest";
            expected_index := eligible[1];
            batch_assignment[expected_index] := candidate_count;
            fprintf "/dev/stderr",
                    "PAX_BATCH_MAP|%o|%o|%o|%o|%o|%o\n",
                    candidate_count, expected_index,
                    PAXExpectedLabels[expected_index], label,
                    allowed_order, stabilizer_order;
            if PAXChunkFirst le expected_index and
               expected_index le PAXChunkLast then
                PAXEmitEngine(G, stabilizer_order,
                              PAXExpectedLabels[expected_index]);
                selected_count +:= 1;
            end if;
        elif candidate_count eq PAXSelectedCandidate then
            PAXEmitEngine(G, stabilizer_order, label);
            selected_count +:= 1;
        end if;
    end for;
end for;

if PAXMode eq "metadata" then
    printf "PA_LAYER_COMPLETE|%o|allowed_orders=%o|subgroup_classes_checked=%o|candidates=%o\n",
           task_label, #allowed_orders, subgroup_classes_checked,
           candidate_count;
elif PAXBatchEngineMode then
    require candidate_count eq PAXExpectedCandidates:
        "batch PA candidate count disagrees with metadata ledger";
    require selected_count eq PAXChunkLast - PAXChunkFirst + 1:
        "batch PA exporter did not emit the exact canonical chunk";

    require &and[index gt 0 : index in batch_assignment]:
        "manifest contains an unmatched PA candidate signature";
    fprintf "/dev/stderr", "PAX_BATCH_CHUNK_COMPLETE label=%o candidates=%o first=%o last=%o emitted=%o\n",
            task_label, candidate_count, PAXChunkFirst, PAXChunkLast,
            selected_count;
else
    require selected_count eq 1: "selected PA candidate was not found exactly once";
    fprintf "/dev/stderr", "PAX_SELECTED_COMPLETE label=%o\n",
            Sprintf("%o_sub%o", task_label, PAXSelectedCandidate);
end if;
quit;
