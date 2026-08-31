///////////////////////////////////////////////////////////////////////////
// Almost-simple frontier generator for 10^6 < degree <= 10^7.
//
// The metadata pass constructs no large coset actions.  It dispositions every
// action as an order obstruction, a published theorem exclusion, or an exact
// computation target.  A selected computation target can then be exported in
// PRIMITIVE_SAXL_V1 format for the independent C++ engine.
//
// Required globals:
//   ASXSimpleId
// Optional globals:
//   ASXMinDegree      (default 1000001)
//   ASXMaxDegree      (default 10000000)
//   ASXMode           ("metadata" or "engine", default "metadata")
//   ASXEngineMode     (boolean/numeric CLI-safe alias for engine mode)
//   ASXSelectedLabel  (one engine selector)
//   ASXSelectedExtension, ASXSelectedMaximal (CLI-safe engine selectors)
//   ASXBatchEngineMode (opt-in: export every compute action for this socle)
//   ASXExpectedActions (required positive batch action count)
///////////////////////////////////////////////////////////////////////////

require assigned ASXSimpleId: "set ASXSimpleId";
if Type(ASXSimpleId) eq MonStgElt then
    ASXSimpleId := StringToInteger(ASXSimpleId);
end if;
if not assigned ASXMinDegree then ASXMinDegree := 1000001; end if;
if not assigned ASXMaxDegree then ASXMaxDegree := 10000000; end if;
if not assigned ASXMode then ASXMode := "metadata"; end if;
if assigned ASXEngineMode then
    if Type(ASXEngineMode) eq MonStgElt then
        ASXEngineMode := ASXEngineMode in {"1", "true", "True"};
    end if;
    if Type(ASXEngineMode) eq RngIntElt then
        ASXEngineMode := ASXEngineMode ne 0;
    end if;
    if ASXEngineMode then ASXMode := "engine"; end if;
end if;
if not assigned ASXBatchEngineMode then ASXBatchEngineMode := false; end if;
if Type(ASXBatchEngineMode) eq MonStgElt then
    ASXBatchEngineMode := ASXBatchEngineMode in {"1", "true", "True"};
end if;
if Type(ASXBatchEngineMode) eq RngIntElt then
    ASXBatchEngineMode := ASXBatchEngineMode ne 0;
end if;
if ASXBatchEngineMode then ASXMode := "engine"; end if;
if not assigned ASXExpectedActions then ASXExpectedActions := 0; end if;
if Type(ASXExpectedActions) eq MonStgElt then
    ASXExpectedActions := StringToInteger(ASXExpectedActions);
end if;
if Type(ASXMinDegree) eq MonStgElt then
    ASXMinDegree := StringToInteger(ASXMinDegree);
end if;
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
require 1 le ASXSimpleId and ASXSimpleId le NumberOfSimpleGroups():
    "bad simple-group ID";
require 8192 le ASXMinDegree and ASXMinDegree le ASXMaxDegree and
        ASXMaxDegree le 10000000: "bad degree interval";
require (ASXMode eq "metadata") or (ASXMode eq "engine"): "bad ASXMode";
if not assigned ASXSelectedLabel then ASXSelectedLabel := ""; end if;
if not assigned ASXSelectedExtension then ASXSelectedExtension := 0; end if;
if not assigned ASXSelectedMaximal then ASXSelectedMaximal := 0; end if;
if Type(ASXSelectedExtension) eq MonStgElt then
    ASXSelectedExtension := StringToInteger(ASXSelectedExtension);
end if;
if Type(ASXSelectedMaximal) eq MonStgElt then
    ASXSelectedMaximal := StringToInteger(ASXSelectedMaximal);
end if;
if ASXMode eq "engine" then
    if ASXBatchEngineMode then
        require ASXSelectedLabel eq "" and ASXSelectedExtension eq 0 and
                ASXSelectedMaximal eq 0 and ASXExpectedActions ge 1:
            "batch engine mode needs only a positive ASXExpectedActions";
    else
        require ASXSelectedLabel ne "" or
                (ASXSelectedExtension ge 1 and ASXSelectedMaximal ge 1):
            "engine mode needs an action selector";
    end if;
end if;
SetColumns(0);

ASXPrintPermutation := procedure(g, d)
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

ASXTransitiveActionGenerators := function(P, G, action_map)
    d := Degree(P);
    for length in [2, 3] do
        for attempt in [1 .. 3] do
            source_trial := [Random(G) : i in [1 .. length]];
            trial := [g @ action_map : g in source_trial];
            if #Orbit(sub<P | trial>, 1) eq d then return trial; end if;
        end for;
    end for;
    chosen := [P | ];
    for g in Setseq(Generators(G)) do
        Append(~chosen, g @ action_map);
        if #Orbit(sub<P | chosen>, 1) eq d then return chosen; end if;
    end for;
    assert false;
    return chosen;
end function;

ASXDeterministicTransitiveActionGenerators := function(P, G, action_map)
    chosen := [P | ];
    for g in Setseq(Generators(G)) do
        Append(~chosen, g @ action_map);
        if #Orbit(sub<P | chosen>, 1) eq Degree(P) then return chosen; end if;
    end for;
    assert false;
    return chosen;
end function;

ASXEmitEngine := procedure(G, H, label)
    d := Index(G, H);
    h_order := Order(H);
    action_map, P := CosetAction(G, H);
    assert Degree(P) eq d and IsPrimitive(P);
    assert Order(P) eq Order(G) and Order(Kernel(action_map)) eq 1;
    HP := Stabilizer(P, 1);
    assert Order(HP) eq h_order;

    print "action";
    printf "label %o\n", label;
    printf "degree %o\n", d;
    printf "stabilizer_order %o\n", h_order;
    print "classification compute";
    print "regular_orbits 0";
    print "regular_count 0";
    hgens := Setseq(Generators(HP));
    assert #hgens gt 0;
    printf "hgens %o\n", #hgens;
    for g in hgens do ASXPrintPermutation(g, d); end for;
    if ASXBatchEngineMode then
        gens := ASXDeterministicTransitiveActionGenerators(P, G, action_map);
    else
        gens := ASXTransitiveActionGenerators(P, G, action_map);
    end if;
    assert #Orbit(sub<P | gens>, 1) eq d;
    printf "gens %o\n", #gens;
    for g in gens do ASXPrintPermutation(g, d); end for;
    print "end";
end procedure;

simple_tuple, simple_order := SimpleGroupId(ASXSimpleId);
simple_name := SimpleGroupName(ASXSimpleId);
require simple_order le ASXMaxDegree * (ASXMaxDegree - 1):
    "simple-group ID lies beyond the base-two order bound";

if ASXMode eq "engine" then print "PRIMITIVE_SAXL_V1"; end if;

// Burness--Huang, Algebraic Combinatorics 5 (2022), Theorem 4.22:
// every base-two almost-simple primitive group with socle PSL(2,q)
// satisfies the Burness--Giudici conjecture, including outer extensions.
if simple_tuple[1] eq 1 and simple_tuple[2] eq 1 then
    if ASXMode eq "metadata" then
        printf "SOCLE|%o|%o|theorem_bh_psl2|0|0|%o\n",
               ASXSimpleId, simple_order, simple_name;
        printf "SOCLE_COMPLETE|%o|theorem_bh_psl2|0|0\n", ASXSimpleId;
    end if;
    quit;
end if;

// Chen--Du, arXiv:2512.22461v2, Theorem 1.2 (preprint): every
// base-two primitive action with Suzuki or Ree socle satisfies the
// conjecture, including outer field extensions.
if (simple_tuple[1] eq 11 and simple_tuple[2] eq 2) or
   (simple_tuple[1] eq 14 and simple_tuple[2] eq 2) then
    if ASXMode eq "metadata" then
        printf "SOCLE|%o|%o|theorem_preprint_sz_ree|0|0|%o\n",
               ASXSimpleId, simple_order, simple_name;
        printf "SOCLE_COMPLETE|%o|theorem_preprint_sz_ree|0|0\n",
               ASXSimpleId;
    end if;
    quit;
end if;

// The exact degree-window audit in audit_alternating_window.m leaves only
// the A15 and S15 actions on partitions into five parts of size three.
// Morris--Spiga, J. Algebra 587 (2021), Theorems 1.1 and 1.2, give base
// size 3 for both, so no alternating-socle base-two target survives.
if simple_tuple[1] eq 17 then
    if ASXMode eq "metadata" then
        printf "SOCLE|%o|%o|theorem_published_alt_window|0|0|%o\n",
               ASXSimpleId, simple_order, simple_name;
        printf "SOCLE_COMPLETE|%o|theorem_published_alt_window|0|0\n",
               ASXSimpleId;
    end if;
    quit;
end if;

// Burness--Giudici, Theorem 6.1, covers precisely these sporadic socles,
// including all outer extensions.  The pinned CTblLib window audit shows
// that every order-eligible sporadic row here has one of these socles.
if simple_tuple[1] eq 18 and
   simple_tuple[2] in {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                       13, 14, 15, 16, 20, 21} then
    if ASXMode eq "metadata" then
        printf "SOCLE|%o|%o|theorem_bg_sporadic|0|0|%o\n",
               ASXSimpleId, simple_order, simple_name;
        printf "SOCLE_COMPLETE|%o|theorem_bg_sporadic|0|0\n", ASXSimpleId;
    end if;
    quit;
end if;

is_psu3 := simple_tuple[1] eq 10 and simple_tuple[2] eq 2;
psu_q := 0;
psu_d := 0;
psu_pso_degree := 0;
psu_pso_order := 0;
if is_psu3 then
    psu_q := simple_tuple[3];
    psu_d := Gcd(3, psu_q + 1);
    psu_pso_degree := psu_q^2 * (psu_q^3 + 1) div psu_d;
    psu_pso_order := psu_q * (psu_q^2 - 1);
end if;

// The two PSU(3,q) preprints together leave, in this exact degree window,
// only the outer PSO(3,q) actions at q=19,23,25,29.  Route every other
// PSU3 socle before constructing Aut(T) or enumerating maximal subgroups.
// This is a preprint-backed window disposition, not a published exclusion.
if is_psu3 and not psu_q in {19, 23, 25, 29} then
    if ASXMode eq "metadata" then
        printf "SOCLE|%o|%o|theorem_preprint_psu3_window_complete|0|0|%o\n",
               ASXSimpleId, simple_order, simple_name;
        printf "SOCLE_COMPLETE|%o|theorem_preprint_psu3_window_complete|0|0\n",
               ASXSimpleId;
    end if;
    quit;
end if;

A := AutomorphismGroupSimpleGroup(simple_tuple);
T := Socle(A);
assert Order(T) eq simple_order;
Q, quotient_map := quo<A | T>;
extension_records := Subgroups(Q);
action_count := 0;
compute_count := 0;
selected_count := 0;
selected_label := "";
extension_number := 0;
batch_records := [* *];

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
        if d lt ASXMinDegree or d gt ASXMaxDegree then continue; end if;

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
        action_count +:= 1;
        label := Sprintf("AS7_sid%o_ext%o_max%o_d%o",
                         ASXSimpleId, extension_number, maximal_number, d);
        h_order := Order(H);
        psu_pso_lane := is_psu3 and d eq psu_pso_degree and
                        Order(H meet T) eq psu_pso_order;
        if h_order gt d - 1 then
            disposition := "order_obstruction";
        elif IsSolvable(H) then
            // Burness--Huang, Algebraic Combinatorics 5 (2022), Theorem 1.1:
            // every primitive base-two action
            // with soluble point stabiliser satisfies the conjecture.
            disposition := "theorem_solvable_stabilizer";
        elif is_psu3 and not psu_pso_lane then
            // Chen--Du--Li, arXiv:2512.22456v2, Theorem 1.2 (preprint).
            disposition := "theorem_preprint_psu3_nonpso";
        elif psu_pso_lane and Order(G) eq Order(T) then
            // Chen--Du--Li, arXiv:2512.22459v2, Theorem 1.2(i) (preprint).
            disposition := "theorem_preprint_psu3_simple_pso";
        elif psu_pso_lane and
             ((psu_d eq 1 and psu_q gt 17^2) or
              (psu_d eq 3 and psu_q gt 45^2)) then
            // Theorem 1.2(ii).  The paper states >= at the lower endpoint
            // but Theorem 1.3 lists <= as residual, so equality is retained.
            disposition := "theorem_preprint_psu3_large_outer_pso";
        else
            disposition := "compute";
            compute_count +:= 1;
        end if;

        if ASXMode eq "metadata" then
            printf "ACTION|%o|%o|%o|%o|%o|%o|%o|%o\n",
                   label, ASXSimpleId, extension_number, maximal_number,
                   d, Order(G), h_order, disposition;
        elif ASXBatchEngineMode then
            if disposition eq "compute" then
                Append(~batch_records, <G, H, label>);
                selected_count +:= 1;
            end if;
        elif label eq ASXSelectedLabel or
             (extension_number eq ASXSelectedExtension and
              maximal_number eq ASXSelectedMaximal) then
            require disposition eq "compute":
                "selected action is already theorem/order disposed";
            ASXEmitEngine(G, H, label);
            selected_count +:= 1;
            selected_label := label;
        end if;
    end for;
end for;

if ASXMode eq "metadata" then
    printf "SOCLE|%o|%o|enumerated|%o|%o|%o\n",
           ASXSimpleId, simple_order, action_count, compute_count, simple_name;
    printf "SOCLE_COMPLETE|%o|enumerated|%o|%o\n",
           ASXSimpleId, action_count, compute_count;
elif ASXBatchEngineMode then
    require selected_count eq compute_count and
            selected_count eq ASXExpectedActions:
        "batch engine action count does not match metadata manifest";
    for record in batch_records do
        ASXEmitEngine(record[1], record[2], record[3]);
    end for;
    fprintf "/dev/stderr",
            "ASX_BATCH_SELECTED_COMPLETE simple_id=%o actions=%o\n",
            ASXSimpleId, selected_count;
else
    require selected_count eq 1: "selected label was not found exactly once";
    fprintf "/dev/stderr", "ASX_SELECTED_COMPLETE label=%o\n", selected_label;
end if;
quit;
