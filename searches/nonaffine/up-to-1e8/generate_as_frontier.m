///////////////////////////////////////////////////////////////////////////
// Almost-simple frontier generator for 10^7 < degree <= 10^8.
//
// The metadata pass constructs no large coset actions.  It dispositions every
// action as an order obstruction, a published theorem exclusion, or an exact
// computation target.  A selected computation target can then be exported in
// PRIMITIVE_SAXL_V1 format for the independent C++ engine.
//
// Required globals:
//   ASXSimpleId
// Optional globals:
//   ASXMinDegree      (default 10000001)
//   ASXMaxDegree      (default 100000000)
//   ASXMode           ("metadata", "fpr", "engine", "compact",
//                      "doublecoset", "densitysearch" or "graphwitness",
//                      default "metadata")
//   ASXEngineMode     (boolean/numeric CLI-safe alias for engine mode)
//   ASXFPRMode        (boolean/numeric CLI-safe alias for FPR mode)
//   ASXSelectedLabel  (one engine selector)
//   ASXSelectedExtension, ASXSelectedMaximal (CLI-safe engine selectors)
//   ASXBatchEngineMode (opt-in: export every compute action for this socle)
//   ASXExpectedActions (required positive batch action count)
//   ASXBatchSelectionMask (optional bit mask on canonical compute sequence)
//   ASXHybridBatchMode (one-pass FPR plus exact residual engine export)
//   ASXDoubleCosetMode (boolean/numeric CLI-safe alias for a compact exact
//                       H\\G/H benchmark on one selected action)
//   ASXCompactMode    (boolean/numeric CLI-safe alias for one-action FPR,
//                      followed when needed by exact compact H\\G/H analysis)
//   ASXSocleObstructionMode (boolean/numeric CLI-safe alias for an exact
//                         one-action (H meet T)\\T/(H meet T) obstruction;
//                         zero regular T-suborbits proves base size > 2 for G)
//   ASXSocleCohortMode (boolean/numeric CLI-safe alias for an exact metadata
//                       pass which groups compute actions by Aut(T)-conjugacy
//                       of K = H meet T; it performs no double-coset work)
//   ASXSocleCohortComputeMode (same exact cohort pass, followed by one
//                              K\\T/K decomposition per cohort)
//   ASXActionCohortMode (selector-free exact Aut(T)-conjugacy classification
//                        of full action stabilizers H; no coset action)
//   ASXCanonicalCompactBatchMode (fresh one-pass FPR plus exact H\\G/H for
//                                 every compute action of the selected socle;
//                                 no numeric action selectors)
//   ASXCanonicalFPRBatchMode (fresh selector-free exact FPR-only pass for
//                             every compute action; never builds a coset
//                             action or computes double cosets)
//   ASXIntrinsicCompactBlockMode (exact FPR plus H\G/H for every action in
//                                 one intrinsic degree/order block;
//                                 never selects by unstable maximal ordinal)
//   ASXIntrinsicPackMode (enumerate one intrinsic block once and serialize
//                         its exact permutation groups G,H for independent
//                         parallel consumers; performs no FPR/double cosets)
//   ASXSelectedDegree, ASXSelectedGroupOrder, ASXSelectedStabilizerOrder
//                                (required intrinsic block selector fields)
//   ASXPackDirectory             (required output directory in pack mode)
//   ASXCanonicalGraphWitnessBatchMode (selector-free exact common-neighbour
//                                      certificates for a canonical bit mask
//                                      of compute actions of one socle)
//   ASXCanonicalGraphWitnessVerifyBatchMode (selector-free independent replay
//                                            of that certificate batch)
//   ASXDensitySearchMode (boolean/numeric CLI-safe alias for a compact exact
//                         regular-double-coset density benchmark)
//   ASXDensityTrials  (maximum random candidates, default 10000)
//   ASXSeed           (positive deterministic random seed for density search)
//   ASXGraphWitnessMode (boolean/numeric CLI-safe alias for an exact compact
//                        common-neighbour certificate on H\G/H)
//   ASXGraphWitnessTrials (maximum global witness candidates, default 10000)
//   ASXGraphWitnessCertificateFile (mandatory output path for production)
//   ASXGraphWitnessCertificateDirectory (mandatory output directory for the
//                                        canonical graph-witness batch mode)
//   ASXGraphWitnessVerifyMode (reconstruct and exactly replay that certificate)
///////////////////////////////////////////////////////////////////////////

require assigned ASXSimpleId: "set ASXSimpleId";
if Type(ASXSimpleId) eq MonStgElt then
    ASXSimpleId := StringToInteger(ASXSimpleId);
end if;
if not assigned ASXMinDegree then ASXMinDegree := 10000001; end if;
if not assigned ASXMaxDegree then ASXMaxDegree := 100000000; end if;
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
if assigned ASXFPRMode then
    if Type(ASXFPRMode) eq MonStgElt then
        ASXFPRMode := ASXFPRMode in {"1", "true", "True"};
    end if;
    if Type(ASXFPRMode) eq RngIntElt then
        ASXFPRMode := ASXFPRMode ne 0;
    end if;
    if ASXFPRMode then ASXMode := "fpr"; end if;
end if;
if assigned ASXDoubleCosetMode then
    if Type(ASXDoubleCosetMode) eq MonStgElt then
        ASXDoubleCosetMode := ASXDoubleCosetMode in {"1", "true", "True"};
    end if;
    if Type(ASXDoubleCosetMode) eq RngIntElt then
        ASXDoubleCosetMode := ASXDoubleCosetMode ne 0;
    end if;
    if ASXDoubleCosetMode then ASXMode := "doublecoset"; end if;
end if;
if assigned ASXCompactMode then
    if Type(ASXCompactMode) eq MonStgElt then
        ASXCompactMode := ASXCompactMode in {"1", "true", "True"};
    end if;
    if Type(ASXCompactMode) eq RngIntElt then
        ASXCompactMode := ASXCompactMode ne 0;
    end if;
    if ASXCompactMode then ASXMode := "compact"; end if;
end if;
if assigned ASXSocleObstructionMode then
    if Type(ASXSocleObstructionMode) eq MonStgElt then
        ASXSocleObstructionMode :=
            ASXSocleObstructionMode in {"1", "true", "True"};
    end if;
    if Type(ASXSocleObstructionMode) eq RngIntElt then
        ASXSocleObstructionMode := ASXSocleObstructionMode ne 0;
    end if;
    if ASXSocleObstructionMode then ASXMode := "socleobstruction"; end if;
end if;
if assigned ASXSocleCohortMode then
    if Type(ASXSocleCohortMode) eq MonStgElt then
        ASXSocleCohortMode :=
            ASXSocleCohortMode in {"1", "true", "True"};
    end if;
    if Type(ASXSocleCohortMode) eq RngIntElt then
        ASXSocleCohortMode := ASXSocleCohortMode ne 0;
    end if;
    if ASXSocleCohortMode then ASXMode := "soclecohorts"; end if;
end if;
if assigned ASXSocleCohortComputeMode then
    if Type(ASXSocleCohortComputeMode) eq MonStgElt then
        ASXSocleCohortComputeMode :=
            ASXSocleCohortComputeMode in {"1", "true", "True"};
    end if;
    if Type(ASXSocleCohortComputeMode) eq RngIntElt then
        ASXSocleCohortComputeMode := ASXSocleCohortComputeMode ne 0;
    end if;
    if ASXSocleCohortComputeMode then
        ASXMode := "soclecohortcompute";
    end if;
end if;
if assigned ASXActionCohortMode then
    if Type(ASXActionCohortMode) eq MonStgElt then
        ASXActionCohortMode :=
            ASXActionCohortMode in {"1", "true", "True"};
    end if;
    if Type(ASXActionCohortMode) eq RngIntElt then
        ASXActionCohortMode := ASXActionCohortMode ne 0;
    end if;
    if ASXActionCohortMode then ASXMode := "actioncohorts"; end if;
end if;
if assigned ASXCanonicalCompactBatchMode then
    if Type(ASXCanonicalCompactBatchMode) eq MonStgElt then
        ASXCanonicalCompactBatchMode :=
            ASXCanonicalCompactBatchMode in {"1", "true", "True"};
    end if;
    if Type(ASXCanonicalCompactBatchMode) eq RngIntElt then
        ASXCanonicalCompactBatchMode := ASXCanonicalCompactBatchMode ne 0;
    end if;
    if ASXCanonicalCompactBatchMode then
        ASXMode := "canonicalcompactbatch";
    end if;
end if;
if assigned ASXCanonicalFPRBatchMode then
    if Type(ASXCanonicalFPRBatchMode) eq MonStgElt then
        ASXCanonicalFPRBatchMode :=
            ASXCanonicalFPRBatchMode in {"1", "true", "True"};
    end if;
    if Type(ASXCanonicalFPRBatchMode) eq RngIntElt then
        ASXCanonicalFPRBatchMode := ASXCanonicalFPRBatchMode ne 0;
    end if;
    if ASXCanonicalFPRBatchMode then
        ASXMode := "canonicalfprbatch";
    end if;
end if;
if assigned ASXIntrinsicCompactBlockMode then
    if Type(ASXIntrinsicCompactBlockMode) eq MonStgElt then
        ASXIntrinsicCompactBlockMode :=
            ASXIntrinsicCompactBlockMode in {"1", "true", "True"};
    end if;
    if Type(ASXIntrinsicCompactBlockMode) eq RngIntElt then
        ASXIntrinsicCompactBlockMode := ASXIntrinsicCompactBlockMode ne 0;
    end if;
    if ASXIntrinsicCompactBlockMode then
        ASXMode := "intrinsiccompactblock";
    end if;
end if;
if assigned ASXIntrinsicPackMode then
    if Type(ASXIntrinsicPackMode) eq MonStgElt then
        ASXIntrinsicPackMode :=
            ASXIntrinsicPackMode in {"1", "true", "True"};
    end if;
    if Type(ASXIntrinsicPackMode) eq RngIntElt then
        ASXIntrinsicPackMode := ASXIntrinsicPackMode ne 0;
    end if;
    if ASXIntrinsicPackMode then ASXMode := "intrinsicpack"; end if;
end if;
if assigned ASXCanonicalGraphWitnessBatchMode then
    if Type(ASXCanonicalGraphWitnessBatchMode) eq MonStgElt then
        ASXCanonicalGraphWitnessBatchMode :=
            ASXCanonicalGraphWitnessBatchMode in {"1", "true", "True"};
    end if;
    if Type(ASXCanonicalGraphWitnessBatchMode) eq RngIntElt then
        ASXCanonicalGraphWitnessBatchMode :=
            ASXCanonicalGraphWitnessBatchMode ne 0;
    end if;
    if ASXCanonicalGraphWitnessBatchMode then
        ASXMode := "canonicalgraphwitnessbatch";
    end if;
end if;
if assigned ASXCanonicalGraphWitnessVerifyBatchMode then
    if Type(ASXCanonicalGraphWitnessVerifyBatchMode) eq MonStgElt then
        ASXCanonicalGraphWitnessVerifyBatchMode :=
            ASXCanonicalGraphWitnessVerifyBatchMode in
                {"1", "true", "True"};
    end if;
    if Type(ASXCanonicalGraphWitnessVerifyBatchMode) eq RngIntElt then
        ASXCanonicalGraphWitnessVerifyBatchMode :=
            ASXCanonicalGraphWitnessVerifyBatchMode ne 0;
    end if;
    if ASXCanonicalGraphWitnessVerifyBatchMode then
        ASXMode := "canonicalgraphwitnessverifybatch";
    end if;
end if;
if assigned ASXDensitySearchMode then
    if Type(ASXDensitySearchMode) eq MonStgElt then
        ASXDensitySearchMode := ASXDensitySearchMode in {"1", "true", "True"};
    end if;
    if Type(ASXDensitySearchMode) eq RngIntElt then
        ASXDensitySearchMode := ASXDensitySearchMode ne 0;
    end if;
    if ASXDensitySearchMode then ASXMode := "densitysearch"; end if;
end if;
if assigned ASXGraphWitnessMode then
    if Type(ASXGraphWitnessMode) eq MonStgElt then
        ASXGraphWitnessMode :=
            ASXGraphWitnessMode in {"1", "true", "True"};
    end if;
    if Type(ASXGraphWitnessMode) eq RngIntElt then
        ASXGraphWitnessMode := ASXGraphWitnessMode ne 0;
    end if;
    if ASXGraphWitnessMode then ASXMode := "graphwitness"; end if;
end if;
if assigned ASXGraphWitnessVerifyMode then
    if Type(ASXGraphWitnessVerifyMode) eq MonStgElt then
        ASXGraphWitnessVerifyMode :=
            ASXGraphWitnessVerifyMode in {"1", "true", "True"};
    end if;
    if Type(ASXGraphWitnessVerifyMode) eq RngIntElt then
        ASXGraphWitnessVerifyMode := ASXGraphWitnessVerifyMode ne 0;
    end if;
    if ASXGraphWitnessVerifyMode then ASXMode := "graphwitnessverify"; end if;
end if;
if not assigned ASXBatchEngineMode then ASXBatchEngineMode := false; end if;
if Type(ASXBatchEngineMode) eq MonStgElt then
    ASXBatchEngineMode := ASXBatchEngineMode in {"1", "true", "True"};
end if;
if Type(ASXBatchEngineMode) eq RngIntElt then
    ASXBatchEngineMode := ASXBatchEngineMode ne 0;
end if;
if ASXBatchEngineMode then ASXMode := "engine"; end if;
if not assigned ASXHybridBatchMode then ASXHybridBatchMode := false; end if;
if Type(ASXHybridBatchMode) eq MonStgElt then
    ASXHybridBatchMode := ASXHybridBatchMode in {"1", "true", "True"};
end if;
if Type(ASXHybridBatchMode) eq RngIntElt then
    ASXHybridBatchMode := ASXHybridBatchMode ne 0;
end if;
require not (ASXBatchEngineMode and ASXHybridBatchMode):
    "batch engine and hybrid batch modes are mutually exclusive";
if ASXHybridBatchMode then ASXMode := "hybrid"; end if;
if not assigned ASXExpectedActions then ASXExpectedActions := 0; end if;
if Type(ASXExpectedActions) eq MonStgElt then
    ASXExpectedActions := StringToInteger(ASXExpectedActions);
end if;
if not assigned ASXBatchSelectionMask then ASXBatchSelectionMask := 0; end if;
if Type(ASXBatchSelectionMask) eq MonStgElt then
    ASXBatchSelectionMask := StringToInteger(ASXBatchSelectionMask);
end if;
require ASXBatchSelectionMask ge 0: "negative batch selection mask";
if not assigned ASXDensityTrials then ASXDensityTrials := 10000; end if;
if Type(ASXDensityTrials) eq MonStgElt then
    ASXDensityTrials := StringToInteger(ASXDensityTrials);
end if;
require ASXDensityTrials ge 1: "density-search trial count must be positive";
if not assigned ASXGraphWitnessTrials then ASXGraphWitnessTrials := 10000; end if;
if Type(ASXGraphWitnessTrials) eq MonStgElt then
    ASXGraphWitnessTrials := StringToInteger(ASXGraphWitnessTrials);
end if;
require ASXGraphWitnessTrials ge 1:
    "graph-witness trial count must be positive";
if not assigned ASXGraphWitnessCertificateFile then
    ASXGraphWitnessCertificateFile := "";
end if;
if not assigned ASXGraphWitnessCertificateDirectory then
    ASXGraphWitnessCertificateDirectory := "";
end if;
if not assigned ASXPackDirectory then ASXPackDirectory := ""; end if;
if not assigned ASXSeed then ASXSeed := 1729; end if;
if Type(ASXSeed) eq MonStgElt then
    ASXSeed := StringToInteger(ASXSeed);
end if;
require ASXSeed ge 1: "density-search seed must be positive";
SetSeed(ASXSeed);
if Type(ASXMinDegree) eq MonStgElt then
    ASXMinDegree := StringToInteger(ASXMinDegree);
end if;
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
require 1 le ASXSimpleId and ASXSimpleId le NumberOfSimpleGroups():
    "bad simple-group ID";
require ASXMinDegree eq 10000001 and ASXMaxDegree eq 100000000:
    "this source is bound to the new degree decade";
require (ASXMode eq "metadata") or (ASXMode eq "fpr") or
        (ASXMode eq "engine") or (ASXMode eq "hybrid") or
        (ASXMode eq "compact") or
        (ASXMode eq "socleobstruction") or
        (ASXMode eq "soclecohorts") or
        (ASXMode eq "soclecohortcompute") or
        (ASXMode eq "actioncohorts") or
        (ASXMode eq "canonicalcompactbatch") or
        (ASXMode eq "canonicalfprbatch") or
        (ASXMode eq "intrinsiccompactblock") or
        (ASXMode eq "intrinsicpack") or
        (ASXMode eq "canonicalgraphwitnessbatch") or
        (ASXMode eq "canonicalgraphwitnessverifybatch") or
        (ASXMode eq "doublecoset") or
        (ASXMode eq "densitysearch") or
        (ASXMode eq "graphwitness") or
        (ASXMode eq "graphwitnessverify"): "bad ASXMode";
if not assigned ASXSelectedLabel then ASXSelectedLabel := ""; end if;
if not assigned ASXSelectedExtension then ASXSelectedExtension := 0; end if;
if not assigned ASXSelectedMaximal then ASXSelectedMaximal := 0; end if;
if Type(ASXSelectedExtension) eq MonStgElt then
    ASXSelectedExtension := StringToInteger(ASXSelectedExtension);
end if;
if Type(ASXSelectedMaximal) eq MonStgElt then
    ASXSelectedMaximal := StringToInteger(ASXSelectedMaximal);
end if;
if not assigned ASXSelectedDegree then ASXSelectedDegree := 0; end if;
if not assigned ASXSelectedGroupOrder then ASXSelectedGroupOrder := 0; end if;
if not assigned ASXSelectedStabilizerOrder then
    ASXSelectedStabilizerOrder := 0;
end if;
if Type(ASXSelectedDegree) eq MonStgElt then
    ASXSelectedDegree := StringToInteger(ASXSelectedDegree);
end if;
if Type(ASXSelectedGroupOrder) eq MonStgElt then
    ASXSelectedGroupOrder := StringToInteger(ASXSelectedGroupOrder);
end if;
if Type(ASXSelectedStabilizerOrder) eq MonStgElt then
    ASXSelectedStabilizerOrder :=
        StringToInteger(ASXSelectedStabilizerOrder);
end if;
if ASXMode in {"intrinsiccompactblock", "intrinsicpack"} then
    require ASXSelectedExtension eq 0 and ASXSelectedDegree ge 1 and
            ASXSelectedGroupOrder ge 1 and ASXSelectedStabilizerOrder ge 1 and
            ASXSelectedMaximal eq 0 and
            ASXExpectedActions ge 1:
        "intrinsic block operation needs its full key and positive multiplicity";
    if ASXMode eq "intrinsiccompactblock" then
        require ASXSelectedLabel eq "" or ASXExpectedActions eq 1:
            "an intrinsic compact label selector must be a singleton";
    else
        require ASXSelectedLabel eq "" and ASXPackDirectory ne "":
            "intrinsic pack needs an output directory and forbids label selectors";
    end if;
end if;
if ASXMode in {"fpr", "engine", "hybrid", "compact", "doublecoset",
               "socleobstruction",
               "densitysearch", "graphwitness", "graphwitnessverify"} then
    if ASXHybridBatchMode then
        require ASXSelectedLabel eq "" and ASXSelectedExtension eq 0 and
                ASXSelectedMaximal eq 0 and
                ASXExpectedActions ge 1:
            "hybrid batch mode needs a payload directory and expected count";
    elif ASXBatchEngineMode then
        require ASXMode eq "engine": "batch mode is available only for engine export";
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

// Exact compact alternative to constructing the degree-sized coset action.
// The H-orbits on G/H are the H-H double cosets.  If D = H g H, then
// |D|/|H| is the corresponding subdegree, so the returned double-coset sizes
// give an independent completeness check without materialising permutations
// on [G:H] points.  A suborbit is regular precisely when |D| = |H|^2.
ASXEmitDoubleCosetBenchmark := function(G, H, label)
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

// Every primitive almost-simple action G/H has G = T H for its transitive
// socle T.  Put K = H meet T.  Every G-double coset has a representative in
// T, and H meet H^t contains K meet K^t.  Consequently, if the exact K\\T/K
// decomposition has no regular K-suborbit, then G/H has no regular
// H-suborbit and therefore has base size greater than two.  A positive
// regular K-suborbit is deliberately only advisory: outer elements can still
// obstruct a regular H-suborbit.
ASXEmitSocleObstruction := function(T, G, H, label)
    K := H meet T;
    degree := Index(G, H);
    assert G eq sub<G | Setseq(Generators(T)) cat Setseq(Generators(H))> and
           Index(T, K) eq degree;
    k_order := Order(K);
    started := Cputime();
    representatives, double_coset_sizes :=
        DoubleCosetRepresentatives(T, K, K);
    elapsed := Cputime(started);
    assert #representatives eq #double_coset_sizes;
    assert #representatives ge 1 and representatives[1] eq Identity(T);

    subdegree_sum := 0;
    regular_orbits := 0;
    printf "AS8_SOCLE_DOUBLE_COSET_V1|%o|degree=%o|T=%o|K=%o|rank=%o|cpu_ms=%o\n",
           label, degree, Order(T), k_order, #representatives,
           Round(1000 * elapsed);
    for i in [1 .. #representatives] do
        double_size := double_coset_sizes[i];
        assert double_size mod k_order eq 0;
        subdegree := double_size div k_order;
        subdegree_sum +:= subdegree;
        if double_size eq k_order^2 then regular_orbits +:= 1; end if;
        printf "SOCLE_DOUBLE_COSET|sequence=%o|double_size=%o|subdegree=%o|regular=%o\n",
               i, double_size, subdegree, double_size eq k_order^2;
    end for;
    assert subdegree_sum eq degree;
    printf "AS8_SOCLE_DOUBLE_COSET_COMPLETE|%o|rank=%o|regular_orbits=%o|subdegree_sum=%o|degree=%o|cpu_ms=%o\n",
           label, #representatives, regular_orbits, subdegree_sum, degree,
           Round(1000 * elapsed);
    if regular_orbits eq 0 then
        printf "AS8_SOCLE_OBSTRUCTION_COMPLETE|%o|disposition=base_gt_2|degree=%o|rank=%o|regular_orbits=0\n",
               label, degree, #representatives;
    else
        printf "AS8_SOCLE_OBSTRUCTION_INCONCLUSIVE|%o|degree=%o|rank=%o|regular_orbits=%o\n",
               label, degree, #representatives, regular_orbits;
    end if;
    return regular_orbits, #representatives;
end function;

// A strict-density proof needs only enough pairwise distinct regular
// H-suborbits to cover more than half the vertices.  A representative g gives
// a regular suborbit exactly when H meet H^g is trivial, and two such
// suborbits are distinct exactly when their representatives lie in distinct
// H-H double cosets.  This produces a positive exact certificate without
// enumerating all double cosets.  Failure to find enough representatives is
// advisory only and never becomes a negative scientific conclusion.
ASXEmitDensitySearchBenchmark := procedure(G, H, label)
    degree := Index(G, H);
    h_order := Order(H);
    required_orbits := (degree div (2 * h_order)) + 1;
    regular_representatives := [G | ];
    regular_canonical_images := [* *];
    canonical_base := [];
    trivial_H := sub<H | >;
    trials := 0;
    started := Cputime();
    while #regular_representatives lt required_orbits and
          trials lt ASXDensityTrials do
        candidate := Random(G);
        trials +:= 1;
        if Order(H meet H^candidate) ne 1 then continue; end if;
        if #regular_representatives eq 0 then
            canonical_image, canonical_base :=
                DoubleCosetCanonical(G, H, candidate, H : M := trivial_H);
        else
            canonical_image, _ :=
                DoubleCosetCanonical(G, H, candidate, H :
                                     B := canonical_base, M := trivial_H);
        end if;
        is_new := true;
        for old_image in regular_canonical_images do
            if canonical_image eq old_image then
                is_new := false;
                break;
            end if;
        end for;
        if is_new then
            Append(~regular_representatives, candidate);
            Append(~regular_canonical_images, canonical_image);
        end if;
    end while;
    elapsed := Cputime(started);
    regular_lower_bound := #regular_representatives * h_order;
    success := 2 * regular_lower_bound gt degree;
    printf "AS8_DENSITY_SEARCH_V1|%o|degree=%o|G=%o|H=%o|required_orbits=%o|found_orbits=%o|trials=%o|regular_lower_bound=%o|success=%o|cpu_ms=%o\n",
           label, degree, Order(G), h_order, required_orbits,
           #regular_representatives, trials, regular_lower_bound, success,
           Round(1000 * elapsed);
    if success then
        // Recheck the exact certificate immediately before promotion.
        for i in [1 .. #regular_representatives] do
            assert Order(H meet H^regular_representatives[i]) eq 1;
            rechecked_image, _ :=
                DoubleCosetCanonical(G, H, regular_representatives[i], H :
                                     B := canonical_base, M := trivial_H);
            assert rechecked_image eq regular_canonical_images[i];
            for j in [1 .. i - 1] do
                assert regular_canonical_images[i] ne
                    regular_canonical_images[j];
            end for;
        end for;
        assert 2 * #regular_representatives * h_order gt degree;
        printf "AS8_DENSITY_SEARCH_COMPLETE|%o|regular_orbits=%o|regular_points_lower_bound=%o|degree=%o|cpu_ms=%o\n",
               label, #regular_representatives, regular_lower_bound, degree,
               Round(1000 * elapsed);
    else
        printf "AS8_DENSITY_SEARCH_INCONCLUSIVE|%o|found_orbits=%o|required_orbits=%o|trials=%o|cpu_ms=%o\n",
               label, #regular_representatives, required_orbits, trials,
               Round(1000 * elapsed);
    end if;
end procedure;

// Let alpha = H and let Hx run through the H-orbits on G/H.  A vertex Hy is
// adjacent in the Saxl graph to both alpha and Hx precisely when
// H meet H^y = 1 and H^x meet H^y = 1.  Thus one explicitly checked witness
// for every exact H-H double-coset representative proves that every pair of
// vertices has a common neighbour, by G-transitivity.  Random search affects
// only whether a positive certificate is found: every promoted certificate
// is rechecked exactly, while exhaustion is advisory.
ASXEmitGraphWitnessBenchmark := procedure(G, H, label, certificate_path)
    degree := Index(G, H);
    h_order := Order(H);
    started := Cputime();
    representatives, double_coset_sizes :=
        DoubleCosetRepresentatives(G, H, H);
    assert #representatives eq #double_coset_sizes;
    assert #representatives ge 1 and representatives[1] eq Identity(G);

    subdegree_sum := 0;
    for double_size in double_coset_sizes do
        assert double_size mod h_order eq 0;
        subdegree_sum +:= double_size div h_order;
    end for;
    assert subdegree_sum eq degree;

    conjugates := [H^representative : representative in representatives];
    assigned_witness := [0 : i in [1 .. #representatives]];
    uncovered := [1 .. #representatives];
    witnesses := [G | ];
    trials := 0;
    while #uncovered gt 0 and trials lt ASXGraphWitnessTrials do
        candidate := Random(G);
        trials +:= 1;
        candidate_conjugate := H^candidate;
        if Order(H meet candidate_conjugate) ne 1 then continue; end if;

        witness_number := #witnesses + 1;
        newly_covered := [];
        still_uncovered := [];
        for i in uncovered do
            if Order(conjugates[i] meet candidate_conjugate) eq 1 then
                Append(~newly_covered, i);
            else
                Append(~still_uncovered, i);
            end if;
        end for;
        if #newly_covered eq 0 then continue; end if;
        Append(~witnesses, candidate);
        for i in newly_covered do
            assigned_witness[i] := witness_number;
        end for;
        uncovered := still_uncovered;
        printf "AS8_GRAPH_WITNESS_PROGRESS|%o|trial=%o|witnesses=%o|newly_covered=%o|remaining=%o\n",
               label, trials, #witnesses, #newly_covered, #uncovered;
    end while;

    elapsed := Cputime(started);
    printf "AS8_GRAPH_WITNESS_V1|%o|degree=%o|G=%o|H=%o|rank=%o|witnesses=%o|trials=%o|remaining=%o|cpu_ms=%o\n",
           label, degree, Order(G), h_order, #representatives, #witnesses,
           trials, #uncovered, Round(1000 * elapsed);
    if #uncovered eq 0 then
        // Exact semantic replay immediately before promotion.
        for witness in witnesses do
            assert Order(H meet H^witness) eq 1;
        end for;
        for i in [1 .. #representatives] do
            witness_number := assigned_witness[i];
            assert 1 le witness_number and witness_number le #witnesses;
            assert Order(conjugates[i] meet H^witnesses[witness_number]) eq 1;
            printf "AS8_GRAPH_WITNESS|sequence=%o|double_size=%o|subdegree=%o|witness=%o\n",
                   i, double_coset_sizes[i],
                   double_coset_sizes[i] div h_order, witness_number;
        end for;
        if certificate_path ne "" then
            certificate_file :=
                Open(certificate_path, "w");
            fprintf certificate_file, "<%o,%o,%o,%o>",
                    [Eltseq(representative) :
                     representative in representatives],
                    double_coset_sizes,
                    [Eltseq(witness) : witness in witnesses],
                    assigned_witness;
            delete certificate_file;
        end if;
        printf "AS8_GRAPH_WITNESS_COMPLETE|%o|disposition=common_neighbour|degree=%o|rank=%o|witnesses=%o|trials=%o|subdegree_sum=%o|cpu_ms=%o\n",
               label, degree, #representatives, #witnesses, trials,
               subdegree_sum, Round(1000 * elapsed);
    else
        printf "AS8_GRAPH_WITNESS_INCONCLUSIVE|%o|degree=%o|rank=%o|witnesses=%o|trials=%o|remaining=%o|cpu_ms=%o\n",
               label, degree, #representatives, #witnesses, trials,
               #uncovered, Round(1000 * elapsed);
    end if;
end procedure;

ASXVerifyGraphWitnessCertificate := procedure(G, H, label, path)
    assert path ne "";
    started := Cputime();
    certificate := eval Read(path);
    assert Type(certificate) eq Tup and #certificate eq 4;
    representative_images := certificate[1];
    certificate_sizes := certificate[2];
    witness_images := certificate[3];
    assigned_witness := certificate[4];
    ambient := Sym(Degree(G));
    representatives := [G |
        G!(ambient!images) : images in representative_images];
    witnesses := [G | G!(ambient!images) : images in witness_images];

    assert #representatives ge 1 and
            representatives[1] eq Identity(G) and
            #certificate_sizes eq #representatives and
            #assigned_witness eq #representatives;
    degree := Index(G, H);
    h_order := Order(H);
    canonical_images := [* *];
    canonical_intersection :=
        H meet H^(representatives[1]^-1);
    canonical_image, canonical_base :=
        DoubleCosetCanonical(G, H, representatives[1], H :
                             M := canonical_intersection);
    Append(~canonical_images, canonical_image);
    subdegree_sum := 0;
    for i in [1 .. #representatives] do
        if i gt 1 then
            canonical_intersection :=
                H meet H^(representatives[i]^-1);
            canonical_image, _ := DoubleCosetCanonical(
                G, H, representatives[i], H : B := canonical_base,
                M := canonical_intersection);
            for old_image in canonical_images do
                assert canonical_image ne old_image;
            end for;
            Append(~canonical_images, canonical_image);
        end if;
        intersection_order := Order(canonical_intersection);
        assert certificate_sizes[i] * intersection_order eq h_order^2;
        assert certificate_sizes[i] mod h_order eq 0;
        subdegree_sum +:= certificate_sizes[i] div h_order;
        witness_number := assigned_witness[i];
        assert 1 le witness_number and witness_number le #witnesses;
        assert Order(H meet H^witnesses[witness_number]) eq 1 and
                Order((H^representatives[i]) meet
                      H^witnesses[witness_number]) eq 1;
    end for;
    assert subdegree_sum eq degree and
            &+certificate_sizes eq Order(G);
    elapsed := Cputime(started);
    printf "AS8_GRAPH_WITNESS_REPLAY_COMPLETE|%o|disposition=common_neighbour|degree=%o|rank=%o|witnesses=%o|subdegree_sum=%o|cpu_ms=%o\n",
           label, degree, #representatives, #witnesses, subdegree_sum,
           Round(1000 * elapsed);
end procedure;

ASXFilePrintPermutation := procedure(file, g, d)
    width := 1;
    capacity := 94;
    while capacity lt d do
        width +:= 1;
        capacity *:= 94;
    end while;
    fprintf file, "packed_gen %o ", width;
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
        fprintf file, "%o", &cat pieces;
    end for;
    fprintf file, "\n";
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

ASXEmitEngineFile := procedure(G, H, label, path)
    d := Index(G, H);
    h_order := Order(H);
    action_map, P := CosetAction(G, H);
    assert Degree(P) eq d and IsPrimitive(P);
    assert Order(P) eq Order(G) and Order(Kernel(action_map)) eq 1;
    HP := Stabilizer(P, 1);
    assert Order(HP) eq h_order;
    file := Open(path, "w");
    fprintf file, "PRIMITIVE_SAXL_V1\n";
    fprintf file, "action\n";
    fprintf file, "label %o\n", label;
    fprintf file, "degree %o\n", d;
    fprintf file, "stabilizer_order %o\n", h_order;
    fprintf file, "classification compute\n";
    fprintf file, "regular_orbits 0\n";
    fprintf file, "regular_count 0\n";
    hgens := Setseq(Generators(HP));
    assert #hgens gt 0;
    fprintf file, "hgens %o\n", #hgens;
    for g in hgens do ASXFilePrintPermutation(file, g, d); end for;
    gens := ASXDeterministicTransitiveActionGenerators(P, G, action_map);
    assert #Orbit(sub<P | gens>, 1) eq d;
    fprintf file, "gens %o\n", #gens;
    for g in gens do ASXFilePrintPermutation(file, g, d); end for;
    fprintf file, "end\n";
    delete file;
end procedure;

// Compute the exact prime-order fixed-point-ratio upper bound without
// constructing the degree-sized coset action.  For the action G/H,
//   Qhat(G,H) = sum_x |x^G cap H|^2 / |x^G|,
// where x runs over prime-order G-classes meeting H.  Qhat < 1/2 implies
// that every two Saxl-neighbourhoods intersect.
ASXEmitFPR := function(G, H, label, extension_number, maximal_number)
    d := Index(G, H);
    printf "AS_ACTION|%o|sid=%o|ext=%o|max=%o|degree=%o|G=%o|H=%o\n",
           label, ASXSimpleId, extension_number, maximal_number, d,
           Order(G), Order(H);

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
    printf "AS8_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
           label, #prime_classes, #intersections, Numerator(bound),
           Denominator(bound), bound lt 1/2;
    return bound;
end function;

simple_tuple, simple_order := SimpleGroupId(ASXSimpleId);
simple_name := SimpleGroupName(ASXSimpleId);
require simple_order le ASXMaxDegree * (ASXMaxDegree - 1):
    "simple-group ID lies beyond the base-two order bound";

if ASXMode eq "engine" then print "PRIMITIVE_SAXL_V1"; end if;
if ASXMode eq "fpr" then print "AS8_FPR_V1"; end if;
if ASXMode eq "hybrid" then print "AS8_HYBRID_V1"; end if;

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

// Burness--Giudici, MPCPS 168 (2020), Theorem 5.1, plus the exact
// degree-window audit and Morris--Spiga base-size formula recorded in
// LITERATURE_LEDGER.md: no alternating-socle base-two target survives.
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
// only the outer PSO(3,q) actions at q=31,32,37.  Route every other
// PSU3 socle before constructing Aut(T) or enumerating maximal subgroups.
// This is a preprint-backed window disposition, not a published exclusion.
if is_psu3 and not psu_q in {31, 32, 37} then
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
hybrid_closed_count := 0;
hybrid_residual_count := 0;
socle_cohort_records := [* *];
socle_cohort_action_count := 0;
action_cohort_records := [* *];
action_cohort_action_count := 0;
canonical_compact_records := [* *];
intrinsic_compact_records := [* *];

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
        label := Sprintf("AS8_sid%o_ext%o_max%o_d%o",
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
        elif ASXMode in {"soclecohorts", "soclecohortcompute"} then
            if disposition eq "compute" then
                K := H meet T;
                require G eq sub<G | Setseq(Generators(T)) cat
                                  Setseq(Generators(H))> and Index(T, K) eq d:
                    "invalid socle action binding";
                cohort_number := 0;
                for i in [1 .. #socle_cohort_records] do
                    old_K := socle_cohort_records[i][1];
                    if Order(K) eq Order(old_K) and
                       IsConjugate(A, K, old_K) then
                        cohort_number := i;
                        break;
                    end if;
                end for;
                if cohort_number eq 0 then
                    Append(~socle_cohort_records, <K, d, label, G, H>);
                    cohort_number := #socle_cohort_records;
                end if;
                socle_cohort_action_count +:= 1;
                printf "AS8_SOCLE_COHORT_V1|label=%o|sid=%o|cohort=%o|degree=%o|K=%o|G=%o|H=%o|extension=%o|maximal=%o|representative=%o\n",
                       label, ASXSimpleId, cohort_number, d, Order(K),
                       Order(G), h_order, extension_number, maximal_number,
                       socle_cohort_records[cohort_number][3];
            end if;
        elif ASXMode eq "actioncohorts" then
            if disposition eq "compute" then
                require G eq sub<G | Setseq(Generators(T)) cat
                                  Setseq(Generators(H))> and
                        Index(T, H meet T) eq d:
                    "invalid full-action binding";
                action_cohort_number := 0;
                for i in [1 .. #action_cohort_records] do
                    old_action_H := action_cohort_records[i][1];
                    if Order(H) eq Order(old_action_H) and
                       IsConjugate(A, H, old_action_H) then
                        action_cohort_number := i;
                        break;
                    end if;
                end for;
                if action_cohort_number eq 0 then
                    Append(~action_cohort_records, <H, d, label, G>);
                    action_cohort_number := #action_cohort_records;
                end if;
                action_cohort_action_count +:= 1;
                printf "AS8_ACTION_COHORT_V1|label=%o|sid=%o|cohort=%o|degree=%o|G=%o|H=%o|extension=%o|maximal=%o|representative=%o\n",
                       label, ASXSimpleId, action_cohort_number, d, Order(G),
                       h_order, extension_number, maximal_number,
                       action_cohort_records[action_cohort_number][3];
            end if;
        elif ASXMode in {"intrinsiccompactblock", "intrinsicpack"} then
            if disposition eq "compute" and
               d eq ASXSelectedDegree and Order(G) eq ASXSelectedGroupOrder and
               h_order eq ASXSelectedStabilizerOrder and
               (ASXMode eq "intrinsicpack" or ASXSelectedLabel eq "" or
                label eq ASXSelectedLabel) then
                Append(~intrinsic_compact_records,
                       <G, H, label, d, extension_number, maximal_number>);
            end if;
        elif ASXMode in {"canonicalcompactbatch", "canonicalfprbatch",
                         "canonicalgraphwitnessbatch",
                         "canonicalgraphwitnessverifybatch"} then
            if disposition eq "compute" then
                Append(~canonical_compact_records,
                       <G, H, label, d, extension_number, maximal_number>);
            end if;
        elif ASXHybridBatchMode then
            if disposition eq "compute" then
                bound := ASXEmitFPR(G, H, label, extension_number,
                                    maximal_number);
                if bound lt 1/2 then
                    hybrid_closed_count +:= 1;
                else
                    ASXEmitEngineFile(G, H, label, label cat ".input");
                    hybrid_residual_count +:= 1;
                end if;
                selected_count +:= 1;
            end if;
        elif ASXBatchEngineMode then
            if disposition eq "compute" then
                selected_by_mask := ASXBatchSelectionMask eq 0 or
                    ((ASXBatchSelectionMask div 2^(compute_count - 1)) mod 2 eq 1);
                if selected_by_mask then
                    Append(~batch_records, <G, H, label>);
                    selected_count +:= 1;
                end if;
            end if;
        elif label eq ASXSelectedLabel or
             (extension_number eq ASXSelectedExtension and
              maximal_number eq ASXSelectedMaximal) then
            require disposition eq "compute":
                "selected action is already theorem/order disposed";
            if ASXMode eq "fpr" then
                _ := ASXEmitFPR(G, H, label, extension_number, maximal_number);
            elif ASXMode eq "doublecoset" then
                benchmark_regular_orbits, benchmark_rank :=
                    ASXEmitDoubleCosetBenchmark(G, H, label);
            elif ASXMode eq "compact" then
                bound := ASXEmitFPR(G, H, label, extension_number,
                                    maximal_number);
                if bound lt 1/2 then
                    printf "AS8_COMPACT_COMPLETE|%o|disposition=fpr_density|degree=%o|regular_orbits=0|regular_points=0\n",
                           label, d;
                else
                    regular_orbits, rank :=
                        ASXEmitDoubleCosetBenchmark(G, H, label);
                    regular_points := regular_orbits * h_order;
                    if regular_orbits eq 0 then
                        compact_disposition := "base_gt_2";
                    elif 2 * regular_points gt d then
                        compact_disposition := "exact_density";
                    else
                        compact_disposition := "graph_residual";
                    end if;
                    printf "AS8_COMPACT_COMPLETE|%o|disposition=%o|degree=%o|rank=%o|regular_orbits=%o|regular_points=%o\n",
                           label, compact_disposition, d, rank,
                           regular_orbits, regular_points;
                end if;
            elif ASXMode eq "socleobstruction" then
                socle_regular_orbits, socle_rank :=
                    ASXEmitSocleObstruction(T, G, H, label);
            elif ASXMode eq "densitysearch" then
                ASXEmitDensitySearchBenchmark(G, H, label);
            elif ASXMode eq "graphwitness" then
                ASXEmitGraphWitnessBenchmark(
                    G, H, label, ASXGraphWitnessCertificateFile);
            elif ASXMode eq "graphwitnessverify" then
                ASXVerifyGraphWitnessCertificate(
                    G, H, label, ASXGraphWitnessCertificateFile);
            else
                ASXEmitEngine(G, H, label);
            end if;
            selected_count +:= 1;
            selected_label := label;
            fprintf "/dev/stderr", "ASX_SELECTED_COMPLETE label=%o\n",
                    selected_label;
            // The audited metadata plan already proves that the canonical
            // selector occurs exactly once.  Do not enumerate unrelated
            // extensions and maximals after the selected action is sealed.
            quit;
        end if;
    end for;
end for;

if ASXMode eq "metadata" then
    printf "SOCLE|%o|%o|enumerated|%o|%o|%o\n",
           ASXSimpleId, simple_order, action_count, compute_count, simple_name;
    printf "SOCLE_COMPLETE|%o|enumerated|%o|%o\n",
           ASXSimpleId, action_count, compute_count;
elif ASXMode eq "soclecohorts" then
    require socle_cohort_action_count eq compute_count:
        "socle cohort action count does not match compute count";
    printf "AS8_SOCLE_COHORT_COMPLETE_V1|sid=%o|actions=%o|cohorts=%o\n",
           ASXSimpleId, socle_cohort_action_count, #socle_cohort_records;
elif ASXMode eq "soclecohortcompute" then
    require socle_cohort_action_count eq compute_count:
        "socle cohort action count does not match compute count";
    socle_closed_count := 0;
    socle_inconclusive_count := 0;
    for i in [1 .. #socle_cohort_records] do
        record := socle_cohort_records[i];
        cohort_label := Sprintf("AS8_sid%o_kcohort%o_d%o_k%o",
                                ASXSimpleId, i, record[2], Order(record[1]));
        socle_regular_orbits, socle_rank :=
            ASXEmitSocleObstruction(T, record[4], record[5], cohort_label);
        if socle_regular_orbits eq 0 then
            socle_closed_count +:= 1;
        else
            socle_inconclusive_count +:= 1;
        end if;
    end for;
    printf "AS8_SOCLE_COHORT_COMPUTE_COMPLETE_V1|sid=%o|actions=%o|cohorts=%o|closed=%o|inconclusive=%o\n",
           ASXSimpleId, socle_cohort_action_count, #socle_cohort_records,
           socle_closed_count, socle_inconclusive_count;
elif ASXMode eq "actioncohorts" then
    require action_cohort_action_count eq compute_count:
        "full-action cohort count does not match compute count";
    printf "AS8_ACTION_COHORT_COMPLETE_V1|sid=%o|actions=%o|cohorts=%o\n",
           ASXSimpleId, action_cohort_action_count, #action_cohort_records;
elif ASXMode eq "canonicalcompactbatch" then
    require #canonical_compact_records eq compute_count:
        "canonical compact action count does not match compute count";
    compact_fpr_count := 0;
    compact_base_count := 0;
    compact_density_count := 0;
    compact_residual_count := 0;
    for record in canonical_compact_records do
        compact_G := record[1];
        compact_H := record[2];
        compact_label := record[3];
        compact_degree := record[4];
        compact_h_order := Order(compact_H);
        compact_bound := ASXEmitFPR(compact_G, compact_H, compact_label,
                                    record[5], record[6]);
        if compact_bound lt 1/2 then
            compact_fpr_count +:= 1;
            printf "AS8_COMPACT_COMPLETE|%o|disposition=fpr_density|degree=%o|regular_orbits=0|regular_points=0\n",
                   compact_label, compact_degree;
        else
            compact_regular_orbits, compact_rank :=
                ASXEmitDoubleCosetBenchmark(compact_G, compact_H,
                                            compact_label);
            compact_regular_points := compact_regular_orbits * compact_h_order;
            if compact_regular_orbits eq 0 then
                compact_disposition := "base_gt_2";
                compact_base_count +:= 1;
            elif 2 * compact_regular_points gt compact_degree then
                compact_disposition := "exact_density";
                compact_density_count +:= 1;
            else
                compact_disposition := "graph_residual";
                compact_residual_count +:= 1;
            end if;
            printf "AS8_COMPACT_COMPLETE|%o|disposition=%o|degree=%o|rank=%o|regular_orbits=%o|regular_points=%o\n",
                   compact_label, compact_disposition, compact_degree,
                   compact_rank, compact_regular_orbits,
                   compact_regular_points;
        end if;
    end for;
    printf "AS8_CANONICAL_COMPACT_BATCH_COMPLETE_V1|sid=%o|actions=%o|fpr_density=%o|base_gt_2=%o|exact_density=%o|graph_residual=%o\n",
           ASXSimpleId, #canonical_compact_records, compact_fpr_count,
           compact_base_count, compact_density_count,
           compact_residual_count;
elif ASXMode eq "canonicalfprbatch" then
    require #canonical_compact_records eq compute_count and
            ASXExpectedActions eq compute_count:
        "canonical FPR action count does not match expected metadata count";
    fpr_closed_count := 0;
    fpr_inconclusive_count := 0;
    for record in canonical_compact_records do
        fpr_bound := ASXEmitFPR(record[1], record[2], record[3],
                               record[5], record[6]);
        if fpr_bound lt 1/2 then
            fpr_disposition := "fpr_density";
            fpr_closed_count +:= 1;
        else
            fpr_disposition := "fpr_inconclusive";
            fpr_inconclusive_count +:= 1;
        end if;
        printf "AS8_CANONICAL_FPR_ACTION_COMPLETE_V1|%o|disposition=%o|degree=%o\n",
               record[3], fpr_disposition, record[4];
    end for;
    printf "AS8_CANONICAL_FPR_BATCH_COMPLETE_V1|sid=%o|actions=%o|fpr_density=%o|fpr_inconclusive=%o\n",
           ASXSimpleId, #canonical_compact_records, fpr_closed_count,
           fpr_inconclusive_count;
elif ASXMode eq "intrinsiccompactblock" then
    require #intrinsic_compact_records eq ASXExpectedActions:
        "intrinsic compact block multiplicity mismatch";
    intrinsic_fpr_count := 0;
    intrinsic_base_count := 0;
    intrinsic_density_count := 0;
    intrinsic_residual_count := 0;
    for record in intrinsic_compact_records do
        intrinsic_G := record[1];
        intrinsic_H := record[2];
        intrinsic_label := record[3];
        intrinsic_degree := record[4];
        intrinsic_h_order := Order(intrinsic_H);
        intrinsic_bound := ASXEmitFPR(intrinsic_G, intrinsic_H,
                                      intrinsic_label, record[5], record[6]);
        if intrinsic_bound lt 1/2 then
            intrinsic_fpr_count +:= 1;
            printf "AS8_COMPACT_COMPLETE|%o|disposition=fpr_density|degree=%o|regular_orbits=0|regular_points=0\n",
                   intrinsic_label, intrinsic_degree;
            printf "AS8_INTRINSIC_COMPACT_ACTION_COMPLETE_V1|%o|disposition=fpr_density|degree=%o|regular_orbits=0|regular_points=0\n",
                   intrinsic_label, intrinsic_degree;
        else
            intrinsic_regular_orbits, intrinsic_rank :=
                ASXEmitDoubleCosetBenchmark(intrinsic_G, intrinsic_H,
                                            intrinsic_label);
            intrinsic_regular_points :=
                intrinsic_regular_orbits * intrinsic_h_order;
            if intrinsic_regular_orbits eq 0 then
                intrinsic_disposition := "base_gt_2";
                intrinsic_base_count +:= 1;
            elif 2 * intrinsic_regular_points gt intrinsic_degree then
                intrinsic_disposition := "exact_density";
                intrinsic_density_count +:= 1;
            else
                intrinsic_disposition := "graph_residual";
                intrinsic_residual_count +:= 1;
            end if;
            printf "AS8_COMPACT_COMPLETE|%o|disposition=%o|degree=%o|rank=%o|regular_orbits=%o|regular_points=%o\n",
                   intrinsic_label, intrinsic_disposition, intrinsic_degree,
                   intrinsic_rank, intrinsic_regular_orbits,
                   intrinsic_regular_points;
            printf "AS8_INTRINSIC_COMPACT_ACTION_COMPLETE_V1|%o|disposition=%o|degree=%o|rank=%o|regular_orbits=%o|regular_points=%o\n",
                   intrinsic_label, intrinsic_disposition, intrinsic_degree,
                   intrinsic_rank, intrinsic_regular_orbits,
                   intrinsic_regular_points;
        end if;
    end for;
    printf "AS8_INTRINSIC_COMPACT_BLOCK_COMPLETE_V1|sid=%o|degree=%o|G=%o|H=%o|actions=%o|fpr_density=%o|base_gt_2=%o|exact_density=%o|graph_residual=%o\n",
           ASXSimpleId, ASXSelectedDegree,
           ASXSelectedGroupOrder, ASXSelectedStabilizerOrder,
           #intrinsic_compact_records, intrinsic_fpr_count,
           intrinsic_base_count, intrinsic_density_count,
           intrinsic_residual_count;
elif ASXMode eq "intrinsicpack" then
    require #intrinsic_compact_records eq ASXExpectedActions:
        "intrinsic pack multiplicity mismatch";
    for i in [1 .. #intrinsic_compact_records] do
        record := intrinsic_compact_records[i];
        packed_G := record[1];
        packed_H := record[2];
        packed_path := Sprintf("%o/action_%o.m", ASXPackDirectory, i);
        packed_file := Open(packed_path, "w");
        fprintf packed_file, "// AS8_INTRINSIC_ACTION_PACK_V1\n";
        fprintf packed_file, "AS8PackedSimpleId := %o;\n", ASXSimpleId;
        fprintf packed_file, "AS8PackedOrdinal := %o;\n", i;
        fprintf packed_file, "AS8PackedMultiplicity := %o;\n",
                #intrinsic_compact_records;
        fprintf packed_file, "AS8PackedEnumerationLabel := \"%o\";\n",
                record[3];
        fprintf packed_file, "AS8PackedDegree := %o;\n", record[4];
        fprintf packed_file, "AS8PackedGroupOrder := %o;\n", Order(packed_G);
        fprintf packed_file, "AS8PackedStabilizerOrder := %o;\n",
                Order(packed_H);
        fprintf packed_file, "AS8PackedAmbientDegree := %o;\n",
                Degree(packed_G);
        fprintf packed_file, "AS8PackedAmbient := Sym(AS8PackedAmbientDegree);\n";
        fprintf packed_file, "AS8PackedG := sub<AS8PackedAmbient | %o>;\n",
                Setseq(Generators(packed_G));
        fprintf packed_file, "AS8PackedH := sub<AS8PackedG | %o>;\n",
                Setseq(Generators(packed_H));
        fprintf packed_file, "assert Order(AS8PackedG) eq AS8PackedGroupOrder;\n";
        fprintf packed_file, "assert Order(AS8PackedH) eq AS8PackedStabilizerOrder;\n";
        fprintf packed_file, "assert Index(AS8PackedG, AS8PackedH) eq AS8PackedDegree;\n";
        delete packed_file;
        printf "AS8_INTRINSIC_PACK_ACTION_V1|sid=%o|ordinal=%o|actions=%o|degree=%o|G=%o|H=%o|ambient_degree=%o|enumeration_label=%o|path=%o\n",
               ASXSimpleId, i, #intrinsic_compact_records, record[4],
               Order(packed_G), Order(packed_H), Degree(packed_G),
               record[3], packed_path;
    end for;
    printf "AS8_INTRINSIC_PACK_COMPLETE_V1|sid=%o|degree=%o|G=%o|H=%o|actions=%o\n",
           ASXSimpleId, ASXSelectedDegree, ASXSelectedGroupOrder,
           ASXSelectedStabilizerOrder, #intrinsic_compact_records;
elif ASXMode eq "canonicalgraphwitnessbatch" then
    require #canonical_compact_records eq compute_count:
        "canonical graph-witness action count does not match compute count";
    require ASXExpectedActions ge 1 and
            ASXBatchSelectionMask gt 0 and
            ASXBatchSelectionMask lt 2^compute_count and
            ASXGraphWitnessCertificateDirectory ne "":
        "canonical graph-witness batch needs count, mask and directory";
    graph_witness_selected_count := 0;
    for i in [1 .. #canonical_compact_records] do
        if (ASXBatchSelectionMask div 2^(i - 1)) mod 2 eq 1 then
            record := canonical_compact_records[i];
            ASXGraphWitnessCertificateFile :=
                ASXGraphWitnessCertificateDirectory cat "/" cat
                record[3] cat ".m";
            ASXEmitGraphWitnessBenchmark(
                record[1], record[2], record[3],
                ASXGraphWitnessCertificateDirectory cat "/" cat
                record[3] cat ".m");
            graph_witness_selected_count +:= 1;
        end if;
    end for;
    require graph_witness_selected_count eq ASXExpectedActions:
        "canonical graph-witness selected count mismatch";
    printf "AS8_CANONICAL_GRAPH_WITNESS_BATCH_COMPLETE_V1|sid=%o|compute_actions=%o|selected=%o\n",
           ASXSimpleId, compute_count, graph_witness_selected_count;
elif ASXMode eq "canonicalgraphwitnessverifybatch" then
    require #canonical_compact_records eq compute_count:
        "canonical graph-witness replay count does not match compute count";
    require ASXExpectedActions ge 1 and
            ASXBatchSelectionMask gt 0 and
            ASXBatchSelectionMask lt 2^compute_count and
            ASXGraphWitnessCertificateDirectory ne "":
        "canonical graph-witness replay needs count, mask and directory";
    graph_witness_replayed_count := 0;
    for i in [1 .. #canonical_compact_records] do
        if (ASXBatchSelectionMask div 2^(i - 1)) mod 2 eq 1 then
            record := canonical_compact_records[i];
            certificate_path :=
                ASXGraphWitnessCertificateDirectory cat "/" cat
                record[3] cat ".m";
            ASXVerifyGraphWitnessCertificate(
                record[1], record[2], record[3], certificate_path);
            graph_witness_replayed_count +:= 1;
        end if;
    end for;
    require graph_witness_replayed_count eq ASXExpectedActions:
        "canonical graph-witness replay selected count mismatch";
    printf "AS8_CANONICAL_GRAPH_WITNESS_REPLAY_BATCH_COMPLETE_V1|sid=%o|compute_actions=%o|selected=%o\n",
           ASXSimpleId, compute_count, graph_witness_replayed_count;
elif ASXHybridBatchMode then
    require selected_count eq compute_count and
            selected_count eq ASXExpectedActions:
        "hybrid batch action count does not match expected metadata count";
    fprintf "/dev/stderr",
            "ASX_HYBRID_COMPLETE simple_id=%o actions=%o closed=%o residual=%o\n",
            ASXSimpleId, selected_count, hybrid_closed_count,
            hybrid_residual_count;
elif ASXBatchEngineMode then
    if ASXBatchSelectionMask eq 0 then
        require selected_count eq compute_count and
                selected_count eq ASXExpectedActions:
            "unfiltered batch action count does not match metadata manifest";
    else
        require ASXBatchSelectionMask lt 2^compute_count and
                selected_count eq ASXExpectedActions:
            "filtered batch action count does not match selection mask";
    end if;
    for record in batch_records do
        ASXEmitEngine(record[1], record[2], record[3]);
    end for;
    fprintf "/dev/stderr",
            "ASX_BATCH_SELECTED_COMPLETE simple_id=%o actions=%o\n",
            ASXSimpleId, selected_count;
else
    require selected_count eq 1: "selected label was not found exactly once";
end if;
quit;
