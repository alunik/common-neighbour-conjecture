///////////////////////////////////////////////////////////////////////////
// Export the first order-eligible nonsoluble irreducible frontier below GL.
// Exact GL-conjugacy fusion prevents the maximal-subgroup tree from
// materialising the same linear group through many different overgroups.
//
// Required command-line globals: ASXP, ASXD.
// Optional globals:
//   ASXMaxNodes  -- hard traversal guard (default 100000)
//   ASXNoFusion  -- if true, tolerate duplicate over-threshold branches
//                   instead of testing GL-conjugacy (default false)
//   ASXFingerprintFusion -- if true, compare the invariant multiset of
//                   line-orbit lengths before exact GL-conjugacy (default
//                   false)
//   ASXShallowBypass -- if true, bypass fusion for |H| <= 2(|V|-1);
//                   every proper child is then already order-eligible
//                   (default false)
//   ASXBypassFactor -- integer factor F >= 1; bypass fusion for
//                   |H| <= F(|V|-1).  Duplicate descent then has depth at
//                   most ceiling(log_2(F)).  Overrides ASXShallowBypass.
//   ASXTraceProgress -- if true, log each over-threshold node before its
//                   expensive operation (default false; diagnostics only)
//   ASXSeedMode  -- when assigned, a file loaded before this script must
//                   define ASXSeedGroups, ASXSeedLabels, and optionally
//                   ASXSkipGroups.  The proper maximal subgroups of the
//                   seeds then form the initial layer.
///////////////////////////////////////////////////////////////////////////

require assigned ASXP and assigned ASXD: "set ASXP and ASXD";
if Type(ASXP) eq MonStgElt then ASXP := StringToInteger(ASXP); end if;
if Type(ASXD) eq MonStgElt then ASXD := StringToInteger(ASXD); end if;
if assigned ASXMaxNodes then
    if Type(ASXMaxNodes) eq MonStgElt then
        ASXMaxNodes := StringToInteger(ASXMaxNodes);
    end if;
else
    ASXMaxNodes := 100000;
end if;
require IsPrime(ASXP): "ASXP must be prime";
require ASXD ge 4 and ASXD le 17: "ASXD must lie in [4,17]";
require ASXMaxNodes ge 1: "ASXMaxNodes must be positive";
if assigned ASXNoFusion and Type(ASXNoFusion) eq MonStgElt then
    ASXNoFusion := ASXNoFusion eq "true";
elif not assigned ASXNoFusion then
    ASXNoFusion := false;
end if;
if assigned ASXFingerprintFusion and
   Type(ASXFingerprintFusion) eq MonStgElt then
    ASXFingerprintFusion := ASXFingerprintFusion eq "true";
elif not assigned ASXFingerprintFusion then
    ASXFingerprintFusion := false;
end if;
if assigned ASXShallowBypass and Type(ASXShallowBypass) eq MonStgElt then
    ASXShallowBypass := ASXShallowBypass eq "true";
elif not assigned ASXShallowBypass then
    ASXShallowBypass := false;
end if;
if assigned ASXBypassFactor then
    if Type(ASXBypassFactor) eq MonStgElt then
        ASXBypassFactor := StringToInteger(ASXBypassFactor);
    end if;
else
    ASXBypassFactor := ASXShallowBypass select 2 else 1;
end if;
require ASXBypassFactor ge 1: "ASXBypassFactor must be positive";
if assigned ASXTraceProgress and Type(ASXTraceProgress) eq MonStgElt then
    ASXTraceProgress := ASXTraceProgress eq "true";
elif not assigned ASXTraceProgress then
    ASXTraceProgress := false;
end if;

ASXNodeFormat := recformat<group, label, parent>;

procedure ASXEmit(ASXH, ASXLabel)
    ASXGens := Setseq(Generators(ASXH));
    print "action";
    printf "label %o\n", ASXLabel;
    printf "p %o\n", ASXP;
    printf "n %o\n", ASXD;
    printf "order %o\n", Order(ASXH);
    print "orientation row";
    printf "gens %o\n", #ASXGens;
    for ASXGenerator in ASXGens do
        ASXMatrix := Matrix(ASXGenerator);
        for ASXRow in [1 .. ASXD] do
            for ASXColumn in [1 .. ASXD] do
                printf "%o ", Integers() ! ASXMatrix[ASXRow, ASXColumn];
            end for;
            print "";
        end for;
    end for;
    print "end";
end procedure;

ASXDegree := ASXP^ASXD;
ASXBound := ASXDegree - 1;
ASXQueue := [];
ASXSeenGroups := [];
ASXSeenOrders := [];
ASXSeenFingerprints := [];
ASXSeenParents := [];
if not assigned ASXSkipGroups then ASXSkipGroups := []; end if;
ASXRootCount := 0;
ASXReducible := 0;
ASXSoluble := 0;
ASXDuplicates := 0;
ASXSkipped := 0;
ASXShallowExpanded := 0;
ASXNextParent := 0;

if assigned ASXSeedMode then
    require assigned ASXSeedGroups and assigned ASXSeedLabels:
        "seed file must define ASXSeedGroups and ASXSeedLabels";
    require #ASXSeedGroups eq #ASXSeedLabels:
        "seed group and label counts differ";
    if not assigned ASXSkipGroups then ASXSkipGroups := []; end if;
    for ASXSeedIndex in [1 .. #ASXSeedGroups] do
        ASXNextParent +:= 1;
        ASXSeed := ASXSeedGroups[ASXSeedIndex];
        ASXMaximals := MaximalSubgroups(ASXSeed);
        for ASXIndex in [1 .. #ASXMaximals] do
            ASXChild := ASXMaximals[ASXIndex]`subgroup;
            Append(~ASXQueue,
                   rec<ASXNodeFormat |
                       group := ASXChild,
                       label := ASXSeedLabels[ASXSeedIndex] cat "m" cat
                                IntegerToString(ASXIndex),
                       parent := ASXNextParent>);
        end for;
    end for;
    fprintf "/dev/stderr",
            "ASX_FRONTIER_SEEDED d=%o p=%o seeds=%o raw_children=%o skip_groups=%o\n",
            ASXD, ASXP, #ASXSeedGroups, #ASXQueue, #ASXSkipGroups;
else
    ASXRoots := ClassicalMaximals("GL", ASXD, ASXP);
    ASXRootCount := #ASXRoots;
    for ASXIndex in [1 .. #ASXRoots] do
        Append(~ASXQueue,
               rec<ASXNodeFormat |
                   group := ASXRoots[ASXIndex],
                   label := "MagmaNonSolFrontier_" cat IntegerToString(ASXD)
                            cat "_" cat IntegerToString(ASXP) cat "_r"
                            cat IntegerToString(ASXIndex),
                   parent := 0>);
    end for;
end if;

print "AFFINE_SAXL_V1";
ASXProcessed := 0;
ASXEnqueued := #ASXQueue;
ASXExpanded := 0;
ASXExported := 0;

// A depth-first stack releases processed matrix groups immediately.  This
// is especially important when duplicate branches are deliberately retained.
while #ASXQueue gt 0 and ASXProcessed lt ASXMaxNodes do
    ASXNode := ASXQueue[#ASXQueue];
    Prune(~ASXQueue);
    ASXProcessed +:= 1;
    ASXH := ASXNode`group;
    ASXLabel := ASXNode`label;
    ASXParent := ASXNode`parent;

    if not IsIrreducible(ASXH) then
        ASXReducible +:= 1;
        continue;
    end if;
    if IsSolvable(ASXH) then
        ASXSoluble +:= 1;
        continue;
    end if;

    ASXOrder := Order(ASXH);
    // Once the regular-orbit order bound has been reached, exporting a
    // repeated GL-conjugate is vastly cheaper than proving the conjugacy.
    // Duplicate graph tests do not affect coverage or correctness.  Exact
    // GL-fusion remains essential only above the bound, where it prevents
    // repeated structural descent through the same large group.
    if ASXOrder le ASXBound then
        ASXEmit(ASXH, ASXLabel);
        ASXExported +:= 1;
        continue;
    end if;
    if ASXTraceProgress then
        fprintf "/dev/stderr", "ASX_NODE position=%o queued=%o order=%o parent=%o label=%o\n",
                ASXProcessed, #ASXQueue, ASXOrder, ASXParent, ASXLabel;
    end if;

    // A proper subgroup has index at least two.  Thus when |H| <= 2B,
    // every proper child has order at most B and is emitted rather than
    // expanded.  Skipping conjugacy here can introduce only one layer of
    // duplicate graph tests, never a duplicate descent tree.
    if ASXBypassFactor gt 1 and ASXOrder le ASXBypassFactor * ASXBound then
        ASXMaximals := MaximalSubgroups(ASXH);
        ASXExpanded +:= 1;
        ASXShallowExpanded +:= 1;
        ASXNextParent +:= 1;
        for ASXIndex in [1 .. #ASXMaximals] do
            Append(~ASXQueue,
                   rec<ASXNodeFormat |
                       group := ASXMaximals[ASXIndex]`subgroup,
                       label := ASXLabel cat "m" cat IntegerToString(ASXIndex),
                       parent := ASXNextParent>);
        end for;
        ASXEnqueued +:= #ASXMaximals;
        continue;
    end if;

    if not ASXNoFusion then
        ASXFingerprint := [];
        ASXDuplicate := false;
        for ASXSeenIndex in [1 .. #ASXSeenGroups] do
            if ASXSeenParents[ASXSeenIndex] ne ASXParent and
               ASXSeenOrders[ASXSeenIndex] eq ASXOrder then
                ASXFingerprintMatch := true;
                if ASXFingerprintFusion then
                    if #ASXFingerprint eq 0 then
                        ASXFingerprint :=
                            [ASXOrbit[1] : ASXOrbit in OrbitsOfSpaces(ASXH, 1)];
                        Sort(~ASXFingerprint);
                    end if;
                    if #ASXSeenFingerprints[ASXSeenIndex] eq 0 then
                        ASXSeenFingerprint :=
                            [ASXOrbit[1] : ASXOrbit in
                                OrbitsOfSpaces(ASXSeenGroups[ASXSeenIndex], 1)];
                        Sort(~ASXSeenFingerprint);
                        ASXSeenFingerprints[ASXSeenIndex] := ASXSeenFingerprint;
                    end if;
                    ASXFingerprintMatch :=
                        ASXSeenFingerprints[ASXSeenIndex] eq ASXFingerprint;
                end if;
                if ASXFingerprintMatch and
                   IsGLConjugate(ASXH, ASXSeenGroups[ASXSeenIndex]) then
                    ASXDuplicate := true;
                    break;
                end if;
            end if;
        end for;
        if ASXDuplicate then
            ASXDuplicates +:= 1;
            continue;
        end if;
        Append(~ASXSeenGroups, ASXH);
        Append(~ASXSeenOrders, ASXOrder);
        Append(~ASXSeenFingerprints, ASXFingerprint);
        Append(~ASXSeenParents, ASXParent);
    end if;

    ASXMaximals := MaximalSubgroups(ASXH);
    ASXExpanded +:= 1;
    ASXNextParent +:= 1;
    for ASXIndex in [1 .. #ASXMaximals] do
        Append(~ASXQueue,
               rec<ASXNodeFormat |
                   group := ASXMaximals[ASXIndex]`subgroup,
                   label := ASXLabel cat "m" cat IntegerToString(ASXIndex),
                   parent := ASXNextParent>);
    end for;
    ASXEnqueued +:= #ASXMaximals;
end while;

ASXComplete := #ASXQueue eq 0;
fprintf "/dev/stderr",
        "ASX_FRONTIER_COMPLETE d=%o p=%o degree=%o roots=%o processed=%o queued=%o expanded=%o shallow_expanded=%o reducible=%o soluble=%o duplicates=%o skipped=%o exported=%o fusion=%o fingerprint_fusion=%o bypass_factor=%o complete=%o\n",
        ASXD, ASXP, ASXDegree, ASXRootCount, ASXProcessed, ASXEnqueued,
        ASXExpanded, ASXShallowExpanded, ASXReducible, ASXSoluble,
        ASXDuplicates, ASXSkipped, ASXExported, not ASXNoFusion,
        ASXFingerprintFusion, ASXBypassFactor, ASXComplete;
require ASXComplete: "ASXMaxNodes guard reached before frontier completion";
quit;
