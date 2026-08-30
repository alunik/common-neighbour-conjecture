///////////////////////////////////////////////////////////////////////////
// Export one exact layer of the nonsoluble irreducible subgroup tree.
//
// Required command-line globals: ASXP, ASXD.
// In root mode, emit the relevant ClassicalMaximals("GL", d, p) roots.
// In seed mode, a file loaded first must define ASXSeedGroups and
// ASXSeedLabels; emit the relevant proper maximal subgroups of every seed.
// The Python driver handles the order bound, Saxl tests, and subsequent
// layers.  Keeping each invocation to one layer allows bounded seed batches
// and returns all Magma heap state to the OS between calls.
///////////////////////////////////////////////////////////////////////////

require assigned ASXP and assigned ASXD: "set ASXP and ASXD";
if Type(ASXP) eq MonStgElt then ASXP := StringToInteger(ASXP); end if;
if Type(ASXD) eq MonStgElt then ASXD := StringToInteger(ASXD); end if;
require IsPrime(ASXP): "ASXP must be prime";
require ASXD ge 4 and ASXD le 19: "ASXD must lie in [4,19]";

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

print "AFFINE_SAXL_V1";
ASXParents := 0;
ASXRawChildren := 0;
ASXReducible := 0;
ASXSoluble := 0;
ASXExported := 0;

if assigned ASXSeedMode then
    require assigned ASXSeedGroups and assigned ASXSeedLabels:
        "seed file must define ASXSeedGroups and ASXSeedLabels";
    require #ASXSeedGroups eq #ASXSeedLabels:
        "seed group and label counts differ";
    ASXParents := #ASXSeedGroups;
    for ASXSeedIndex in [1 .. #ASXSeedGroups] do
        ASXMaximals := MaximalSubgroups(ASXSeedGroups[ASXSeedIndex]);
        ASXRawChildren +:= #ASXMaximals;
        for ASXIndex in [1 .. #ASXMaximals] do
            ASXChild := ASXMaximals[ASXIndex]`subgroup;
            if not IsIrreducible(ASXChild) then
                ASXReducible +:= 1;
                continue;
            end if;
            if IsSolvable(ASXChild) then
                ASXSoluble +:= 1;
                continue;
            end if;
            ASXEmit(ASXChild,
                    ASXSeedLabels[ASXSeedIndex] cat "m" cat
                    IntegerToString(ASXIndex));
            ASXExported +:= 1;
        end for;
    end for;
else
    ASXRoots := ClassicalMaximals("GL", ASXD, ASXP);
    ASXParents := 1;
    ASXRawChildren := #ASXRoots;
    for ASXIndex in [1 .. #ASXRoots] do
        ASXRoot := ASXRoots[ASXIndex];
        if not IsIrreducible(ASXRoot) then
            ASXReducible +:= 1;
            continue;
        end if;
        if IsSolvable(ASXRoot) then
            ASXSoluble +:= 1;
            continue;
        end if;
        ASXEmit(ASXRoot,
                "MagmaNonSolLayer_" cat IntegerToString(ASXD) cat "_" cat
                IntegerToString(ASXP) cat "_r" cat IntegerToString(ASXIndex));
        ASXExported +:= 1;
    end for;
end if;

fprintf "/dev/stderr",
        "ASX_LAYER_COMPLETE d=%o p=%o parents=%o raw_children=%o reducible=%o soluble=%o exported=%o\n",
        ASXD, ASXP, ASXParents, ASXRawChildren, ASXReducible,
        ASXSoluble, ASXExported;
quit;
