///////////////////////////////////////////////////////////////////////////
// Export the non-geometric Aschbacher-C9 modules that survive the necessary
// affine base-two order bound in dimensions 12--19.
//
// Required command-line globals: ASXP, ASXD.
// The exhaustive Brauer-character sieve is implemented independently in
// enumerate_c9_base2_brauer.g.  Its only non-self-dual survivors are the two
// degree-12 modules for 2.L2(23) over F_3 and the two degree-15 modules for
// L2(31) over F_2.  Self-dual survivors preserve a nondegenerate form and
// therefore belong to C8, which the geometric campaign already covers.
///////////////////////////////////////////////////////////////////////////

require assigned ASXP and assigned ASXD: "set ASXP and ASXD";
if Type(ASXP) eq MonStgElt then ASXP := StringToInteger(ASXP); end if;
if Type(ASXD) eq MonStgElt then ASXD := StringToInteger(ASXD); end if;
require <ASXP, ASXD> in {<3, 12>, <2, 13>, <2, 14>, <2, 15>, <2, 16>, <2, 17>,
                          <2, 18>, <2, 19>}:
    "unsupported C9 supplement target";

procedure ASXEmit(ASXH, ASXLabel, ASXP, ASXD)
    assert IsIrreducible(ASXH);
    ASXOrder := Order(ASXH);
    assert ASXOrder le ASXP^ASXD - 1;
    ASXGens := Setseq(Generators(ASXH));
    print "action";
    printf "label %o\n", ASXLabel;
    printf "p %o\n", ASXP;
    printf "n %o\n", ASXD;
    printf "order %o\n", ASXOrder;
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

SetSeed(1);
print "AFFINE_SAXL_V1";
ASXExported := 0;

if ASXP eq 3 and ASXD eq 12 then
    ASXSource := SL(2, 23);
    ASXModules := AbsolutelyIrreducibleModules(ASXSource, GF(3));
    ASXModules := [ASXM : ASXM in ASXModules |
        Dimension(ASXM) eq 12 and #BaseRing(ASXM) eq 3];
    assert #ASXModules eq 2;
    // The two dual Brauer characters are interchanged by the diagonal outer
    // automorphism, so their image groups form one GL-conjugacy class; see
    // audit_c9_outer_fusion.g.
    ASXH := ActionGroup(ASXModules[1]);
    assert Order(ASXH) eq 12144;
    ASXEmit(ASXH, "C9_2L2_23_d12_p3", ASXP, ASXD);
    ASXExported +:= 1;
elif ASXP eq 2 and ASXD eq 15 then
    ASXSource := PSL(2, 31);
    ASXModules := AbsolutelyIrreducibleModules(ASXSource, GF(2));
    ASXModules := [ASXM : ASXM in ASXModules |
        Dimension(ASXM) eq 15 and #BaseRing(ASXM) eq 2];
    assert #ASXModules eq 2;
    ASXH := ActionGroup(ASXModules[1]);
    assert Order(ASXH) eq 14880;
    ASXEmit(ASXH, "C9_L2_31_d15_p2", ASXP, ASXD);
    ASXExported +:= 1;
end if;

fprintf "/dev/stderr",
    "ASX_C9_BASE2_COMPLETE d=%o p=%o degree=%o exported=%o\n",
    ASXD, ASXP, ASXP^ASXD, ASXExported;
quit;
