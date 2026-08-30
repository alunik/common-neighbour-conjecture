///////////////////////////////////////////////////////////////////////////
// Export every order-eligible nonsoluble irreducible subgroup supplied by
// Magma's complete dimension-2/3 irreducible-subgroup constructor.
//
// Required command-line globals: ASXP, ASXD.
// The completeness guarantee used here is valid for d in {2,3} and
// characteristic at least 5.  Affine actions are represented over the prime
// field, so ASXP is required to be prime.
///////////////////////////////////////////////////////////////////////////

require assigned ASXP and assigned ASXD: "set ASXP and ASXD";
if Type(ASXP) eq MonStgElt then ASXP := StringToInteger(ASXP); end if;
if Type(ASXD) eq MonStgElt then ASXD := StringToInteger(ASXD); end if;
require IsPrime(ASXP) and ASXP ge 5: "ASXP must be a prime at least 5";
require ASXD in {2, 3}: "the direct complete constructor is used only for d=2 or d=3";

ASXDegree := ASXP^ASXD;
ASXGroups := IrreducibleSubgroups(ASXD, ASXP : Soluble := false);
ASXEligible := 0;

print "AFFINE_SAXL_V1";
for ASXIndex in [1 .. #ASXGroups] do
    ASXH := ASXGroups[ASXIndex];
    assert IsIrreducible(ASXH);
    assert not IsSolvable(ASXH);
    ASXOrder := Order(ASXH);
    if ASXOrder le ASXDegree - 1 then
        ASXEligible +:= 1;
        ASXGens := Setseq(Generators(ASXH));
        print "action";
        printf "label MagmaIrredNonSol_%o_%o_%o\n", ASXD, ASXP, ASXIndex;
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
    end if;
end for;

fprintf "/dev/stderr",
        "ASX_DIRECT_COMPLETE d=%o p=%o degree=%o nonsolvable_irreducible=%o order_eligible=%o\n",
        ASXD, ASXP, ASXDegree, #ASXGroups, ASXEligible;
quit;
