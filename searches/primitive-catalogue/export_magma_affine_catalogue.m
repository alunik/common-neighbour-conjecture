///////////////////////////////////////////////////////////////////////////
// Export a stable global-index slice of all affine primitive groups of
// degree 2..8191 from Magma's primitive-group database.
//
// Required command-line globals:
//   ASXFirstAction, ASXLastAction
//
// The global index orders actions lexicographically by (degree, database ID).
// Magma numbers the affine primitive groups first at each degree; this script
// additionally verifies every exported socle is elementary abelian and
// regular before constructing its module action.
///////////////////////////////////////////////////////////////////////////

require assigned ASXFirstAction and assigned ASXLastAction:
    "set ASXFirstAction and ASXLastAction";
if Type(ASXFirstAction) eq MonStgElt then
    ASXFirstAction := StringToInteger(ASXFirstAction);
end if;
if Type(ASXLastAction) eq MonStgElt then
    ASXLastAction := StringToInteger(ASXLastAction);
end if;
require ASXFirstAction ge 1 and ASXFirstAction le ASXLastAction:
    "invalid global action-index range";

ASXTotal := 0;
for ASXDegree in [2 .. 8191] do
    ASXTotal +:= NumberOfPrimitiveAffineGroups(ASXDegree);
end for;
require ASXLastAction le ASXTotal:
    "ASXLastAction exceeds the affine catalogue size";

print "AFFINE_SAXL_V1";
ASXPosition := 0;
ASXExported := 0;
for ASXDegree in [2 .. 8191] do
    ASXCount := NumberOfPrimitiveAffineGroups(ASXDegree);
    for ASXIndex in [1 .. ASXCount] do
        ASXPosition +:= 1;
        if ASXPosition ge ASXFirstAction and ASXPosition le ASXLastAction then
            ASXG := PrimitiveGroup(ASXDegree, ASXIndex);
            ASXN := Socle(ASXG);
            assert IsElementaryAbelian(ASXN);
            assert #ASXN eq ASXDegree;

            ASXH := Stabilizer(ASXG, 1);
            ASXM, ASXToModule := GModule(ASXG, ASXN);
            ASXAction := GModuleAction(ASXM);
            ASXGens := Setseq(Generators(ASXH));

            print "action";
            printf "label Primitive_%o_%o\n", ASXDegree, ASXIndex;
            printf "p %o\n", #BaseRing(ASXM);
            printf "n %o\n", Dimension(ASXM);
            printf "order %o\n", #ASXH;
            print "orientation row";
            printf "gens %o\n", #ASXGens;
            for ASXGenerator in ASXGens do
                ASXMatrix := ASXAction(ASXGenerator);
                for ASXRow in [1 .. Nrows(ASXMatrix)] do
                    for ASXColumn in [1 .. Ncols(ASXMatrix)] do
                        printf "%o ", Integers() ! ASXMatrix[ASXRow, ASXColumn];
                    end for;
                    print "";
                end for;
            end for;
            print "end";
            ASXExported +:= 1;
        end if;
    end for;
end for;
assert ASXPosition eq ASXTotal;
assert ASXExported eq ASXLastAction - ASXFirstAction + 1;
quit;
