#############################################################################
## Count the IRREDSOL actions, order-eligible actions, and immediate order
## obstructions in every guardian root covered by build_irredsol_roots.g.
#############################################################################

if LoadPackage("irredsol") = fail then Error("IRREDSOL is required"); fi;

ASX_Lower := 8192;
ASX_Upper := 2^24 - 1;
ASX_Specs := [];
for ASX_p in Filtered([2 .. RootInt(ASX_Upper, 2)], IsPrimeInt) do
    ASX_n := 2;
    while ASX_p^ASX_n <= ASX_Upper do
        if ASX_p^ASX_n >= ASX_Lower then
            Add(ASX_Specs, [ASX_p^ASX_n, ASX_n, ASX_p]);
        fi;
        ASX_n := ASX_n + 1;
    od;
od;
Sort(ASX_Specs);

Print("root_index\tdegree\tn\tp\td\tguardian\tactions\torder_eligible\torder_obstructed\n");
ASX_RootIndex := 0;
for ASX_Spec in ASX_Specs do
    ASX_Degree := ASX_Spec[1];
    ASX_n := ASX_Spec[2];
    ASX_p := ASX_Spec[3];
    for ASX_d in DivisorsInt(ASX_n) do
        ASX_Indices := IndicesIrreducibleSolubleMatrixGroups(ASX_n, ASX_p, ASX_d);
        if Length(ASX_Indices) > 0 then
            ASX_m := ASX_n / ASX_d;
            ASX_q := ASX_p^ASX_d;
            LoadAbsolutelyIrreducibleSolubleGroupData(ASX_m, ASX_q);
            if ASX_m = 1 then
                ASX_Guardians := [1];
            else
                ASX_Guardians := Set(List(ASX_Indices,
                    ASX_k -> IRREDSOL_DATA.GROUPS[ASX_m][ASX_q][ASX_k][1]));
            fi;
            for ASX_Guardian in ASX_Guardians do
                ASX_RootIndex := ASX_RootIndex + 1;
                if ASX_m = 1 then
                    ASX_RootIndices := ASX_Indices;
                    ASX_Orders := List(ASX_RootIndices,
                        ASX_k -> IRREDSOL_DATA.GROUPS_DIM1[ASX_q][ASX_k][1]);
                else
                    ASX_RootIndices := Filtered(ASX_Indices,
                        ASX_k -> IRREDSOL_DATA.GROUPS[ASX_m][ASX_q][ASX_k][1]
                                 = ASX_Guardian);
                    ASX_Orders := List(ASX_RootIndices,
                        ASX_k -> IRREDSOL_DATA.GROUPS[ASX_m][ASX_q][ASX_k][4]);
                fi;
                ASX_Eligible := Number(ASX_Orders, ASX_o -> ASX_o <= ASX_Degree - 1);
                Print(ASX_RootIndex, "\t", ASX_Degree, "\t", ASX_n, "\t",
                      ASX_p, "\t", ASX_d, "\t", ASX_Guardian, "\t",
                      Length(ASX_RootIndices), "\t", ASX_Eligible, "\t",
                      Length(ASX_RootIndices) - ASX_Eligible, "\n");
            od;
        fi;
    od;
od;
quit;
