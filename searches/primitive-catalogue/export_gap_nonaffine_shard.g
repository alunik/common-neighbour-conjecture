#############################################################################
# Export the manifest interval ASX_FIRST_LINE..ASX_LAST_LINE from GAP's
# extended primitive catalogue to PRIMITIVE_SAXL_V1.  GAP constructs groups
# and stabilisers; C++ computes regular suborbits and the Saxl graph.
#############################################################################

if LoadPackage("primgrp") = fail then Error("PrimGrp is required"); fi;
SetPrintFormattingStatus("*stdout*", false);

if not IsBound(ASX_Manifest) then Error("set ASX_Manifest"); fi;
if not IsBound(ASX_First) or not IsBound(ASX_Last) or
   ASX_First < 1 or ASX_Last < ASX_First then
    Error("bad manifest-line interval");
fi;

ASX_PrintPermutation := function(ASX_Generator, ASX_Degree)
    local ASX_Point;
    for ASX_Point in [1 .. ASX_Degree] do
        Print(ASX_Point ^ ASX_Generator, " ");
    od;
    Print("\n");
end;

ASX_Stream := InputTextFile(ASX_Manifest);
if ASX_Stream = fail then Error("cannot open manifest"); fi;
Print("PRIMITIVE_SAXL_V1\n");
ASX_LineNumber := 0;
ASX_Exported := 0;
while true do
    ASX_Line := ReadLine(ASX_Stream);
    if ASX_Line = fail then break; fi;
    ASX_LineNumber := ASX_LineNumber + 1;
    if ASX_LineNumber < ASX_First then continue; fi;
    if ASX_LineNumber > ASX_Last then break; fi;
    ASX_Fields := SplitString(Chomp(ASX_Line), "\t");
    if Length(ASX_Fields) < 4 then Error("bad manifest row"); fi;
    ASX_Ordinal := Int(ASX_Fields[1]);
    ASX_Degree := Int(ASX_Fields[2]);
    ASX_Id := Int(ASX_Fields[3]);
    if ASX_Ordinal <> ASX_LineNumber then Error("manifest ordinal mismatch"); fi;

    ASX_GroupOrder := Int(ASX_Fields[4]);
    if ASX_GroupOrder = fail or ASX_GroupOrder mod ASX_Degree <> 0 then
        Error("bad group order in manifest");
    fi;
    ASX_StabilizerOrder := ASX_GroupOrder / ASX_Degree;
    if ASX_StabilizerOrder = 1 then
        ASX_Classification := "base1";
    elif ASX_StabilizerOrder > ASX_Degree - 1 then
        ASX_Classification := "order_obstruction";
    else
        ASX_Classification := "compute";
    fi;

    Print("action\n");
    Print("label GAPPrimitive_", ASX_Degree, "_", ASX_Id, "\n");
    Print("degree ", ASX_Degree, "\n");
    Print("stabilizer_order ", ASX_StabilizerOrder, "\n");
    Print("classification ", ASX_Classification, "\n");
    Print("regular_orbits 0\nregular_count 0\n");
    if ASX_Classification = "compute" then
        # The manifest-generating pass already certified every non-affine
        # action.  Reconstruct only the cases that survive the order bound.
        ASX_Group := PrimitiveGroup(ASX_Degree, ASX_Id);
        ASX_Type := ONanScottType(ASX_Group);
        if ASX_Type = 1 or ASX_Type = "1" then
            Error("affine action in non-affine manifest");
        fi;
        if Size(ASX_Group) <> ASX_GroupOrder then
            Error("group-order manifest mismatch");
        fi;
        ASX_Stabilizer := Stabilizer(ASX_Group, 1);
        if Size(ASX_Stabilizer) <> ASX_StabilizerOrder then
            Error("stabilizer-order mismatch");
        fi;
        ASX_HGenerators := GeneratorsOfGroup(ASX_Stabilizer);
        Print("hgens ", Length(ASX_HGenerators), "\n");
        for ASX_Generator in ASX_HGenerators do
            ASX_PrintPermutation(ASX_Generator, ASX_Degree);
        od;
        ASX_Generators := GeneratorsOfGroup(ASX_Group);
        Print("gens ", Length(ASX_Generators), "\n");
        for ASX_Generator in ASX_Generators do
            ASX_PrintPermutation(ASX_Generator, ASX_Degree);
        od;
    fi;
    Print("end\n");
    ASX_Exported := ASX_Exported + 1;
od;
CloseStream(ASX_Stream);
if ASX_Exported <> ASX_Last - ASX_First + 1 then
    Error("wrong number of exported actions");
fi;
PrintTo("/dev/stderr", "GAP_EXPORT_COMPLETE first=", ASX_First,
        " last=", ASX_Last, " actions=", ASX_Exported, "\n");
QUIT;
