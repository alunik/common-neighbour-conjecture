#############################################################################
# Build the exact list of non-affine primitive actions in the extended GAP
# PrimGrp catalogue.  Output is TSV: ordinal, degree, primitive ID, order,
# O'Nan--Scott type.  Degrees 4096..8191 require the official extended data.
#############################################################################

if LoadPackage("primgrp") = fail then Error("PrimGrp is required"); fi;
SetPrintFormattingStatus("*stdout*", false);

if not IsBound(ASX_MinDegree) then ASX_MinDegree := 4096; fi;
if not IsBound(ASX_MaxDegree) then ASX_MaxDegree := 8191; fi;
if ASX_MinDegree < 2 or ASX_MaxDegree > 8191 or ASX_MinDegree > ASX_MaxDegree then
    Error("invalid degree interval");
fi;

ASX_Ordinal := 0;
for ASX_Degree in [ASX_MinDegree .. ASX_MaxDegree] do
    if not PrimitiveGroupsAvailable(ASX_Degree) then
        Error("primitive groups unavailable at degree ", ASX_Degree);
    fi;
    for ASX_Id in [1 .. NrPrimitiveGroups(ASX_Degree)] do
        if ASX_Degree > 4095 then
            # Extended-data records put the O'Nan--Scott type in their header,
            # before the (potentially enormous) generator list.  Skip affine
            # records after reading only that header; parse the full official
            # GAP record only for the non-affine actions retained below.
            ASX_Filename := PrimGrpArtifactFilename(ASX_Degree, ASX_Id);
            ASX_Stream := InputTextFile(ASX_Filename);
            if ASX_Stream = fail then Error("cannot open primitive-group record"); fi;
            ASX_TypeLine := ReadLine(ASX_Stream);
            while ASX_TypeLine <> fail and
                  PositionSublist(ASX_TypeLine, "ONanScottType") = fail do
                ASX_TypeLine := ReadLine(ASX_Stream);
            od;
            CloseStream(ASX_Stream);
            if ASX_TypeLine = fail then Error("missing O'Nan--Scott header"); fi;
            ASX_IsAffine := PositionSublist(
                ASX_TypeLine, "ONanScottType := \"1\""
            ) <> fail;
            if not ASX_IsAffine then
                ASX_Data := PRIMGrp(ASX_Degree, ASX_Id);
                if ASX_Data[1] <> ASX_Id then Error("primitive ID mismatch"); fi;
                ASX_Type := ASX_Data[4];
                ASX_Order := ASX_Data[2];
            fi;
        else
            ASX_Group := PrimitiveGroup(ASX_Degree, ASX_Id);
            ASX_Type := ONanScottType(ASX_Group);
            ASX_Order := Size(ASX_Group);
            ASX_IsAffine := ASX_Type = 1 or ASX_Type = "1";
        fi;
        if not ASX_IsAffine then
            if ASX_Type = 1 or ASX_Type = "1" then Error("affine header mismatch"); fi;
            ASX_Ordinal := ASX_Ordinal + 1;
            Print(ASX_Ordinal, "\t", ASX_Degree, "\t", ASX_Id, "\t",
                  ASX_Order, "\t", ASX_Type, "\n");
        fi;
    od;
od;
PrintTo("/dev/stderr", "NONAFFINE_MANIFEST_COMPLETE actions=", ASX_Ordinal,
        " degrees=", ASX_MinDegree, "..", ASX_MaxDegree, "\n");
QUIT;
