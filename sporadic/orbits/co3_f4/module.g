# Export the standard ATLAS Co3 generators on the irreducible
# 22-dimensional GF(2)-module as packed row masks.
#
# Usage:
#   gap -q module.g > module.txt

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;
if LoadPackage("ctbllib") <> true then
  Error("CTblLib is required");
fi;

F := GF(2);;
d := 22;;
expected_order := 495766656000;;
info := OneAtlasGeneratingSetInfo(
  "Co3", Characteristic, 2, Dimension, d
);;
if info = fail then
  Error("ATLAS Co3 GF(2), dimension-22 representation is unavailable");
fi;
atlas := AtlasGenerators(info.identifier);;
if atlas.repname <> "Co3G1-f2r22B0" then
  Error("unexpected ATLAS representation");
fi;
gens := atlas.generators;;
G := Group(gens);;
if Size(G) <> expected_order then
  Error("unexpected Co3 order");
fi;
module := GModuleByMats(gens, F);;
if not MTX.IsIrreducible(module) then
  Error("the ATLAS module is not irreducible");
fi;
if not MTX.IsAbsolutelyIrreducible(module) then
  Error("the ATLAS module is not absolutely irreducible");
fi;
modular_table := CharacterTable("Co3") mod 2;;
if Length(Filtered(Irr(modular_table), chi -> chi[1] = d)) <> 1 then
  Error("the degree-22 irreducible Brauer character is not unique");
fi;

RowMask := function(row)
  local mask, i;
  mask := 0;
  for i in [1 .. d] do
    if not IsZero(row[i]) then
      mask := mask + 2^(i - 1);
    fi;
  od;
  return mask;
end;;

Print("CO3_F4_INITIAL_V1\n");
Print("ATLAS_REP ", atlas.repname, "\n");
Print("CHARACTER ", atlas.charactername, "\n");
Print("FIELD 2\n");
Print("DIM ", d, "\n");
Print("GROUP_ORDER ", Size(G), "\n");
Print("GENERATOR_ORDERS ",
  JoinStringsWithSeparator(List(gens, g -> String(Order(g))), " "), "\n");
Print("GENERATOR_PRODUCT_ORDER ", Order(gens[1] * gens[2]), "\n");
Print("ABSOLUTELY_IRREDUCIBLE true\n");
Print("UNIQUE_DEGREE_22_BRAUER_IRREDUCIBLE true\n");
Print("GENERATORS ", Length(gens), "\n");
for g in gens do
  Print("MATRIX");
  for row in g do
    Print(" ", RowMask(row));
  od;
  Print("\n");
od;
Print("END\n");

QUIT_GAP(0);
