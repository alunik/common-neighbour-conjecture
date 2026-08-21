# Export the standard ATLAS generators for 3.Suz.2 on its 24-dimensional
# GF(2)-module as packed row masks.

if LoadPackage("atlasrep") <> true then
  Error("AtlasRep is required");
fi;

F := GF(2);;
d := 24;;
expectedOrder := 2690072985600;;
info := OneAtlasGeneratingSetInfo(
  "3.Suz.2", Characteristic, 2, Dimension, d
);;
if info = fail then
  Error("ATLAS 3.Suz.2 GF(2), dimension-24 representation unavailable");
fi;
atlas := AtlasGenerators(info.identifier);;
if atlas.repname <> "3Suzd2G1-f2r24B0" then
  Error("unexpected ATLAS representation");
fi;
generators := atlas.generators;;
G := Group(generators);;
SetSize(G, expectedOrder);;
if Size(G) <> expectedOrder then
  Error("unexpected 3.Suz.2 order");
fi;
module := GModuleByMats(generators, F);;
if not MTX.IsIrreducible(module) then
  Error("the ATLAS module is reducible");
fi;

# The Table 1 row also allows external GF(4)^*-scalars.  The central element
# already inside 3.Suz.2 is not one of those scalar matrices on this module,
# so adjoining GF(4)^* really multiplies the group order by three.
kernelProgram := AtlasProgram(
  "3.Suz.2", atlas.standardization, "kernel", "Suz.2"
);;
centralGenerator := ResultOfStraightLineProgram(
  kernelProgram.program, generators
)[1];;
if Order(centralGenerator) <> 3 then
  Error("unexpected internal central-generator order");
fi;
F4 := GF(4);;
omega := Z(4);;
centralGenerator4 := ImmutableMatrix(F4, centralGenerator);;
externalScalar := omega * IdentityMat(d, F4);;
module4 := GModuleByMats(
  List(generators, generator -> ImmutableMatrix(F4, generator)), F4
);;
if not MTX.IsIrreducible(module4) or
    not MTX.IsAbsolutelyIrreducible(module4) then
  Error("the extended ATLAS module is not absolutely irreducible over GF(4)");
fi;
if centralGenerator4 = externalScalar or
    centralGenerator4 = externalScalar^2 then
  Error("internal centre was confused with an external scalar");
fi;

RowMask := function(row)
  local mask, position;
  mask := 0;
  for position in [1..d] do
    if not IsZero(row[position]) then
      mask := mask + 2^(position - 1);
    fi;
  od;
  return mask;
end;;

Print("SUZ3D2_F4_INITIAL_V1\n");
Print("ATLAS_REP ", atlas.repname, "\n");
Print("FIELD 2\n");
Print("DIM ", d, "\n");
Print("GROUP_ORDER ", Size(G), "\n");
Print("INTERNAL_CENTER_ORDER ", Order(centralGenerator), "\n");
Print("INTERNAL_CENTER_NOT_EXTERNAL_SCALAR true\n");
Print("HMAX_ORDER ", 3 * Size(G), "\n");
Print("GENERATOR_ORDERS ",
  JoinStringsWithSeparator(
    List(generators, generator -> String(Order(generator))), " "),
  "\n");
Print("GENERATOR_PRODUCT_ORDER ",
  Order(generators[1] * generators[2]), "\n");
Print("IRREDUCIBLE_OVER_GF2 true\n");
Print("ABSOLUTELY_IRREDUCIBLE_OVER_GF4 true\n");
Print("GENERATORS ", Length(generators), "\n");
for generator in generators do
  Print("MATRIX");
  for row in generator do
    Print(" ", RowMask(row));
  od;
  Print("\n");
od;
Print("END\n");
QUIT_GAP(0);
