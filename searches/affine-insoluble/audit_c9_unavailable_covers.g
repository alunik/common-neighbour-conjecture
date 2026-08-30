# Complete the universal-cover gaps from enumerate_c9_base2_brauer.g without
# constructing large matrix groups.  For odd q >= 37, the cross-characteristic
# projective-degree bound for L2(q) is (q-1)/2, already larger than 12.  Covers
# of order greater than 3^12-1 fail the regular-orbit order bound instead.

LoadPackage("ctbllib");
SizeScreen( [ 10000, 10000 ] );

bound := 3^12 - 1;
crossCharacteristic := [ 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79 ];
orderObstructed := [ 83, 89, 97, 101 ];

for q in crossCharacteristic do
    simpleOrder := q * (q^2 - 1) / 2;
    coverOrder := 2 * simpleOrder;
    lowerBound := (q - 1) / 2;
    if q mod 3 = 0 or coverOrder > bound or lowerBound <= 12 then
        Error( "invalid L2 cross-characteristic bound case" );
    fi;
    Print( "C9_PROJECTIVE_DEGREE_BOUND|group=2.L2(", q,
           ")|p=3|d=12|cover_order=", coverOrder,
           "|lower_bound=", lowerBound, "|excluded=true\n" );
od;

for q in orderObstructed do
    simpleOrder := q * (q^2 - 1) / 2;
    coverOrder := 2 * simpleOrder;
    if simpleOrder > bound or coverOrder <= bound then
        Error( "invalid L2 cover-order obstruction" );
    fi;
    Print( "C9_COVER_ORDER_BOUND|group=2.L2(", q,
           ")|p=3|d=12|cover_order=", coverOrder,
           "|vector_bound=", bound, "|excluded=true\n" );
od;

# The exceptional universal cover of L3(4) is not stored, but every possible
# prime-field scalar image has centre of order at most two.  The relevant
# double-cover table is stored and its 3-Brauer table has no degree-12
# prime-field character.
tbl := CharacterTable( "2.L3(4)" );
if tbl = fail then
    Error( "missing scalar-compatible L3(4) cover table" );
fi;
modtbl := BrauerTable( tbl, 3 );
if modtbl = fail then
    Error( "missing scalar-compatible L3(4) cover table" );
fi;
hits := Filtered( Irr( modtbl ), chi ->
    chi[1] = 12 and
    List( chi, value -> GaloisCyc( value, 3 ) ) = chi
);
if not IsEmpty( hits ) then
    Error( "unexpected degree-12 module for a cover of L3(4)" );
fi;
Print( "C9_EXCEPTIONAL_COVER_AUDIT|group=2.L3(4)|p=3|d=12",
       "|brauer_table=true|prime_field_hits=0\n" );
Print( "C9_UNAVAILABLE_COVERS_COMPLETE\n" );
QUIT;
