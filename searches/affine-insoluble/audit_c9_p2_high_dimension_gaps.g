# Complete the missing 2-Brauer-table cases for dimensions 18 and 19.  For
# odd q >= 41, the cross-characteristic projective-degree lower bound for
# L2(q) is (q-1)/2 > 19.  The boundary q=37 and the defining-characteristic
# group L2(64) are constructed exactly in audit_c9_missing_brauer.m.

SizeScreen( [ 10000, 10000 ] );
qs := [ 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 ];
for q in qs do
    lowerBound := (q - 1) / 2;
    if q mod 2 = 0 or lowerBound <= 19 then
        Error( "invalid high-dimensional L2 bound case" );
    fi;
    Print( "C9_P2_PROJECTIVE_DEGREE_BOUND|group=L2(", q,
           ")|p=2|dimensions=[18,19]|lower_bound=", lowerBound,
           "|excluded=true\n" );
od;
Print( "C9_P2_HIGH_DIMENSION_GAPS_COMPLETE\n" );
QUIT;
