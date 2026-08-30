// Exact Magma completion of the small CTblLib Brauer-table gaps reported by
// enumerate_c9_base2_brauer.g.  Each relevant absolutely irreducible module
// is constructed, and the requested prime-field dimensions are asserted
// absent.

SetSeed(1);

procedure AuditMissing(ASXQ, ASXP, ASXDimensions, ASXUseCover)
    ASXG := ASXUseCover select SL(2, ASXQ) else PSL(2, ASXQ);
    ASXModules := AbsolutelyIrreducibleModules(ASXG, GF(ASXP));
    ASXPrimeField := [ASXM : ASXM in ASXModules |
        #BaseRing(ASXM) eq ASXP and Dimension(ASXM) gt 1];
    ASXHits := [Dimension(ASXM) : ASXM in ASXPrimeField |
        Dimension(ASXM) in ASXDimensions];
    assert IsEmpty(ASXHits);
    ASXMinimum := IsEmpty(ASXPrimeField) select 0 else
        Min([Dimension(ASXM) : ASXM in ASXPrimeField]);
    printf "C9_MISSING_BRAUER_AUDIT|q=%o|p=%o|cover=%o|order=%o|dimensions=%o|prime_field_modules=%o|min_nontrivial_dimension=%o|hits=0\n",
           ASXQ, ASXP, ASXUseCover, Order(ASXG), ASXDimensions,
           #ASXPrimeField, ASXMinimum;
end procedure;

// Missing characteristic-2 Brauer tables.  A central involution acts
// trivially on an irreducible module, so the simple quotient is sufficient.
AuditMissing(37, 2, {15, 16, 17, 18, 19}, false);
AuditMissing(41, 2, {16, 17}, false);
AuditMissing(43, 2, {16, 17}, false);
AuditMissing(47, 2, {16, 17}, false);
AuditMissing(53, 2, {17}, false);
AuditMissing(59, 2, {17}, false);
AuditMissing(61, 2, {17}, false);
AuditMissing(64, 2, {18, 19}, false);

// The 2.L2(81) table exists but its 3-Brauer table does not.  This is the
// only defining-characteristic gap; unavailable cross-characteristic covers
// are excluded uniformly in audit_c9_unavailable_covers.g.
AuditMissing(81, 3, {12}, true);
printf "C9_MISSING_BRAUER_AUDIT_COMPLETE\n";
