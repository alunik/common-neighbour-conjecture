# Check the outer-automorphism fusion for the two non-self-dual C9 pairs.
# Each pair induces to one irreducible character of twice the degree, so the
# two modules form one outer orbit and hence one GL-conjugacy class of image
# subgroups.  The absence of a same-degree character for the extension also
# proves that the linear normaliser does not contain the outer automorphism.

LoadPackage("ctbllib");
SizeScreen( [ 10000, 10000 ] );

AuditPair := function( subgroupName, extensionName, p, d )
    local sub, ext, positions, induced, extensionSameDegree, pos;
    sub := BrauerTable( CharacterTable( subgroupName ), p );
    ext := BrauerTable( CharacterTable( extensionName ), p );
    positions := Filtered( [ 1 .. Length( Irr( sub ) ) ], pos ->
        Irr( sub )[pos][1] = d and
        List( Irr( sub )[pos], value -> GaloisCyc( value, p ) ) =
            Irr( sub )[pos]
    );
    induced := List( positions,
        pos -> InducedClassFunction( Irr( sub )[pos], ext ) );
    extensionSameDegree := Filtered( Irr( ext ), chi -> chi[1] = d );
    if Length( positions ) <> 2 or
       induced[1] <> induced[2] or
       induced[1][1] <> 2*d or
       Position( Irr( ext ), induced[1] ) = fail or
       Length( extensionSameDegree ) <> 0 then
        Error( "unexpected C9 outer fusion" );
    fi;
    Print( "C9_OUTER_FUSION|subgroup=", subgroupName,
           "|extension=", extensionName,
           "|p=", p, "|d=", d,
           "|subgroup_positions=", positions,
           "|induced_position=", Position( Irr( ext ), induced[1] ),
           "|induced_degree=", induced[1][1],
           "|extension_degree_d_characters=", Length( extensionSameDegree ),
           "|gl_classes=1|normaliser_index=1\n" );
end;

AuditPair( "2.L2(23)", "2.L2(23).2", 3, 12 );
AuditPair( "L2(31)", "L2(31).2", 2, 15 );
Print( "C9_OUTER_FUSION_COMPLETE\n" );
QUIT;
