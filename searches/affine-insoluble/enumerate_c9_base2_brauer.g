# Necessary-condition sieve for missing Aschbacher-C9 affine stabilisers.
#
# A base-two affine action V:H requires a regular H-orbit on V, hence
# |H| <= |V|-1.  For a C9 group the quasisimple layer E is faithful modulo
# scalars, so |E| has the same bound.  We enumerate the relevant modular
# irreducible characters before constructing any matrix groups.

LoadPackage("ctbllib");
LoadPackage("AtlasRep");
SizeScreen( [ 10000, 10000 ] );

IsPrimeFieldCharacter := function( chi, p )
    return List( chi, value -> GaloisCyc( value, p ) ) = chi;
end;

CentreSizeFromTable := function( tbl )
    local size;
    size := Size( tbl );
    return Number( SizesCentralizers( tbl ), value -> value = size );
end;

CharacterKernelSize := function( tbl, chi )
    local positions, classes;
    positions := ClassPositionsOfKernel( chi );
    classes := SizesConjugacyClasses( tbl );
    return Sum( positions, pos -> classes[pos] );
end;

ScanTarget := function( p, d )
    local bound, names, name, tbl, centre, modtbl, irr, pos, chi, kernel,
          atlas, atlasCount, selfDual, eligibleTables, missingBrauer,
          candidates, simpleOrder, imageOrder, imageCentre, simpleGroups,
          simpleOrders, simpleNames, libraryOrders, simpleName, covers,
          unavailableCovers;
    bound := p^d - 1;
    # This is the completeness gate.  AllSmallNonabelianSimpleGroups is GAP's
    # classification-based enumerator, independent of CTblLib.  Compare its
    # order multiset with the simple tables, then scan the universal cover of
    # each simple factor.  Every representation of every quasisimple central
    # cover is a representation of the universal cover with central kernel.
    simpleGroups := AllSmallNonabelianSimpleGroups( [ 1 .. bound ] );
    simpleOrders := SortedList( List( simpleGroups, Size ) );
    simpleNames := AllCharacterTableNames(
        IsSimple, true, IsAbelian, false, IsDuplicateTable, false,
        Size, [ 1 .. bound ]
    );
    libraryOrders := SortedList(
        List( simpleNames, name -> Size( CharacterTable( name ) ) ) );
    if simpleOrders <> libraryOrders then
        Error( "CTblLib does not cover the classified simple factors" );
    fi;
    SortBy( simpleNames,
        name -> [ Size( CharacterTable( name ) ), name ] );
    unavailableCovers := 0;
    for simpleName in simpleNames do
        covers := AllCharacterTableNames(
            Identifier, simpleName, OfThose, SchurCover );
        if Length( covers ) <> 1 then
            Error( "ambiguous universal-cover identifier" );
        fi;
        if CharacterTable( covers[1] ) = fail then
            unavailableCovers := unavailableCovers + 1;
            Print( "UNAVAILABLE_UNIVERSAL_COVER|p=", p, "|d=", d,
                   "|simple_group=", simpleName,
                   "|simple_order=", Size( CharacterTable( simpleName ) ),
                   "|cover=", covers[1], "\n" );
        fi;
    od;
    # Scan every available quasisimple table, not just one cover per simple
    # factor.  This includes all stored central quotients and isoclinism
    # types.  The unavailable scalar-compatible covers are completed by the
    # exact Magma audit in audit_c9_missing_brauer.m.
    names := AllCharacterTableNames(
        IsQuasisimple, true, IsDuplicateTable, false );
    names := Filtered( names, function( name )
        local table, centreSize;
        table := CharacterTable( name );
        centreSize := CentreSizeFromTable( table );
        return Size( table ) / centreSize <= bound;
    end );
    SortBy( names, name -> [ Size( CharacterTable( name ) ), name ] );
    Print( "SIMPLE_SCOPE|p=", p, "|d=", d, "|bound=", bound,
           "|classified_simple_groups=", Length( simpleGroups ),
           "|simple_tables=", Length( simpleNames ),
           "|unavailable_universal_covers=", unavailableCovers, "\n" );
    Print( "TARGET|p=", p, "|d=", d, "|bound=", bound,
           "|available_quasisimple_tables=", Length( names ), "\n" );
    eligibleTables := 0;
    missingBrauer := 0;
    candidates := 0;
    for name in names do
        tbl := CharacterTable( name );
        centre := CentreSizeFromTable( tbl );
        simpleOrder := Size( tbl ) / centre;
        eligibleTables := eligibleTables + 1;
        modtbl := BrauerTable( tbl, p );
        if modtbl = fail then
            missingBrauer := missingBrauer + 1;
            Print( "MISSING_BRAUER|p=", p, "|d=", d,
                   "|group=", name, "|order=", Size( tbl ),
                   "|centre=", centre, "\n" );
            continue;
        fi;
        irr := Irr( modtbl );
        for pos in [ 1 .. Length( irr ) ] do
            chi := irr[pos];
            if chi[1] <> d or not IsPrimeFieldCharacter( chi, p ) then
                continue;
            fi;
            kernel := CharacterKernelSize( modtbl, chi );
            # A nontrivial irreducible representation of a quasisimple group
            # has central kernel.  Work with its faithful image: this also
            # catches a 2.S image represented via a larger universal cover.
            if centre mod kernel <> 0 then
                Error( "noncentral kernel in quasisimple character" );
            fi;
            imageOrder := Size( tbl ) / kernel;
            imageCentre := centre / kernel;
            if imageOrder > bound or imageCentre > p - 1 then
                continue;
            fi;
            selfDual := List( chi, ComplexConjugate ) = chi;
            atlas := AllAtlasGeneratingSetInfos(
                name, IsMatrixGroup, true,
                Ring, GF(p), Dimension, d
            );
            atlasCount := Length( atlas );
            candidates := candidates + 1;
            Print( "CANDIDATE|p=", p, "|d=", d,
                   "|group=", name, "|order=", Size( tbl ),
                   "|centre=", centre, "|kernel=", kernel,
                   "|image_order=", imageOrder,
                   "|image_centre=", imageCentre,
                   "|brauer_position=", pos,
                   "|self_dual=", selfDual,
                   "|atlas_prime_field_representations=", atlasCount,
                   "\n" );
        od;
    od;
    Print( "TARGET_COMPLETE|p=", p, "|d=", d,
           "|eligible_quasisimple_tables=", eligibleTables,
           "|missing_brauer_tables=", missingBrauer,
           "|candidate_characters=", candidates, "\n" );
end;

if IsBound( C9OnlyP ) and IsBound( C9OnlyD ) then
    ScanTarget( C9OnlyP, C9OnlyD );
else
    ScanTarget( 3, 12 );
    for d in [ 13 .. 19 ] do
        ScanTarget( 2, d );
    od;
fi;

QUIT;
