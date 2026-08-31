///////////////////////////////////////////////////////////////////////////
// Exact catalogue-side PA component/layer inventory for
// 10^7 < d <= 10^8.  Components 8192..10000 are constructed separately by
// export_pa8_seam_structured.m; this source is fail-closed at the catalogue
// boundary.
///////////////////////////////////////////////////////////////////////////

if not assigned PAXMinDegree then PAXMinDegree := 10000001; end if;
if not assigned PAXMaxDegree then PAXMaxDegree := 100000000; end if;
if not assigned PAXFirstComponentDegree then PAXFirstComponentDegree := 5; end if;
if not assigned PAXLastComponentDegree then PAXLastComponentDegree := 8191; end if;
if Type(PAXMinDegree) eq MonStgElt then
    PAXMinDegree := StringToInteger(PAXMinDegree);
end if;
if Type(PAXMaxDegree) eq MonStgElt then
    PAXMaxDegree := StringToInteger(PAXMaxDegree);
end if;
if Type(PAXFirstComponentDegree) eq MonStgElt then
    PAXFirstComponentDegree := StringToInteger(PAXFirstComponentDegree);
end if;
if Type(PAXLastComponentDegree) eq MonStgElt then
    PAXLastComponentDegree := StringToInteger(PAXLastComponentDegree);
end if;
require PAXMinDegree eq 10000001 and PAXMaxDegree eq 100000000:
    "this source is bound to the new degree decade";
require PrimitiveGroupDatabaseLimit() ge 8191:
    "primitive catalogue does not cover its pinned side of the seam";
require 5 le PAXFirstComponentDegree and
        PAXFirstComponentDegree le PAXLastComponentDegree and
        PAXLastComponentDegree le 8191:
    "bad catalogue component-degree interval";
SetColumns(0);

print "PA8_CATALOGUE_COMPONENT_INVENTORY_V1";
component_count := 0;
layer_count := 0;

// Reject component degrees for which no product power n^k, k >= 2, lies in
// the target decade.  This arithmetic test must precede every catalogue call:
// in particular 1029 <= n <= 3162 has n^2 < 10^7 < 10^8 < n^3, so the
// primitive-group database cannot contribute a PA layer there.
PAXHasPowerInWindow := function(n)
    degree := n^2;
    while degree le PAXMaxDegree do
        if degree ge PAXMinDegree then return true; end if;
        if degree gt PAXMaxDegree div n then break; end if;
        degree *:= n;
    end while;
    return false;
end function;

for n in [PAXFirstComponentDegree .. PAXLastComponentDegree] do
    if not PAXHasPowerInWindow(n) then continue; end if;
    if NumberOfPrimitiveAlmostSimpleGroups(n) eq 0 then continue; end if;
    candidates := [* P : P in PrimitiveGroups(n : Filter := "AlmostSimple") |
                         Order(Socle(P)) le n * (n - 1) *];
    while #candidates gt 0 do
        component_socle := Socle(candidates[1]);
        cohort_positions := [i : i in [1 .. #candidates] |
            Socle(candidates[i]) eq component_socle or
            IsConjugate(Sym(n), Socle(candidates[i]), component_socle)];
        largest_position := cohort_positions[1];
        for position in cohort_positions do
            if Order(candidates[position]) gt
               Order(candidates[largest_position]) then
                largest_position := position;
            end if;
        end for;
        component_normalizer := candidates[largest_position];
        component_id, identified_degree :=
            PrimitiveGroupIdentification(component_normalizer);
        assert identified_degree eq n;
        socle_order := Order(Socle(component_normalizer));
        outer_order := Order(component_normalizer) div socle_order;
        metadata_method := IsPrimitive(Socle(component_normalizer)) select
            "quotient_first_primitive_socle" else
            "quotient_first_component_and_top";
        assert socle_order le n * (n - 1);
        component_count +:= 1;

        degree := n^2;
        exponent := 2;
        while degree le PAXMaxDegree do
            if degree ge PAXMinDegree then
                socle_power := socle_order^exponent;
                maximum_quotient_order :=
                    degree * (degree - 1) div socle_power;
                if maximum_quotient_order ge exponent then
                    printf "LAYER|PA8_comp%o_%o_k%o_d%o|%o|%o|%o|%o|%o|%o|%o|%o\n",
                           n, component_id, exponent, degree,
                           n, component_id, exponent, degree, socle_order,
                           outer_order, maximum_quotient_order,
                           metadata_method;
                    layer_count +:= 1;
                end if;
            end if;
            if degree gt PAXMaxDegree div n then break; end if;
            degree *:= n;
            exponent +:= 1;
        end while;

        for position in Reverse(Sort(cohort_positions)) do
            Remove(~candidates, position);
        end for;
    end while;
end for;

printf "PA8_CATALOGUE_COMPONENT_INVENTORY_COMPLETE|min_degree=%o|max_degree=%o|first_component_degree=%o|last_component_degree=%o|components=%o|layers=%o|catalogue_boundary=8191\n",
       PAXMinDegree, PAXMaxDegree, PAXFirstComponentDegree,
       PAXLastComponentDegree, component_count, layer_count;
quit;
