///////////////////////////////////////////////////////////////////////////
// Exact product-action component/layer inventory for 10^6 < d <= 10^7.
//
// Compound-diagonal components are handled by the separate exact CD degree
// exclusion.  Every PA component here is almost simple, has degree n <=
// floor(sqrt(10^7)) < 8192, and is therefore inside Magma's complete
// primitive-group catalogue.
///////////////////////////////////////////////////////////////////////////

if not assigned PAXMinDegree then PAXMinDegree := 1000001; end if;
if not assigned PAXMaxDegree then PAXMaxDegree := 10000000; end if;
if not assigned PAXFirstComponentDegree then PAXFirstComponentDegree := 5; end if;
if not assigned PAXLastComponentDegree then
    PAXLastComponentDegree := Floor(Sqrt(PAXMaxDegree));
end if;
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
require 5 le PAXMinDegree and PAXMinDegree le PAXMaxDegree and
        PAXMaxDegree le 10000000: "bad degree interval";
require PrimitiveGroupDatabaseLimit() ge Floor(Sqrt(PAXMaxDegree)):
    "primitive catalogue does not cover every component degree";
require 5 le PAXFirstComponentDegree and
        PAXFirstComponentDegree le PAXLastComponentDegree and
        PAXLastComponentDegree le Floor(Sqrt(PAXMaxDegree)):
    "bad component-degree interval";
SetColumns(0);

print "PA_COMPONENT_INVENTORY_V1";
component_count := 0;
layer_count := 0;
maximum_component_degree := Floor(Sqrt(PAXMaxDegree));

for n in [PAXFirstComponentDegree .. PAXLastComponentDegree] do
    if NumberOfPrimitiveAlmostSimpleGroups(n) eq 0 then continue; end if;
    candidates := [* P : P in PrimitiveGroups(n : Filter := "AlmostSimple") |
                         Order(Socle(P)) le n * (n - 1) *];
    while #candidates gt 0 do
        component_socle := Socle(candidates[1]);
        cohort_positions := [i : i in [1 .. #candidates] |
                                  Socle(candidates[i]) eq component_socle or
                                  IsConjugate(Sym(n), Socle(candidates[i]),
                                                       component_socle)];
        largest_position := cohort_positions[1];
        for position in cohort_positions do
            if Order(candidates[position]) gt Order(candidates[largest_position]) then
                largest_position := position;
            end if;
        end for;
        component_normalizer := candidates[largest_position];
        component_id, identified_degree :=
            PrimitiveGroupIdentification(component_normalizer);
        assert identified_degree eq n;
        socle_order := Order(Socle(component_normalizer));
        outer_order := Order(component_normalizer) div socle_order;
        assert socle_order le n * (n - 1);
        component_count +:= 1;

        degree := n^2;
        exponent := 2;
        while degree le PAXMaxDegree do
            if degree ge PAXMinDegree then
                socle_power := socle_order^exponent;
                maximum_quotient_order := degree * (degree - 1) div socle_power;
                if maximum_quotient_order ge exponent then
                    printf "LAYER|PA7_comp%o_%o_k%o_d%o|%o|%o|%o|%o|%o|%o|%o\n",
                           n, component_id, exponent, degree,
                           n, component_id, exponent, degree, socle_order,
                           outer_order, maximum_quotient_order;
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

printf "PA_COMPONENT_INVENTORY_COMPLETE|min_degree=%o|max_degree=%o|first_component_degree=%o|last_component_degree=%o|components=%o|layers=%o|max_component_degree=%o\n",
       PAXMinDegree, PAXMaxDegree,
       PAXFirstComponentDegree, PAXLastComponentDegree,
       component_count, layer_count, maximum_component_degree;
quit;
