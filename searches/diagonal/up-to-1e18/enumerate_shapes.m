///////////////////////////////////////////////////////////////////////////
// Exact arithmetic SD shapes for 10^8 < degree <= 10^18.
//
// The SimpleGroupId list is ordered by group order.  An SD action with
// socle T^k has degree |T|^(k-1).  These are shapes, not quotient/action
// counts.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);

lower := 100000000;
upper := 1000000000000000000;
simple_order_limit := 1000000000; // sqrt(upper), since k >= 3

sd_count := 0;
simple_count := 0;
max_k := 0;

sid := 1;
while true do
    tuple, t_order := SimpleGroupId(sid);
    if t_order gt simple_order_limit then break; end if;
    simple_count +:= 1;
    name := SimpleGroupName(sid);
    printf "SIMPLE|sid=%o|name=%o|order=%o\n", sid, name, t_order;

    degree := t_order^2;
    k := 3;
    while degree le upper do
        if degree gt lower then
            sd_count +:= 1;
            max_k := Maximum(max_k, k);
            printf "SD_SHAPE|sid=%o|name=%o|order=%o|k=%o|degree=%o\n",
                   sid, name, t_order, k, degree;
        end if;

        if degree gt upper div t_order then break; end if;
        degree *:= t_order;
        k +:= 1;
    end while;
    sid +:= 1;
end while;

printf "SHAPE_SUMMARY|simple_groups=%o|sd_shapes=%o|max_k=%o\n",
       simple_count, sd_count, max_k;
assert simple_count eq 277;
assert sd_count eq 357;
assert max_k eq 11;
quit;
