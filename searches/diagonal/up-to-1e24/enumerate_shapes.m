///////////////////////////////////////////////////////////////////////////
// Exact new simple-diagonal (SD) arithmetic shapes for
// 10^18 < degree <= 10^24.
//
// SimpleGroupId is ordered by group order.  Since k >= 3, only simple
// groups of order at most sqrt(10^24)=10^12 can occur.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);

lower := 1000000000000000000;
upper := 1000000000000000000000000;
simple_order_limit := 1000000000000;

sd_count := 0;
simple_count := 0;
max_k := 0;

sid := 1;
while true do
    tuple, t_order := SimpleGroupId(sid);
    if t_order gt simple_order_limit then
        printf "CATALOGUE_CUTOFF|next_sid=%o|next_order=%o|limit=%o\n",
               sid, t_order, simple_order_limit;
        break;
    end if;
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

assert simple_count eq 1650;
assert sd_count eq 1556;
assert max_k eq 14;
printf "SHAPE_SUMMARY|simple_groups=%o|sd_shapes=%o|max_k=%o\n",
       simple_count, sd_count, max_k;
quit;
