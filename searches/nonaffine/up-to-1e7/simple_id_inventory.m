///////////////////////////////////////////////////////////////////////////
// Exact simple-group inventory under the degree-10^7 base-two order bound.
///////////////////////////////////////////////////////////////////////////

if not assigned ASXMaxDegree then ASXMaxDegree := 10000000; end if;
if Type(ASXMaxDegree) eq MonStgElt then
    ASXMaxDegree := StringToInteger(ASXMaxDegree);
end if;
require 5 le ASXMaxDegree and ASXMaxDegree le 10000000: "bad cap";
SetColumns(0);

bound := ASXMaxDegree * (ASXMaxDegree - 1);
count := 0;
psl2_count := 0;
sz_ree_count := 0;
alt_window_count := 0;
sporadic_theorem_count := 0;
print "SIMPLE_ID_INVENTORY_V1";
for simple_id in [1 .. NumberOfSimpleGroups()] do
    simple_tuple, simple_order := SimpleGroupId(simple_id);
    if simple_order gt bound then break; end if;
    count +:= 1;
    is_psl2 := simple_tuple[1] eq 1 and simple_tuple[2] eq 1;
    is_sz_ree := (simple_tuple[1] eq 11 and simple_tuple[2] eq 2) or
                 (simple_tuple[1] eq 14 and simple_tuple[2] eq 2);
    is_alt_window := simple_tuple[1] eq 17;
    is_sporadic_theorem := simple_tuple[1] eq 18 and
        simple_tuple[2] in {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                            13, 14, 15, 16, 20, 21};
    if is_psl2 then psl2_count +:= 1; end if;
    if is_sz_ree then sz_ree_count +:= 1; end if;
    if is_alt_window then alt_window_count +:= 1; end if;
    if is_sporadic_theorem then sporadic_theorem_count +:= 1; end if;
    disposition := is_psl2 select "theorem_bh_psl2" else
                   (is_sz_ree select "theorem_preprint_sz_ree" else
                   (is_alt_window select "theorem_published_alt_window" else
                   (is_sporadic_theorem select "theorem_bg_sporadic" else
                                                        "enumerate")));
    printf "SIMPLE|%o|%o|%o|%o\n",
           simple_id, simple_order, disposition, SimpleGroupName(simple_id);
end for;
printf "SIMPLE_ID_INVENTORY_COMPLETE|max_degree=%o|order_bound=%o|simple_ids=%o|psl2=%o|sz_ree=%o|alt_window=%o|sporadic_theorem=%o|enumerate=%o\n",
       ASXMaxDegree, bound, count, psl2_count, sz_ree_count,
       alt_window_count, sporadic_theorem_count,
       count - psl2_count - sz_ree_count - alt_window_count -
       sporadic_theorem_count;
quit;
