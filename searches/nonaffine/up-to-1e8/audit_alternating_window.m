///////////////////////////////////////////////////////////////////////////
// Exact necessary-order audit for alternating socles in
// 10^7 < degree <= 10^8.
//
// The eight surviving maximal-action rows all have H meet A_n primitive on
// the natural n points.  Burness--Giudici, Theorem 5.1, therefore supplies
// the full common-neighbour conclusion conditionally on base size two.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);
lower := 10^7;
upper := 10^8;
assert Order(AlternatingGroup(19)) gt upper * (upper - 1);

expected := [
    <13, 1, 39916800, 78, true, true>,
    <13, 2, 39916800, 156, true, true>,
    <14, 1, 39916800, 1092, true, true>,
    <14, 2, 39916800, 2184, true, true>,
    <15, 1, 32432400, 20160, true, true>,
    <15, 1, 32432400, 20160, true, true>,
    <16, 1, 32432400, 322560, true, true>,
    <16, 1, 32432400, 322560, true, true>
];
observed := [];

for n in [5 .. 18] do
    for kind in [1 .. 2] do
        G := kind eq 1 select AlternatingGroup(n) else SymmetricGroup(n);
        for record in MaximalSubgroups(G : IndexLimit := upper) do
            H := record`subgroup;
            degree := Index(G, H);
            if degree le lower or degree gt upper or Order(H) gt degree - 1 then
                continue;
            end if;
            row := <n, kind, degree, Order(H), IsTransitive(H), IsPrimitive(H)>;
            Append(~observed, row);
            printf "ALT8_WINDOW|n=%o|group=%o|degree=%o|H_order=%o|natural_transitive=%o|natural_primitive=%o\n",
                   n, kind eq 1 select "A" else "S", degree, Order(H),
                   IsTransitive(H), IsPrimitive(H);
        end for;
    end for;
end for;

Sort(~observed);
Sort(~expected);
assert observed eq expected;
assert forall{row : row in observed | row[5] and row[6]};
printf "ALT8_WINDOW_AUDIT_PASS|rows=%o|necessary_condition=H_order_le_degree_minus_1|all_natural_stabilizers_primitive=true\n",
       #observed;
quit;
