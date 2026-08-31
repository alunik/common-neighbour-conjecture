///////////////////////////////////////////////////////////////////////////
// Exact metadata audit for almost-simple groups with alternating socle in
// the new degree window 10^6 < d <= 10^7.
//
// This does not calculate a base.  It applies only the necessary condition
// |H| <= d-1 for a base-two coset action G/H, then identifies the natural
// block system in every surviving row.  Morris--Spiga, J. Algebra 587
// (2021), 569--593, Theorems 1.1 and 1.2, give base size 3 for both
// surviving actions (partitions into five parts of size three).
///////////////////////////////////////////////////////////////////////////

SetColumns(0);
degree_lower := 10^6;
degree_upper := 10^7;
version, release, patch := GetVersion();
printf "MAGMA_VERSION|%o.%o-%o\n", version, release, patch;

// If b(G,H)=2, then |G|=d|H| <= d(d-1).  Since |A_17| already
// exceeds the upper-window bound, only A_n with 5 <= n <= 16 can occur.
assert Order(AlternatingGroup(17)) gt degree_upper * (degree_upper - 1);

FindFiveTripleBlocks := function(H)
    for block in Subsets({1 .. 15}, 3) do
        if 1 notin block then
            continue;
        end if;
        block_orbit := Orbit(H, block);
        if #block_orbit ne 5 then
            continue;
        end if;
        blocks := Setseq(block_orbit);
        if #(&join blocks) ne 15 then
            continue;
        end if;
        if exists{i : i in [1 .. 5] |
                  exists{j : j in [i + 1 .. 5] |
                         #(blocks[i] meet blocks[j]) ne 0}} then
            continue;
        end if;
        return true;
    end for;
    return false;
end function;

rows := 0;
for n in [5 .. 16] do
    for kind in [1 .. 2] do
        G := kind eq 1 select AlternatingGroup(n) else SymmetricGroup(n);
        group_label := kind eq 1 select Sprintf("A%o", n)
                                        else Sprintf("S%o", n);
        for maximal_record in MaximalSubgroups(G : IndexLimit := degree_upper) do
            H := maximal_record`subgroup;
            d := Index(G, H);
            if d le degree_lower or d gt degree_upper or Order(H) gt d - 1 then
                continue;
            end if;

            // The audit is expected to leave precisely the five-triples
            // partition stabilizer in A_15 and S_15.
            assert n eq 15;
            assert d eq Factorial(15) div (Factorial(3)^5 * Factorial(5));
            expected_h_order := Factorial(3)^5 * Factorial(5);
            if kind eq 1 then
                expected_h_order div:= 2;
            end if;
            assert Order(H) eq expected_h_order;
            assert IsTransitive(H) and not IsPrimitive(H);
            assert FindFiveTripleBlocks(H);

            printf "ALT_WINDOW|G=%o|degree=%o|H_order=%o|natural_transitive=%o|natural_primitive=%o|blocks=5|block_size=3\n",
                   group_label, d, Order(H), IsTransitive(H), IsPrimitive(H);
            rows +:= 1;
        end for;
    end for;
end for;

assert rows eq 2;
print "ALT_WINDOW_AUDIT_PASS|rows=2|necessary_condition=H_order_le_degree_minus_1";
quit;
