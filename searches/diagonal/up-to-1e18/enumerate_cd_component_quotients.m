///////////////////////////////////////////////////////////////////////////
// Exact quotient-class inventory for every SD component occurring in a CD
// shape with 10^8 < degree <= 10^18.
//
// For component socle T^k the diagonal normalizer quotient is
// Out(T) x S_k.  Conjugacy-class representatives B with transitive S_k
// projection give the induced primitive component actions.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);
lower := 100000000;
upper := 1000000000000000000;
simple_order_limit := 1000000000;

pairs := {@ @};
sid := 1;
while true do
    tuple, t_order := SimpleGroupId(sid);
    if t_order gt simple_order_limit then break; end if;
    component_degree := t_order^2;
    k := 3;
    while component_degree^2 le upper do
        compound_degree := component_degree^2;
        occurs := false;
        while compound_degree le upper do
            if compound_degree gt lower then occurs := true; end if;
            if compound_degree gt upper div component_degree then break; end if;
            compound_degree *:= component_degree;
        end while;
        if occurs then Include(~pairs, <sid, k>); end if;
        if component_degree gt upper div t_order then break; end if;
        component_degree *:= t_order;
        k +:= 1;
    end while;
    sid +:= 1;
end while;

total_classes := 0;
CanonicalSubgroupKey := function(B)
    rows := Sort([Eltseq(g) : g in Set(B)]);
    return Sprint(rows);
end function;
for pair in pairs do
    sid := pair[1];
    k := pair[2];
    tuple, t_order := SimpleGroupId(sid);
    A := AutomorphismGroupSimpleGroup(tuple);
    T := Socle(A);
    Q, qmap := quo<A | T>;
    S := Sym(k);
    D, dins, dprojs := DirectProduct(Q, S);
    component_degree := t_order^(k - 1);
    candidates := [];
    candidate_keys := [];
    for record in Subgroups(D) do
        B := record`subgroup;
        top := sub<S | [g @ dprojs[2] : g in Generators(B)]>;
        if not IsTransitive(top) then continue; end if;
        key := CanonicalSubgroupKey(B);
        position := 1;
        while position le #candidate_keys and candidate_keys[position] lt key do
            position +:= 1;
        end while;
        require position gt #candidate_keys or candidate_keys[position] ne key:
            "duplicate canonical subgroup key";
        Insert(~candidate_keys, position, key);
        Insert(~candidates, position, B);
    end for;
    count := 0;
    for B in candidates do
        top := sub<S | [g @ dprojs[2] : g in Generators(B)]>;
        outer_image := sub<Q | [g @ dprojs[1] : g in Generators(B)]>;
        top_kernel := B meet Kernel(dprojs[2]);
        // Huang's theorem gives base size two for every retained component;
        // retain the elementary necessary-order assertion independently.
        assert t_order * Order(B) le component_degree - 1;
        count +:= 1;
        printf "CD_COMPONENT_CLASS|sid=%o|name=%o|order=%o|k=%o|component_degree=%o|case=%o|B=%o|top=%o|outer_image=%o|top_kernel=%o|Out=%o\n",
               sid, SimpleGroupName(sid), t_order, k, component_degree,
               count, Order(B), Order(top), Order(outer_image),
               Order(top_kernel), Order(Q);
    end for;
    total_classes +:= count;
    printf "CD_COMPONENT_PAIR|sid=%o|name=%o|order=%o|k=%o|component_degree=%o|classes=%o|Out=%o\n",
           sid, SimpleGroupName(sid), t_order, k, component_degree, count,
           Order(Q);
end for;
printf "CD_COMPONENT_SUMMARY|pairs=%o|classes=%o\n", #pairs, total_classes;
quit;
