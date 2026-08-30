///////////////////////////////////////////////////////////////////////////
// Exact prime-order FPR bounds for every component quotient class at one
// fixed pair (T,k).  Keeping the complete Subgroups(D) enumeration in one
// process avoids any dependence on representative ordering across Magma
// processes.
///////////////////////////////////////////////////////////////////////////

require assigned CD18SimpleId and assigned CD18K and
        assigned CD18ExpectedCases and assigned CD18ExpectedOrder and
        assigned CD18ExpectedDegree: "set all CD18 pair task parameters";
if Type(CD18SimpleId) eq MonStgElt then CD18SimpleId := StringToInteger(CD18SimpleId); end if;
if Type(CD18K) eq MonStgElt then CD18K := StringToInteger(CD18K); end if;
if Type(CD18ExpectedCases) eq MonStgElt then CD18ExpectedCases := StringToInteger(CD18ExpectedCases); end if;
if Type(CD18ExpectedOrder) eq MonStgElt then CD18ExpectedOrder := StringToInteger(CD18ExpectedOrder); end if;
if Type(CD18ExpectedDegree) eq MonStgElt then CD18ExpectedDegree := StringToInteger(CD18ExpectedDegree); end if;
SetColumns(0);

tuple, t_order := SimpleGroupId(CD18SimpleId);
require t_order eq CD18ExpectedOrder: "simple order mismatch";
A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle mismatch";
Q, qmap := quo<A | T>;
S := Sym(CD18K);
D, dins, dprojs := DirectProduct(Q, S);
degree := t_order^(CD18K - 1);
require degree eq CD18ExpectedDegree: "component degree mismatch";

candidates := [];
for record in Subgroups(D) do
    B := record`subgroup;
    top := sub<S | [g @ dprojs[2] : g in Generators(B)]>;
    if IsTransitive(top) then Append(~candidates, B); end if;
end for;
require #candidates eq CD18ExpectedCases: "component quotient census mismatch";

W, bins, topin, topproj := WreathProduct(A, S);
LiftOuterTop := function(b)
    q := b @ dprojs[1];
    a := q @@ qmap;
    sigma := b @ dprojs[2];
    diagonal_a := &*[W | a @ bins[i] : i in [1 .. CD18K]];
    return diagonal_a * (sigma @ topin);
end function;

label := Sprintf("CD18CP_sid%o_k%o", CD18SimpleId, CD18K);
print "CD18_COMPONENT_PAIR_FPR_V1";
printf "CD18_COMPONENT_PAIR|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|cases=%o|Out=%o\n",
       label, CD18SimpleId, SimpleGroupName(CD18SimpleId), t_order,
       CD18K, degree, #candidates, Order(Q);

for case_number in [1 .. #candidates] do
    B := candidates[case_number];
    top_B := sub<S | [g @ dprojs[2] : g in Generators(B)]>;
    outer_image := sub<Q | [g @ dprojs[1] : g in Generators(B)]>;
    top_kernel := B meet Kernel(dprojs[2]);
    b_lifts := [LiftOuterTop(b) : b in Generators(B)];
    independent_t := [t @ bins[i] : i in [1 .. CD18K], t in Generators(T)];
    diagonal_t := [&*[W | t @ bins[i] : i in [1 .. CD18K]] : t in Generators(T)];
    G := sub<W | independent_t cat b_lifts>;
    H := sub<W | diagonal_t cat b_lifts>;
    require Order(G) eq t_order^CD18K * Order(B): "bad component G order";
    require Order(H) eq t_order * Order(B): "bad component H order";

    printf "CD18_COMPONENT_ACTION|%o|case=%o|B=%o|top=%o|outer_image=%o|top_kernel=%o|G=%o|H=%o\n",
           label, case_number, Order(B), Order(top_B), Order(outer_image),
           Order(top_kernel), Order(G), Order(H);
    hclasses := Classes(H);
    prime_classes := [c : c in hclasses | IsPrime(c[1])];
    g_representatives := [W | ];
    intersections := [];
    g_class_sizes := [];
    for c in prime_classes do
        x := W!c[3];
        position := 0;
        for i in [1 .. #g_representatives] do
            if IsConjugate(G, x, g_representatives[i]) then
                position := i;
                break;
            end if;
        end for;
        if position eq 0 then
            Append(~g_representatives, x);
            Append(~intersections, c[2]);
            Append(~g_class_sizes, Order(G) div Order(Centralizer(G, x)));
        else
            intersections[position] +:= c[2];
        end if;
    end for;
    bound := &+[Rationals() | intersections[i]^2 / g_class_sizes[i]
                : i in [1 .. #intersections]];
    for sequence in [1 .. #intersections] do
        term := Rationals()!(intersections[sequence]^2) / g_class_sizes[sequence];
        printf "FPR_CLASS|case=%o|sequence=%o|intersection=%o|G_class=%o|term_num=%o|term_den=%o\n",
               case_number, sequence, intersections[sequence], g_class_sizes[sequence],
               Numerator(term), Denominator(term);
    end for;
    printf "CD18_COMPONENT_CASE_COMPLETE|%o|case=%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
           label, case_number, #prime_classes, #intersections, Numerator(bound),
           Denominator(bound), bound lt 1/2;
end for;
printf "CD18_COMPONENT_PAIR_COMPLETE|%o|cases=%o\n", label, #candidates;
quit;
