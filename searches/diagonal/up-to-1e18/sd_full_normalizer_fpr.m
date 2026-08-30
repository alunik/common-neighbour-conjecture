///////////////////////////////////////////////////////////////////////////
// Exact prime-order FPR union bound for the full SD normalizer
//   N = T^k.(Out(T) x S_k)
// on Diag(T)\T^k, without constructing its degree-sized action.
//
// If Qhat < 1/2, every two Saxl vertices for N have a common neighbour.
// Since every primitive diagonal group with the same (T,k) is a subgroup
// of N on the same point set, the closure is inherited by all of them.
///////////////////////////////////////////////////////////////////////////

require assigned SD18SimpleId and assigned SD18K:
    "set SD18SimpleId and SD18K";
if Type(SD18SimpleId) eq MonStgElt then
    SD18SimpleId := StringToInteger(SD18SimpleId);
end if;
if Type(SD18K) eq MonStgElt then
    SD18K := StringToInteger(SD18K);
end if;
SetColumns(0);

lower := 100000000;
upper := 1000000000000000000;
tuple, t_order := SimpleGroupId(SD18SimpleId);
degree := t_order^(SD18K - 1);
if assigned SD18ComponentMode then
    if Type(SD18ComponentMode) eq MonStgElt then
        SD18ComponentMode := SD18ComponentMode eq "true";
    end if;
else
    SD18ComponentMode := false;
end if;
if SD18ComponentMode then
    require degree^2 le upper and degree^2 gt lower:
        "component does not occur in the exact CD window";
else
    require lower lt degree and degree le upper: "shape outside exact SD window";
end if;
require 3 le SD18K and SD18K le 11: "invalid SD18 k";

A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle order mismatch";
Q, qmap := quo<A | T>;
S := Sym(SD18K);

W, bins, topin, topproj := WreathProduct(A, S);
outer_lifts := [];
for q in Generators(Q) do
    a := q @@ qmap;
    Append(~outer_lifts,
           &*[W | a @ bins[i] : i in [1 .. SD18K]]);
end for;
top_lifts := [W | s @ topin : s in Generators(S)];
independent_t := [t @ bins[i] : i in [1 .. SD18K],
                  t in Generators(T)];
diagonal_t := [&*[W | t @ bins[i] : i in [1 .. SD18K]]
               : t in Generators(T)];
G := sub<W | independent_t cat outer_lifts cat top_lifts>;
H := sub<W | diagonal_t cat outer_lifts cat top_lifts>;

expected_g := t_order^SD18K * Order(Q) * Factorial(SD18K);
expected_h := t_order * Order(Q) * Factorial(SD18K);
require Order(G) eq expected_g: "bad full-normalizer order";
require Order(H) eq expected_h: "bad point-stabilizer order";

label := Sprintf("SD18_sid%o_k%o", SD18SimpleId, SD18K);
print "SD18_FULL_FPR_V1";
printf "SD18_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|Out=%o|G=%o|H=%o\n",
       label, SD18SimpleId, SimpleGroupName(SD18SimpleId), t_order,
       SD18K, degree, Order(Q), Order(G), Order(H);

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
for i in [1 .. #intersections] do
    term := Rationals()!(intersections[i]^2) / g_class_sizes[i];
    printf "FPR_CLASS|sequence=%o|intersection=%o|G_class=%o|term_num=%o|term_den=%o\n",
           i, intersections[i], g_class_sizes[i], Numerator(term),
           Denominator(term);
end for;
printf "SD18_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
       label, #prime_classes, #intersections, Numerator(bound),
       Denominator(bound), bound lt 1/2;
quit;
