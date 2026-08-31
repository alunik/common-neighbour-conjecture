///////////////////////////////////////////////////////////////////////////
// Independent exact SD fixed-point-ratio preflight for the full normalizer.
//
// This factorized Aut(T) x S_k calculation still constructs the ambient group
// and asks Magma for every G-centralizer.  It removes only impossible
// conjugacy tests: when
// the top permutation fixes a coordinate, that coordinate pins the Aut(T)
// class, so its H-class cannot fuse with a previously listed H-class.  When
// there is no fixed coordinate, fusion is still decided by IsConjugate in G.
///////////////////////////////////////////////////////////////////////////

require assigned SD24SimpleId and assigned SD24K:
    "set SD24SimpleId and SD24K";
if Type(SD24SimpleId) eq MonStgElt then
    SD24SimpleId := StringToInteger(SD24SimpleId);
end if;
if Type(SD24K) eq MonStgElt then
    SD24K := StringToInteger(SD24K);
end if;
SetColumns(0);

lower := 1000000000000000000;
upper := 1000000000000000000000000;
tuple, t_order := SimpleGroupId(SD24SimpleId);
degree := t_order^(SD24K - 1);
require lower lt degree and degree le upper: "shape outside exact SD window";
require 3 le SD24K and SD24K le 14: "invalid SD24 k";

A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle order mismatch";
Q, qmap := quo<A | T>;
S := Sym(SD24K);

W, bins, topin, topproj := WreathProduct(A, S);
outer_lifts := [];
for q in Generators(Q) do
    a := q @@ qmap;
    Append(~outer_lifts, &*[W | a @ bins[i] : i in [1 .. SD24K]]);
end for;
top_lifts := [W | s @ topin : s in Generators(S)];
independent_t := [t @ bins[i] : i in [1 .. SD24K],
                  t in Generators(T)];
G := sub<W | independent_t cat outer_lifts cat top_lifts>;

expected_g := t_order^SD24K * Order(Q) * Factorial(SD24K);
expected_h := Order(A) * Factorial(SD24K);
require Order(G) eq expected_g: "bad full-normalizer order";
require expected_h eq t_order * Order(Q) * Factorial(SD24K):
    "bad point-stabilizer order";

label := Sprintf("SD24F_sid%o_k%o", SD24SimpleId, SD24K);
print "SD24_FACTORIZED_FPR_V1";
printf "SD24F_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|Out=%o|G=%o|H=%o\n",
       label, SD24SimpleId, SimpleGroupName(SD24SimpleId), t_order,
       SD24K, degree, Order(Q), Order(G), expected_h;

aclasses := Classes(A);
sclasses := Classes(S);
g_representatives := [W | ];
representative_top_class := [Integers() | ];
intersections := [];
g_class_sizes := [];
h_prime_classes := 0;
for ca in aclasses do
    for cs_index in [1 .. #sclasses] do
        cs := sclasses[cs_index];
        element_order := LCM(ca[1], cs[1]);
        if not IsPrime(element_order) then
            continue;
        end if;
        h_prime_classes +:= 1;
        diagonal_a := &*[W | ca[3] @ bins[i] : i in [1 .. SD24K]];
        x := diagonal_a * (W!(cs[3] @ topin));
        fixed_coordinates := #[i : i in [1 .. SD24K] | i^cs[3] eq i];
        position := 0;
        if fixed_coordinates eq 0 then
            for i in [1 .. #g_representatives] do
                if representative_top_class[i] eq cs_index and
                   IsConjugate(G, x, g_representatives[i]) then
                    position := i;
                    break;
                end if;
            end for;
        end if;
        h_class_size := ca[2] * cs[2];
        if position eq 0 then
            Append(~g_representatives, x);
            Append(~representative_top_class, cs_index);
            Append(~intersections, h_class_size);
            Append(~g_class_sizes, Order(G) div Order(Centralizer(G, x)));
        else
            intersections[position] +:= h_class_size;
        end if;
    end for;
end for;

bound := &+[Rationals() | intersections[i]^2 / g_class_sizes[i]
            : i in [1 .. #intersections]];
for i in [1 .. #intersections] do
    term := Rationals()!(intersections[i]^2) / g_class_sizes[i];
    printf "FPR_CLASS|sequence=%o|intersection=%o|G_class=%o|term_num=%o|term_den=%o\n",
           i, intersections[i], g_class_sizes[i], Numerator(term),
           Denominator(term);
end for;
printf "SD24F_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
       label, h_prime_classes, #intersections, Numerator(bound),
       Denominator(bound), bound lt 1/2;
quit;
