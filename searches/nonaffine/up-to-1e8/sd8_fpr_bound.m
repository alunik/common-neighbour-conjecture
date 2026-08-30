///////////////////////////////////////////////////////////////////////////
// Exact prime-order fixed-point-ratio upper bound for one SD8 action.
//
// For G/H, the non-base probability from a fixed point is at most
//   Qhat(G,H) = sum_x |x^G cap H|^2 / |x^G|,
// over prime-order G-classes.  Qhat < 1/2 makes the Saxl neighbourhood
// larger than half the vertices, so every two vertices have a common
// neighbour.  No degree-sized permutation representation is constructed.
///////////////////////////////////////////////////////////////////////////

require assigned SD8SimpleId and assigned SD8K and assigned SD8Case:
    "set SD8SimpleId, SD8K and SD8Case";
if Type(SD8SimpleId) eq MonStgElt then SD8SimpleId := StringToInteger(SD8SimpleId); end if;
if Type(SD8K) eq MonStgElt then SD8K := StringToInteger(SD8K); end if;
if Type(SD8Case) eq MonStgElt then SD8Case := StringToInteger(SD8Case); end if;
SetColumns(0);

// simple ID, k, exact selected quotient-class count
rows := [
    <9,3,5>, <10,3,8>, <11,3,5>, <12,3,5>, <13,3,5>,
    <14,3,16>, <15,3,2>, <16,3,12>, <3,4,16>, <1,5,13>
];
matching := [row : row in rows | row[1] eq SD8SimpleId and row[2] eq SD8K];
require #matching eq 1: "selected socle/k is not in the exact SD8 window";
expected_cases := matching[1][3];

tuple, t_order := SimpleGroupId(SD8SimpleId);
A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle order mismatch";
Q, qmap := quo<A | T>;
S := Sym(SD8K);
D, dins, dprojs := DirectProduct(Q, S);
degree := t_order^(SD8K - 1);
require 10000000 lt degree and degree le 100000000: "degree outside window";

candidates := [];
for record in Subgroups(D) do
    B := record`subgroup;
    top := sub<S | [g @ dprojs[2] : g in Generators(B)]>;
    if not IsTransitive(top) then continue; end if;
    if t_order * Order(B) gt degree - 1 then continue; end if;
    is_residual := Order(top) in {Factorial(SD8K) div 2, Factorial(SD8K)};
    selected := SD8K eq 3 or (SD8SimpleId eq 1 and SD8K eq 5) or is_residual;
    if selected then Append(~candidates, B); end if;
end for;
require #candidates eq expected_cases: "exact SD8 quotient count mismatch";
require 1 le SD8Case and SD8Case le #candidates: "bad SD8Case";
B := candidates[SD8Case];
top_B := sub<S | [g @ dprojs[2] : g in Generators(B)]>;

W, bins, topin, topproj := WreathProduct(A, S);
LiftOuterTop := function(b)
    q := b @ dprojs[1];
    a := q @@ qmap;
    sigma := b @ dprojs[2];
    diagonal_a := &*[W | a @ bins[i] : i in [1 .. SD8K]];
    return diagonal_a * (sigma @ topin);
end function;

b_lifts := [LiftOuterTop(b) : b in Generators(B)];
independent_t := [t @ bins[i] : i in [1 .. SD8K], t in Generators(T)];
diagonal_t := [&*[W | t @ bins[i] : i in [1 .. SD8K]] : t in Generators(T)];
G := sub<W | independent_t cat b_lifts>;
H := sub<W | diagonal_t cat b_lifts>;
require Order(G) eq t_order^SD8K * Order(B): "bad SD8 G order";
require Order(H) eq t_order * Order(B): "bad SD8 H order";

label := Sprintf("SD8_sid%o_k%o_case%o", SD8SimpleId, SD8K, SD8Case);
print "SD8_FPR_V1";
printf "SD_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|case=%o|cases=%o|B=%o|top=%o|G=%o|H=%o\n",
       label, SD8SimpleId, SimpleGroupName(SD8SimpleId), t_order, SD8K,
       degree, SD8Case, expected_cases, Order(B), Order(top_B), Order(G), Order(H);

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
           i, intersections[i], g_class_sizes[i], Numerator(term), Denominator(term);
end for;
printf "SD8_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
       label, #prime_classes, #intersections, Numerator(bound),
       Denominator(bound), bound lt 1/2;
quit;
