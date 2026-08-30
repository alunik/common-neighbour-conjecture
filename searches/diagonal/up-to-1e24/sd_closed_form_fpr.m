///////////////////////////////////////////////////////////////////////////
// Exact prime-order FPR union bound for the full SD normalizer, using the
// closed centralizer and fusion formula for constant-diagonal elements.
//
// Let x=(a,...,a)pi in D=Aut(T) x S_k have prime order p.  Write pi as r
// p-cycles and f fixed points.  If f>0, then
//
// |C_G(x)| = |C_Sk(pi)| |T|^r |C_T(a)|^f
//              |C_A(a)T/T|.
//
// If f=0, the fixed-coordinate constraint disappears and
//
// |C_G(x)| = |C_Sk(pi)| |T|^r |C_Out(T)(aT)|.
//
// With f>0 the D-class does not fuse with a different A-class.  With f=0,
// precisely the A-classes having conjugate images in Out(T) fuse.  Thus no
// wreath-product construction or large permutation-group centralizer is
// required.
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
require lower lt degree and degree le upper: "shape outside exact SD24 window";
require 3 le SD24K and SD24K le 14: "invalid SD24 k";

A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle order mismatch";
Q, qmap := quo<A | T>;
S := Sym(SD24K);
expected_g := t_order^SD24K * Order(Q) * Factorial(SD24K);
expected_h := Order(A) * Factorial(SD24K);

label := Sprintf("SD24C_sid%o_k%o", SD24SimpleId, SD24K);
print "SD24_CLOSED_FORM_FPR_V1";
printf "SD24C_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|Out=%o|G=%o|H=%o\n",
       label, SD24SimpleId, SimpleGroupName(SD24SimpleId), t_order,
       SD24K, degree, Order(Q), expected_g, expected_h;

aclasses := Classes(A);
sclasses := Classes(S);
qclasses := Classes(Q);
intersections := [Integers() | ];
g_class_sizes := [Integers() | ];
keys := [PowerSequence(Integers()) | ];
h_prime_classes := 0;

for ai in [1 .. #aclasses] do
    ca := aclasses[ai];
    a := ca[3];
    for si in [1 .. #sclasses] do
        cs := sclasses[si];
        element_order := LCM(ca[1], cs[1]);
        if not IsPrime(element_order) then continue; end if;
        h_prime_classes +:= 1;
        f := 0;
        r := 0;
        for cycle in CycleStructure(cs[3]) do
            if cycle[1] eq 1 then
                f := cycle[2];
            elif cycle[1] eq element_order then
                r := cycle[2];
            else
                require false: "unexpected prime-order cycle structure";
            end if;
        end for;
        require f + element_order*r eq SD24K: "cycle census mismatch";
        cs_centralizer := Factorial(SD24K) div cs[2];
        if f gt 0 then
            ct := Order(Centralizer(T, a));
            ca_centralizer := Order(A) div ca[2];
            require ca_centralizer mod ct eq 0: "centralizer image mismatch";
            image_centralizer := ca_centralizer div ct;
            cg := cs_centralizer * t_order^r * ct^f * image_centralizer;
            key := [1, si, ai];
        else
            alpha := a @ qmap;
            qi := 0;
            for j in [1 .. #qclasses] do
                if IsConjugate(Q, alpha, qclasses[j][3]) then qi := j; break; end if;
            end for;
            require qi ne 0: "outer fusion class missing";
            outer_centralizer := Order(Q) div qclasses[qi][2];
            cg := cs_centralizer * t_order^r * outer_centralizer;
            key := [0, si, qi];
        end if;
        require expected_g mod cg eq 0: "centralizer does not divide G";
        g_class_size := expected_g div cg;
        position := Index(keys, key);
        h_class_size := ca[2] * cs[2];
        if position eq 0 then
            Append(~keys, key);
            Append(~intersections, h_class_size);
            Append(~g_class_sizes, g_class_size);
        else
            require g_class_sizes[position] eq g_class_size:
                    "fused class-size mismatch";
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
printf "SD24C_FPR_COMPLETE|%o|prime_H_classes=%o|prime_G_classes=%o|bound_num=%o|bound_den=%o|lt_half=%o\n",
       label, h_prime_classes, #intersections, Numerator(bound),
       Denominator(bound), bound lt 1/2;
quit;
