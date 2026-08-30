///////////////////////////////////////////////////////////////////////////
// Exact prime-order union bound for rigid k-subsets in Hol(T), for the
// incremental SD24 degree window.
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

tuple, t_order := SimpleGroupId(SD24SimpleId);
degree := t_order^(SD24K - 1);
require 1000000000000000000 lt degree and
        degree le 1000000000000000000000000:
        "shape outside exact SD24 window";
require t_order le 1000000:
        "rigid-subset engine restricted to small residual factors";
require 3 le SD24K and SD24K lt t_order: "bad subset size";
A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
require Order(T) eq t_order: "simple socle mismatch";
Hol := Holomorph(T);
require Degree(Hol) eq t_order: "bad holomorph degree";
require Order(Hol) eq t_order * Order(A): "bad holomorph order";

FixedKSubsets := function(g, k)
    coefficients := [Integers() | 1] cat [Integers() | 0 : i in [1 .. k]];
    for pair in CycleStructure(g) do
        length := pair[1];
        multiplicity := pair[2];
        updated := [Integers() | 0 : i in [0 .. k]];
        for subset_degree in [0 .. k] do
            if coefficients[subset_degree + 1] eq 0 then continue; end if;
            for chosen in [0 .. multiplicity] do
                target := subset_degree + chosen * length;
                if target gt k then break; end if;
                updated[target + 1] +:= coefficients[subset_degree + 1] *
                                        Binomial(multiplicity, chosen);
            end for;
        end for;
        coefficients := updated;
    end for;
    return coefficients[k + 1];
end function;

classes := Classes(Hol);
prime_classes := [c : c in classes | IsPrime(c[1])];
bad_upper := 0;
print "SD24_RIGID_SUBSET_BOUND_V1";
label := Sprintf("SD24R_sid%o_k%o", SD24SimpleId, SD24K);
printf "SD24_RIGID_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|degree=%o|A=%o|Hol=%o\n",
       label, SD24SimpleId, SimpleGroupName(SD24SimpleId), t_order, SD24K,
       degree, Order(A), Order(Hol);
for sequence in [1 .. #prime_classes] do
    c := prime_classes[sequence];
    fixed := FixedKSubsets(c[3], SD24K);
    contribution := c[2] * fixed;
    bad_upper +:= contribution;
    printf "SUBSET_CLASS|sequence=%o|prime=%o|class_size=%o|fixed_subsets=%o|contribution=%o|cycles=",
           sequence, c[1], c[2], fixed, contribution;
    cycle_structure := CycleStructure(c[3]);
    for i in [1 .. #cycle_structure] do
        if i gt 1 then printf ","; end if;
        printf "%o:%o", cycle_structure[i][1], cycle_structure[i][2];
    end for;
    print "";
end for;
total_subsets := Binomial(t_order, SD24K);
rigid_lower := Maximum(0, total_subsets - bad_upper);
density_lower := Rationals()!(rigid_lower * Factorial(SD24K)) /
                 t_order^SD24K;
matching_lower := &*[Integers() | t_order - 2*j : j in [0 .. SD24K - 1]];
bad_ordered_upper := bad_upper * Factorial(SD24K);
rainbow_common := matching_lower gt 2 * bad_ordered_upper;
printf "SD24_RIGID_SUBSET_COMPLETE|%o|sid=%o|k=%o|prime_classes=%o|total_subsets=%o|bad_upper=%o|rigid_lower=%o|density_num=%o|density_den=%o|gt_half=%o|matching_lower=%o|bad_ordered_upper=%o|rainbow_common=%o\n",
       label, SD24SimpleId, SD24K, #prime_classes, total_subsets, bad_upper,
       rigid_lower, Numerator(density_lower), Denominator(density_lower),
       density_lower gt 1/2, matching_lower, bad_ordered_upper,
       rainbow_common;
quit;
