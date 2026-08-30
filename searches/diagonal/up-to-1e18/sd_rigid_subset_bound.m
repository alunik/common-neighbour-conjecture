///////////////////////////////////////////////////////////////////////////
// Exact prime-order union bound for rigid k-subsets in Hol(T).
//
// For the full SD normalizer with S_k top, a regular point corresponds to
// an ordered k-subset whose underlying set has trivial stabilizer in Hol(T).
// Every nontrivial set stabilizer contains a prime-order element.  Summing
// fixed k-subsets over prime-order Hol(T)-classes therefore gives a rigorous
// lower bound on the regular-neighbour density.
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

tuple, t_order := SimpleGroupId(SD18SimpleId);
require t_order le 1000000: "rigid-subset engine restricted to small residual factors";
require 3 le SD18K and SD18K lt t_order: "bad subset size";
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
        for degree in [0 .. k] do
            if coefficients[degree + 1] eq 0 then continue; end if;
            for chosen in [0 .. multiplicity] do
                target := degree + chosen * length;
                if target gt k then break; end if;
                updated[target + 1] +:= coefficients[degree + 1] *
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
print "SD18_RIGID_SUBSET_BOUND_V1";
label := Sprintf("SD18R_sid%o_k%o", SD18SimpleId, SD18K);
printf "SD18_RIGID_ACTION|%o|sid=%o|T=%o|T_order=%o|k=%o|A=%o|Hol=%o\n",
       label, SD18SimpleId, SimpleGroupName(SD18SimpleId), t_order, SD18K,
       Order(A), Order(Hol);
for sequence in [1 .. #prime_classes] do
    c := prime_classes[sequence];
    fixed := FixedKSubsets(c[3], SD18K);
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
total_subsets := Binomial(t_order, SD18K);
rigid_lower := Maximum(0, total_subsets - bad_upper);
density_lower := Rationals()!(rigid_lower * Factorial(SD18K)) /
                 t_order^SD18K;
matching_lower := &*[Integers() | t_order - 2*j : j in [0 .. SD18K - 1]];
bad_ordered_upper := bad_upper * Factorial(SD18K);
rainbow_common := matching_lower gt 2 * bad_ordered_upper;
printf "SD18_RIGID_SUBSET_COMPLETE|%o|sid=%o|k=%o|prime_classes=%o|total_subsets=%o|bad_upper=%o|rigid_lower=%o|density_num=%o|density_den=%o|gt_half=%o|matching_lower=%o|bad_ordered_upper=%o|rainbow_common=%o\n",
       label, SD18SimpleId, SD18K, #prime_classes, total_subsets, bad_upper,
       rigid_lower, Numerator(density_lower), Denominator(density_lower),
       density_lower gt 1/2, matching_lower, bad_ordered_upper,
       rainbow_common;
quit;
