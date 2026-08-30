///////////////////////////////////////////////////////////////////////////
// Export the exact regular-source and target-orbit representatives for one
// already-certified A5,k=3,m=6 recursive-frontier group.
//
// A generated case descriptor is evaluated immediately before this file.
///////////////////////////////////////////////////////////////////////////

require assigned CaseMode: "missing CaseMode";
require assigned CaseNumber: "missing CaseNumber";
require assigned CaseExpectedOrder: "missing CaseExpectedOrder";
require assigned CaseExpectedRegular: "missing CaseExpectedRegular";
require assigned CaseTops: "missing CaseTops";
require assigned CaseOuters: "missing CaseOuters";
require assigned CasePerms: "missing CasePerms";
assert (CaseMode eq "full") or (CaseMode eq "s3");
assert #CaseTops gt 0 and #CaseOuters eq #CaseTops and
       #CasePerms eq #CaseTops;
SetColumns(0);

X := Sym(55);
top3 := X![20,2,11,22,9,34,52,6,46,41,21,12,1,27,47,54,26,42,36,13,3,25,28,31,4,37,48,40,19,53,50,15,43,8,7,29,17,30,18,23,45,39,51,44,10,5,32,14,16,24,33,35,38,49,55];
top2 := X![1,12,21,22,46,6,41,34,9,52,11,2,20,47,27,54,36,42,26,13,3,4,40,50,25,19,15,28,37,43,31,48,53,8,45,17,29,51,39,23,7,18,30,44,35,5,14,32,49,24,38,10,33,16,55];
outer_action := X![11,12,13,4,8,9,10,5,6,7,1,2,3,17,19,18,14,16,15,21,20,22,24,23,25,27,26,31,32,33,28,29,30,46,45,47,48,51,49,50,52,54,53,55,35,34,36,37,39,40,38,41,43,42,44];

Y := Sym(77);
y_top3 := Y![1,7,31,54,71,2,6,8,29,10,19,34,17,47,70,14,64,55,30,20,9,39,65,73,38,56,49,3,21,11,28,32,33,37,40,43,12,50,66,53,27,72,68,23,57,4,16,15,41,25,42,26,35,46,63,52,69,58,77,76,59,5,18,13,44,22,24,36,45,48,62,51,67,74,60,75,61];
y_top2 := Y![1,7,31,54,71,6,2,8,9,20,30,34,64,14,55,47,17,70,19,10,29,65,39,73,49,56,38,28,21,11,3,33,32,12,53,68,37,27,23,40,50,57,43,66,72,46,16,63,25,41,69,52,35,4,15,26,42,58,59,76,77,62,48,13,22,44,67,36,51,18,5,45,24,74,75,60,61];
y_outer := Y![1,2,3,5,4,6,7,8,19,20,21,12,16,17,18,13,14,15,9,10,11,25,27,26,22,24,23,28,30,29,31,33,32,34,36,35,37,39,38,43,44,45,40,41,42,62,64,63,65,66,69,67,68,71,70,73,72,74,75,77,76,46,48,47,49,50,52,53,51,55,54,57,56,58,59,61,60];

outer := CyclicGroup(2);
s3 := Sym(3);
s3_action := hom<s3 -> X | [top3, top2]>;
y_s3_action := hom<s3 -> Y | [y_top3, y_top2]>;
Rfull, rins, rprojs := DirectProduct(outer, s3);
R := CaseMode eq "full" select Rfull else
     sub<Rfull | [g @ rins[2] : g in Generators(s3)]>;
RImage := function(g)
    a := IsIdentity(g @ rprojs[1]) select X!1 else outer_action;
    return a * ((g @ rprojs[2]) @ s3_action);
end function;
RImageY := function(g)
    a := IsIdentity(g @ rprojs[1]) select Y!1 else y_outer;
    return a * ((g @ rprojs[2]) @ y_s3_action);
end function;

m := 6;
compound := Sym(m);
D, base_inclusions, top_inclusion, top_projection := WreathProduct(R, compound);
width := Degree(R);

ComponentElement := function(bit, permutation)
    assert bit in {0,1};
    a := bit eq 0 select outer!1 else outer.1;
    return (a @ rins[1]) * ((s3!permutation) @ rins[2]);
end function;

input_generators := [D | ];
for number in [1 .. #CaseTops] do
    assert #CaseTops[number] eq m and #CaseOuters[number] eq m and
           #CasePerms[number] eq m;
    base := D!1;
    for coordinate in [1 .. m] do
        component := ComponentElement(CaseOuters[number][coordinate],
                                      CasePerms[number][coordinate]);
        assert component in R;
        base *:= component @ base_inclusions[coordinate];
    end for;
    sigma := compound!CaseTops[number];
    Append(~input_generators, base * (sigma @ top_inclusion));
end for;
B := sub<D | input_generators>;
assert Order(B) eq CaseExpectedOrder;

rx := sub<X | [RImage(g) : g in Generators(R)]>;
r_orbits := Sort([Setseq(o) : o in Orbits(rx)], func<a,b | a[1]-b[1]>);
representatives := [o[1] : o in r_orbits];
stabilisers := [sub<R | [g : g in R | representatives[i]^RImage(g)
                                  eq representatives[i]]>
               : i in [1 .. #representatives]];
profiles := [PowerSequence(Integers()) | ];
BuildProfiles := procedure(~output, partial, minimum)
    if #partial eq m then Append(~output, partial); return; end if;
    for value in [minimum .. #representatives] do
        $$(~output, partial cat [value], value);
    end for;
end procedure;
BuildProfiles(~profiles, [], 1);

ProfileStabiliserWith := function(profile, alphabet_stabilisers)
    generators := [D | ];
    for coordinate in [1 .. m] do
        for g in Generators(alphabet_stabilisers[profile[coordinate]]) do
            Append(~generators, g @ base_inclusions[coordinate]);
        end for;
    end for;
    for coordinate in [1 .. m-1] do
        if profile[coordinate] eq profile[coordinate+1] then
            Append(~generators, (compound!(coordinate,coordinate+1)) @ top_inclusion);
        end if;
    end for;
    return sub<D | generators>;
end function;
profile_stabilisers := [ProfileStabiliserWith(p, stabilisers) : p in profiles];

ry := sub<Y | [RImageY(g) : g in Generators(R)]>;
y_orbits := Sort([Setseq(o) : o in Orbits(ry)], func<a,b | a[1]-b[1]>);
y_representatives := [o[1] : o in y_orbits];
y_stabilisers := [sub<R | [g : g in R |
    y_representatives[i]^RImageY(g) eq y_representatives[i]]>
    : i in [1 .. #y_representatives]];
y_profiles := [PowerSequence(Integers()) | ];
BuildYProfiles := procedure(~output, partial, minimum)
    if #partial eq m then Append(~output, partial); return; end if;
    for value in [minimum .. #y_representatives] do
        $$(~output, partial cat [value], value);
    end for;
end procedure;
BuildYProfiles(~y_profiles, [], 1);
y_profile_stabilisers := [ProfileStabiliserWith(p, y_stabilisers)
                          : p in y_profiles];

Decompose := function(element)
    sigma := element @ top_projection;
    base := element * (sigma @ top_inclusion)^-1;
    components := [R | ];
    for coordinate in [1 .. m] do
        offset := (coordinate-1)*width;
        Append(~components, R![(offset+j)^base-offset : j in [1..width]]);
    end for;
    return components, sigma;
end function;

ActProfile := function(profile, alphabet_representatives, image, element)
    components, sigma := Decompose(element);
    input := [alphabet_representatives[i] : i in profile];
    output := [Integers() | 0 : i in [1 .. m]];
    for i in [1 .. m] do
        output[i^sigma] := input[i]^image(components[i]);
    end for;
    return output;
end function;

PrintTuple := procedure(kind, orbit_number, tuple)
    printf "%o|case=%o|orbit=%o|tuple=[", kind, CaseNumber, orbit_number;
    for i in [1 .. #tuple] do
        if i gt 1 then printf ","; end if;
        printf "%o", tuple[i];
    end for;
    print "]";
end procedure;

printf "ORBIT_REP_SHAPE|component=%o|m=%o|R=%o|D=%o|x_orbits=%o|x_profiles=%o|y_orbits=%o|y_profiles=%o\n",
       CaseMode, m, Order(R), Order(D), #r_orbits, #profiles,
       #y_orbits, #y_profiles;
printf "ORBIT_REP_CASE|case=%o|M=%o|regular=%o|regular_orbits=%o\n",
       CaseNumber, Order(B), CaseExpectedRegular, CaseExpectedRegular div Order(B);
for number in [1 .. #input_generators] do
    components, sigma := Decompose(input_generators[number]);
    printf "ORBIT_REP_GENERATOR|case=%o|generator=%o|top=%o|components=[",
           CaseNumber, number, Eltseq(sigma);
    for coordinate in [1 .. m] do
        if coordinate gt 1 then printf ";"; end if;
        outer_bit := IsIdentity(components[coordinate] @ rprojs[1]) select 0 else 1;
        printf "%o:%o", outer_bit, Eltseq(components[coordinate] @ rprojs[2]);
    end for;
    print "]";
end for;

x_count := 0;
for i in [1 .. #profiles] do
    K := profile_stabilisers[i];
    for d in DoubleCosetRepresentatives(D,B,K) do
        if Order(B meet K^(d^-1)) eq 1 then
            x_count +:= 1;
            PrintTuple("X_REGULAR_REP", x_count,
                ActProfile(profiles[i], representatives, RImage, d^-1));
        end if;
    end for;
end for;
assert x_count * Order(B) eq CaseExpectedRegular;
y_count := 0;
for i in [1 .. #y_profiles] do
    K := y_profile_stabilisers[i];
    for d in DoubleCosetRepresentatives(D,B,K) do
        y_count +:= 1;
        PrintTuple("Y_TARGET_REP", y_count,
            ActProfile(y_profiles[i], y_representatives, RImageY, d^-1));
    end for;
end for;
printf "ORBIT_REP_CASE_COMPLETE|case=%o|x=%o|y=%o\n",
       CaseNumber, x_count, y_count;
quit;
