///////////////////////////////////////////////////////////////////////////
// Exact recursive frontier for full-component A5,k=3 compound quotients.
//
// Starting at D=(C2 x S3) wr S_m, descend only through maximal subgroups
// whose compound top is transitive and whose induced component is still the
// full C2 x S3.  A branch stops at its first base-two node.  Branches with a
// proper component are routed to the separately certified component-wreath
// cases.  D-conjugate nodes are deduplicated exactly.
///////////////////////////////////////////////////////////////////////////

require assigned CD18M: "set CD18M";
if Type(CD18M) eq MonStgElt then CD18M := StringToInteger(CD18M); end if;
if not assigned CD18RMode then CD18RMode := "full"; end if;
if not assigned CD18ExportReps then CD18ExportReps := false; end if;
if Type(CD18ExportReps) eq MonStgElt then
    CD18ExportReps := CD18ExportReps eq "true";
end if;
require CD18M in {3,4,5}: "m must be 3, 4 or 5";
require CD18RMode in {"full", "s3"}: "component mode must be full or s3";
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
R := CD18RMode eq "full" select Rfull else
     sub<Rfull | [g @ rins[2] : g in Generators(s3)]>;
RImage := function(g)
    a := IsIdentity(g @ rprojs[1]) select X!1 else outer_action;
    return a * ((g @ rprojs[2]) @ s3_action);
end function;
RImageY := function(g)
    a := IsIdentity(g @ rprojs[1]) select Y!1 else y_outer;
    return a * ((g @ rprojs[2]) @ y_s3_action);
end function;
compound := Sym(CD18M);
D, base_inclusions, top_inclusion, top_projection := WreathProduct(R, compound);
width := Degree(R);
first_block := {1 .. width};

rx := sub<X | [RImage(g) : g in Generators(R)]>;
r_orbits := Sort([Setseq(o) : o in Orbits(rx)], func<a,b | a[1]-b[1]>);
representatives := [o[1] : o in r_orbits];
stabilisers := [sub<R | [g : g in R | representatives[i]^RImage(g)
                                  eq representatives[i]]>
               : i in [1 .. #representatives]];
profiles := [PowerSequence(Integers()) | ];
BuildProfiles := procedure(~output, partial, minimum)
    if #partial eq CD18M then Append(~output, partial); return; end if;
    for value in [minimum .. #representatives] do
        $$(~output, partial cat [value], value);
    end for;
end procedure;
BuildProfiles(~profiles, [], 1);

ProfileStabiliserWith := function(profile, alphabet_stabilisers)
    generators := [D | ];
    for coordinate in [1 .. CD18M] do
        for g in Generators(alphabet_stabilisers[profile[coordinate]]) do
            Append(~generators, g @ base_inclusions[coordinate]);
        end for;
    end for;
    for coordinate in [1 .. CD18M-1] do
        if profile[coordinate] eq profile[coordinate+1] then
            Append(~generators, (compound!(coordinate,coordinate+1)) @ top_inclusion);
        end if;
    end for;
    return sub<D | generators>;
end function;
profile_stabilisers := [ProfileStabiliserWith(p, stabilisers) : p in profiles];
y_orbits := [];
y_representatives := [];
y_stabilisers := [];
y_profiles := [];
y_profile_stabilisers := [];
if CD18ExportReps then
    ry := sub<Y | [RImageY(g) : g in Generators(R)]>;
    y_orbits := Sort([Setseq(o) : o in Orbits(ry)], func<a,b | a[1]-b[1]>);
    y_representatives := [o[1] : o in y_orbits];
    y_stabilisers := [sub<R | [g : g in R |
        y_representatives[i]^RImageY(g) eq y_representatives[i]]>
        : i in [1 .. #y_representatives]];
    y_profiles := [PowerSequence(Integers()) | ];
    BuildYProfiles := procedure(~output, partial, minimum)
        if #partial eq CD18M then Append(~output, partial); return; end if;
        for value in [minimum .. #y_representatives] do
            $$(~output, partial cat [value], value);
        end for;
    end procedure;
    BuildYProfiles(~y_profiles, [], 1);
    y_profile_stabilisers := [ProfileStabiliserWith(p, y_stabilisers)
                              : p in y_profiles];
end if;

Component := function(B)
    block_stabilizer := Stabilizer(B, first_block);
    return sub<R | [R![j^g : j in [1..width]]
                    : g in Generators(block_stabilizer)]>;
end function;

HasRegularPoint := function(B)
    if Order(B) gt 55^CD18M then return false, 0, D!1; end if;
    for i in [1 .. #profile_stabilisers] do
        K := profile_stabilisers[i];
        for d in DoubleCosetRepresentatives(D,B,K) do
            if Order(B meet K^(d^-1)) eq 1 then
                return true, i, d^-1;
            end if;
        end for;
    end for;
    return false, 0, D!1;
end function;

RegularPointCount := function(B)
    regular := 0;
    for K in profile_stabilisers do
        regular +:= #[d : d in DoubleCosetRepresentatives(D,B,K) |
                      Order(B meet K^(d^-1)) eq 1] * Order(B);
    end for;
    return regular;
end function;

Decompose := function(element)
    sigma := element @ top_projection;
    base := element * (sigma @ top_inclusion)^-1;
    components := [R | ];
    for coordinate in [1 .. CD18M] do
        offset := (coordinate-1)*width;
        Append(~components, R![(offset+j)^base-offset : j in [1..width]]);
    end for;
    return components, sigma;
end function;

ActProfile := function(profile, alphabet_representatives, image, element)
    components, sigma := Decompose(element);
    input := [alphabet_representatives[i] : i in profile];
    output := [Integers() | 0 : i in [1 .. CD18M]];
    for i in [1 .. CD18M] do
        output[i^sigma] := input[i]^image(components[i]);
    end for;
    return output;
end function;

PrintTuple := procedure(kind, node, orbit_number, tuple)
    printf "%o|case=%o|orbit=%o|tuple=[", kind, node, orbit_number;
    for i in [1 .. #tuple] do
        if i gt 1 then printf ","; end if;
        printf "%o", tuple[i];
    end for;
    print "]";
end procedure;

PrintGenerators := procedure(node, B)
    generators := [g : g in Generators(B)];
    for generator_number in [1 .. #generators] do
        components, sigma := Decompose(generators[generator_number]);
        printf "FRONTIER_GENERATOR|node=%o|generator=%o|top=%o|components=[",
               node, generator_number, Eltseq(sigma);
        for coordinate in [1 .. CD18M] do
            if coordinate gt 1 then printf ";"; end if;
            outer_bit := IsIdentity(components[coordinate] @ rprojs[1])
                select 0 else 1;
            printf "%o:%o", outer_bit, Eltseq(components[coordinate] @ rprojs[2]);
        end for;
        print "]";
        if CD18ExportReps then
            printf "ORBIT_REP_GENERATOR|case=%o|generator=%o|top=%o|components=[",
                   node, generator_number, Eltseq(sigma);
            for coordinate in [1 .. CD18M] do
                if coordinate gt 1 then printf ";"; end if;
                outer_bit := IsIdentity(components[coordinate] @ rprojs[1])
                    select 0 else 1;
                printf "%o:%o", outer_bit,
                       Eltseq(components[coordinate] @ rprojs[2]);
            end for;
            print "]";
        end if;
    end for;
end procedure;

seen := [D];
queue := [D];
depths := [0];
head := 1;
frontier := 0;
proper_routes := 0;
inherited_routes := 0;
zero_nodes := 1;
closed_conjugates := [];
printf "FRONTIER_SHAPE|component=%o|m=%o|R=%o|D=%o|profiles=%o\n",
       CD18RMode, CD18M, Order(R), Order(D), #profiles;
if CD18ExportReps then
    printf "ORBIT_REP_SHAPE|component=%o|m=%o|R=%o|D=%o|x_orbits=%o|x_profiles=%o|y_orbits=%o|y_profiles=%o\n",
           CD18RMode, CD18M, Order(R), Order(D), #r_orbits, #profiles,
           #y_orbits, #y_profiles;
end if;
while head le #queue do
    B := queue[head];
    depth := depths[head];
    head +:= 1;
    for record in MaximalSubgroups(B) do
        C := record`subgroup;
        top := C @ top_projection;
        if not IsTransitive(top) then continue; end if;
        component := Component(C);
        if Order(component) ne Order(R) then
            proper_routes +:= 1;
            continue;
        end if;
        duplicate := false;
        for prior in seen do
            if Order(prior) eq Order(C) and IsConjugate(D, prior, C) then
                duplicate := true;
                break;
            end if;
        end for;
        if duplicate then continue; end if;
        Append(~seen, C);
        has_regular, profile_number, transporter := HasRegularPoint(C);
        if has_regular then
            inherited := false;
            for overgroup in closed_conjugates do
                if C subset overgroup then inherited := true; break; end if;
            end for;
            if inherited then
                inherited_routes +:= 1;
                printf "FRONTIER_INHERITED|serial=%o|depth=%o|order=%o|index=%o|top=%o\n",
                       inherited_routes, depth+1, Order(C), Index(D,C), Order(top);
                continue;
            end if;
            frontier +:= 1;
            regular := RegularPointCount(C);
            printf "FRONTIER_BASE2|node=%o|depth=%o|order=%o|index=%o|top=%o|profile=%o|regular=%o|gt_half=%o|generators=%o\n",
                   frontier, depth+1, Order(C), Index(D,C), Order(top),
                   profile_number, regular, 2*regular gt 55^CD18M,
                   #Generators(C);
            if not CD18ExportReps then
                PrintGenerators(frontier, C);
            elif 2*regular le 55^CD18M then
                printf "ORBIT_REP_CASE|case=%o|M=%o|regular=%o|regular_orbits=%o\n",
                       frontier, Order(C), regular, regular div Order(C);
                PrintGenerators(frontier, C);
                x_count := 0;
                for i in [1 .. #profiles] do
                    K := profile_stabilisers[i];
                    for d in DoubleCosetRepresentatives(D,C,K) do
                        if Order(C meet K^(d^-1)) eq 1 then
                            x_count +:= 1;
                            PrintTuple("X_REGULAR_REP", frontier, x_count,
                                ActProfile(profiles[i], representatives,
                                           RImage, d^-1));
                        end if;
                    end for;
                end for;
                y_count := 0;
                for i in [1 .. #y_profiles] do
                    K := y_profile_stabilisers[i];
                    for d in DoubleCosetRepresentatives(D,C,K) do
                        y_count +:= 1;
                        PrintTuple("Y_TARGET_REP", frontier, y_count,
                            ActProfile(y_profiles[i], y_representatives,
                                       RImageY, d^-1));
                    end for;
                end for;
                printf "ORBIT_REP_CASE_COMPLETE|case=%o|x=%o|y=%o\n",
                       frontier, x_count, y_count;
            end if;
            normalizer := Normalizer(D,C);
            for d in Transversal(D,normalizer) do
                Append(~closed_conjugates, C^d);
            end for;
        else
            zero_nodes +:= 1;
            Append(~queue, C);
            Append(~depths, depth+1);
            printf "FRONTIER_ZERO|serial=%o|depth=%o|order=%o|index=%o|top=%o\n",
                   zero_nodes, depth+1, Order(C), Index(D,C), Order(top);
        end if;
    end for;
end while;
printf "FRONTIER_COMPLETE|component=%o|m=%o|seen=%o|zero_nodes=%o|base2_frontier=%o|proper_routes=%o|inherited_routes=%o|max_depth=%o\n",
       CD18RMode, CD18M, #seen, zero_nodes, frontier, proper_routes, inherited_routes,
       Maximum(depths);
quit;
