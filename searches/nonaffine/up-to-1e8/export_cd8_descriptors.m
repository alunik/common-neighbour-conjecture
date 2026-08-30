///////////////////////////////////////////////////////////////////////////
// Exact quotient/action descriptors for the 39 CD8 classes.
//
// Each generator is encoded as a signed permutation of the six A5 socle
// coordinates.  The bit records the standard outer involution of A5.
///////////////////////////////////////////////////////////////////////////

SetColumns(0);
tuple, t_order := SimpleGroupId(1);
require t_order eq 60: "expected A5";
A := AutomorphismGroupSimpleGroup(tuple);
T := Socle(A);
Q, qmap := quo<A | T>;
require Order(Q) eq 2: "expected Out(A5)=C2";
S3 := Sym(3);
R, rins, rprojs := DirectProduct(Q, S3);
S2 := Sym(2);
D, dbins, dtopin, dtopproj := WreathProduct(R, S2);

candidates := [];
first_block := {1 .. Degree(R)};
for record in Subgroups(D) do
    B := record`subgroup;
    if not IsTransitive(B @ dtopproj) then continue; end if;
    block_stabilizer := Stabilizer(B, first_block);
    coordinate_generators := [R |
        R![j^g : j in [1 .. Degree(R)]]
        : g in Generators(block_stabilizer)];
    coordinate := sub<R | coordinate_generators>;
    coordinate_top := sub<S3 |
        [g @ rprojs[2] : g in Generators(coordinate)]>;
    if not IsTransitive(coordinate_top) then continue; end if;
    if 3600 * Order(B) gt 12960000 - 1 then continue; end if;
    Append(~candidates, B);
end for;
require #candidates eq 39: "bad CD8 candidate count";

print "CD8_DESCRIPTORS_V1";
histogram := AssociativeArray(Integers());
for case_id in [1 .. #candidates] do
    B := candidates[case_id];
    b_order := Order(B);
    if not IsDefined(histogram, b_order) then histogram[b_order] := 0; end if;
    histogram[b_order] +:= 1;
    signed_generators := [Sym(12) | ];
    generators := Setseq(Generators(B));
    printf "CD_CASE|case=%o|B=%o|H=%o|degree=12960000|generators=%o\n",
           case_id, b_order, 3600 * b_order, #generators;
    for sequence -> b in generators do
        top := b @ dtopproj;
        base := b * (top @ dtopin)^-1;
        targets := [];
        bits := [];
        width := Degree(R);
        for component in [1 .. 2] do
            offset := (component - 1) * width;
            r := R![(offset + j)^base - offset : j in [1 .. width]];
            q := r @ rprojs[1];
            internal := r @ rprojs[2];
            bit := q ne Q!1 select 1 else 0;
            for coordinate in [1 .. 3] do
                target_component := component ^ top;
                target_coordinate := coordinate ^ internal;
                Append(~targets, (target_component - 1) * 3 + target_coordinate);
                Append(~bits, bit);
            end for;
        end for;
        images := [];
        for coordinate in [1 .. 6] do
            for sign in [0, 1] do
                Append(~images, 2 * (targets[coordinate] - 1) +
                                (sign + bits[coordinate]) mod 2 + 1);
            end for;
        end for;
        Append(~signed_generators, Sym(12)!images);
        printf "CD_GEN|case=%o|sequence=%o|perm=", case_id, sequence;
        for i in [1 .. 6] do
            if i gt 1 then printf ","; end if;
            printf "%o", targets[i];
        end for;
        printf "|outer=";
        for i in [1 .. 6] do
            if i gt 1 then printf ","; end if;
            printf "%o", bits[i];
        end for;
        print "";
    end for;
    signed_group := sub<Sym(12) | signed_generators>;
    require Order(signed_group) eq b_order: "signed descriptor loses quotient order";
    printf "CD_CASE_COMPLETE|case=%o|B=%o|signed_closure=%o\n",
           case_id, b_order, Order(signed_group);
end for;

expected := AssociativeArray(Integers());
for pair in [<6,2>, <12,6>, <18,1>, <24,6>, <36,6>, <48,1>,
             <72,10>, <144,6>, <288,1>] do
    expected[pair[1]] := pair[2];
end for;
require Sort(Setseq(Keys(histogram))) eq Sort(Setseq(Keys(expected))):
    "CD8 order histogram support mismatch";
require &and[histogram[key] eq expected[key] : key in Keys(expected)]:
    "CD8 order histogram multiplicity mismatch";
print "CD8_DESCRIPTORS_COMPLETE|cases=39|degree=12960000";
quit;
