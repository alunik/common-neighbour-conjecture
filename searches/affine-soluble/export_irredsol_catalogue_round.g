#############################################################################
## Export one whole maximal-uncovered layer of each requested IRREDSOL
## guardian.  This is the scalable version of catalogue-frontier descent.
##
## ASX_TASKS records have fields n, p, d, guardian, tested, safe.  Both lists
## contain IRREDSOL catalogue indices for this guardian.  Safe entries cover
## all their catalogue subgroups; tested-but-unsafe entries are removed so
## the next maximal layer becomes visible.
#############################################################################

if LoadPackage("irredsol") = fail then Error("IRREDSOL is required"); fi;
if not IsBound(ASX_TASKS) then Error("ASX_TASKS must be assigned"); fi;
SizeScreen([1000000, 1000000]);

ASX_Context := function(task)
    local m, q, guardianData, P, hom, c, matrix, matrixGroup;
    m := task.n / task.d;
    q := task.p^task.d;
    LoadAbsolutelyIrreducibleSolubleGroupData(m, q);
    if m = 1 then
        P := IRREDSOL_DATA.GUARDIANS[1][q][1];
        c := MinimalGeneratingSet(P)[1];
        matrix := [[Z(q)]];
        matrixGroup := GroupWithGenerators([matrix], IdentityMat(1, GF(q)));
        hom := GroupHomomorphismByImagesNC(P, matrixGroup, [c], [matrix]);
        SetIsBijective(hom, true);
        return rec(m := m, q := q, P := P, hom := hom,
                   fpcgs := FamilyPcgs(P), groups := []);
    fi;
    guardianData := IRREDSOL_DATA.GUARDIANS[m][q][task.guardian];
    P := Source(guardianData[3]);
    return rec(m := m, q := q, P := P, hom := guardianData[3],
               fpcgs := FamilyPcgs(P), groups := []);
end;

ASX_Order := function(context, k)
    if context.m = 1 then
        return IRREDSOL_DATA.GROUPS_DIM1[context.q][k][1];
    fi;
    return IRREDSOL_DATA.GROUPS[context.m][context.q][k][4];
end;

ASX_GuardianIndex := function(context, k)
    if context.m = 1 then return 1; fi;
    return IRREDSOL_DATA.GROUPS[context.m][context.q][k][1];
end;

ASX_Group := function(context, k)
    local order, generator, code;
    if not IsBound(context.groups[k]) then
        if context.m = 1 then
            order := ASX_Order(context, k);
            generator := context.fpcgs[1]^((context.q - 1) / order);
            context.groups[k] := SubgroupNC(context.P, [generator]);
        else
            code := IRREDSOL_DATA.GROUPS[context.m][context.q][k][2];
            context.groups[k] := GroupOfPcgs(
                CanonicalPcgsByNumber(context.fpcgs, code));
        fi;
    fi;
    return context.groups[k];
end;

ASX_PrimeFieldMatrices := function(task, context, K)
    local matrices, basis;
    matrices := List(SmallGeneratingSet(K), g -> ImageElm(context.hom, g));
    if task.d > 1 then
        basis := CanonicalBasis(AsVectorSpace(GF(task.p), GF(context.q)));
        matrices := List(matrices, matrix -> BlownUpMat(basis, matrix));
    fi;
    return matrices;
end;

ASX_Emit := function(task, context, k)
    local K, matrices, matrix, row, label;
    K := ASX_Group(context, k);
    matrices := ASX_PrimeFieldMatrices(task, context, K);
    label := Concatenation("ICR_", String(task.n), "_", String(task.p), "_",
                           String(task.d), "_", String(task.guardian), "_",
                           String(k));
    Print("action\n");
    Print("label ", label, "\n");
    Print("p ", task.p, "\n");
    Print("n ", task.n, "\n");
    Print("order ", Size(K), "\n");
    Print("orientation row\n");
    Print("gens ", Length(matrices), "\n");
    for matrix in matrices do
        for row in matrix do
            Print(JoinStringsWithSeparator(
                List(row, entry -> String(IntFFE(entry))), " "), "\n");
        od;
    od;
    Print("end\n");
end;

ASX_Round := function(task, context)
    local threshold, indices, safeGroups, candidates, maximalIndices,
          maximalGroups, i, j, H, covered, position;
    threshold := task.p^task.n - 1;
    indices := IndicesIrreducibleSolubleMatrixGroups(task.n, task.p, task.d);
    indices := Filtered(indices, i ->
        ASX_GuardianIndex(context, i) = task.guardian and
        ASX_Order(context, i) <= threshold and not i in task.tested);

    safeGroups := List(task.safe, i -> ASX_Group(context, i));
    candidates := [];
    for i in indices do
        H := ASX_Group(context, i);
        covered := false;
        for position in [1 .. Length(task.safe)] do
            if ASX_Order(context, task.safe[position]) mod ASX_Order(context, i) = 0
                    and IsSubgroup(safeGroups[position], H) then
                covered := true;
                break;
            fi;
        od;
        if not covered then Add(candidates, i); fi;
    od;
    Sort(candidates, function(i, j)
        if ASX_Order(context, i) <> ASX_Order(context, j) then
            return ASX_Order(context, i) > ASX_Order(context, j);
        fi;
        return i < j;
    end);

    maximalIndices := [];
    maximalGroups := [];
    for i in candidates do
        H := ASX_Group(context, i);
        covered := false;
        for position in [1 .. Length(maximalIndices)] do
            if ASX_Order(context, maximalIndices[position]) mod
                    ASX_Order(context, i) = 0 and
                    IsSubgroup(maximalGroups[position], H) then
                covered := true;
                break;
            fi;
        od;
        if not covered then
            Add(maximalIndices, i);
            Add(maximalGroups, H);
        fi;
    od;
    for i in maximalIndices do ASX_Emit(task, context, i); od;
end;

Print("AFFINE_SAXL_V1\n");
for ASX_Task in ASX_TASKS do
    ASX_ContextRecord := ASX_Context(ASX_Task);
    ASX_Round(ASX_Task, ASX_ContextRecord);
od;
quit;
