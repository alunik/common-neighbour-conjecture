# Exact check of the O'Nan order obstruction used in
# generate_almost_simple_candidates.m.
#
# GAP's Character Table Library stores the complete maximal-subgroup table
# names for ON and ON.2.  A primitive almost-simple coset action has a
# core-free maximal point stabilizer.  The assertions below prove that every
# such action of degree at most 10^6 is excluded by |H| > degree - 1, the
# necessary order condition for base size two.

degreeCap := 10^6;;

CheckAlmostSimpleGroup := function(name, socleOrder)
    local table, tableOrder, maximalNames, corefree, maximalName,
          maximalTable, subgroupOrder, index, belowCap;
    table := CharacterTable(name);
    if table = fail then
        Error("missing character table for ", name);
    fi;
    tableOrder := Size(table);
    maximalNames := Maxes(table);
    if maximalNames = fail then
        Error("missing maximal-subgroup list for ", name);
    fi;

    corefree := [];
    for maximalName in maximalNames do
        maximalTable := CharacterTable(maximalName);
        if maximalTable = fail then
            Error("missing maximal-subgroup character table ", maximalName);
        fi;
        subgroupOrder := Size(maximalTable);
        index := tableOrder / subgroupOrder;
        Print("ON_MAXIMAL group=", name,
              " subgroup=", maximalName,
              " order=", subgroupOrder,
              " index=", index, "\n");
        # ON is the unique non-core-free maximal subgroup of ON.2.
        if subgroupOrder <> socleOrder then
            Add(corefree, [subgroupOrder, index]);
        fi;
    od;

    belowCap := Filtered(corefree, pair -> pair[2] <= degreeCap);
    if not ForAll(belowCap, pair -> pair[1] > pair[2] - 1) then
        Error("unobstructed O'Nan action found for ", name);
    fi;
    Print("ON_GROUP_CHECKED group=", name,
          " corefree_below_cap=", Length(belowCap), "\n");
    return belowCap;
end;

onOrder := 460815505920;;
onActions := CheckAlmostSimpleGroup("ON", onOrder);;
on2Actions := CheckAlmostSimpleGroup("ON.2", onOrder);;

if Set(onActions) <> [[3753792, 122760]] then
    Error("unexpected ON actions below the cap");
fi;
if on2Actions <> [] then
    Error("unexpected ON.2 actions below the cap");
fi;
Print("ON_ORDER_OBSTRUCTION_CHECK_PASS degree_cap=", degreeCap, "\n");
QUIT;
