# Exact CTblLib metadata audit for almost-simple sporadic groups in the new
# degree window 10^6 < d <= 10^7.
#
# This is not a base-size computation.  It lists every maximal coset action
# satisfying the necessary condition |H| <= d-1.  Burness--Giudici,
# Theorem 6.1, is then applied conditionally to the base-two rows.

degreeLower := 10^6;;
degreeUpper := 10^7;;
orderBound := degreeUpper * (degreeUpper - 1);;
SizeScreen([1000, 1000]);;
Print("GAP_VERSION|", GAPInfo.Version,
      "|CTBLLIB_VERSION|", InstalledPackageVersion("ctbllib"), "\n");

# The 26 sporadic simple groups and every nontrivial outer extension.
groupNames := [
    "M11", "M12", "M12.2", "M22", "M22.2", "M23", "M24",
    "J1", "J2", "J2.2", "J3", "J3.2", "J4",
    "HS", "HS.2", "McL", "McL.2", "Suz", "Suz.2", "He", "He.2",
    "Ru", "ON", "ON.2", "Co1", "Co2", "Co3", "Fi22", "Fi22.2",
    "Fi23", "Fi24'", "Fi24'.2", "HN", "HN.2", "Ly", "Th", "B", "M"
];;

expected := [
    ["M24", "L3(2)", 1457280, 168],
    ["Suz", "2^2+8(a5xs3)", 1216215, 368640],
    ["Suz", "M12.2", 2358720, 190080],
    ["Suz", "3^2+4:2(2^2xa4)2", 3203200, 139968],
    ["Suz.2", "2^(2+8):(S5xS3)", 1216215, 737280],
    ["Suz.2", "M12.2x2", 2358720, 380160],
    ["Suz.2", "3^(2+4):2(S4xD8)", 3203200, 279936],
    ["He", "7:3xpsl(3,2)", 1142400, 3528],
    ["He", "5^2:4A4", 3358656, 1200],
    ["He.2", "7:6xL3(2)", 1142400, 7056],
    ["He.2", "Fi22N5", 3358656, 2400],
    ["Ru", "L2(25).2^2", 4677120, 31200],
    ["Ru", "A8", 7238400, 20160],
    ["ON", "J1", 2624832, 175560],
    ["ON", "4_2.L3(4).2_1", 2857239, 161280],
    ["ON.2", "J1x2", 2624832, 351120],
    ["ON.2", "4_2.L3(4).(2^2)_{12*3}", 2857239, 322560],
    ["Co3", "2^4.a8", 1536975, 322560],
    ["Co3", "L3(4).D12", 2049300, 241920],
    ["Co3", "2xm12", 2608200, 190080]
];;

observed := [];;
for groupName in groupNames do
    table := CharacterTable(groupName);
    if table = fail then
        Error("missing character table for ", groupName);
    fi;
    if Size(table) > orderBound then
        continue;
    fi;
    maximalNames := Maxes(table);
    if maximalNames = fail then
        Error("missing CTblLib maximal-subgroup list for ", groupName);
    fi;
    for maximalName in maximalNames do
        maximalTable := CharacterTable(maximalName);
        if maximalTable = fail then
            Error("missing maximal-subgroup table ", maximalName);
        fi;
        subgroupOrder := Size(maximalTable);
        index := Size(table) / subgroupOrder;
        if index > degreeLower and index <= degreeUpper and
           subgroupOrder <= index - 1 then
            Add(observed, [groupName, maximalName, index, subgroupOrder]);
        fi;
    od;
od;

Sort(observed);
Sort(expected);
if observed <> expected then
    Error("sporadic degree-window inventory differs from pinned expectation");
fi;

for row in observed do
    Print("SPORADIC_WINDOW|G=", row[1], "|H=", row[2],
          "|degree=", row[3], "|H_order=", row[4], "\n");
od;
Print("SPORADIC_WINDOW_AUDIT_PASS|rows=", Length(observed),
      "|necessary_condition=H_order_le_degree_minus_1\n");
QUIT;
