# Exact CTblLib necessary-order audit for sporadic socles in
# 10^7 < degree <= 10^8.  This is not a base-size computation.

degreeLower := 10^7;;
degreeUpper := 10^8;;
orderBound := degreeUpper * (degreeUpper - 1);;
SizeScreen([1000, 1000]);;

groupNames := [
    "M11", "M12", "M12.2", "M22", "M22.2", "M23", "M24",
    "J1", "J2", "J2.2", "J3", "J3.2", "J4",
    "HS", "HS.2", "McL", "McL.2", "Suz", "Suz.2", "He", "He.2",
    "Ru", "ON", "ON.2", "Co1", "Co2", "Co3", "Fi22", "Fi22.2",
    "Fi23", "Fi24'", "Fi24'.2", "HN", "HN.2", "Ly", "Th", "B", "M"
];;

# Multiplicity is intentional: CTblLib lists distinct maximal classes with
# the same name, order and index for Fi22 and Suz.
expected := [
    ["Co2", "3^1+4:2^1+4.s5", 45337600, 933120],
    ["Co3", "2^2.[2^7*3^2].S3", 17931375, 27648],
    ["Co3", "s3xpsl(2,8).3", 54648000, 9072],
    ["Fi22", "3^(1+6):2^(3+4):3^2:2", 12812800, 5038848],
    ["Fi22", "A10.2", 17791488, 3628800],
    ["Fi22", "A10.2", 17791488, 3628800],
    ["Fi22.2", "3^(1+6)_+:2^(3+4):(S3xS3)", 12812800, 10077696],
    ["Fi22.2", "G2(3).2", 15206400, 8491392],
    ["HN", "2^(1+8).(A5xA5).2", 74064375, 3686400],
    ["HN.2", "2^(1+8)_+.(A5xA5).2^2", 74064375, 7372800],
    ["ON", "3^4:2^(1+4)D10", 17778376, 25920],
    ["ON", "4^3.L3(2)", 42858585, 10752],
    ["ON", "L2(31)", 30968784, 14880],
    ["ON", "M11", 58183776, 7920],
    ["ON", "ONM11", 58183776, 7920],
    ["ON", "ONM5", 17778376, 25920],
    ["ON", "ONM8", 30968784, 14880],
    ["ON.2", "(3^2:4xA6).2^2", 17778376, 51840],
    ["ON.2", "3^4:2^(1+4).(5:4)", 17778376, 51840],
    ["ON.2", "4^3.(L3(2)x2)", 42858585, 21504],
    ["ON.2", "7^(1+2)_+:(3xD16)", 55978560, 16464],
    ["Ru", "3.A6.2^2", 33779200, 4320],
    ["Ru", "5^1+2:(2^5)", 36481536, 4000],
    ["Ru", "5^2:4s5", 12160512, 12000],
    ["Ru", "L2(13).2", 66816000, 2184],
    ["Ru", "L2(29)", 11980800, 12180],
    ["Suz", "(3^2:4xa6).2", 17297280, 25920],
    ["Suz", "(a6xa5).2", 10378368, 43200],
    ["Suz", "L2(25)", 57480192, 7800],
    ["Suz", "L3(3).2", 39916800, 11232],
    ["Suz", "L3(3).2", 39916800, 11232],
    ["Suz.2", "(3^2:8xA6).2", 17297280, 51840],
    ["Suz.2", "(A6:2_2xA5).2", 10378368, 86400],
    ["Suz.2", "L2(25).2_2", 57480192, 15600]
];;

observed := [];;
for groupName in groupNames do
    table := CharacterTable(groupName);
    if table = fail then Error("missing character table for ", groupName); fi;
    if Size(table) > orderBound then continue; fi;
    maximalNames := Maxes(table);
    if maximalNames = fail then Error("missing maximal list for ", groupName); fi;
    for maximalName in maximalNames do
        maximalTable := CharacterTable(maximalName);
        if maximalTable = fail then Error("missing maximal table ", maximalName); fi;
        subgroupOrder := Size(maximalTable);
        index := Size(table) / subgroupOrder;
        if index > degreeLower and index <= degreeUpper and
           subgroupOrder <= index - 1 then
            Add(observed, [groupName, maximalName, index, subgroupOrder]);
        fi;
    od;
od;

Sort(observed);;
Sort(expected);;
if observed <> expected then
    Error("sporadic 10^8 window inventory differs from pinned expectation");
fi;
for row in observed do
    Print("SPORADIC8_WINDOW|G=", row[1], "|H=", row[2],
          "|degree=", row[3], "|H_order=", row[4], "\n");
od;
Print("SPORADIC8_WINDOW_AUDIT_PASS|rows=", Length(observed),
      "|necessary_condition=H_order_le_degree_minus_1\n");
QUIT;
