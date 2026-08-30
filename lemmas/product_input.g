# The product-action input (Lemma 4.1 of the paper).
#
# K = PSL(2,11) acts on the 55 cosets of a dihedral subgroup of order 12.
# The script checks that the two regular suborbits give self-paired
# orbitals, exhibits a pair (alpha, beta) and a regular orbital O such that
# no point gamma has both (alpha, gamma) and (gamma, beta) in O, and checks
# that for each regular orbital any two points are joined by a three-step
# walk whose three ordered pairs lie in that orbital.
#
# It then constructs the top group Q = {(a,b) in S4 x S4 : a*b^-1 in A4}
# in its transitive action of degree 12 and enumerates all 2^12 binary
# vectors: exactly 576 are distinguishing, all of weight six, forming
# two Q-orbits of length 288 interchanged by complementation.

Require := function(condition, message)
  if not condition then
    Error(message);
  fi;
end;;

K := PSL(2, 11);;
J := Normalizer(K, SylowSubgroup(K, 3));;
Require(Size(J) = 12 and IsDihedralGroup(J),
        "J should be dihedral of order 12");
G := Image(FactorCosetAction(K, J));;
Require(NrMovedPoints(G) = 55 and Size(G) = 660,
        "the coset action should be faithful of degree 55");
Print("PSL(2,11) constructed on 55 points\n");

H1 := Stabilizer(G, 1);;
regs := Filtered(Orbits(H1, [1..55]), o -> Length(o) = 12);;
Require(Length(regs) = 2, "there should be exactly two regular suborbits");

orbitalLabel := NullMat(55, 55);;
for a in [1, 2] do
  for pr in Orbit(G, [1, regs[a][1]], OnPairs) do
    orbitalLabel[pr[1]][pr[2]] := a;
  od;
od;

Require(ForAll([1..55], u -> ForAll([1..55],
        v -> orbitalLabel[u][v] = orbitalLabel[v][u])),
        "both regular orbitals should be self-paired");
Print("two self-paired regular orbitals found\n");

# (i): the regular-orbital obstruction.
found := fail;;
for alpha in [1..55] do
  for beta in [1..55] do
    if found = fail and alpha <> beta then
      for c in [1, 2] do
        if found = fail and ForAll([1..55],
             gamma -> not (orbitalLabel[alpha][gamma] = c
                           and orbitalLabel[gamma][beta] = c)) then
          found := [alpha, beta, c];
        fi;
      od;
    fi;
  od;
od;
Require(found <> fail, "(i) fails: no witness pair found");
Print("(i) verified: witness (alpha, beta, orbital label) = ", found, "\n");

# (ii): three-step walks inside each regular orbital.
for a in [1, 2] do
  A := NullMat(55, 55);
  for u in [1..55] do
    for v in [1..55] do
      if orbitalLabel[u][v] = a then A[u][v] := 1; fi;
    od;
  od;
  M := A * A * A;
  Require(ForAll(M, row -> ForAll(row, e -> e > 0)),
          "(ii) fails: a three-step walk is missing");
od;
Print("(ii) verified: three-step walks in each orbital\n");

# (iii): the distinguishing binary words of the top group Q.
S4 := SymmetricGroup(4);;
DP := DirectProduct(S4, S4);;
e1 := Embedding(DP, 1);;  e2 := Embedding(DP, 2);;
p1 := Projection(DP, 1);; p2 := Projection(DP, 2);;
Qgens := List(GeneratorsOfGroup(S4), g -> Image(e1, g) * Image(e2, g));;
Add(Qgens, Image(e1, (1,2,3)));;
QQ := Subgroup(DP, Qgens);;
Require(Size(QQ) = 288, "Q should have order 288");
Require(ForAll(QQ, q -> Image(p1, q) * Image(p2, q)^-1
                        in AlternatingGroup(4)),
        "Q should consist of the pairs (a,b) with a*b^-1 in A4");
dom := AsSSortedList(AlternatingGroup(4));;
Q12 := Image(ActionHomomorphism(QQ, dom,
         function(x, q) return Image(p1, q)^-1 * x * Image(p2, q); end));;
Require(Size(Q12) = 288 and IsTransitive(Q12, [1..12]),
        "Q should act faithfully and transitively on 12 points");
Print("top group Q of order 288 constructed on 12 points\n");

dws := Filtered(Combinations([1..12]),
                s -> IsTrivial(Stabilizer(Q12, s, OnSets)));;
Require(Length(dws) = 576 and ForAll(dws, s -> Length(s) = 6),
        "(iii) fails: there should be 576 distinguishing words, all of weight 6");
worbs := Orbits(Q12, dws, OnSets);;
Require(Length(worbs) = 2 and ForAll(worbs, o -> Length(o) = 288),
        "(iii) fails: the words should form two orbits of length 288");
Require(Set(worbs[1], s -> Difference([1..12], s)) = Set(worbs[2]),
        "(iii) fails: complementation should interchange the two orbits");
Print("(iii) verified: 576 words of weight 6, two orbits swapped by complementation\n");

Print("product-action input: all checks passed\n");
QUIT;
