# The twisted wreath input (Lemma 5.1 of the paper).
#
# In the coordinate model of the twisted wreath product A5 twr S6, with
# transversal a_1 = (), a_i = (1 i) and q_{p,i} = a_i^-1 * p * a_{i^p},
# an element z of the socle N0 = T^6 (T = A5 on the points {2,...,6})
# satisfies z^p * w = z if and only if
#     z_{i^p} = z_i^{q_{p,i}} * w_{i^p}   for all i.
# On each cycle of p one entry of z determines all the others, so the
# equation can be decided by testing at most 60 initial values per cycle.
#
# The script performs this check for all 720 elements p of S6 and the
# witness w from the paper (whose entries are written there as
# permutations of {1,...,5} and are shifted here to {2,...,6}), and
# verifies that z^p * w = z never has a solution.

Require := function(condition, message)
  if not condition then
    Error(message);
  fi;
end;;

S6 := SymmetricGroup(6);;
a := [(), (1,2), (1,3), (1,4), (1,5), (1,6)];;

shift := MappingPermListList([1, 2, 3, 4, 5], [2, 3, 4, 5, 6]);;
wpaper := [(1,2,5), (1,3,4,2,5), (1,5,4), (1,2,4,3,5), (2,3,4), (3,4,5)];;
w := List(wpaper, x -> x^shift);;

T := AlternatingGroup([2..6]);;
Require(ForAll(w, x -> x in T), "the entries of w should lie in T");
Telts := AsSSortedList(T);;

for p in Elements(S6) do
  solvable := true;
  for cyc in Cycles(p, [1..6]) do
    len := Length(cyc);
    ok := false;
    for z0 in Telts do
      z := z0;
      i := cyc[1];
      good := true;
      for step in [1..len] do
        j := i^p;
        q := a[i]^-1 * p * a[j];
        Require(1^q = 1, "q_{p,i} should fix the point 1");
        znext := (z^q) * w[j];
        if step < len then
          z := znext;
          i := j;
        elif znext <> z0 then
          good := false;
        fi;
      od;
      if good then
        ok := true;
        break;
      fi;
    od;
    if not ok then
      solvable := false;
      break;
    fi;
  od;
  if solvable then
    Error("the equation z^p * w = z has a solution for p = ", p);
  fi;
od;

Print("twisted wreath input: all 720 cycle checks passed\n");
QUIT;
