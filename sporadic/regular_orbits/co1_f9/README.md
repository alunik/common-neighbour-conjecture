# `2.Co1` on `GF(9)^24`

Let `K = 2.Co1` in its standard 24-dimensional representation over
`GF(3)`, extended to `GF(9)`, and let

```text
Khat = <K, GF(9)^*>.
```

This program finds 267 distinct regular `Khat`-orbits.  Since

```text
|Khat|                       = 33,262,214,452,346,880,000
267 |Khat|                   =  8,881,011,258,776,616,960,000
9^23 - 1                    =  8,862,938,119,652,501,095,928
```

the set of regular vectors used by the program has size greater than
`9^23 - 1`.  The projective density argument then gives
`R_Khat + R_Khat = GF(9)^24`.  The same vectors work for every scalar
subgroup occurring in the corresponding row of Lee--Popiel's Table 1.

## What the program does

The file `representation.txt` contains two matrices for the standard
`GF(3)` representation, a vector fixed by an `M24` complement, and a dual
seed.  The file `m24_action.g` contains the signed permutation action of
that complement.

`search.cpp` constructs the dual orbit of length 196560 and generates a
deterministic stream of vectors `w`.  It retains `w` when a zero-count
invariant distinguishes the fixed line `<v>` from the other three lines in
`<v,w>`.  For the vectors `v + omega*w`, this reduces the remaining tests
to the signed `M24` action.

`check.g` uses GAP to compute every stabiliser and every canonical orbit
representative exactly.  It stops after finding 267 pairwise distinct
regular orbits and prints the density calculation.  There are no saved
runtime files.

## Run

Use a C++17 compiler and GAP with the Images package.  The full calculation
should be run on CREATE:

```bash
GAP_BIN=/path/to/gap CXX=c++ ./run.sh
```

The final lines include

```text
regular_orbits 267
regular_vectors 8881011258776616960000
projective_threshold 8862938119652501095928
strict_margin 18073139124115864072
regular_plus_regular_is_whole_space true
```
