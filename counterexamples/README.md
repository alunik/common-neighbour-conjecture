# The five affine counterexamples

One self-contained GAP file per counterexample listed in the computational
section of the paper.  Each file constructs the point stabiliser
`H <= GL(d, 3)`, and the shared verifier [`verify.g`](verify.g) then checks
that `H` is irreducible with exactly one regular orbit `R` on `V = GF(3)^d`,
that exactly the stated number of vectors lie outside `R + R`, and that all
of them lie in `R + R + R`, so that the Saxl graph of `V:H` has diameter
exactly 3.

| file | group | outside `R + R` |
|---|---|---|
| [`degree_19683.g`](degree_19683.g) | `F_3^9 : (C_2^6 : D18)` | 96 |
| [`degree_59049.g`](degree_59049.g) | `F_3^10 : (C_2^5 : S_5)` | 64 |
| [`degree_531441_12T35.g`](degree_531441_12T35.g) | `F_3^12 : (C_2^8 : (12T35))` | 1600 |
| [`degree_531441_12T38.g`](degree_531441_12T38.g) | `F_3^12 : (C_2^8 : (12T38))` | 1600 |
| [`degree_14348907.g`](degree_14348907.g) | `F_3^15 : AGL_4(2)` | 32 |

Run from this directory with, for example,

    gap -q degree_19683.g

Each file needs only core GAP, prints the checks it performs, stops with an
error if any check fails, and finishes in seconds (about twenty seconds for
the last one).
