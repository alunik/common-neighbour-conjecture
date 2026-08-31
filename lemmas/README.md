# Finite computations

Self-contained GAP files for the computational inputs in Sections 3--5.

| file | lemma | contents |
|---|---|---|
| [`affine_input.g`](affine_input.g) | Lemma 3.3 | the monomial group `L0 = D:P` of order 1152 in `GL(9,3)`: irreducibility, the two regular orbits, their self-sum data, and the two-sum pattern of the all-ones vector |
| [`product_input.g`](product_input.g) | Lemma 4.1 | `PSL(2,11)` on 55 points: its two self-paired regular orbitals, the obstruction in part (i), the three-step walks in part (ii), and the 576 distinguishing binary vectors for the degree-12 top group `Q` |
| [`twisted_wreath_census.g`](twisted_wreath_census.g) | Lemma 5.1 | a table-of-marks census giving the exact number `64,790,243` of regular suborbits of `A5 twr S6` |
| [`twisted_wreath_input.g`](twisted_wreath_input.g) | Lemma 5.1 | the 720 cycle-by-cycle recurrence checks showing `z^p * w <> z` for the twisted wreath witness `w` in `A5^6` |

Run with, for example,

    gap -q affine_input.g

Each file prints the checks it performs and stops with an error if any check
fails. The census uses the TomLib package, which is distributed with GAP.
The other files need only core GAP.
