# Seed lemmas

One self-contained GAP file per computational check in the paper.

| file | lemma | contents |
|---|---|---|
| [`affine_seed.g`](affine_seed.g) | Lemma 3.1 | the monomial group `L0 = D:P` of order 1152 in `GL(9,3)`: irreducibility, the two regular orbits, their self-sum data, and the two-sum pattern of the all-ones vector |
| [`product_seed.g`](product_seed.g) | Lemma 4.1 | `PSL(2,11)` on 55 points: the two self-paired regular orbitals, the switching witness, the monochromatic three-step walks, and the 576 distinguishing binary words of the degree-12 top group `Q` |
| [`twisted_wreath_seed.g`](twisted_wreath_seed.g) | Lemma 5.1 | the 720 cycle-by-cycle fibre checks showing `z^p * w <> z` for the twisted wreath witness `w` in `A5^6` |

Run with, for example,

    gap -q affine_seed.g

Each file needs only core GAP, prints the checks it performs, stops with an
error if any check fails, and finishes in a few seconds.
