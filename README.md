# Common-neighbour conjecture

This repository contains computational source code and a Lean formalization
accompanying our paper on the common-neighbour conjecture for Saxl graphs.
It is intentionally source-focused. Generated outputs, run logs, and
computational certificates are not included.

## Contents

- [`lemmas/`](lemmas/): the finite GAP computations used in Sections 3--5.
- [`counterexamples/`](counterexamples/): the four affine counterexamples in
  Section 7 and the perfect-stabiliser example in Section 3.
- [`searches/`](searches/): source code for the searches reported in Section 7.
- [`sporadic/`](sporadic/): computations for the ten cases in Table 1 of
  Lee--Popiel.
- [`lean/`](lean/): a self-contained Lean 4 formalization of Theorems 1.2 and
  1.3.

The search code uses GAP and Magma to construct groups and export actions,
Python for orchestration and exact bookkeeping, and C++ for the larger exact
orbit and graph calculations. Individual directories list the packages and
language versions they require.
