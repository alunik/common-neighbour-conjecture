# Common-neighbour conjecture

This repository contains the computational code and Lean formalization
accompanying our work on the common-neighbour conjecture for Saxl graphs.

The repository is intentionally small and source-focused. It contains only
the standalone verification material used by the paper.

## Contents

- [`lemmas/`](lemmas/): the seed computations supporting the computational
  checks in the paper, one short GAP file per lemma.
- [`counterexamples/`](counterexamples/): construction and verification of
  the five affine counterexamples listed in the computational section.
- [`sporadic/`](sporadic/): computations for the ten cases in Table 1 of
  Lee--Popiel.
- [`lean/`](lean/): a self-contained Lean 4 formalization of Theorems 1.2 and
  1.3.
