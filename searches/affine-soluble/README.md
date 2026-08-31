# Soluble affine stabilisers

These files implement the descent through H\"ofling's IRREDSOL catalogue.

- `build_irredsol_roots.g` and `build_irredsol_root_counts.g` enumerate the
  guardian groups covering the required prime-field actions.
- `export_irredsol_catalogue_round.g` exports one inclusion layer to the
  affine graph engine.
- `run_irredsol_catalogue_rounds.py` performs the downward traversal and
  applies the overgroup pruning described in the paper.
- `count_irredsol_dimension1.py` handles the one-dimensional cases and
  `summarize_irredsol_census.py` combines the source-generated records.

The code requires GAP 4.15.1 with IRREDSOL 1.4.4, Python 3.10 or later, and a
C++17 compiler. Manifests and results are generated locally and are not
stored in this repository. With this IRREDSOL version the complete root range
is selected by `--root-first 1 --root-last 2864`.

The census reported in the paper was completed on 5 August 2026. For
`8192 <= p^n <= 2^24 - 1` with `n >= 2`, all 2,864 guardian roots were
exhausted in the three disjoint ranges `1--1905`, `1906--2402`, and
`2403--2864`. The terminal states had empty queues and no resource-guard
failures. Together with the one-dimensional calculation above, the census
found exactly three soluble affine actions of diameter at least three: the
degree-19,683 action and the two degree-531,441 actions listed in
[`counterexamples/`](../../counterexamples/), where their diameter is
verified to be exactly three.
