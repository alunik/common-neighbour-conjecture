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
