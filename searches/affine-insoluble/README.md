# Insoluble affine stabilisers

`nonsolvable_targets.py` enumerates the prime-field degree and dimension
range. The three `export_magma_nonsolvable_*.m` files give the core group
constructors and generic maximal-subgroup descent, exporting survivors to the
affine graph engine.

The separate Aschbacher class C9 supplement consists of
`enumerate_c9_base2_brauer.g`, `export_magma_c9_base2_candidates.m`, and the
four `audit_c9_*` coverage checks. These are the
mathematical constructors and exporters used by the search. The production
calculation also used deterministic classical-group fast paths and
machine-specific orchestration for difficult branches. Those controllers,
along with generated frontier files, are not included.

The code requires Magma, GAP with CTblLib and AtlasRep, Python 3.10 or later,
and the affine C++ engine in [`../engines/`](../engines/).

The projective-degree lower bound instantiated by two of the C9 checks is an
external mathematical input from Landazuri and Seitz, [*On the minimal
degrees of projective representations of the finite Chevalley
groups*](https://doi.org/10.1016/0021-8693(74)90150-1). The outer-fusion check
uses Brauer-character induction to infer fusion and the absence of a linear
normaliser extension. It does not calculate matrix normalisers directly.
Finally, the C9 Magma exporter writes an action only for the two surviving
parameter pairs. Its header-only output for the other pairs is not input for
the affine engine.
