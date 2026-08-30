# Diagonal-type groups

This directory contains separate source for simple-diagonal and
compound-diagonal groups. The Magma programs enumerate the admissible shapes
and quotient groups and apply exact fixed-point-ratio and
partial-transversal bounds. The Python programs perform the full-wreath
partition and Hall reductions. The C++ programs resolve the residual
`A5`-component orbit-incidence cases.

- [`up-to-1e18/`](up-to-1e18/) contains the first range of simple- and
  compound-diagonal calculations.
- [`up-to-1e24/`](up-to-1e24/) extends the shape enumeration and bounds and
  contains the additional six-component exact calculations.
- [`shared/`](shared/) contains the finite-group action tables and coupled
  orbit engines reused by the exact C++ programs.

The source directories contain no enumerated quotient lists, saved outputs,
or certificates. Required inputs are produced by the preceding Magma source
in the same directory.

`up-to-1e24/sd_factorized_fpr_fixed_coordinate.m` is retained as an
independent reference calculation. The current
`a5_k3_cd_correlated_m6_exact.cpp` contains a corrected base-six encoding and
is not byte-identical to the earlier run source. It is published as the
corrected current implementation, not as an archival source snapshot. This
revision has not been rerun as part of the present source-only release.
