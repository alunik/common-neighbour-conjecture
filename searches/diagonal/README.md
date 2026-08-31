# Diagonal-type groups

This directory uses the five-family O'Nan--Scott convention used by Huang and
by Liebeck, Praeger and Saxl. Thus diagonal type here is the family often
called simple diagonal (SD) in the refined terminology. Product actions with
a diagonal-type component, often called compound diagonal (CD), belong to
product type in this convention and are not included here. The paper already
provides product-type counterexamples, so no positive claim about that family
is intended in this directory.

The Magma programs enumerate the admissible SD shapes and apply exact
fixed-point-ratio bounds. The seven residual cases, with `T = A5` and
`8 <= k <= 14`, are resolved by exact rigid-subset and partial-transversal
bounds.

- [`up-to-1e18/`](up-to-1e18/) contains the first window,
  `10^8 < n <= 10^18`.
- [`up-to-1e24/`](up-to-1e24/) contains the extension window,
  `10^18 < n <= 10^24`.
- [`shared/`](shared/) contains exact small diagonal-action tables and
  common-neighbour engines, some of which are also reused by the low-degree
  non-affine search.

The source directories contain no generated shape inventories, saved
outputs, or certificates. `up-to-1e24/sd_factorized_fpr_fixed_coordinate.m`
is retained as an independent reference calculation.
