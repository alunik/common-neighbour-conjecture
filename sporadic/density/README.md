# Density computations

This directory contains two short density calculations.  Two further cases
use regular-orbit counts followed by the same projective-density argument;
they are under [`regular_orbits/co1_f9`](../regular_orbits/co1_f9/) and
[`regular_orbits/suz6_f19`](../regular_orbits/suz6_f19/).  Thus four of the
ten sporadic cases are settled by density in total.

In every density case, if `R` is the set of nonzero regular vectors in
`GF(q)^d`, the calculation proves

```text
|R| > q^(d-1) - 1.
```

The elementary projective-density lemma then gives `R + R = V`, which is the
common-neighbour condition for the corresponding affine group.

The two cases housed in this directory are:

- `Co3` on its 23-dimensional module over `GF(5)`;
- `3.Suz` on its 12-dimensional module, extended from `GF(4)` to `GF(16)`.

The C++ programs contain the short integer calculations and print the final
conclusion.  Run them with:

```bash
./run.sh
```

The GAP programs recompute the representation-theoretic inputs from
AtlasRep.  They require GAP with the AtlasRep and CTblLib packages and take
considerably longer:

```bash
./run.sh --full
```

`GAP_BIN` and `CXX` may be used to select a GAP executable or C++ compiler.
