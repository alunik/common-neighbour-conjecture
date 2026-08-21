# Density computations

This directory contains the two sporadic cases settled by a density
argument.  In both cases, if `R` is the set of nonzero regular vectors in
`GF(q)^d`, the calculation proves

```text
|R| > q^(d-1) - 1.
```

The elementary projective-density lemma then gives `R + R = V`, which is the
common-neighbour condition for the corresponding affine group.

The two cases are:

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
