# `Co3` on `GF(4)^22`

The program writes a vector over `GF(4)` as an ordered pair of binary
22-vectors.  It enumerates the five `Co3`-orbits on the binary module, the
320 orbits on the full `GF(4)`-module, and a regular-plus-regular expression
for every target orbit for `Co3` itself.  Thus `R + R = V` when the scalar
subgroup is trivial.  The same run also finds that `GF(4)^* x Co3` has no
regular vector, so that larger group is not a base-two case.

Run on CREATE with GAP (including AtlasRep and CTblLib) and a C++17 compiler:

```bash
GAP_BIN=/path/to/gap CXX=c++ ./run.sh
```

The final output includes

```text
total_pair_orbits 320
total_regular_pair_orbits 2
regular_full_scalar_vectors 0
all_target_orbits_have_verified_decomposition true
sumset_equals_GF4_22 true
```
