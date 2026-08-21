# `3.Suz.2` on `GF(4)^24`

The program writes a vector over `GF(4)` as a pair of binary 24-vectors.  It
enumerates all 643 vector orbits, finds the regular orbits for both
`3.Suz.2` and its full scalar extension, and computes a regular-plus-regular
expression for every target orbit.

Run on CREATE with GAP (including AtlasRep) and a C++20 compiler:

```bash
GAP_BIN=/path/to/gap CXX=c++ ./run.sh
```

The final output includes

```text
GF4_vector_orbits 643
regular_G_vector_orbits 7
regular_Hmax_vector_orbits 2
G_sumset_equals_GF4_24 true
Hmax_sumset_equals_GF4_24 true
```
