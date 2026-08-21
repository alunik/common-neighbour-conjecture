# The two 12-dimensional modules for 6.Suz over F_19

This program handles the `Z o (6.Suz), d = 12, r = 19` row of Table 1 in
Lee and Popiel, *Saxl graphs of primitive affine groups with sporadic point
stabilisers*.  There are two contragredient absolutely irreducible modules,
and the program treats both of them.

For each module, the program constructs a 32760-point orbit of projective
functionals.  A projective stabiliser of a vector must preserve the subset of
these functionals which vanish on it.  GAP computes the set stabiliser for
each of fifteen displayed representatives and finds it to be trivial.  The
fifteen zero-set sizes are different, so the representatives lie in fifteen
different projective orbits of regular vectors.

The centre of `6.Suz` has order 6.  After adjoining all scalars from
`F_19^*`, each regular projective orbit therefore contains

```
|6.Suz| * 18 / 6 = 8070218956800
```

regular nonzero vectors.  The program finishes by checking

```
15 * 8070218956800 = 121053284352000 > 19^11 - 1.
```

The projective density argument then gives `R + R = V`.  Regularity for the
full scalar group also implies regularity for every smaller allowed scalar
subgroup, so this handles every `Z` in the Table 1 row.

## Running the computation

The only requirements are a C++17 compiler and GAP with the AtlasRep package.
From this directory, run

```bash
./run.sh
```

The C++ program builds the projective actions.  GAP checks those actions,
computes the thirty set stabilisers, and prints the density conclusion.  The
intermediate orbit data are placed in a temporary directory and removed when
the program exits.

The matrices are in `matrices.g`, the thirty representatives are in
`vectors.g`, and the mathematical checks are in `compute.g`.
