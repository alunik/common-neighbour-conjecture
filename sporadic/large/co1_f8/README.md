# `Co1` on `F8^24`

This program handles the action of

```text
H = F_8^* x Co1  on  V = F_8^24.
```

Write the 24-dimensional binary `Co1` module as `W`.  Relative to the basis
`1, eta, eta^2` of `F_8` over `F_2`, a vector of `V` is a triple in `W^3` and
multiplication by `eta` is

```text
(a,b,c) -> (c,a+c,b).
```

The program first enumerates the `Co1`-orbits on `W^3`.  It does this one
coordinate at a time: there are four orbits on `W`, 46 on `W^2`, and each
pair stabilizer acts on the final copy of `W`.  The resulting list has 9,511
orbits, 119 of which have trivial `Co1` stabilizer.  Following multiplication
by `eta` identifies 17 regular `H`-orbits.

For each of the 9,511 target orbits, the final stage chooses regular vectors
`y`, classifies `x-y`, and checks that both summands are regular.  It finishes
by printing

```text
R + R = V
```

where `R` is the set of regular vectors for `H`.  Since `F_8^*` has prime
order 7, the same computation also handles every scalar subgroup in this
case.

## Running the computation

The requirements are GAP with the `atlasrep` package, a C++17 compiler, and
Python 3.  No precomputed orbit data is used.  Intermediate files are made in
a fresh temporary directory and removed at the end.

This is a large computation and should be submitted on CREATE:

```sh
cd sporadic/large/co1_f8
sbatch run_create.slurm
```

The shell driver is `run.sh`.  `GAP`, `CXX`, `JOBS`, and `TRIALS` may be set
in the environment; the defaults are `gap`, `g++`, four workers, and 128
summand trials per target.
