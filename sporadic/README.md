# Sporadic point stabilisers

Let $H\leqslant \operatorname{GL}(V)$, and let

$$
R_H=\{v\in V:H_v=1\}.
$$

For an affine group $V{:}H$ with base size two, the common-neighbour
condition is equivalent to $R_H+R_H=V$.  The programs in this directory
perform the calculations for the ten cases left open in Table 1 of
Lee--Popiel.

The ten cases split into six settled by explicit orbit-and-sumset
calculations and four settled by the projective-density bound.  The four
density cases occur in two computational forms:

- a prime-order eigenspace count followed by a projective-density bound;
- a regular-orbit count followed by the same density bound;

The remaining calculations use:

- a complete orbit enumeration followed by a sumset calculation;
- larger orbit calculations for the `Co1` and `Fi22` modules.

The two shorter density calculations are under `density/`.  The other two
are under `regular_orbits/`, since their main step is to certify enough
distinct regular projective orbits before applying the same density bound.

Each directory contains commented source code and brief running instructions.
The programs print the relevant orbit counts, density bounds, or sumset result.
The exhaustive jobs should be run on CREATE.

## Cases

| program or directory | group and module | method |
| --- | --- | --- |
| `density/co3_f5.g` | `Co3` on `GF(5)^23` | prime-order eigenspaces, then density |
| `density/suz3_f16.g` | `3.Suz` on `GF(16)^12` | projective regular orbits, then density |
| `orbits/co3_f4` | `Co3` on `GF(4)^22` | complete orbit and sumset enumeration |
| `orbits/suz3d2_f4` | `3.Suz.2` on `GF(4)^24` | complete orbit and sumset enumeration |
| `large/co1_f8` | `Co1` on `GF(8)^24` | hierarchical orbit enumeration |
| `regular_orbits/co1_f9` | `2.Co1` on `GF(9)^24` | projective regular-orbit count, then density |
| `large/fi22` | the three `Fi22` rows | orbit enumeration and exact stabilisers |
| `regular_orbits/suz6_f19` | `6.Suz` on `GF(19)^12` | projective regular-orbit count, then density |
