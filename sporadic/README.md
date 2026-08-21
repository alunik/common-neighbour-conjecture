# Sporadic point stabilisers

Let $H\leqslant \operatorname{GL}(V)$, and let

$$
R_H=\{v\in V:H_v=1\}.
$$

For an affine group $V{:}H$ with base size two, the common-neighbour
condition is equivalent to $R_H+R_H=V$.  The programs in this directory
perform the calculations for the ten cases left open in Table 1 of
Lee--Popiel.

There are two short calculations and two larger computational patterns:

- a prime-order eigenspace count followed by a projective-density bound;
- a complete orbit enumeration followed by a sumset calculation;
- a regular-orbit count followed by the same density bound;
- larger orbit calculations for the `Co1` and `Fi22` modules.

Each directory contains commented source code and brief running instructions.
The programs print the relevant orbit counts, density bounds, or sumset result.
The exhaustive jobs should be run on CREATE.

## Cases

| program or directory | group and module | method |
| --- | --- | --- |
| `density/co3_f5.g` | `Co3` on `GF(5)^23` | prime-order eigenspaces |
| `density/suz3_f16.g` | `3.Suz` on `GF(16)^12` | projective regular orbits |
| `orbits/co3_f4` | `Co3` on `GF(4)^22` | complete orbit and sumset enumeration |
| `orbits/suz3d2_f4` | `3.Suz.2` on `GF(4)^24` | complete orbit and sumset enumeration |
| `large/co1_f8` | `Co1` on `GF(8)^24` | hierarchical orbit enumeration |
| `regular_orbits/co1_f9` | `2.Co1` on `GF(9)^24` | projective regular-orbit count |
| `large/fi22` | the three `Fi22` rows | orbit enumeration and exact stabilisers |
| `regular_orbits/suz6_f19` | `6.Suz` on `GF(19)^12` | projective regular-orbit count |
