# Computational searches

This directory contains the source code for the searches described in
Section 7 of the paper. It contains no generated actions, saved results, run
logs, or computational certificates.

The files are cleaned current implementations of the search algorithms. They
are not an archival copy of every controller or of the exact source snapshot
used on each compute node.

| range in the paper | source |
|---|---|
| primitive groups of degree below 8192 | [`primitive-catalogue/`](primitive-catalogue/) |
| affine groups with soluble stabiliser | [`affine-soluble/`](affine-soluble/) |
| affine groups with insoluble stabiliser | [`affine-insoluble/`](affine-insoluble/) |
| non-affine groups through degree `10^8` | [`nonaffine/`](nonaffine/) |
| simple- and compound-diagonal groups through degree `10^24` | [`diagonal/`](diagonal/) |

The shared exact graph engines are in [`engines/`](engines/). GAP and Magma
construct the relevant groups and export exact actions. The C++ programs
perform the larger orbit, sumset, and common-neighbour calculations. Python
programs enumerate degree ranges and coordinate multi-stage searches.

The source is arranged by mathematical role rather than by the machines on
which the calculations were performed. Cluster controllers, diagnostic
experiments, superseded implementations, and all generated material have
therefore been omitted.
