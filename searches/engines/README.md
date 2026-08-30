# Exact graph engines

`affine_saxl_engine.cpp` accepts prime-field matrix generators and computes
regular vector orbits, their two-sum, and graph distance by an
orbit-compressed breadth-first search.

`primitive_saxl_engine.cpp` accepts a finite transitive permutation action,
constructs its regular suborbits, and tests the common-neighbour condition on
point-stabiliser orbit representatives.

Both input formats and standalone C++17 build commands are documented at the
top of the corresponding source file.
