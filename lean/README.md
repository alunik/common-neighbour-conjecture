# Lean formalization

This directory contains the Lean 4 formalization of the main theorem in
Aluna Rizzoli and Adam R. Thomas, *Counterexamples to common-neighbour
conjectures for Saxl graphs*.

The reader-facing file is
[`Examples/MainTheorems.lean`](Examples/MainTheorems.lean). It displays only
the four non-Mathlib definitions needed to read the result and its single
theorem statement; its proof term is kept behind a one-line alias. The
remaining `.lean` files are the minimal transitive local module dependency
closure.

Certificate data, exploratory developments, computational audits, and
unrelated families are not included.

The project pins Lean `v4.34.0-rc2` and Mathlib commit
`87f6d5ec4c780581c9a78b06a9c5f1cf86dc5a70`. From this directory, run:

```sh
lake exe cache get
lake build
```

The Mathlib cache command avoids compiling Mathlib from source. The default
build target is only `Examples.MainTheorems` and its proof closure.
