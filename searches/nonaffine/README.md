# Non-affine primitive groups

The three subdirectories correspond to successive degree windows. Their
Magma programs generate almost simple, product-action, simple-diagonal, and
compound-diagonal candidates and their exact degree inventories. C++
programs handle the larger structured actions and orbit-incidence
calculations.

- [`up-to-1e6/`](up-to-1e6/) contains the initial type-by-type generators and
  the separate O'Nan order-obstruction check.
- [`up-to-1e7/`](up-to-1e7/) adds the structured product-action exporter and
  the exceptional `PSL(2,7)` diagonal calculation.
- [`up-to-1e8/`](up-to-1e8/) contains the almost-simple packed-action code,
  the product-action exporters, and the simple- and compound-diagonal
  calculations needed in the final window.

The common-neighbour tests use the primitive engine in [`../engines/`](../engines/).
Some exporters accept global arrays prepared by the preceding inventory
script, as documented in their source headers. Generated action streams and
result files are deliberately omitted.
