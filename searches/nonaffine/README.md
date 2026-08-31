# Non-affine primitive groups

The three subdirectories correspond to successive degree windows. Their
Magma programs generate almost simple, diagonal-type, and product-action
candidates and their exact degree inventories. The product-action lane
includes actions often called compound diagonal (CD) in the refined
terminology. C++ programs handle the larger structured actions and
orbit-incidence calculations.

- [`up-to-1e6/`](up-to-1e6/) contains the initial type-by-type generators and
  the separate O'Nan order-obstruction check.
- [`up-to-1e7/`](up-to-1e7/) adds the structured product-action exporter and
  the exceptional `PSL(2,7)` diagonal calculation.
- [`up-to-1e8/`](up-to-1e8/) contains the almost-simple packed-action code,
  the product-action exporters, and the diagonal-type and product-action
  calculations with diagonal components needed in the final window.

Files whose names contain `cd` implement these product actions with diagonal
components; they are retained here under product type.

The common-neighbour tests use the primitive engine in [`../engines/`](../engines/).
Some exporters accept global arrays prepared by the preceding inventory
script, as documented in their source headers. Generated action streams and
result files are deliberately omitted.
