# Primitive catalogue below degree 8192

`build_gap_nonaffine_manifest.g` enumerates the non-affine actions in GAP's
PrimGrp catalogue. `export_gap_nonaffine_shard.g` exports those actions for
the primitive graph engine. `export_magma_affine_catalogue.m` exports the
affine catalogue actions for the affine engine.

The GAP scripts require PrimGrp, including the extended catalogue data for
degrees 4096 through 8191. The exporters write the text formats documented in
[`../engines/`](../engines/). The calculations used GAP 4.15.1 and PrimGrp
4.0.1.
