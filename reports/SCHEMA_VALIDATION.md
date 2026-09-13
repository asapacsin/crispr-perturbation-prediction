# Schema validation record

Validated on 2026-09-12 after Cursor Grok Fast implementation and independent orchestrator review.

- `python -m unittest discover -s tests -q`: **37 tests passed**. Covers label canonicalization, missing labels, categorical HDF5 layouts and invalid codes, unknown partitions, target leakage, and bridge exclusions.
- `python scripts/validate_dataset_schema.py`: **SCHEMA VALIDATION PASSED** against the actual H5AD, including file SHA-256, metadata, target IDs, source annotation checks, and proposed split isolation.
- `python scripts/write_verified_guide_sequences.py`: generated 105 rows; 103 exact matches, one documented dataset-local RHOXF2 mapping, one unresolved SLC38A2 sequence. All shared source sequences agree and supplied sequences contain 20 ACGT bases. Pandas emitted a non-failing `verify_integrity` deprecation warning.
- Independent checks: all 105 schema target-to-Ensembl mappings agree with the target table; split and joined-metadata barcodes are unique and match exactly (91,205); all four proposals' cell/control totals agree with their manifests.
- Dataset SHA-256 remains `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`.

The initial validator failed on the legacy categorical layout; the corrected reader passes the real file and regression fixtures. An additional check rejects controls incorrectly marked as excluded validation/test gene bridges.

This validates the documented prepared dataset and proposed manifests; it does not reproduce the full upstream count matrix or establish a leakage-free prospective feature-selection pipeline. No model was trained. Final protocol and training approval remain false. Worker handoffs predate these final root checks.
