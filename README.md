# AIT33 Region-Specific Single-Cell Reference Retrieval

Retrieve a balanced, region-specific single-cell RNA-seq reference from the Allen Institute **Whole Mouse Brain Consensus taxonomy (AIT33)** as an AnnData `.h5ad` file.

Expression comes from **Consensus-WMB-AIBS-10X** (raw counts reused from WMB-10Xv2/v3). Cell-type labels come from the integrated AIT33 taxonomy. Macosko snRNA-seq data is not included.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

List selectable anatomical regions:

```bash
python retrieve_ait33_reference.py --list-regions --cache-dir ./abc_atlas_cache
```

Build a cerebellum reference (~1000 cells per abundant subclass):

```bash
python retrieve_ait33_reference.py \
    --region cerebellum \
    --cells-per-subclass 1000 \
    --min-cells-per-subclass 500 \
    --seed 0 \
    --cache-dir ./abc_atlas_cache \
    --output cerebellum_ait33_reference.h5ad
```

Metadata-only dry run (no expression downloads):

```bash
python retrieve_ait33_reference.py \
    --region cerebellum \
    --metadata-only \
    --cache-dir ./abc_atlas_cache \
    --output cerebellum_ait33_reference.h5ad
```

## Outputs

```text
cerebellum_ait33_reference.h5ad
cerebellum_ait33_reference_sampling_summary.csv
```

`adata.X` contains raw counts (sparse). Taxonomy and anatomy columns live in `adata.obs`. Provenance is stored in `adata.uns["reference_info"]`.

## Notes

- First runs download Allen metadata and only the expression matrices needed for selected cells.
- Region matching is case-insensitive (`cerebellum` → `CB`).
- `--require-full-sample` forces every retained subclass to contribute exactly `--cells-per-subclass` cells.
