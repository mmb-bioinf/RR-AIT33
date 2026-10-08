# AIT33 Region-Specific Single-Cell Reference Retrieval

## Goal

Implement a Python script that retrieves a **region-specific, balanced single-cell reference dataset** from the Allen Institute **Whole Mouse Brain Consensus taxonomy (AIT33)**.

The intended downstream use is reference-informed / weighted SNMF benchmarking.

The user should be able to select an anatomical brain region, for example:

```text
cerebellum
```

The script should then create an `.h5ad` file containing approximately:

```text
1000 single-cell expression profiles per AIT33 subclass
```

for subclasses that are sufficiently abundant in the selected region.

Rare subclasses should be excluded.

---

# Data source

Use the current Allen Brain Cell Atlas Python infrastructure:

```python
abc_atlas_access
```

with:

```python
from abc_atlas_access.abc_project_cache import AbcProjectCache
```

Use the **AIT33 / 2025 Whole Mouse Brain Consensus taxonomy**:

```text
20251030_BICCN-Adult_Mouse_WholeBrain_consensus
AIT33
```

The relevant Allen data packages are:

```text
Consensus-WMB-AIBS-10X
Consensus-WMB-integrated-taxonomy
```

For this script, use **only `Consensus-WMB-AIBS-10X` expression data**.

Do not mix the Macosko single-nucleus data into this reference dataset.

The reason is that we want a relatively homogeneous **single-cell RNA-seq reference** for downstream deconvolution.

The integrated AIT33 taxonomy should nevertheless be used for the cell-type labels.

---

# High-level workflow

The script should implement:

```text
AIT33 metadata
    │
    ├── obtain anatomical region annotations
    ├── obtain AIT33 taxonomy assignments
    │
    ▼
user selects anatomical region
    │
    ▼
identify AIBS cells belonging to region
    │
    ▼
determine AIT33 subclasses represented in region
    │
    ▼
count cells per subclass
    │
    ├── remove rare subclasses
    │
    ▼
sample ~1000 cells per retained subclass
    │
    ▼
retrieve raw gene-expression profiles
    │
    ▼
construct AnnData
    │
    ▼
region_reference.h5ad
```

---

# 1. Command-line interface

Implement the script as a command-line program.

Example:

```bash
python retrieve_ait33_reference.py \
    --region cerebellum \
    --cells-per-subclass 1000 \
    --min-cells-per-subclass 500 \
    --output cerebellum_ait33_reference.h5ad
```

Important arguments:

```text
--region
    Anatomical region selected by the user.

--cells-per-subclass
    Target number of cells sampled for every retained subclass.
    Default: 1000

--min-cells-per-subclass
    Minimum number of cells from that subclass that must occur
    in the selected region for the subclass to be retained.
    Default: 500

--seed
    Random seed for reproducible sampling.
    Default: 42

--output
    Output .h5ad filename.

--cache-dir
    Directory used by AbcProjectCache.
```

Optional useful arguments:

```text
--list-regions
--overwrite
--metadata-only
```

`--list-regions` should print valid selectable regions and exit.

---

# 2. Region selection

Do **not** hard-code a small list such as:

```python
["cerebellum", "cortex", "thalamus"]
```

Instead, determine the available anatomical region labels from the Allen metadata and expose them to the user.

For example:

```bash
python retrieve_ait33_reference.py --list-regions
```

might return something resembling:

```text
Cerebellum
Hippocampal formation
Hypothalamus
Isocortex
Midbrain
Olfactory areas
Pallidum
Pons
Striatum
Thalamus
...
```

The exact names must come from the downloaded metadata rather than being assumed.

Region matching should be case-insensitive.

Prefer:

```python
"cerebellum"
```

matching:

```text
"Cerebellum"
```

If several anatomical metadata levels exist, retain both the coarse and fine anatomical annotations in the output.

---

# 3. Distinguish anatomical filtering from storage regions

Allen expression matrices are divided into coarse anatomical packages.

These file divisions should **not automatically define biological inclusion criteria**.

For example:

```text
CB
HPF
HY
Isocortex
MB
...
```

may determine which expression file contains a cell.

The actual decision:

```text
Does this cell belong to the requested anatomical region?
```

should preferably be made from the appropriate anatomical metadata.

Therefore distinguish:

```text
expression_matrix_region
```

from:

```text
biological_anatomical_region
```

The former determines which `.h5ad` file needs to be accessed.

The latter determines whether the cell should enter the output dataset.

---

# 4. AIT33 taxonomy mapping

Use:

```text
Consensus-WMB-integrated-taxonomy
```

to map AIBS cells onto the consensus taxonomy.

Relevant taxonomy levels are:

```text
neighborhood
class
subclass
supertype
cluster
```

The primary target label for this dataset is:

```text
subclass
```

Each selected cell should therefore ultimately have at least:

```text
AIT33 subclass
AIT33 class
AIT33 supertype
AIT33 cluster
```

stored in `adata.obs`.

Use the official taxonomy membership tables rather than trying to infer taxonomy membership from label strings.

Relevant metadata currently include objects such as:

```text
cell_to_cluster_membership
cluster
cluster_annotation_term
cluster_annotation_term_set
cluster_to_cluster_annotation_membership
```

Construct an explicit cell-to-taxonomy mapping from these tables.

Do not rely on taxonomy row order.

Do not infer parent-child relationships from human-readable names.

---

# 5. Determine subclasses present in the selected region

After selecting cells belonging to the requested anatomical region, calculate:

```python
region_counts = metadata.groupby("subclass").size()
```

Conceptually:

\[
N_{k,r}
=
\left|
\left\{
c :
\operatorname{subclass}(c)=k
\land
\operatorname{region}(c)=r
\right\}
\right|
\]

where:

- \(k\) = AIT33 subclass
- \(r\) = selected anatomical region.

---

# 6. Filter low-abundance subclasses

Default:

```python
MIN_CELLS_PER_SUBCLASS = 500
```

Only retain subclasses satisfying:

\[
N_{k,r} \geq N_{\min}.
\]

For example:

```python
valid_subclasses = region_counts[
    region_counts >= min_cells_per_subclass
].index
```

The threshold must be configurable.

This filtering is intentional.

The goal is **not** to represent every extremely rare subclass. The goal is to construct a stable reference dataset containing cell populations that are meaningfully represented in the chosen tissue region.

Print a summary such as:

```text
Selected region: Cerebellum

AIT33 subclasses observed:       47
Subclasses with >=500 cells:     31
Subclasses excluded:             16
```

Also print the excluded subclasses and their counts.

---

# 7. Balanced sampling

For every retained subclass, randomly sample approximately:

```text
1000 cells
```

without replacement.

Default:

```python
CELLS_PER_SUBCLASS = 1000
```

Sampling rule:

\[
n_k
=
\min(N_{k,r},1000).
\]

More generally:

```python
n_select = min(
    cells_per_subclass,
    number_of_available_cells
)
```

Because the default minimum abundance is 500, some subclasses may contribute between 500 and 999 cells.

Subclasses with >=1000 available cells should contribute exactly 1000 cells.

Use:

```python
random_state=seed
```

for reproducibility.

---

# 8. Avoid unnecessary expression-data downloads

The complete AIT33 dataset is extremely large.

Do **not** download the full AIT33 expression matrix.

The program should:

1. download/load metadata first;
2. identify the exact selected cell IDs;
3. determine which expression matrix files contain those cells;
4. access/download only the necessary AIBS expression files;
5. subset those matrices to the selected cells;
6. concatenate the resulting subsets.

This is an important implementation requirement.

The script should remain practical on a workstation.

---

# 9. Expression values

Retrieve **raw count expression** wherever possible.

Do not log-transform the data.

Do not normalize expression during the retrieval step.

The resulting:

```python
adata.X
```

should contain raw counts.

Before writing the file, verify that values are compatible with counts:

```python
X.min() >= 0
```

and sampled non-zero values should be integer-valued or numerically extremely close to integers.

Avoid accidentally retrieving:

```text
log2(CPM + 1)
```

instead of counts.

The Allen data release explicitly distinguishes `raw` and `log2` expression matrices, so select the `raw` matrix.

---

# 10. Gene representation

Use the complete available gene set from the selected AIBS expression matrices.

Genes should be consistent across all concatenated anatomical expression files.

`adata.var` should preserve at least:

```text
gene_identifier
gene_symbol
```

if both fields are available.

Gene identifiers must be unique.

Prefer stable gene IDs as `var_names` if duplicate gene symbols occur.

Do not silently remove genes simply because their expression is zero in a particular sampled subclass.

---

# 11. Required `adata.obs`

The final AnnData should contain useful metadata for every selected cell.

At minimum:

```text
cell_id
class
subclass
supertype
cluster
anatomical_region
```

If available, also retain:

```text
neighborhood
neurotransmitter
cluster_alias
donor
sex
library
dataset
feature_matrix_label
coarse_anatomical_region
fine_anatomical_region
```

The exact Allen column names may differ.

Do not invent missing metadata fields.

Instead, map available official columns to clear output column names where appropriate.

---

# 12. Preserve the original cell IDs

The original Allen cell identifier must be retained.

Prefer:

```python
adata.obs_names = original_cell_ids
```

and verify:

```python
adata.obs_names.is_unique
```

Do not replace cell IDs with sequential numbers.

---

# 13. `adata.uns` provenance

Store retrieval information in:

```python
adata.uns["reference_info"]
```

Example:

```python
adata.uns["reference_info"] = {
    "taxonomy": "AIT33",
    "taxonomy_name": "20251030_BICCN-Adult_Mouse_WholeBrain_consensus",
    "expression_source": "Consensus-WMB-AIBS-10X",
    "selected_region": region,
    "taxonomy_level": "subclass",
    "cells_per_subclass_target": 1000,
    "min_cells_per_subclass": 500,
    "random_seed": 42,
    "expression_scale": "raw_counts",
}
```

Also store:

```text
retrieval date
abc_atlas_access version
Allen manifest/release version
```

when available.

This provenance is important because the Allen data releases can change.

---

# 14. Sampling summary

Store a table with statistics for every subclass.

For example:

| subclass | cells_in_region | retained | cells_sampled |
|---|---:|---:|---:|
| subclass_A | 12,450 | True | 1,000 |
| subclass_B | 834 | True | 834 |
| subclass_C | 112 | False | 0 |

This can either be:

```python
adata.uns["subclass_sampling"]
```

or written as an additional:

```text
<output_prefix>_sampling_summary.csv
```

Prefer writing the CSV in addition to storing summary information in AnnData.

---

# 15. Output structure

Example:

```text
output/
├── cerebellum_ait33_reference.h5ad
└── cerebellum_ait33_reference_sampling_summary.csv
```

Expected AnnData dimensions might look like:

```text
AnnData object with n_obs × n_vars = 28000 × 21205
```

if 28 sufficiently abundant subclasses were identified and approximately 1,000 cells were sampled from each.

The number of subclasses must be determined from the data rather than assumed.

---

# 16. Output summary

At completion, print something like:

```text
AIT33 reference successfully generated

Region:
    Cerebellum

Dataset:
    Consensus-WMB-AIBS-10X

Taxonomy:
    AIT33

Taxonomy level:
    subclass

Subclasses retained:
    31

Cells:
    29,842

Genes:
    21,205

Median cells/subclass:
    1,000

Minimum cells/subclass:
    612

Maximum cells/subclass:
    1,000

Expression:
    raw counts

Output:
    cerebellum_ait33_reference.h5ad
```

---

# 17. Validation checks

Before writing the `.h5ad`, automatically verify:

## Cell IDs

```python
assert adata.obs_names.is_unique
```

## Gene IDs

```python
assert adata.var_names.is_unique
```

## No missing subclass labels

```python
assert not adata.obs["subclass"].isna().any()
```

## Correct region

Verify that all selected cells satisfy the intended region-selection rule.

## Minimum subclass abundance

Every retained subclass must have passed the original regional abundance filter.

## Maximum sampled cells

```python
assert (
    adata.obs["subclass"].value_counts()
    <= cells_per_subclass
).all()
```

## Count scale

Verify that the matrix is non-negative and count-like.

## Sparse representation

The final expression matrix should remain sparse.

Do not inadvertently convert a large sparse matrix to a dense NumPy array.

---

# 18. Important memory requirement

Never execute operations such as:

```python
adata.X.toarray()
```

on the entire expression matrix.

Keep the expression matrix as:

```python
scipy.sparse
```

throughout retrieval, subsetting, concatenation, and output.

The script must be capable of processing tens of thousands of cells × ~20,000 genes without excessive RAM consumption.

---

# 19. Region-name handling

Make region selection user-friendly.

Support:

```bash
--list-regions
```

and case-insensitive exact matching.

Optionally allow aliases such as:

```text
cerebellum -> CB
hippocampus -> HPF
cortex -> Isocortex
```

but aliases must map onto official metadata values.

Do not use an alias to silently broaden the biological region.

When a region name is ambiguous, print the matching official regions and exit rather than selecting one arbitrarily.

For example:

```text
Region 'hippocampus' matched multiple Allen regions:

1. Hippocampal formation
2. Hippocampal region
3. Dentate gyrus

Please select one explicitly.
```

---

# 20. Separate metadata selection from expression loading

Structure the implementation into reusable functions.

Suggested architecture:

```python
def initialize_cache(cache_dir):
    ...

def load_ait33_taxonomy(cache):
    ...

def load_aibs_cell_metadata(cache):
    ...

def build_cell_taxonomy_table(...):
    ...

def get_available_regions(metadata):
    ...

def resolve_region(region_query, available_regions):
    ...

def select_region_cells(metadata, region):
    ...

def filter_subclasses(
    metadata,
    min_cells_per_subclass
):
    ...

def sample_cells_per_subclass(
    metadata,
    cells_per_subclass,
    seed
):
    ...

def determine_required_expression_files(
    selected_cells
):
    ...

def load_selected_expression(
    cache,
    selected_cells,
    required_files
):
    ...

def build_anndata(
    expression,
    metadata,
    gene_metadata
):
    ...

def validate_reference(adata, config):
    ...

def write_outputs(adata, summary, output):
    ...
```

Avoid putting the entire workflow in a single function.

---

# 21. Important interpretation of "subclass contained in a region"

A subclass should be considered present in the region based on **cell-level anatomical metadata**, not based merely on its general anatomical annotation in the taxonomy description.

In other words, use:

\[
\text{cells observed in region}
\]

rather than:

\[
\text{subclass annotated as generally associated with region}.
\]

This makes the selection criterion directly data-driven.

---

# 22. Recommended default filtering strategy

Use:

```text
target cells/subclass:        1000
minimum cells/subclass:        500
```

Therefore:

```text
N < 500:
    discard subclass

500 <= N < 1000:
    keep all N cells

N >= 1000:
    randomly sample 1000 cells
```

This produces a relatively balanced dataset while avoiding exclusion of moderately represented populations.

Make both thresholds user-configurable.

---

# 23. Potential extension: stricter balanced mode

Optionally support:

```bash
--require-full-sample
```

In this mode:

```text
min_cells_per_subclass = cells_per_subclass
```

Thus every retained subclass contributes exactly 1000 cells.

This mode may be useful for machine-learning experiments requiring completely balanced classes.

It should not be the default.

---

# 24. Potential extension: donor-aware sampling

If donor information is available, avoid drawing all 1000 cells for a subclass from a single donor.

A future or optional strategy is:

```text
subclass
    ├── donor 1
    ├── donor 2
    └── donor 3
```

followed by approximately balanced sampling across donors.

For the initial implementation, ordinary seeded random sampling is sufficient, but retain donor metadata so this can be implemented later.

---

# 25. Do not perform downstream preprocessing

This script is responsible for **retrieval and reference construction only**.

Do not perform:

```text
TP10K normalization
log transformation
HVG selection
PCA
UMAP
batch correction
reference averaging
pseudobulk construction
marker-gene selection
SNMF
```

These should happen in separate preprocessing/analysis scripts.

The output should preserve the most reusable representation possible: selected single-cell raw counts plus high-quality metadata.

---

# 26. Dependencies

Expected main dependencies:

```text
abc_atlas_access
anndata
scanpy
pandas
numpy
scipy
```

Use the current supported `abc_atlas_access` API rather than hardcoding AWS URLs whenever possible.

The Allen examples currently use:

```python
abc_cache = AbcProjectCache.from_cache_dir(
    cache_dir=download_base
)
```

and metadata access patterns such as:

```python
abc_cache.get_metadata_dataframe(...)
```

and:

```python
abc_cache.list_expression_matrix_files(...)
```

Use these APIs where appropriate.

---

# 27. Robustness to Allen release changes

Do not assume that:

- row order remains stable;
- specific metadata columns always occur at fixed positions;
- a particular `.h5ad` filename remains unchanged;
- the current manifest is permanent.

Use IDs and metadata joins.

When necessary, inspect:

```python
abc_cache.current_manifest
```

and discover available files using the cache API.

Fail with a useful error message if expected AIT33 resources are missing.

---

# 28. Acceptance criteria

The implementation is complete when all of the following work.

### Region discovery

```bash
python retrieve_ait33_reference.py --list-regions
```

returns valid anatomical regions.

### Cerebellum retrieval

```bash
python retrieve_ait33_reference.py \
    --region cerebellum \
    --cells-per-subclass 1000 \
    --min-cells-per-subclass 500 \
    --output cerebellum_ait33_reference.h5ad
```

produces a valid AnnData file.

### Dataset properties

The resulting AnnData:

- contains only AIBS single-cell profiles;
- uses AIT33 consensus taxonomy labels;
- contains cells associated with the selected region;
- uses `subclass` as the primary cell-type label;
- removes subclasses with fewer than the configured number of regional cells;
- contains at most 1000 cells per retained subclass;
- preserves raw counts;
- remains sparse;
- retains original cell IDs;
- includes hierarchical taxonomy metadata;
- includes anatomical metadata;
- records retrieval provenance.

### Reproducibility

Running the script twice with the same:

```text
region
thresholds
random seed
Allen release
```

must select the same cells.

---

# 29. First development target

Use:

```text
Cerebellum
```

as the initial test case.

Development order:

```text
1. Connect to current Allen manifest.
2. Load AIBS AIT33 cell metadata.
3. Join AIT33 taxonomy.
4. Inspect anatomical columns.
5. Implement --list-regions.
6. Select cerebellar cells.
7. Count subclasses.
8. Apply abundance threshold.
9. Sample cells.
10. Determine required expression matrices.
11. retrieve raw expression.
12. construct AnnData.
13. validate.
14. write .h5ad and summary CSV.
```

Do not optimize prematurely.

First obtain a correct cerebellum reference, inspect it manually, and only then generalize the code to arbitrary regions.
