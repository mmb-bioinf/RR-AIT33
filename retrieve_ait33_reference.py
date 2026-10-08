#!/usr/bin/env python3
"""Retrieve a region-specific, balanced AIT33 single-cell reference dataset.

Uses Consensus-WMB-AIBS-10X raw expression (via reused WMB-10X matrices) and
AIT33 consensus taxonomy labels from Consensus-WMB-integrated-taxonomy.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "abc_atlas_access is required. Install with:\n"
        "  pip install git+https://github.com/alleninstitute/abc_atlas_access.git"
    ) from exc


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AIBS_DIRECTORY = "Consensus-WMB-AIBS-10X"
TAXONOMY_DIRECTORY = "Consensus-WMB-integrated-taxonomy"
GENE_DIRECTORY = "WMB-10X"
TAXONOMY_NAME = "20251030_BICCN-Adult_Mouse_WholeBrain_consensus"
TAXONOMY_ID = "AIT33"

DEFAULT_CELLS_PER_SUBCLASS = 1000
DEFAULT_MIN_CELLS_PER_SUBCLASS = 500
DEFAULT_SEED = 42

# Friendly aliases mapped onto official anatomical_division_label values.
# Only used when the target label exists in the loaded metadata.
REGION_ALIASES: dict[str, str] = {
    "cerebellum": "CB",
    "cb": "CB",
    "hippocampus": "HPF",
    "hippocampal formation": "HPF",
    "hpf": "HPF",
    "cortex": "Isocortex",
    "isocortex": "Isocortex",
    "thalamus": "TH",
    "th": "TH",
    "hypothalamus": "HY",
    "hy": "HY",
    "midbrain": "MB",
    "mb": "MB",
    "striatum": "STR",
    "str": "STR",
    "pallidum": "PAL",
    "pal": "PAL",
    "olfactory areas": "OLF",
    "olfactory": "OLF",
    "olf": "OLF",
    "cortical subplate": "CTXsp",
    "ctxsp": "CTXsp",
    "pons": "MY-Pons-BS",
    "medulla": "MY-Pons-BS",
    "hindbrain": "MY-Pons-BS",
    "my-pons-bs": "MY-Pons-BS",
}

OBS_COLUMN_MAP = {
    "anatomical_division_label": "anatomical_region",
    "region_of_interest_acronym": "fine_anatomical_region",
    "donor_label": "donor",
    "donor_sex": "sex",
    "library_label": "library",
    "dataset_label": "dataset",
    "feature_matrix_label": "feature_matrix_label",
    "cluster_alias": "cluster_alias",
    "neighborhood": "neighborhood",
    "class": "class",
    "subclass": "subclass",
    "supertype": "supertype",
    "cluster": "cluster",
    "neurotransmitter": "neurotransmitter",
}


# ---------------------------------------------------------------------------
# Cache / metadata loading
# ---------------------------------------------------------------------------


def initialize_cache(cache_dir: str | Path) -> AbcProjectCache:
    """Create an AbcProjectCache rooted at *cache_dir*."""
    cache_path = Path(cache_dir).expanduser().resolve()
    cache_path.mkdir(parents=True, exist_ok=True)
    cache = AbcProjectCache.from_cache_dir(cache_dir=cache_path)
    print(f"ABC Atlas manifest: {cache.current_manifest}")
    return cache


def _require_directory(cache: AbcProjectCache, directory: str) -> None:
    # list_directories is a property (list) in current abc_atlas_access.
    directories = cache.list_directories
    if callable(directories):
        directories = directories()
    available = set(directories)
    if directory not in available:
        raise RuntimeError(
            f"Required Allen directory '{directory}' is missing from the "
            f"current manifest ({cache.current_manifest}). "
            f"Available directories include: {sorted(available)[:20]}..."
        )


def load_aibs_cell_metadata(cache: AbcProjectCache) -> pd.DataFrame:
    """Load AIBS consensus cell metadata joined with library and donor tables."""
    _require_directory(cache, AIBS_DIRECTORY)

    cells = cache.get_metadata_dataframe(
        directory=AIBS_DIRECTORY,
        file_name="cell_metadata",
        dtype={"cell_label": str},
        low_memory=False,
    ).set_index("cell_label")

    library = cache.get_metadata_dataframe(
        directory=AIBS_DIRECTORY,
        file_name="library",
    )
    if "library_label" in library.columns:
        library = library.set_index("library_label")
    cells = cells.join(library, on="library_label", how="left", rsuffix="_library")

    donor = cache.get_metadata_dataframe(
        directory=AIBS_DIRECTORY,
        file_name="donor",
    )
    if "donor_label" in donor.columns:
        donor = donor.set_index("donor_label")
    cells = cells.join(donor, on="donor_label", how="left", rsuffix="_donor")

    required = [
        "anatomical_division_label",
        "feature_matrix_label",
        "dataset_label",
    ]
    missing = [c for c in required if c not in cells.columns]
    if missing:
        raise RuntimeError(
            f"AIBS cell metadata is missing required columns: {missing}"
        )

    print(f"Loaded AIBS cell metadata: {len(cells):,} cells")
    return cells


def load_ait33_taxonomy(cache: AbcProjectCache) -> dict[str, pd.DataFrame]:
    """Load official AIT33 taxonomy membership tables."""
    _require_directory(cache, TAXONOMY_DIRECTORY)

    term_set = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIRECTORY,
        file_name="cluster_annotation_term_set",
        skip_hash_check=True,
    ).set_index("label")

    membership = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIRECTORY,
        file_name="cluster_to_cluster_annotation_membership",
    )

    cell_to_cluster = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIRECTORY,
        file_name="cell_to_cluster_membership",
        dtype={"cell_label": str},
    ).set_index("cell_label")

    return {
        "term_set": term_set,
        "membership": membership,
        "cell_to_cluster": cell_to_cluster,
    }


def build_cell_taxonomy_table(
    cells: pd.DataFrame,
    taxonomy: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Join AIBS cells onto AIT33 taxonomy labels via cluster_alias."""
    membership = taxonomy["membership"]
    term_set = taxonomy["term_set"]
    cell_to_cluster = taxonomy["cell_to_cluster"]

    if "cluster_alias" not in membership.columns:
        raise RuntimeError(
            "cluster_to_cluster_annotation_membership lacks 'cluster_alias'"
        )
    if "cluster_annotation_term_set_name" not in membership.columns:
        raise RuntimeError(
            "cluster_to_cluster_annotation_membership lacks "
            "'cluster_annotation_term_set_name'"
        )
    if "cluster_annotation_term_name" not in membership.columns:
        raise RuntimeError(
            "cluster_to_cluster_annotation_membership lacks "
            "'cluster_annotation_term_name'"
        )

    cluster_details = (
        membership.groupby(
            ["cluster_alias", "cluster_annotation_term_set_name"]
        )["cluster_annotation_term_name"]
        .first()
        .unstack()
    )

    # Order columns by official term-set order when available.
    if "name" in term_set.columns:
        ordered = [n for n in term_set["name"].tolist() if n in cluster_details.columns]
        remaining = [c for c in cluster_details.columns if c not in ordered]
        cluster_details = cluster_details[ordered + remaining]

    expected_levels = ["subclass", "class", "supertype", "cluster"]
    missing_levels = [lvl for lvl in expected_levels if lvl not in cluster_details.columns]
    if missing_levels:
        raise RuntimeError(
            f"AIT33 taxonomy pivot is missing expected levels: {missing_levels}. "
            f"Available: {list(cluster_details.columns)}"
        )

    annotated = cells.join(cell_to_cluster, how="inner")
    if "cluster_alias" not in annotated.columns:
        raise RuntimeError("cell_to_cluster_membership did not provide cluster_alias")

    annotated = annotated.join(cluster_details, on="cluster_alias", how="left")
    n_missing = int(annotated["subclass"].isna().sum())
    if n_missing:
        warnings.warn(
            f"Dropping {n_missing:,} cells without AIT33 subclass labels",
            stacklevel=2,
        )
        annotated = annotated.loc[~annotated["subclass"].isna()].copy()

    print(
        f"Joined AIT33 taxonomy: {len(annotated):,} cells, "
        f"{annotated['subclass'].nunique()} subclasses"
    )
    return annotated


def load_region_value_sets(cache: AbcProjectCache) -> pd.DataFrame:
    """Load value_sets rows describing anatomical_division_label terms."""
    try:
        value_sets = cache.get_metadata_dataframe(
            directory=AIBS_DIRECTORY,
            file_name="value_sets",
        )
    except Exception:
        return pd.DataFrame()

    if "field" not in value_sets.columns:
        return pd.DataFrame()
    region_sets = value_sets.loc[
        value_sets["field"] == "anatomical_division_label"
    ].copy()
    return region_sets


# ---------------------------------------------------------------------------
# Region discovery / resolution
# ---------------------------------------------------------------------------


def get_available_regions(
    metadata: pd.DataFrame,
    value_sets: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Return official anatomical regions present in *metadata*.

    Columns: label, description, n_cells
    """
    if "anatomical_division_label" not in metadata.columns:
        raise RuntimeError("Metadata lacks anatomical_division_label")

    counts = (
        metadata.groupby("anatomical_division_label", dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={"anatomical_division_label": "label"})
    )

    desc_map: dict[str, str] = {}
    if value_sets is not None and not value_sets.empty:
        label_col = "label" if "label" in value_sets.columns else None
        desc_col = "description" if "description" in value_sets.columns else None
        if label_col and desc_col:
            for _, row in value_sets.iterrows():
                desc_map[str(row[label_col])] = str(row[desc_col])

    counts["description"] = counts["label"].map(
        lambda x: desc_map.get(str(x), str(x))
    )
    counts = counts.sort_values("label").reset_index(drop=True)
    return counts


def resolve_region(
    region_query: str,
    available_regions: pd.DataFrame,
) -> str:
    """Resolve a user region query to one official anatomical_division_label.

    Raises SystemExit on zero or ambiguous matches.
    """
    query = region_query.strip()
    if not query:
        raise SystemExit("Empty --region value")

    labels = available_regions["label"].astype(str).tolist()
    descriptions = available_regions["description"].astype(str).tolist()
    label_by_lower = {lab.lower(): lab for lab in labels}
    desc_by_lower = {
        desc.lower(): lab for lab, desc in zip(labels, descriptions)
    }

    matches: list[str] = []

    # Alias first (maps onto official labels when present).
    alias_target = REGION_ALIASES.get(query.lower())
    if alias_target is not None and alias_target in labels:
        matches.append(alias_target)

    if query.lower() in label_by_lower:
        matches.append(label_by_lower[query.lower()])
    if query.lower() in desc_by_lower:
        matches.append(desc_by_lower[query.lower()])

    # Deduplicate while preserving order.
    unique_matches: list[str] = []
    for m in matches:
        if m not in unique_matches:
            unique_matches.append(m)

    if len(unique_matches) == 1:
        return unique_matches[0]

    if len(unique_matches) > 1:
        lines = [
            f"Region '{region_query}' matched multiple Allen regions:",
            "",
        ]
        for i, lab in enumerate(unique_matches, start=1):
            desc = available_regions.loc[
                available_regions["label"] == lab, "description"
            ].iloc[0]
            lines.append(f"{i}. {lab} ({desc})")
        lines.append("")
        lines.append("Please select one explicitly.")
        raise SystemExit("\n".join(lines))

    # Partial containment search for helpful error messages only.
    candidates = []
    q = query.lower()
    for lab, desc in zip(labels, descriptions):
        if q in lab.lower() or q in desc.lower():
            candidates.append((lab, desc))

    if candidates:
        lines = [
            f"Region '{region_query}' matched multiple Allen regions:"
            if len(candidates) > 1
            else f"Region '{region_query}' did not exactly match; closest:",
            "",
        ]
        for i, (lab, desc) in enumerate(candidates, start=1):
            lines.append(f"{i}. {lab} ({desc})")
        lines.append("")
        lines.append("Please select one explicitly (use --list-regions).")
        raise SystemExit("\n".join(lines))

    raise SystemExit(
        f"Region '{region_query}' not found. Use --list-regions to see options."
    )


def select_region_cells(metadata: pd.DataFrame, region: str) -> pd.DataFrame:
    """Return cells whose biological anatomical_division_label equals *region*."""
    mask = metadata["anatomical_division_label"].astype(str) == str(region)
    selected = metadata.loc[mask].copy()
    if selected.empty:
        raise RuntimeError(f"No cells found for anatomical region '{region}'")
    print(f"Selected region: {region} ({len(selected):,} cells)")
    return selected


# ---------------------------------------------------------------------------
# Subclass filtering / sampling
# ---------------------------------------------------------------------------


def filter_subclasses(
    metadata: pd.DataFrame,
    min_cells_per_subclass: int,
) -> tuple[pd.Index, pd.Series]:
    """Return retained subclass labels and full regional subclass counts."""
    region_counts = metadata.groupby("subclass", observed=True).size().sort_values(
        ascending=False
    )
    valid = region_counts[region_counts >= min_cells_per_subclass].index
    excluded = region_counts[region_counts < min_cells_per_subclass]

    print()
    print(f"Selected region subclasses")
    print(f"  AIT33 subclasses observed:       {len(region_counts)}")
    print(
        f"  Subclasses with >={min_cells_per_subclass} cells: "
        f"{len(valid)}"
    )
    print(f"  Subclasses excluded:             {len(excluded)}")
    if len(excluded):
        print()
        print("Excluded subclasses (count):")
        for subclass, count in excluded.items():
            print(f"  {subclass}: {int(count)}")
    print()
    return valid, region_counts


def sample_cells_per_subclass(
    metadata: pd.DataFrame,
    valid_subclasses: Iterable[str],
    region_counts: pd.Series,
    cells_per_subclass: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Balanced sample without replacement; return sampled cells + summary."""
    rng = np.random.default_rng(seed)
    valid_set = set(valid_subclasses)
    sampled_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    for subclass, n_available in region_counts.items():
        retained = subclass in valid_set
        if retained:
            n_select = min(int(cells_per_subclass), int(n_available))
            pool = metadata.loc[metadata["subclass"] == subclass]
            if n_select < len(pool):
                # Deterministic permutation via numpy Generator.
                chosen_idx = rng.choice(
                    pool.index.to_numpy(),
                    size=n_select,
                    replace=False,
                )
                sampled_parts.append(pool.loc[chosen_idx])
            else:
                sampled_parts.append(pool)
            n_sampled = n_select
        else:
            n_sampled = 0

        summary_rows.append(
            {
                "subclass": subclass,
                "cells_in_region": int(n_available),
                "retained": bool(retained),
                "cells_sampled": int(n_sampled),
            }
        )

    if not sampled_parts:
        raise RuntimeError(
            "No subclasses passed the abundance filter; "
            "lower --min-cells-per-subclass or choose another region."
        )

    sampled = pd.concat(sampled_parts, axis=0)
    # Stable order for reproducibility of downstream concat.
    sampled = sampled.sort_index()
    summary = pd.DataFrame(summary_rows).sort_values(
        ["retained", "cells_in_region"], ascending=[False, False]
    )
    print(
        f"Sampled {len(sampled):,} cells across "
        f"{sampled['subclass'].nunique()} subclasses "
        f"(seed={seed})"
    )
    return sampled, summary


# ---------------------------------------------------------------------------
# Expression loading
# ---------------------------------------------------------------------------


def determine_required_expression_files(
    selected_cells: pd.DataFrame,
) -> pd.DataFrame:
    """Return unique (dataset_label, feature_matrix_label) pairs."""
    cols = ["dataset_label", "feature_matrix_label"]
    missing = [c for c in cols if c not in selected_cells.columns]
    if missing:
        raise RuntimeError(f"Selected cells missing columns: {missing}")

    files = (
        selected_cells.groupby(cols, observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(cols)
    )
    print("Required expression matrices:")
    for _, row in files.iterrows():
        print(
            f"  {row['dataset_label']} / {row['feature_matrix_label']}/raw "
            f"({int(row['n_cells']):,} cells)"
        )
    return files


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def load_gene_metadata(cache: AbcProjectCache) -> pd.DataFrame:
    """Load WMB-10X gene metadata indexed by gene_identifier."""
    genes = cache.get_metadata_dataframe(
        directory=GENE_DIRECTORY,
        file_name="gene",
    )
    if "gene_identifier" not in genes.columns:
        raise RuntimeError("WMB-10X gene table lacks gene_identifier")
    genes = genes.set_index("gene_identifier")
    return genes


def load_selected_expression(
    cache: AbcProjectCache,
    selected_cells: pd.DataFrame,
    required_files: pd.DataFrame,
) -> tuple[ad.AnnData, pd.Index]:
    """Load raw expression for selected cells only; keep X sparse.

    Returns AnnData of concatenated subsets and the index of cells that were
    successfully retrieved (may be a subset if some cells are missing from
    reused WMB expression files).
    """
    adatas: list[ad.AnnData] = []
    found_ids: list[str] = []
    missing_total = 0

    for _, row in required_files.iterrows():
        directory = str(row["dataset_label"])
        matrix = str(row["feature_matrix_label"])
        file_name = f"{matrix}/raw"

        cell_ids = selected_cells.loc[
            (selected_cells["dataset_label"] == directory)
            & (selected_cells["feature_matrix_label"] == matrix)
        ].index

        print(f"Loading {directory}:{file_name} ...")
        try:
            path = cache.get_file_path(directory=directory, file_name=file_name)
        except Exception as exc:
            raise RuntimeError(
                f"Failed to locate expression file {directory}/{file_name}: {exc}"
            ) from exc

        backed = ad.read_h5ad(path, backed="r")
        try:
            common = backed.obs_names.intersection(cell_ids)
            n_missing = len(cell_ids) - len(common)
            if n_missing:
                missing_total += n_missing
                warnings.warn(
                    f"{n_missing} selected cells not found in {directory}/{file_name}",
                    stacklevel=2,
                )
            if len(common) == 0:
                continue

            # Load only the selected rows into memory; keep sparse X.
            subset = backed[common].to_memory()
            if sparse.issparse(subset.X):
                pass
            elif subset.X is not None:
                # Ensure sparse representation without densifying unnecessarily.
                subset.X = sparse.csr_matrix(subset.X)
            adatas.append(subset)
            found_ids.extend(list(common))
        finally:
            if hasattr(backed, "file") and backed.file is not None:
                backed.file.close()
            del backed

    if not adatas:
        raise RuntimeError(
            "No expression data could be loaded for the selected cells. "
            "Consensus metadata cells may be absent from reused WMB matrices."
        )

    if missing_total:
        print(
            f"Warning: {missing_total:,} selected cells were absent from "
            "expression matrices and will be dropped."
        )

    combined = ad.concat(adatas, join="outer", merge="same")
    if not sparse.issparse(combined.X):
        combined.X = sparse.csr_matrix(combined.X)
    else:
        combined.X = combined.X.tocsr()

    # Preserve original cell ID order from sampling when possible.
    ordered = [cid for cid in selected_cells.index if cid in combined.obs_names]
    combined = combined[ordered].copy()
    return combined, pd.Index(ordered)


# ---------------------------------------------------------------------------
# AnnData construction / validation / I/O
# ---------------------------------------------------------------------------


def _map_obs_columns(metadata: pd.DataFrame) -> pd.DataFrame:
    """Map available official columns onto clear output names."""
    out = pd.DataFrame(index=metadata.index)
    out["cell_id"] = metadata.index.astype(str)

    for src, dst in OBS_COLUMN_MAP.items():
        if src in metadata.columns:
            out[dst] = metadata[src].values

    # Coarse alias when anatomical_region present.
    if "anatomical_region" in out.columns:
        out["coarse_anatomical_region"] = out["anatomical_region"]

    return out


def build_anndata(
    expression: ad.AnnData,
    metadata: pd.DataFrame,
    gene_metadata: pd.DataFrame,
    reference_info: dict[str, Any],
    sampling_summary: pd.DataFrame,
) -> ad.AnnData:
    """Assemble final AnnData with mapped obs, gene var, and provenance."""
    # Align metadata to expression obs order.
    meta = metadata.loc[expression.obs_names]
    obs = _map_obs_columns(meta)

    adata = ad.AnnData(
        X=expression.X,
        obs=obs,
        var=expression.var.copy() if expression.var is not None else pd.DataFrame(
            index=expression.var_names
        ),
    )
    adata.obs_names = expression.obs_names.astype(str)
    adata.obs_names.name = None

    # Prefer stable gene IDs as var_names.
    if adata.var_names.isna().any() or not adata.var_names.is_unique:
        # Fall back to expression var as-is; attempt gene_identifier later.
        pass

    # Attach gene metadata where identifiers match.
    gene_cols = [
        c
        for c in ["gene_symbol", "gene_identifier", "biotype", "name"]
        if c in gene_metadata.columns or c == gene_metadata.index.name
    ]
    # gene_identifier is the index.
    aligned_genes = gene_metadata.reindex(adata.var_names)
    if "gene_symbol" in aligned_genes.columns:
        adata.var["gene_symbol"] = aligned_genes["gene_symbol"].values
    adata.var["gene_identifier"] = adata.var_names.astype(str)

    # If expression used gene symbols as var_names and IDs are available, swap.
    if (
        "gene_symbol" in expression.var.columns
        and expression.var["gene_symbol"].nunique() == expression.n_vars
    ):
        # Keep current; already unique symbols.
        pass

    if not adata.var_names.is_unique:
        # Force unique IDs.
        adata.var_names = pd.Index(
            [
                f"{name}_{i}" if adata.var_names.duplicated()[i] else name
                for i, name in enumerate(adata.var_names.astype(str))
            ]
        )

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    else:
        adata.X = adata.X.tocsr()

    adata.uns["reference_info"] = reference_info
    adata.uns["subclass_sampling"] = sampling_summary.to_dict(orient="list")
    return adata


def _is_count_like(X: sparse.spmatrix, n_samples: int = 1000, tol: float = 1e-6) -> bool:
    """Check that sampled non-zero entries are nearly integers and non-negative."""
    if X.min() < -tol:
        return False
    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return True
        if data.size > n_samples:
            rng = np.random.default_rng(0)
            sample = rng.choice(data, size=n_samples, replace=False)
        else:
            sample = data
    else:
        flat = np.asarray(X).ravel()
        nonzero = flat[flat != 0]
        if nonzero.size == 0:
            return True
        if nonzero.size > n_samples:
            rng = np.random.default_rng(0)
            sample = rng.choice(nonzero, size=n_samples, replace=False)
        else:
            sample = nonzero
    return bool(np.all(np.abs(sample - np.round(sample)) <= tol))


def validate_reference(
    adata: ad.AnnData,
    config: dict[str, Any],
    region_label: str,
) -> None:
    """Run acceptance checks before writing outputs."""
    assert adata.obs_names.is_unique, "obs_names must be unique"
    assert adata.var_names.is_unique, "var_names must be unique"
    assert "subclass" in adata.obs.columns, "obs must contain subclass"
    assert not adata.obs["subclass"].isna().any(), "subclass has missing values"

    if "anatomical_region" in adata.obs.columns:
        bad = adata.obs["anatomical_region"].astype(str) != str(region_label)
        assert not bad.any(), (
            f"{int(bad.sum())} cells do not match selected region '{region_label}'"
        )

    counts = adata.obs["subclass"].value_counts()
    assert (counts <= config["cells_per_subclass"]).all(), (
        "Some subclasses exceed cells_per_subclass"
    )
    # Retained subclasses in the written object must meet the abundance floor
    # for cells that were available in-region before expression dropouts.
    # After expression loading, counts may dip below min if matrices lack cells;
    # enforce that no subclass is empty.
    assert (counts > 0).all(), "Empty subclass present in output"

    assert sparse.issparse(adata.X), "Expression matrix must remain sparse"
    assert adata.X.min() >= 0, "Expression contains negative values"
    assert _is_count_like(adata.X), (
        "Expression values do not look like raw counts "
        "(expected non-negative near-integers)"
    )
    print("Validation checks passed.")


def write_outputs(
    adata: Optional[ad.AnnData],
    summary: pd.DataFrame,
    output: str | Path,
    overwrite: bool = False,
    metadata_only: bool = False,
    metadata: Optional[pd.DataFrame] = None,
) -> Path:
    """Write .h5ad (unless metadata-only) and sampling summary CSV."""
    output_path = Path(output)
    if output_path.suffix.lower() != ".h5ad":
        output_path = output_path.with_suffix(".h5ad")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    summary_path = output_path.with_name(
        f"{output_path.stem}_sampling_summary.csv"
    )

    if output_path.exists() and not overwrite and not metadata_only:
        raise SystemExit(
            f"Output exists: {output_path}. Pass --overwrite to replace it."
        )

    summary.to_csv(summary_path, index=False)
    print(f"Wrote sampling summary: {summary_path}")

    if metadata_only:
        meta_path = output_path.with_name(f"{output_path.stem}_metadata.csv")
        if metadata is not None:
            mapped = _map_obs_columns(metadata)
            mapped.to_csv(meta_path, index=True)
            print(f"Wrote metadata table: {meta_path}")
        return summary_path

    assert adata is not None
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Wrote AnnData: {output_path}")
    return output_path


def print_completion_summary(
    adata: ad.AnnData,
    region: str,
    output_path: Path,
) -> None:
    """Print the Task §16 completion block."""
    counts = adata.obs["subclass"].value_counts()
    print()
    print("AIT33 reference successfully generated")
    print()
    print("Region:")
    print(f"    {region}")
    print()
    print("Dataset:")
    print(f"    {AIBS_DIRECTORY}")
    print()
    print("Taxonomy:")
    print(f"    {TAXONOMY_ID}")
    print()
    print("Taxonomy level:")
    print("    subclass")
    print()
    print("Subclasses retained:")
    print(f"    {adata.obs['subclass'].nunique()}")
    print()
    print("Cells:")
    print(f"    {adata.n_obs:,}")
    print()
    print("Genes:")
    print(f"    {adata.n_vars:,}")
    print()
    print("Median cells/subclass:")
    print(f"    {int(counts.median()):,}")
    print()
    print("Minimum cells/subclass:")
    print(f"    {int(counts.min()):,}")
    print()
    print("Maximum cells/subclass:")
    print(f"    {int(counts.max()):,}")
    print()
    print("Expression:")
    print("    raw counts")
    print()
    print("Output:")
    print(f"    {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve a region-specific balanced AIT33 single-cell "
            "reference AnnData from Consensus-WMB-AIBS-10X raw counts."
        )
    )
    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help="Anatomical region (case-insensitive; see --list-regions).",
    )
    parser.add_argument(
        "--cells-per-subclass",
        type=int,
        default=DEFAULT_CELLS_PER_SUBCLASS,
        help=f"Target cells sampled per retained subclass "
        f"(default: {DEFAULT_CELLS_PER_SUBCLASS}).",
    )
    parser.add_argument(
        "--min-cells-per-subclass",
        type=int,
        default=DEFAULT_MIN_CELLS_PER_SUBCLASS,
        help=f"Minimum regional cells required to retain a subclass "
        f"(default: {DEFAULT_MIN_CELLS_PER_SUBCLASS}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Random seed for sampling (default: {DEFAULT_SEED}).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output .h5ad path.",
    )
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="./abc_atlas_cache",
        help="Directory used by AbcProjectCache (default: ./abc_atlas_cache).",
    )
    parser.add_argument(
        "--list-regions",
        action="store_true",
        help="Print valid selectable anatomical regions and exit.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output .h5ad.",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Select/sample cells and write summary without downloading expression.",
    )
    parser.add_argument(
        "--require-full-sample",
        action="store_true",
        help=(
            "Require every retained subclass to have at least "
            "--cells-per-subclass regional cells (fully balanced)."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)

    if args.require_full_sample:
        args.min_cells_per_subclass = args.cells_per_subclass

    if args.cells_per_subclass < 1:
        raise SystemExit("--cells-per-subclass must be >= 1")
    if args.min_cells_per_subclass < 1:
        raise SystemExit("--min-cells-per-subclass must be >= 1")

    if not args.list_regions and not args.region:
        raise SystemExit("Provide --region or --list-regions")
    if not args.list_regions and not args.output:
        raise SystemExit("Provide --output (or use --list-regions)")

    cache = initialize_cache(args.cache_dir)
    cells = load_aibs_cell_metadata(cache)
    taxonomy = load_ait33_taxonomy(cache)
    metadata = build_cell_taxonomy_table(cells, taxonomy)
    value_sets = load_region_value_sets(cache)
    available_regions = get_available_regions(metadata, value_sets)

    if args.list_regions:
        print("Available anatomical regions (anatomical_division_label):")
        print()
        for _, row in available_regions.iterrows():
            desc = row["description"]
            label = row["label"]
            if str(desc) != str(label):
                print(f"  {label:16s}  {desc}  ({int(row['n_cells']):,} cells)")
            else:
                print(f"  {label:16s}  ({int(row['n_cells']):,} cells)")
        return 0

    region = resolve_region(args.region, available_regions)
    region_cells = select_region_cells(metadata, region)
    valid_subclasses, region_counts = filter_subclasses(
        region_cells, args.min_cells_per_subclass
    )
    sampled, summary = sample_cells_per_subclass(
        region_cells,
        valid_subclasses,
        region_counts,
        args.cells_per_subclass,
        args.seed,
    )

    if args.metadata_only:
        write_outputs(
            adata=None,
            summary=summary,
            output=args.output,
            overwrite=args.overwrite,
            metadata_only=True,
            metadata=sampled,
        )
        print()
        print("Metadata-only run complete (expression not downloaded).")
        print(f"Region: {region}")
        print(f"Subclasses retained: {len(valid_subclasses)}")
        print(f"Cells sampled: {len(sampled):,}")
        return 0

    required_files = determine_required_expression_files(sampled)
    gene_metadata = load_gene_metadata(cache)
    expression, found_ids = load_selected_expression(
        cache, sampled, required_files
    )
    sampled_found = sampled.loc[found_ids]

    # Drop subclasses that lost all cells during expression retrieval.
    post_counts = sampled_found.groupby("subclass").size()
    empty = [s for s in valid_subclasses if s not in post_counts.index]
    if empty:
        warnings.warn(
            f"{len(empty)} retained subclasses lost all cells during "
            f"expression loading and will be omitted: {empty[:10]}",
            stacklevel=2,
        )
        sampled_found = sampled_found.loc[
            ~sampled_found["subclass"].isin(empty)
        ]
        expression = expression[sampled_found.index].copy()
        summary.loc[summary["subclass"].isin(empty), "cells_sampled"] = 0
        summary.loc[summary["subclass"].isin(empty), "retained"] = False

    # Update sampled counts after expression dropouts.
    actual_counts = sampled_found.groupby("subclass").size()
    for subclass, n in actual_counts.items():
        summary.loc[summary["subclass"] == subclass, "cells_sampled"] = int(n)

    reference_info = {
        "taxonomy": TAXONOMY_ID,
        "taxonomy_name": TAXONOMY_NAME,
        "expression_source": AIBS_DIRECTORY,
        "selected_region": region,
        "taxonomy_level": "subclass",
        "cells_per_subclass_target": args.cells_per_subclass,
        "min_cells_per_subclass": args.min_cells_per_subclass,
        "random_seed": args.seed,
        "expression_scale": "raw_counts",
        "retrieval_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "abc_atlas_access_version": _package_version("abc_atlas_access"),
        "allen_manifest": str(cache.current_manifest),
        "require_full_sample": bool(args.require_full_sample),
    }

    adata = build_anndata(
        expression=expression,
        metadata=sampled_found,
        gene_metadata=gene_metadata,
        reference_info=reference_info,
        sampling_summary=summary,
    )

    validate_reference(
        adata,
        config={
            "cells_per_subclass": args.cells_per_subclass,
            "min_cells_per_subclass": args.min_cells_per_subclass,
        },
        region_label=region,
    )

    output_path = write_outputs(
        adata=adata,
        summary=summary,
        output=args.output,
        overwrite=args.overwrite,
        metadata_only=False,
    )
    print_completion_summary(adata, region, Path(output_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
