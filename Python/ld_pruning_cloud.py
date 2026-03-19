#!/usr/bin/env python3
"""
LD Pruning Script for Gene-Level Sparse LD Matrices

This script performs LD pruning on variants within genes using sparse LD matrices (.npz files).
It supports both local filesystem and Google Cloud Storage (GCS) for LD matrix access.

The pruning algorithm uses a greedy approach: starting from variants with highest priority
(e.g., lowest p-value or highest effect size), it removes variants in LD above the threshold.
"""

import subprocess
import os
import sys
from pathlib import Path
import numpy as np
import polars as pl
from scipy import sparse
from typing import Optional, Dict, List, Tuple, Set
import glob
from tqdm import tqdm
import argparse
from io import BytesIO, StringIO
import fsspec
import pickle

# Optional GCS support
try:
    from google.cloud import storage
except ImportError:
    storage = None

# Optional Hail Batch support
try:
    import hailtop.batch as hb
except ImportError:
    hb = None

script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.append(str(script_dir))

# --- Default Paths ---
TOP_DIR = "/Users/wlu/Dropbox (Partners HealthCare)"
GNOMAD_LOF_METRICS_FILE = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/utility/gnomad.v4.constraint.by_gene.tsv"


def summarize_ld_distribution(ld_values: np.ndarray) -> dict:
    """Summarize the distribution of LD r² values for threshold selection."""
    # Remove diagonal (self-correlations = 1)
    ld_values = ld_values[ld_values < 0.9999]
    
    if len(ld_values) == 0:
        return {}
    
    summary = {
        'n_pairs': len(ld_values),
        'mean_r2': np.mean(ld_values),
        'median_r2': np.median(ld_values),
        'std_r2': np.std(ld_values),
        'pct_above_0.1': 100 * np.mean(ld_values > 0.1),
        'pct_above_0.2': 100 * np.mean(ld_values > 0.2),
        'pct_above_0.5': 100 * np.mean(ld_values > 0.5),
        'pct_above_0.8': 100 * np.mean(ld_values > 0.8),
    }
    return summary


def compute_r2_from_xtx(xtx_by_n: sparse.spmatrix, allele_frequency: np.ndarray) -> sparse.spmatrix:
    """
    Compute the LD r² matrix from X'X/n and allele frequencies.
    
    The correlation between SNPs i and j is:
        r_ij = (XtX_ij/n - 2*p_i*p_j) / sqrt((XtX_ii/n - 2*p_i^2) * (XtX_jj/n - 2*p_j^2))
    
    Args:
        xtx_by_n: Sparse matrix X'X/n where X is the 0-1-2 genotype matrix
        allele_frequency: Array of allele frequencies for each variant
        
    Returns:
        Sparse r² matrix
    """
    # Compute variances: Var = XtX_ii/n - 2*p^2
    diag = np.asarray(xtx_by_n.diagonal()).flatten()
    variances = diag - 2 * allele_frequency**2
    
    # Handle zero or negative variance (monomorphic sites)
    variances = np.maximum(variances, 1e-10)
    std_devs = np.sqrt(variances)
    
    # Convert to COO for efficient element-wise operations
    xtx_coo = xtx_by_n.tocoo()
    rows, cols, data = xtx_coo.row, xtx_coo.col, xtx_coo.data
    
    # r_ij = (XtX_ij/n - 2*p_i*p_j) / (std_i * std_j)
    mean_product = 2 * allele_frequency[rows] * allele_frequency[cols]
    covariances = data - mean_product
    correlations = covariances / (std_devs[rows] * std_devs[cols])
    
    # Square to get r²
    r2_values = correlations ** 2
    
    # Create r² matrix
    r2_matrix = sparse.coo_matrix((r2_values, (rows, cols)), shape=xtx_by_n.shape)
    
    return r2_matrix.tocsr()


def load_gnomad_mapping(gnomad_file: Path) -> Dict[str, str]:
    """Loads the gnomAD file and creates a gene symbol -> gene ID mapping."""
    print(f"Loading gnomAD mapping from: {gnomad_file}")
    gnomad_df = pl.read_csv(
        gnomad_file,
        separator='\t',
        comment_prefix='#',
        null_values="NA",
        columns=['gene', 'gene_id']
    )
    gene_map = dict(zip(gnomad_df['gene'], gnomad_df['gene_id']))
    return gene_map


def ld_prune_gene(
    ld_matrix: sparse.spmatrix,
    snplist_df: pl.DataFrame,
    variant_df: pl.DataFrame,
    r2_threshold: float = 0.2,
    priority_col: Optional[str] = None,
    ascending: bool = True,
    return_removed: bool = False,
) -> pl.DataFrame:
    """
    Perform LD pruning for variants in a single gene.
    
    Uses a greedy algorithm:
    1. Sort variants by priority (e.g., p-value ascending or effect size descending)
    2. Starting from highest priority, keep the variant
    3. Remove all variants in LD (r² > threshold) with the kept variant
    4. Repeat until all variants are processed
    
    Args:
        ld_matrix: Sparse X'X/n matrix for the gene (will be converted to r²)
        snplist_df: DataFrame with SNP info from LD matrix (must have 'SNP' and 'AF' columns)
        variant_df: DataFrame with variant-level data to prune
        r2_threshold: LD r² threshold for pruning (variants with r² > threshold are removed)
        priority_col: Column name to use for prioritization (lower values = higher priority if ascending=True)
        ascending: If True, lower values of priority_col have higher priority
        return_removed: If True, return the removed variants instead of the kept variants
        
    Returns:
        DataFrame of variants that passed LD pruning (or were removed if return_removed=True)
    """
    # Merge variant data with LD matrix SNP list
    merged_df = snplist_df.with_row_index('matrix_index') \
        .join(variant_df, on='SNP', how='inner')
    
    if merged_df.height == 0:
        return variant_df.head(0)  # Return empty DataFrame with same schema
    
    # Sort by priority FIRST (before computing r² matrix indices)
    if priority_col and priority_col in merged_df.columns:
        merged_df = merged_df.sort(priority_col, descending=not ascending)
    
    # Get matrix indices for variants in our data (AFTER sorting)
    matrix_indices = merged_df['matrix_index'].to_numpy()
    
    # Get allele frequencies in the same order as matrix_indices
    af_lookup = dict(zip(
        snplist_df.with_row_index('idx')['idx'].to_list(),
        snplist_df['AF'].to_list()
    ))
    af_ordered = np.array([af_lookup[idx] for idx in matrix_indices])
    
    # Subset LD matrix to variants in our data (in sorted order)
    xtx_subset = ld_matrix[matrix_indices, :][:, matrix_indices]
    
    # Convert X'X/n to r² matrix
    r2_matrix = compute_r2_from_xtx(xtx_subset, af_ordered)
    
    # Greedy LD pruning - directly filter on r² matrix entries
    n_variants = merged_df.height
    kept_mask = np.ones(n_variants, dtype=bool)
    pruned_indices: Set[int] = set()
    
    for i in range(n_variants):
        if i in pruned_indices:
            continue
            
        # Keep this variant, prune variants in LD above threshold
        r2_row = r2_matrix.getrow(i)
        
        # Get indices where r² > threshold (directly from sparse matrix)
        _, cols, vals = sparse.find(r2_row)
        for j, val in zip(cols, vals):
            if j > i and j not in pruned_indices and val > r2_threshold:
                pruned_indices.add(j)
                kept_mask[j] = False
    
    # Filter to kept or removed variants based on return_removed flag
    if return_removed:
        removed_snps = merged_df.filter(~pl.Series(kept_mask))['SNP'].to_list()
        result_df = variant_df.filter(pl.col('SNP').is_in(removed_snps))
    else:
        kept_snps = merged_df.filter(pl.Series(kept_mask))['SNP'].to_list()
        result_df = variant_df.filter(pl.col('SNP').is_in(kept_snps))
    
    return result_df


def prune_single_gene_batch(
    gene_symbol: str,
    gene_id: str,
    variant_data_gcs: str,
    gcs_bucket_name: str,
    gcs_prefix: str,
    output_dir: str,
    r2_threshold: float = 0.5,
    priority_col: Optional[str] = None,
    priority_ascending: bool = True,
    add_chr_prefix: bool = False,
    output_removed: bool = False,
    min_af: float = 0.0,
    max_af: float = 1.0,
):
    """
    LD prune variants for a single gene. Called by Hail Batch jobs.
    Uses pandas instead of polars for compatibility with standard Docker images.
    """
    import pandas as pd
    import numpy as np
    from scipy import sparse
    from google.cloud import storage
    import fsspec
    from io import BytesIO, StringIO

    print(f"Processing gene {gene_symbol} ({gene_id})")

    # Initialize GCS
    client = storage.Client()
    gcs_bucket = client.bucket(gcs_bucket_name)
    fs = fsspec.filesystem("gcs")

    # Load variant data in chunks, keeping only rows for target gene
    # This avoids loading the entire file into memory
    print(f"Loading variant data from: {variant_data_gcs}")
    print(f"Reading in chunks, filtering for gene: {gene_symbol}")

    # Handle gzip/bgz compressed files
    compression = 'gzip' if variant_data_gcs.endswith(('.gz', '.bgz')) else None

    chunks = []
    rows_scanned = 0
    with fs.open(variant_data_gcs, 'rb') as f:
        for chunk in pd.read_csv(f, sep='\t', comment='#', na_values="NA",
                                  compression=compression, chunksize=100000):
            rows_scanned += len(chunk)
            if 'gene' not in chunk.columns:
                raise ValueError("Variant data must have a 'gene' column")
            gene_chunk = chunk[chunk['gene'] == gene_symbol]
            if len(gene_chunk) > 0:
                chunks.append(gene_chunk)

    print(f"Scanned {rows_scanned} total rows")

    if chunks:
        gene_variants = pd.concat(chunks, ignore_index=True)
    else:
        print(f"No variants found for gene {gene_symbol}")
        return

    # Apply AF filter
    if 'AF' in gene_variants.columns:
        gene_variants = gene_variants[
            (gene_variants['AF'] >= min_af) & (gene_variants['AF'] <= max_af)
        ]

    # Deduplicate by SNP
    gene_variants = gene_variants.drop_duplicates(subset='SNP', keep='first')

    if len(gene_variants) == 0:
        print(f"No variants found for gene {gene_symbol}")
        return

    print(f"Found {len(gene_variants)} variants for gene {gene_symbol}")

    # Load LD matrix and SNP list from GCS
    npz_blob = gcs_bucket.blob(f"{gcs_prefix}{gene_id}.npz")
    snp_blob = gcs_bucket.blob(f"{gcs_prefix}{gene_id}.snplist")

    if not (npz_blob.exists() and snp_blob.exists()):
        print(f"No LD matrix found for gene {gene_symbol} ({gene_id}), keeping all variants")
        result_df = gene_variants if not output_removed else gene_variants.head(0)
    else:
        # Load SNP list
        snp_text = snp_blob.download_as_text()
        snplist_df = pd.read_csv(StringIO(snp_text), sep='\t')

        if add_chr_prefix:
            snplist_df['SNP'] = 'chr' + snplist_df['SNP'].astype(str)

        # Load sparse LD matrix
        npz_bytes = npz_blob.download_as_bytes()
        ld_matrix = sparse.load_npz(BytesIO(npz_bytes))

        # Perform LD pruning using pandas version
        result_df = _ld_prune_gene_pandas(
            ld_matrix=ld_matrix,
            snplist_df=snplist_df,
            variant_df=gene_variants,
            r2_threshold=r2_threshold,
            priority_col=priority_col,
            ascending=priority_ascending,
            return_removed=output_removed,
        )

    # Write result to GCS
    output_path = f"{output_dir}/{gene_symbol}.tsv"
    print(f"Writing {len(result_df)} variants to: {output_path}")

    with fs.open(output_path, "w") as f:
        result_df.to_csv(f, sep='\t', index=False)

    print(f"Finished processing gene {gene_symbol}")


def _ld_prune_gene_pandas(
    ld_matrix: sparse.spmatrix,
    snplist_df,  # pandas DataFrame
    variant_df,  # pandas DataFrame
    r2_threshold: float = 0.2,
    priority_col: Optional[str] = None,
    ascending: bool = True,
    return_removed: bool = False,
):
    """
    Pandas version of LD pruning for batch jobs.
    """
    import pandas as pd
    import numpy as np
    from scipy import sparse

    # Merge variant data with LD matrix SNP list
    snplist_df = snplist_df.reset_index(drop=True)
    snplist_df['matrix_index'] = snplist_df.index

    merged_df = snplist_df.merge(variant_df, on='SNP', how='inner')

    if len(merged_df) == 0:
        return variant_df.head(0)

    # Sort by priority
    if priority_col and priority_col in merged_df.columns:
        merged_df = merged_df.sort_values(priority_col, ascending=ascending)

    # Get matrix indices
    matrix_indices = merged_df['matrix_index'].values

    # Get allele frequencies
    af_ordered = snplist_df.set_index(snplist_df.index)['AF'].loc[matrix_indices].values

    # Subset LD matrix
    xtx_subset = ld_matrix[matrix_indices, :][:, matrix_indices]

    # Compute r² matrix
    diag = np.asarray(xtx_subset.diagonal()).flatten()
    variances = diag - 2 * af_ordered**2
    variances = np.maximum(variances, 1e-10)
    std_devs = np.sqrt(variances)

    xtx_coo = xtx_subset.tocoo()
    rows, cols, data = xtx_coo.row, xtx_coo.col, xtx_coo.data

    mean_product = 2 * af_ordered[rows] * af_ordered[cols]
    covariances = data - mean_product
    correlations = covariances / (std_devs[rows] * std_devs[cols])
    r2_values = correlations ** 2

    r2_matrix = sparse.coo_matrix((r2_values, (rows, cols)), shape=xtx_subset.shape).tocsr()

    # Greedy LD pruning
    n_variants = len(merged_df)
    kept_mask = np.ones(n_variants, dtype=bool)
    pruned_indices = set()

    for i in range(n_variants):
        if i in pruned_indices:
            continue
        r2_row = r2_matrix.getrow(i)
        _, cols_idx, vals = sparse.find(r2_row)
        for j, val in zip(cols_idx, vals):
            if j > i and j not in pruned_indices and val > r2_threshold:
                pruned_indices.add(j)
                kept_mask[j] = False

    # Filter results
    if return_removed:
        removed_snps = merged_df[~kept_mask]['SNP'].tolist()
        result_df = variant_df[variant_df['SNP'].isin(removed_snps)]
    else:
        kept_snps = merged_df[kept_mask]['SNP'].tolist()
        result_df = variant_df[variant_df['SNP'].isin(kept_snps)]

    return result_df


TMP_BUCKET = "gs://aou_tmp/burdenEM"


def concatenate_batch_output(
    input_dir: str,
    output_file: str,
    file_pattern: str = "*.tsv",
):
    """
    Concatenate per-gene output files from batch jobs into a single file.

    Args:
        input_dir: GCS directory containing per-gene TSV files (e.g., gs://bucket/path/)
        output_file: Output path for concatenated file (GCS or local)
        file_pattern: Glob pattern for input files (default: *.tsv)
    """
    import pandas as pd

    fs = fsspec.filesystem("gcs")

    # List all matching files
    if not input_dir.endswith("/"):
        input_dir += "/"

    # Remove gs:// prefix for glob
    gcs_path = input_dir.replace("gs://", "")
    all_files = fs.glob(f"{gcs_path}{file_pattern}")

    if not all_files:
        print(f"No files found matching {input_dir}{file_pattern}")
        return

    print(f"Found {len(all_files)} files to concatenate")

    # Read and concatenate all files
    dfs = []
    for i, filepath in enumerate(tqdm(all_files, desc="Reading files")):
        try:
            with fs.open(f"gs://{filepath}", 'r') as f:
                df = pd.read_csv(f, sep='\t')
                if len(df) > 0:
                    dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {filepath}: {e}")

    if not dfs:
        print("No data to concatenate")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined {len(dfs)} files with {len(combined_df)} total rows")

    # Write output
    print(f"Writing to: {output_file}")
    if output_file.startswith("gs://"):
        with fs.open(output_file, 'w') as f:
            combined_df.to_csv(f, sep='\t', index=False)
    else:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(output_file, sep='\t', index=False)

    print(f"Done! Wrote {len(combined_df)} rows to {output_file}")


def main():
    print("Starting LD pruning script...")
    
    parser = argparse.ArgumentParser(
        description="Perform LD pruning on gene-level sparse LD matrices.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments (for non-concatenate modes)
    parser.add_argument(
        "--variant_data", "-v",
        default=None,
        help="Path to variant data file (TSV/CSV with columns: SNP, gene, and optionally a priority column)"
    )
    
    # LD matrix source (local or GCS)
    parser.add_argument(
        "--ld_matrices",
        help="Path to local directory containing LD matrices (.npz and .snplist files)"
    )
    parser.add_argument(
        "--gcs_bucket",
        type=str,
        default=None,
        help="Name of the GCS bucket containing LD .npz and .snplist files"
    )
    parser.add_argument(
        "--gcs_prefix",
        type=str,
        default='',
        help="Prefix path inside the GCS bucket (e.g., 'within_gene_ld_ukbb/')"
    )
    
    # Pruning parameters
    parser.add_argument(
        "--r2_threshold", "-r",
        type=float,
        default=0.5,
        help="LD r² threshold for pruning. Variants with r² > threshold will be pruned."
    )
    parser.add_argument(
        "--priority_col", "-p",
        type=str,
        default=None,
        help="Column name to use for variant prioritization (e.g., 'pvalue' or 'effect_size')"
    )
    parser.add_argument(
        "--priority_ascending",
        action="store_true",
        default=True,
        help="If set, lower values of priority_col have higher priority (use for p-values)"
    )
    parser.add_argument(
        "--priority_descending",
        action="store_false",
        dest="priority_ascending",
        help="If set, higher values of priority_col have higher priority (use for effect sizes)"
    )
    
    # Gene mapping
    parser.add_argument(
        "--gene_map",
        default=GNOMAD_LOF_METRICS_FILE,
        help="Path to gene symbol to Ensembl ID mapping file"
    )
    
    # Filtering
    parser.add_argument(
        "--min_af",
        type=float,
        default=0.0,
        help="Minimum allele frequency threshold"
    )
    parser.add_argument(
        "--max_af",
        type=float,
        default=1.0,
        help="Maximum allele frequency threshold"
    )
    
    # Output
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output file path for pruned variants"
    )
    parser.add_argument(
        "--num_genes", "-n",
        type=int,
        default=None,
        help="Optional: Limit processing to the first N genes"
    )
    parser.add_argument(
        "--add_chr_prefix",
        action="store_true",
        help="Add 'chr' prefix to SNP IDs when matching with LD matrix"
    )
    parser.add_argument(
        "--diagnose",
        action="store_true",
        help="Run in diagnostic mode: summarize LD distribution across genes to help choose threshold. No pruning performed."
    )
    parser.add_argument(
        "--diagnose_genes",
        type=int,
        default=100,
        help="Number of genes to sample for diagnostic mode"
    )
    parser.add_argument(
        "--output_removed",
        action="store_true",
        help="If set, output the removed (pruned) variants instead of the kept variants"
    )

    # Batch mode arguments
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Run in batch mode using Hail Batch (one job per gene)"
    )
    parser.add_argument(
        "--gene_list",
        type=str,
        default=None,
        help="Path to file with list of genes to process (one gene per line). If not provided, will load from gene_pickle_prefix."
    )
    parser.add_argument(
        "--gene_pickle_prefix",
        type=str,
        default="gs://aou_wlu/v8_within_gene_LD/gene_variant_map/gene_map_chrom",
        help="GCS prefix for gene pickle files (e.g., gs://bucket/path/gene_map_chrom). Will load {prefix}{chrom}.pickle for chrom 1-22."
    )
    parser.add_argument(
        "--billing_project",
        type=str,
        default="all-by-aou",
        help="Billing project for Hail Batch"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="8G",
        help="Memory per batch job"
    )
    parser.add_argument(
        "--image",
        type=str,
        default="us-central1-docker.pkg.dev/aou-neale-gwas/hail-batch/within_gene_ld:0.2",
        help="Docker image for batch jobs"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test mode: process only first 10 genes"
    )

    # Concatenate mode arguments
    parser.add_argument(
        "--concatenate-per-gene-outputs",
        action="store_true",
        help="Concatenate per-gene batch output files into a single file"
    )
    parser.add_argument(
        "--concat-input-dir",
        type=str,
        default=None,
        help="GCS directory containing per-gene TSV files to concatenate (e.g., gs://bucket/path/)"
    )
    parser.add_argument(
        "--concat-output",
        type=str,
        default=None,
        help="Output path for concatenated file"
    )
    parser.add_argument(
        "--concat-pattern",
        type=str,
        default="*.tsv",
        help="Glob pattern for input files (default: *.tsv)"
    )

    args = parser.parse_args()

    # Concatenate mode: combine per-gene output files
    if args.concatenate_per_gene_outputs:
        if not args.concat_input_dir:
            parser.error("--concatenate requires --concat_input_dir")
        if not args.concat_output:
            parser.error("--concatenate requires --concat_output")
        concatenate_batch_output(
            input_dir=args.concat_input_dir,
            output_file=args.concat_output,
            file_pattern=args.concat_pattern,
        )
        return

    # Validate arguments for non-concatenate modes
    if not args.variant_data:
        parser.error("--variant_data is required (except in --concatenate mode)")
    if not args.output:
        parser.error("--output is required (except in --concatenate mode)")
    if not args.ld_matrices and not args.gcs_bucket:
        parser.error("Either --ld_matrices or --gcs_bucket must be specified")
    
    # Initialize GCS client if requested
    gcs_client = None
    gcs_bucket = None
    if args.gcs_bucket:
        if storage is None:
            raise ImportError(
                "google-cloud-storage is required to read from GCS. "
                "Please install it with `pip install google-cloud-storage`."
            )
        gcs_client = storage.Client()
        gcs_bucket = gcs_client.bucket(args.gcs_bucket)
    
    # Load gene mapping
    gene_to_id_map = load_gnomad_mapping(Path(args.gene_map))
    # Create reverse mapping (gene_id -> gene_symbol)
    gene_id_to_symbol = {v: k for k, v in gene_to_id_map.items()}

    # Batch mode: submit Hail Batch jobs (handles gene list separately)
    if args.batch:
        if hb is None:
            raise ImportError(
                "hailtop.batch is required for batch mode. "
                "Please install it with `pip install hailtop`."
            )
        if not args.gcs_bucket:
            raise ValueError("Batch mode requires --gcs_bucket to be specified")
        if not args.output.startswith("gs://"):
            raise ValueError("Batch mode requires --output to be a GCS path (gs://...)")

        # Load gene list
        if args.gene_list:
            # Load from text file
            print(f"Loading gene list from: {args.gene_list}")
            with open(args.gene_list, 'r') as f:
                genes_to_process = [line.strip() for line in f if line.strip()]
        else:
            # Load from pickle files on GCS
            print(f"Loading gene list from pickle files: {args.gene_pickle_prefix}*.pickle")
            fs = fsspec.filesystem("gcs")
            genes_to_process = []
            for chrom in range(1, 23):
                pickle_path = f"{args.gene_pickle_prefix}{chrom}.pickle"
                try:
                    with fs.open(pickle_path, 'rb') as f:
                        chrom_genes = pickle.load(f)
                        genes_to_process.extend(chrom_genes)
                        print(f"  Loaded {len(chrom_genes)} genes from chr{chrom}")
                except Exception as e:
                    print(f"  Warning: Could not load {pickle_path}: {e}")

        if args.num_genes is not None and args.num_genes > 0:
            genes_to_process = genes_to_process[:args.num_genes]

        if args.test:
            genes_to_process = genes_to_process[:10]
            print(f"Test mode: limiting to {len(genes_to_process)} genes")

        print(f"Processing {len(genes_to_process)} genes with r² threshold = {args.r2_threshold}")

        backend = hb.ServiceBackend(
            billing_project=args.billing_project,
            remote_tmpdir=TMP_BUCKET,
        )

        batch_name = f"LD_pruning_r2_{args.r2_threshold}"
        if args.output_removed:
            batch_name += "_removed"

        b = hb.Batch(
            name=batch_name,
            backend=backend,
            default_storage="10Gi",
        )

        print(f"\nSubmitting {len(genes_to_process)} batch jobs...")
        jobs_submitted = 0
        genes_skipped = 0

        # Import hail for hadoop_exists
        import hail as hl
        output_dir = args.output if args.output.endswith("/") else args.output + "/"

        for gene in tqdm(genes_to_process, desc="Submitting jobs"):
            # Pickle files contain gene_ids (ENSG...), need to get symbol
            if gene.startswith("ENSG"):
                gene_id = gene
                gene_symbol = gene_id_to_symbol.get(gene_id)
                if not gene_symbol:
                    continue
            else:
                # Assume it's a gene symbol
                gene_symbol = gene
                gene_id = gene_to_id_map.get(gene_symbol)
                if not gene_id:
                    continue

            # Check if output file already exists
            output_file_path = f"{output_dir}{gene_symbol}.tsv"
            if hl.hadoop_exists(output_file_path):
                genes_skipped += 1
                continue

            j = b.new_python_job(name=f"LD_prune_{gene_symbol}")
            j.memory(args.memory)
            j.image(args.image)

            j.call(
                prune_single_gene_batch,
                gene_symbol=gene_symbol,
                gene_id=gene_id,
                variant_data_gcs=args.variant_data,
                gcs_bucket_name=args.gcs_bucket,
                gcs_prefix=args.gcs_prefix,
                output_dir=args.output,
                r2_threshold=args.r2_threshold,
                priority_col=args.priority_col,
                priority_ascending=args.priority_ascending,
                add_chr_prefix=args.add_chr_prefix,
                output_removed=args.output_removed,
                min_af=args.min_af,
                max_af=args.max_af,
            )
            jobs_submitted += 1

        print(f"\nSubmitted {jobs_submitted} jobs, skipped {genes_skipped} (already completed)")
        b.run()

    else:
        # Non-batch mode: load variant data locally
        print(f"Loading variant data from: {args.variant_data}")
        variants_df = pl.read_csv(
            args.variant_data,
            separator='\t',
            has_header=True,
            comment_prefix='#',
            null_values="NA",
        )

        # Apply AF filter if AF column exists
        if 'AF' in variants_df.columns:
            variants_df = variants_df.filter(
                (pl.col('AF') >= args.min_af) & (pl.col('AF') <= args.max_af)
            )

        # Deduplicate by SNP
        variants_df = variants_df.group_by('SNP').first()

        print(f"Loaded {variants_df.height} variants")
        print(f"Columns: {variants_df.columns}")

        # Get list of genes to process
        if 'gene' not in variants_df.columns:
            raise ValueError("Variant data must have a 'gene' column")

        genes_to_process = sorted([g for g in variants_df['gene'].unique() if g is not None])
        if args.num_genes is not None and args.num_genes > 0:
            genes_to_process = genes_to_process[:args.num_genes]
            print(f"\nLimiting processing to the first {len(genes_to_process)} genes.")

        print(f"\nProcessing {len(genes_to_process)} genes with r² threshold = {args.r2_threshold}")

        # Test mode: limit to first 10 genes
        if args.test:
            genes_to_process = genes_to_process[:10]
            print(f"\nTest mode: limiting to {len(genes_to_process)} genes")

        # Diagnostic mode: summarize LD distribution
        if args.diagnose:
            print(f"\n--- DIAGNOSTIC MODE: Sampling LD distribution from {args.diagnose_genes} genes ---")
            print("Note: Matrices store X'X/n, converting to r² using allele frequencies...")
            all_r2_values = []
            sample_genes = genes_to_process[:args.diagnose_genes]
            
            for gene_symbol in tqdm(sample_genes, desc="Sampling LD"):
                gene_id = gene_to_id_map.get(gene_symbol)
                if not gene_id:
                    continue
                
                try:
                    if args.gcs_bucket:
                        npz_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.npz")
                        snp_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.snplist")
                        if not (npz_blob.exists() and snp_blob.exists()):
                            continue
                        npz_bytes = npz_blob.download_as_bytes()
                        ld_matrix = sparse.load_npz(BytesIO(npz_bytes))
                        snp_text = snp_blob.download_as_text()
                        snplist_df = pl.read_csv(StringIO(snp_text), separator='\t', has_header=True)
                    else:
                        ld_matrix_path = Path(args.ld_matrices) / f"{gene_id}.npz"
                        snplist_path = Path(args.ld_matrices) / f"{gene_id}.snplist"
                        if not (ld_matrix_path.exists() and snplist_path.exists()):
                            continue
                        ld_matrix = sparse.load_npz(ld_matrix_path)
                        snplist_df = pl.read_csv(snplist_path, separator='\t', has_header=True)
                    
                    # Get allele frequencies and compute r²
                    af = snplist_df['AF'].to_numpy()
                    r2_matrix = compute_r2_from_xtx(ld_matrix, af)
                    
                    # Get upper triangle values (avoid double-counting)
                    r2_upper = sparse.triu(r2_matrix, k=1)
                    values = r2_upper.data
                    if len(values) > 0:
                        all_r2_values.extend(values)
                except Exception as e:
                    continue
            
            if all_r2_values:
                all_r2_values = np.array(all_r2_values)
                summary = summarize_ld_distribution(all_r2_values)
                
                print(f"\n--- LD r² Distribution Summary ({len(all_r2_values):,} variant pairs) ---")
                print(f"  Mean r²:   {summary['mean_r2']:.4f}")
                print(f"  Median r²: {summary['median_r2']:.4f}")
                print(f"  Std r²:    {summary['std_r2']:.4f}")
                print(f"\n  % pairs with r² > 0.1: {summary['pct_above_0.1']:.1f}%")
                print(f"  % pairs with r² > 0.2: {summary['pct_above_0.2']:.1f}%")
                print(f"  % pairs with r² > 0.5: {summary['pct_above_0.5']:.1f}%")
                print(f"  % pairs with r² > 0.8: {summary['pct_above_0.8']:.1f}%")
                print(f"\nRecommendation: If most pairs have r² near 0 or 1, use threshold 0.5-0.8")
                print("                If r² is more spread out, use threshold 0.1-0.2")
            else:
                print("No LD values found in sampled genes.")
            return
        
        # Process each gene
        pruned_results = []
        genes_without_ld = 0
        genes_without_id = 0
        
        for gene_symbol in tqdm(genes_to_process, desc="LD Pruning genes"):
            gene_id = gene_to_id_map.get(gene_symbol)
            if not gene_id:
                genes_without_id += 1
                continue
            
            # Get variants for this gene
            gene_variants = variants_df.filter(pl.col('gene') == gene_symbol)
            if gene_variants.height == 0:
                continue
            
            # Load LD matrix and SNP list
            try:
                if args.gcs_bucket:
                    # Load from GCS
                    npz_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.npz")
                    snp_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.snplist")
                    
                    if not (npz_blob.exists() and snp_blob.exists()):
                        genes_without_ld += 1
                        # Keep all variants if no LD matrix
                        pruned_results.append(gene_variants)
                        continue
                    
                    # Load SNP list
                    snp_text = snp_blob.download_as_text()
                    snplist_df = pl.read_csv(StringIO(snp_text), separator='\t', has_header=True)
                    
                    if args.add_chr_prefix:
                        snplist_df = snplist_df.with_columns(
                            ('chr' + pl.col('SNP')).alias('SNP')
                        )
                    
                    # Load sparse LD matrix
                    npz_bytes = npz_blob.download_as_bytes()
                    ld_matrix = sparse.load_npz(BytesIO(npz_bytes))
                    
                else:
                    # Load from local filesystem
                    ld_matrix_path = Path(args.ld_matrices) / f"{gene_id}.npz"
                    ld_snplist_path = Path(args.ld_matrices) / f"{gene_id}.snplist"
                    
                    if not (ld_matrix_path.exists() and ld_snplist_path.exists()):
                        genes_without_ld += 1
                        # Keep all variants if no LD matrix
                        pruned_results.append(gene_variants)
                        continue
                    
                    snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True)
                    
                    if args.add_chr_prefix:
                        snplist_df = snplist_df.with_columns(
                            ('chr' + pl.col('SNP')).alias('SNP')
                        )
                    
                    ld_matrix = sparse.load_npz(ld_matrix_path)
                
                # Perform LD pruning
                pruned_gene_variants = ld_prune_gene(
                    ld_matrix=ld_matrix,
                    snplist_df=snplist_df,
                    variant_df=gene_variants,
                    r2_threshold=args.r2_threshold,
                    priority_col=args.priority_col,
                    ascending=args.priority_ascending,
                    return_removed=args.output_removed,
                )
                
                pruned_results.append(pruned_gene_variants)
                
            except Exception as e:
                print(f"\nError processing gene {gene_symbol} ({gene_id}): {e}")
                # Keep all variants on error
                pruned_results.append(gene_variants)
                continue
        
        # Combine results
        if pruned_results:
            final_df = pl.concat(pruned_results, how='diagonal')
        else:
            final_df = variants_df.head(0)  # Empty DataFrame with same schema
        
        # Summary statistics
        original_count = variants_df.height
        output_count = final_df.height
        
        if args.output_removed:
            print(f"\n--- LD Pruning Summary (outputting REMOVED variants) ---")
            print(f"Original variants: {original_count}")
            print(f"Variants removed (in output): {output_count} ({100*output_count/max(original_count,1):.1f}%)")
            print(f"Variants kept (not in output): {original_count - output_count}")
        else:
            print(f"\n--- LD Pruning Summary ---")
            print(f"Original variants: {original_count}")
            print(f"Variants after pruning: {output_count}")
            print(f"Variants removed: {original_count - output_count} ({100*(original_count - output_count)/max(original_count,1):.1f}%)")
        print(f"Genes without Ensembl ID: {genes_without_id}")
        print(f"Genes without LD matrix: {genes_without_ld}")
        
        # Save results - support both local and GCS paths
        output_path = args.output
        print(f"\nWriting {final_df.height} pruned variants to: {output_path}")
        
        if output_path.startswith("gs://"):
            # Write to GCS
            fs = fsspec.filesystem("gcs")
            with fs.open(output_path, "w") as f:
                f.write(final_df.write_csv(separator='\t'))
        else:
            # Write to local filesystem
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            final_df.write_csv(output_path, separator='\t')
        
        print("\nLD pruning complete!")

if __name__ == "__main__":
    main()
