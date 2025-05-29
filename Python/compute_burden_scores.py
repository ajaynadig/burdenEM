import os
import sys
from pathlib import Path
import numpy as np
import polars as pl
from scipy import sparse
from typing import Optional, Dict, List, Tuple
import glob # For finding files
from tqdm import tqdm # For progress bar
import matplotlib.pyplot as plt
import argparse

script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.append(str(script_dir))

# --- Import custom functions/classes --- 
from linkage_disequilibrium import get_burden_score #, XtXOpt, DiagOpt

# --- Define Paths --- 

TOP_DIR = "/Users/lukeoconnor/Dropbox"

# Directory containing the within-gene LD matrices (e.g., .npz files)
LD_MATRIX_DIR = Path(TOP_DIR) / "within_gene_ld_ukbb"

# Directory containing variant-level data
VARIANT_DATA_DIR = Path(TOP_DIR) / "burdenEM_results/data/utility/test_ukbb_exome_vep.txt.bgz"

# Path to the gnomAD LoF metrics file (for gene ID mapping)
GNOMAD_LOF_METRICS_FILE = Path(TOP_DIR) / "burdenEM_results/data/utility/gnomad.v4.constraint.by_gene.tsv"

# Annotations to get from the gnomAD file
GENE_ANNOT_NAMES = ['lof.oe', 'mis_pphen.oe', 'lof.oe_ci.upper', 'cds_length']

# Output directory for results (keep for later)
OUTPUT_DIR = Path("./burden_scores_output")
OUTPUT_DIR.mkdir(exist_ok=True) # Create output directory if it doesn't exist

# >>> Define Output File for R <<<
R_OUTPUT_FILE = Path(TOP_DIR) / "burdenEM_results/data/utility/ukbb_ld_corrected_burden_scores"

# --- Define Constants ---
TRAIT = "50_NA" # Hard-coded trait name based on user request

FUNCTIONAL_CATEGORIES = [
    # 'missense_benign',
    # 'missense_damaging',
    'pLoF', 
    'synonymous'
]
# FUNCTIONAL_CATEGORIES = ['synonymous']
ANNOT_NAMES = FUNCTIONAL_CATEGORIES  # Annotation columns to process
MIN_AF = 0.0
MAX_AF = 0.0001

AF_RTOL = 1

# --- Helper Functions ---

def load_gnomad_mapping(gnomad_file: Path, feature_names: List[str]=[]) -> Dict[str, str]:
    """Loads the gnomAD file and creates a gene symbol -> gene ID mapping."""
    print(f"Loading gnomAD mapping from: {gnomad_file}")
    columns = ['gene', 'gene_id'] + feature_names
    gnomad_df = pl.read_csv(
        gnomad_file,
        separator='\t',
        comment_prefix='#',
        null_values="NA",
        columns=columns # Only load necessary columns
    )
    
    # Create dictionary mapping
    gene_map = dict(zip(gnomad_df['gene'], gnomad_df['gene_id']))
    return gene_map, gnomad_df

def load_sumstats_for_category(base_dir: Path, trait: str, category: str) -> Optional[pl.DataFrame]:
    """Loads and concatenates all sumstat files for a given trait and category."""
    pattern = f"genebass_{trait}_{category}_*.txt.bgz"
    file_pattern = str(base_dir / pattern)
    print(f"Scanning for files matching: {file_pattern}")
    
    # Find files matching the pattern
    files = glob.glob(file_pattern)
    if not files:
        print(f"  No files found for {category}.")
        return None

    print(f"  Found {len(files)} files for {category}. Loading...")
    # Use scan_csv for potentially large numbers of files, adding category column
    df = pl.scan_csv(file_pattern, separator='\t', comment_prefix='#', infer_schema_length=10000)\
            .with_columns(pl.lit(category).alias("functional_category"))\
            .collect()
    print(f"  Loaded {df.height} variants for {category}.")
    print(f"  Columns in {category} sumstats: {df.columns}")
    return df

def calculate_uncorrected_score(
    annot_df: pl.DataFrame, 
    af_col: str, 
    annot_names: List[str]
) -> Dict[str, float]:
    """
    Calculates the uncorrected burden score sum(2*AF*(1-AF)) for each annotation.

    Args:
        annot_df: DataFrame containing variants and one-hot annotation columns.
        af_col: Name of the allele frequency column.
        annot_names: List of annotation column names to calculate scores for.

    Returns:
        A list of uncorrected scores, one for each annotation name. Returns 0.0 for 
        annotations not present or with no variants.
    """
    scores = {}

    for name in annot_names:
        # Filter for variants belonging to this annotation
        annot_variants = annot_df.filter(pl.col("annotation") == name)
        if annot_variants.is_empty() or af_col not in annot_variants.columns:
            scores[name] = 0.0
            continue
            
        # Calculate score for this annotation's variants
        # Ensure AF column is float and handle potential nulls
        af_series = annot_variants[af_col].cast(pl.Float64)
        score_series = 2.0 * af_series * (1.0 - af_series)
        scores[name] = score_series.sum() # Handle potential null sum if all AFs were null

    
    # Ensure the order matches the input annot_names
    return scores

# --- Main Execution --- 

def main():
    print("Starting burden score computation script...")
    
    # --- Argument Parsing --- 
    parser = argparse.ArgumentParser(description="Compute LD-corrected burden scores from variant summary statistics.")
    parser.add_argument("--sumstats", default=VARIANT_DATA_DIR, help="Path to the *directory* containing summary statistics files (e.g., where genebass_TRAIT_CATEGORY_*.txt.bgz files are located)")
    parser.add_argument("--gene_map", default=GNOMAD_LOF_METRICS_FILE, help="Path to the gene symbol to Ensembl ID mapping file.")
    parser.add_argument("--ld_matrices", default=LD_MATRIX_DIR, help="Path to the directory containing LD matrices.")
    parser.add_argument("--output_prefix", "-o", default=R_OUTPUT_FILE, help="Prefix for the output files (e.g., Results/ukbb_ld_corrected_burden_scores)")
    parser.add_argument("--num_genes", "-n", type=int, default=None, help="Optional: Limit processing to the first N genes.")
    parser.add_argument("--no_save", default=False, action="store_true")
    # Optional: Add other parameters like AF limits if needed as arguments
    # parser.add_argument("--min_af", type=float, default=0.0, help="Minimum allele frequency threshold.")
    # parser.add_argument("--max_af", type=float, default=0.001, help="Maximum allele frequency threshold.")

    args = parser.parse_args()

    OUTPUT_DIR = Path(args.output_prefix).parent

    gene_to_id_map, gnomad_df = load_gnomad_mapping(GNOMAD_LOF_METRICS_FILE, GENE_ANNOT_NAMES)

    # 2. Load sumstats for each functional category
    all_sumstats_df = pl.read_csv(
        args.sumstats,
        separator='\t',
        has_header=True,
        comment_prefix='#',
        null_values="NA",
    ).filter((pl.col('AF') >= MIN_AF), (pl.col('AF') <= MAX_AF))

    # Determine genes to process based on --num_genes
    genes_to_process = sorted(list(all_sumstats_df['gene'].unique()))
    if args.num_genes is not None and args.num_genes > 0:
        genes_to_process = genes_to_process[:args.num_genes]
        print(f"\nLimiting processing to the first {len(genes_to_process)} genes.")

    genes_without_id = 0
    gene_results_list = []
    for gene_symbol in tqdm(genes_to_process, desc="Processing genes"): # Updated desc
        gene_id = gene_to_id_map.get(gene_symbol)
        if not gene_id:
            genes_without_id += 1
            continue
        ld_matrix_path = LD_MATRIX_DIR / f"{gene_id}.npz"
        ld_snplist_path = LD_MATRIX_DIR / f"{gene_id}.snplist"
        if not (ld_matrix_path.exists() and ld_snplist_path.exists()):
            genes_without_id += 1
            continue

        gene_sumstats_df = all_sumstats_df\
                        .filter(pl.col("gene") == gene_symbol)\
                        .rename({"AF": "AF_annot"}) # avoid name conflict

        snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True)\
                        .with_columns(('chr' + pl.col('SNP')).name.keep())

        ld_matrix = sparse.load_npz(ld_matrix_path)

        corrected_scores_dict = get_burden_score(
            matrices=[ld_matrix],
            matrix_snplists=[snplist_df],
            annot_snplists=[gene_sumstats_df],
            annot_af_name="AF_annot",
            annot_names=ANNOT_NAMES,
            merge_fields=['SNP'],
            AF_rtol=AF_RTOL,
        )[0] # Assuming it returns a list/tuple with the dict as the first element

        uncorrected_scores_dict = calculate_uncorrected_score(
            gene_sumstats_df,
            "AF_annot",
            ANNOT_NAMES
        )

        for category in ANNOT_NAMES:
            gene_results_list.append({
                "gene": gene_symbol,
                "gene_id": gene_id,
                "functional_category": category,
                "burden_score": corrected_scores_dict[category],
                "burden_score_no_ld": uncorrected_scores_dict[category],
            })

    print("\nConverting results to DataFrame...")
    gnomad_join_df = gnomad_df.select(['gene'] + GENE_ANNOT_NAMES)
    final_results_df = pl.DataFrame(gene_results_list) \
        .join(gnomad_join_df, on='gene', how='inner') \
        .group_by('gene', 'functional_category').first()

    print(f"Final results DataFrame shape: {final_results_df.shape}")
    print(f"Length before merging: {len(gene_results_list)}")
    print(f"Final results columns: {final_results_df.columns}")

    for category in FUNCTIONAL_CATEGORIES:
        # Define final columns including gnomAD features
        output_columns = [
            'gene',
            'gene_id',
            'functional_category',
            'burden_score',
            'burden_score_no_ld',
        ] + GENE_ANNOT_NAMES

        # Filter results for the current category
        category_df = final_results_df.filter(pl.col('functional_category') == category)

        # Construct output filename using the prefix and output directory
        output_filename = f"{args.output_prefix}_{category}_{MIN_AF}_{MAX_AF}.tsv"
        output_path = OUTPUT_DIR / Path(output_filename).name # Ensure it's placed in OUTPUT_DIR
        if not args.no_save:
            print(f"Writing {category_df.height} rows for category '{category}' to {output_path}")
            category_df.select(output_columns).write_csv(output_path, separator='\t')

if __name__ == "__main__":
    main()
