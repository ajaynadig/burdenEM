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

import os, logging

script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.append(str(script_dir))

# --- Import custom functions/classes --- 
from linkage_disequilibrium import get_burden_score #, XtXOpt, DiagOpt

# --- Define Paths --- 

TOP_DIR = "/Users/wlu/Dropbox (Partners HealthCare)"

# Directory containing the within-gene LD matrices (e.g., .npz files)
LD_MATRIX_DIR = Path(TOP_DIR) / "within_gene_ld_ukbb"

# Directory containing variant-level data
VARIANT_DATA_DIR = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/genebass/var_txt" # Corrected path

# Path to the gnomAD LoF metrics file (for gene ID mapping)
GNOMAD_LOF_METRICS_FILE = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/utility/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
GNOMAD_LOF_METRICS_FILE = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/utility/gnomad.v4.constraint.by_gene.tsv"

# Annotations to get from the gnomAD file
GENE_ANNOT_NAMES = ['oe_lof', 'oe_mis_pphen', 'oe_lof_upper', 'cds_length'] # v2
GENE_ANNOT_NAMES = ['lof.oe', 'mis_pphen.oe', 'lof.oe_ci.upper', 'cds_length'] # v4

# Output directory for results (keep for later)
OUTPUT_DIR = Path("./burden_scores_output")
OUTPUT_DIR.mkdir(exist_ok=True) # Create output directory if it doesn't exist


# --- Define Constants ---
# TRAIT = "50_NA" # Hard-coded trait name based on user request
TRAITS = ['30060_NA', '30220_NA', '21001_NA', '100017_NA', '20414_NA', '30740_NA', '30780_NA', '20022_NA', '20127_NA', '4079_NA', '20023_NA', '30650_NA', '30680_NA', '30750_NA', '21002_NA', '30190_NA', '30620_NA', 'WHR_custom_NA', '30770_NA', '78_NA', '4080_NA', '49_NA', '30020_NA', '189_NA', '30600_NA', '3062_NA', '50_NA', '30670_NA', '30730_NA', '30610_NA', '100022_NA', '48_NA', '30760_NA']

# >>> Define Output File for R <<<
R_OUTPUT_FILE = Path(TOP_DIR) / f"burdenEM/burdenEM_results/data/utility/burdenscore/ukbb_ld_corrected_burden_scores"

FUNCTIONAL_CATEGORIES = [
    'missense_benign', 
    'missense_damaging', 
    'pLoF', 
    'synonymous'
]
# FUNCTIONAL_CATEGORIES = ['synonymous']
ANNOT_AF_NAME = "AF"    # Allele frequency column in sumstats
RENAMED_ANNOT_AF_NAME = "AF_annot" # Name to use after renaming within the script
ANNOT_NAMES = FUNCTIONAL_CATEGORIES  # Annotation columns to process
MIN_AF = 0.0
MAX_AF = 0.001

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
    ).drop_nulls()
    
    # Create dictionary mapping
    gene_map = dict(zip(gnomad_df['gene'], gnomad_df['gene_id']))
    print(f"Created mapping for {len(gene_map)} genes.")
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
) -> List[float]:
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
        annot_variants = annot_df.filter(pl.col(name) == 1)
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

class Tee:
    """Duplicate writes to multiple streams."""
    def __init__(self, *streams):
        self.streams = streams
    def write(self, data):
        for s in self.streams:
            s.write(data)
    def flush(self):
        for s in self.streams:
            s.flush()

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
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files.")
    # Optional: Add other parameters like AF limits if needed as arguments
    # parser.add_argument("--min_af", type=float, default=0.0, help="Minimum allele frequency threshold.")
    # parser.add_argument("--max_af", type=float, default=0.001, help="Maximum allele frequency threshold.")

    args = parser.parse_args()

    # Define OUTPUT_DIR based on the output prefix
    OUTPUT_DIR = Path(args.output_prefix).parent

    # 1. Load gnomAD gene mapping
    gene_to_id_map, gnomad_df = load_gnomad_mapping(Path(args.gene_map), GENE_ANNOT_NAMES)

    # 2. Load sumstats for each functional category
    all_category_dfs = []
    sumstats_dir = Path(args.sumstats)

    orig_stdout = sys.stdout
    orig_stderr = sys.stderr

    for TRAIT in TRAITS:
        if Path(f'{args.output_prefix}_{TRAIT}_pLoF.tsv').exists() and not args.overwrite:
            continue
        logname = f'{args.output_prefix}_{TRAIT}.log'
        with open(logname, 'w') as logf:
            # Redirect both stdout & stderr through Tee
            sys.stdout = Tee(orig_stdout, logf)
            sys.stderr = Tee(orig_stderr, logf)
            print(f'Processing trait: {TRAIT}')
            for category in FUNCTIONAL_CATEGORIES:
                glob_pattern = sumstats_dir / f"genebass_{TRAIT}_{category}_*.txt.bgz"
                file_list = sorted(list(glob_pattern.parent.glob(glob_pattern.name)))

                if not file_list:
                    print(f"Warning: No files found matching pattern for category '{category}'. Skipping.")
                    continue

                print(f"Scanning for files matching: {glob_pattern}")
                print(f"  Found {len(file_list)} files for {category}. Loading...")
                
                category_dfs = []
                # Define expected dtypes to handle potential parsing issues, especially for AC_cases
                expected_dtypes = {
                    "gene": pl.Utf8,
                    "beta": pl.Float64,
                    "N": pl.Int64,
                    "variant_variance": pl.Float64,
                    "phenotype_key": pl.Utf8,
                    "description": pl.Utf8,
                    "CHR": pl.Utf8, # Load as string initially
                    "POS": pl.Int64,
                    "AC_cases": pl.Float64, # Load as float first due to potential 'NA'
                    "AF": pl.Float64,
                    "AF_total": pl.Float64,
                    "prevalence": pl.Float64,
                    "trait_type": pl.Utf8,
                    "functional_category": pl.Utf8
                }
                
                for file_path in file_list:
                    try:
                        df = pl.read_csv(
                            file_path,
                            separator='\t',
                            has_header=True,
                            comment_prefix='#',
                            null_values="NA",
                            dtypes=expected_dtypes
                        )
                        # Add the one-hot column for this category
                        df = df.with_columns(pl.lit(1).cast(pl.Int8).alias(category))
                        category_dfs.append(df)
                    except Exception as e:
                        print(f"Error loading file {file_path}: {e}")
                        # Decide if you want to skip the file or exit
                        # continue # Skip this file
                        # sys.exit(1) # Exit script

                if category_dfs:
                    concatenated_category_df = pl.concat(category_dfs)
                    print(f"  Loaded {concatenated_category_df.height} variants for {category}.")
                    print(f"  Columns in {category} sumstats: {concatenated_category_df.columns}")
                    all_category_dfs.append(concatenated_category_df)
                else:
                    print(f"Warning: No data loaded for category '{category}'. Skipping.")

            if not all_category_dfs:
                print("Error: No summary statistics data loaded for any category. Exiting.")
                sys.exit(1)

            # 3. Concatenate all categories and sort
            print("Concatenating data from all categories...")
            all_sumstats_df = pl.concat(all_category_dfs, how='diagonal') # Use diagonal concatenation
            total_variant_count = all_sumstats_df.height
            print(f"Total variants loaded across all categories: {total_variant_count}")

            # Fill NaNs in the newly added category columns with 0
            all_sumstats_df = all_sumstats_df.with_columns([
                pl.col(cat).fill_null(0) for cat in FUNCTIONAL_CATEGORIES if cat in all_sumstats_df.columns
            ])

            # Ensure CHR is string and remove 'chr' prefix, and POS is integer for joining
            all_sumstats_df = all_sumstats_df.with_columns([
                pl.col("CHR").cast(pl.Utf8).str.replace("chr", ""),
                pl.col("POS").cast(pl.Int64)
            ])
            # Also ensure the new category columns are Int8
            for category in FUNCTIONAL_CATEGORIES:
                if category in all_sumstats_df.columns:
                    all_sumstats_df = all_sumstats_df.with_columns(pl.col(category).cast(pl.Int8))

            print(f"Sumstats dtypes after type casting and one-hot encoding: {all_sumstats_df.dtypes}")

            # 3a. Filter by Allele Frequency
            print(f"Filtering variants by Allele Frequency ({MIN_AF} <= AF <= {MAX_AF})...")
            filtered_variant_count_before_unique = all_sumstats_df.height # Store count before filtering
            all_sumstats_df = all_sumstats_df.filter(
                (pl.col(ANNOT_AF_NAME) >= MIN_AF) & (pl.col(ANNOT_AF_NAME) <= MAX_AF)
            )
            filtered_variant_count = all_sumstats_df.height
            print(f"Variants remaining after AF filtering: {filtered_variant_count} (out of {filtered_variant_count_before_unique})")

            # --- Add step to remove duplicate rows based on CHR, POS, gene ---
            print("Removing duplicate variant entries (based on CHR, POS, gene)...")
            unique_check_cols = ['CHR', 'POS', 'gene']
            # Ensure required columns exist before attempting unique
            if all(col in all_sumstats_df.columns for col in unique_check_cols):
                n_before_unique = all_sumstats_df.height
                all_sumstats_df = all_sumstats_df.unique(subset=unique_check_cols, keep='first', maintain_order=True)
                n_after_unique = all_sumstats_df.height
                print(f"Removed {n_before_unique - n_after_unique} duplicate variant entries. {n_after_unique} variants remaining.")
                filtered_variant_count = n_after_unique # Update count after deduplication
            else:
                print(f"Warning: Could not perform duplicate check. Missing one or more columns: {unique_check_cols}")
                # Decide how to handle this - maybe raise error if these are critical?
                # raise ValueError(f"Missing essential columns for duplicate check: {unique_check_cols}")

            print("Sorting by gene...")
            # Sort after potential duplicate removal
            all_sumstats_df = all_sumstats_df.sort("gene")
            unique_genes = all_sumstats_df["gene"].unique().to_list()
            print(f"Found {len(unique_genes)} unique genes after AF filtering.")

            # uniqueness check
            original_length = len(all_sumstats_df)
            filtered_length = len(all_sumstats_df.group_by(['gene','POS']).first())
            print(f"Found {original_length - filtered_length} with the same gene and position.")

            # 4. Check for LD file existence per gene
            print("Checking for corresponding LD matrix and snplist files...")
            total_merged_variants = 0 # Counter for successfully merged variants
            genes_without_id = 0
            genes_without_ld = 0
            gene_results_list = [] # List to store results in long format
            genes_processed_count = 0 # Count successful corrected score calculations
            genes_other_error_count = 0 # Count other errors during score calculation

            discordant_genes_info = [] # List to store info on genes with score ratio < 0.5

            # Determine genes to process based on --num_genes
            unique_genes = sorted(list(all_sumstats_df['gene'].unique()))
            if args.num_genes is not None and args.num_genes > 0:
                genes_to_process = unique_genes[:args.num_genes]
                print(f"\nLimiting processing to the first {len(genes_to_process)} genes.")
            else:
                genes_to_process = unique_genes
                print(f"\nProcessing all {len(genes_to_process)} unique genes found.")

            # genes_to_process = unique_genes # Process all genes
            # print("\nLIMITING TO FIRST 100 GENES FOR DEBUGGING")
            # genes_to_process = genes_to_process[:100]

            for gene_symbol in tqdm(genes_to_process, desc="Processing genes"): # Updated desc
                # Reset flags/variables for each gene
                error_message = None
                gene_id_for_output = 'NA' # Default to NA
                ld_data_available = True  # Assume LD is available initially
                corrected_scores_dict = None
                uncorrected_scores_dict = None

                # Get gene ID
                gene_id = gene_to_id_map.get(gene_symbol)
                if gene_id:
                    gene_id_for_output = gene_id # Store valid ID for output
                else:
                    genes_without_id += 1

                # Check for LD files ONLY if we have a valid gene ID
                if gene_id_for_output != 'NA':
                    ld_matrix_path = LD_MATRIX_DIR / f"{gene_id_for_output}.npz"
                    ld_snplist_path = LD_MATRIX_DIR / f"{gene_id_for_output}.snplist"
                    if not (ld_matrix_path.exists() and ld_snplist_path.exists()):
                        genes_without_ld += 1
                        ld_data_available = False
                else:
                    # If gene ID is NA, LD is not available
                    ld_data_available = False
                    # Avoid double counting if ID was already missing
                    # genes_without_ld += 1 # Counted implicitly by ID lookup failure

                # Get the sumstats slice for this gene
                gene_sumstats_df = all_sumstats_df.filter(pl.col("gene") == gene_symbol)

                # --- Attempt LD-Corrected Score Calculation ---
                if ld_data_available:
                    # Load snplist (needed for shape check and input)
                    snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True)

                    # Prepare snplist (assuming this is needed for get_burden_score)
                    first_snp_id = snplist_df['SNP'][0]
                    chrom = first_snp_id.split(':')[0]
                    snplist_df = snplist_df.with_columns([
                        pl.lit(chrom).cast(pl.Utf8).alias("CHR"),
                        pl.col("POS").cast(pl.Int64)
                    ]).select(["CHR", "POS", "REF", "ALT", "SNP", "AF"])

                    # Load LD Matrix
                    ld_matrix = sparse.load_npz(ld_matrix_path)

                    # Prepare sumstats (renaming AF)
                    gene_sumstats_df_renamed = gene_sumstats_df.rename({ANNOT_AF_NAME: RENAMED_ANNOT_AF_NAME})
                    # Call get_burden_score
                    corrected_scores_dict = get_burden_score(
                        matrices=[ld_matrix],
                        matrix_snplists=[snplist_df],
                        annot_snplists=[gene_sumstats_df_renamed],
                        annot_af_name=RENAMED_ANNOT_AF_NAME,
                        annot_names=ANNOT_NAMES
                    )[0] # Assuming it returns a list/tuple with the dict as the first element
                    genes_processed_count += 1 # Count success here


                # --- Attempt Uncorrected Score Calculation (Always happens) ---
                if gene_sumstats_df.height > 0: # Only calculate if there are variants for the gene
                    uncorrected_scores_dict = calculate_uncorrected_score(
                        gene_sumstats_df,
                        ANNOT_AF_NAME,
                        ANNOT_NAMES
                    )

                # --- Append results for each category --- #
                for cat in ANNOT_NAMES:
                    # Use .get() which returns None if key missing or dict is None
                    corrected_score = corrected_scores_dict.get(cat) if corrected_scores_dict is not None else None
                    uncorrected_score = uncorrected_scores_dict.get(cat) if uncorrected_scores_dict is not None else None

                    # --- Check for discordant scores (ratio < 0.5) ---
                    if corrected_score is not None and uncorrected_score is not None and uncorrected_score != 0:
                        ratio = corrected_score / uncorrected_score
                        if ratio < 0.5:
                            discordant_genes_info.append({
                                "gene": gene_symbol,
                                "category": cat,
                                "burden_score_ld": corrected_score,
                                "burden_score_no_ld": uncorrected_score,
                                "ratio": ratio
                            })
                    # -----------------------------------------------------

                    gene_results_list.append({
                        "gene": gene_symbol,
                        "gene_id": gene_id_for_output, # Will be 'NA' if lookup failed
                        "min_af": MIN_AF,
                        "max_af": MAX_AF,
                        "functional_category": cat,
                        "burden_score": corrected_score, # Will be None if LD/calc failed
                        "burden_score_no_ld": uncorrected_score # Will be None if calc failed
                    })

            print("\n--- Final Merge Results ---")
            print(f"Total variants loaded from summary statistics: {total_variant_count}")
            print(f"Total variants after AF filter ({MIN_AF} <= AF <= {MAX_AF}): {filtered_variant_count}")
            print(f"Total variants successfully merged with LD snplists (on CHR & POS): {total_merged_variants}")
            print(f"Number of unique genes without an Ensembl ID mapping: {genes_without_id}")
            print(f"Number of unique genes missing/failing LD/merge/shape check/score calc: {len(unique_genes) - genes_without_id - len(gene_results_list)}")
            print(f"Number of genes for which burden scores were computed: {genes_processed_count}")

            # Print results for the first few genes as a sample
            print("\n--- Calculating Mean Scores per Category ---")
            corrected_scores_by_cat = {cat: [] for cat in ANNOT_NAMES}
            uncorrected_scores_by_cat = {cat: [] for cat in ANNOT_NAMES}

            for result in gene_results_list:
                cat = result['functional_category']
                corrected_score = result['burden_score']
                uncorrected_score = result['burden_score_no_ld']

                if corrected_score is not None: # Append only if the score exists for this gene/cat
                    corrected_scores_by_cat[cat].append(corrected_score)
                if uncorrected_score is not None:
                    uncorrected_scores_by_cat[cat].append(uncorrected_score)

            # Calculate and print means
            for cat in ANNOT_NAMES:
                mean_corrected = np.mean(corrected_scores_by_cat[cat]) if corrected_scores_by_cat[cat] else 0.0
                mean_uncorrected = np.mean(uncorrected_scores_by_cat[cat]) if uncorrected_scores_by_cat[cat] else 0.0
                count_corrected = len(corrected_scores_by_cat[cat])
                count_uncorrected = len(uncorrected_scores_by_cat[cat])
                print(f"Category: {cat}")
                print(f"  Mean Corrected Score ({count_corrected} genes): {mean_corrected:.4e}")
                print(f"  Mean Uncorrected Score ({count_uncorrected} genes): {mean_uncorrected:.4e}")

            # --- Print Discordant Gene Info --- 
            if discordant_genes_info:
                print("\n--- Genes with Discordant Scores (LD-Corrected / Uncorrected < 0.5) ---")
                print(f"Found {len(discordant_genes_info)} gene-category pairs with potential issues:")
                # Convert to DataFrame for nice printing
                discordant_df = pl.DataFrame(discordant_genes_info)
                with pl.Config(tbl_rows=-1, tbl_cols=-1): # Print all rows/cols
                    print(discordant_df)
            else:
                print("\nNo gene-category pairs found with highly discordant scores (ratio < 0.5).")
            # -------------------------------------

            # 6. Convert results list to DataFrame and save
            print("\nConverting results to DataFrame...")
            if gene_results_list:
                gene_results_df = pl.DataFrame(gene_results_list)

                # --- Join with gnomAD features --- #
                print("\nJoining burden scores with gnomAD gene features...")
                # Select only the necessary columns from gnomad_df for the join
                gnomad_join_df = gnomad_df.select(['gene'] + GENE_ANNOT_NAMES)
                # Perform a left join
                final_results_df = gene_results_df.join(gnomad_join_df, on='gene', how='left')

                print(f"Final results DataFrame shape: {final_results_df.shape}")
                print(f"Final results columns: {final_results_df.columns}")
                # print(final_results_df.head())

                # Write output files per functional category
                print("\nWriting output files per functional category...")
                unique_categories_in_results = final_results_df['functional_category'].unique().to_list()

                for category in unique_categories_in_results:
                    # Define final columns including gnomAD features
                    output_columns = [
                        'gene',
                        'gene_id',
                        'functional_category',
                        'burden_score',
                        'burden_score_no_ld',
                        'min_af',
                        'max_af'
                    ] + GENE_ANNOT_NAMES

                    # Filter results for the current category
                    category_df = final_results_df.filter(pl.col('functional_category') == category)

                    # Construct output filename using the prefix and output directory
                    output_filename = f"{args.output_prefix}_{TRAIT}_{category}.tsv"
                    output_path = OUTPUT_DIR / Path(output_filename).name # Ensure it's placed in OUTPUT_DIR

                    print(f"Writing {category_df.height} rows for category '{category}' to {output_path}")
                    category_df.select(output_columns).write_csv(output_path, separator='\t')
            else:
                print("No results were generated to write to file.")
            
            sys.stdout.flush()
            sys.stderr.flush()
               # Restore originals before next iteration
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr

        # (Optional) notify on console that the log file is closed
        print(f"[logged to {logname}]") 

if __name__ == "__main__":
    main()


# python3 compute_burden_scores.py \
#   --sumstats <path_to_sumstats_directory> \
#   --gene_map <path_to_gnomad_file> \
#   --ld_matrices <path_to_ld_matrix_directory> \
#   -o <output_file_prefix> \
#   # Optional: -n <number_of_genes_to_process> \
#   > <output_file_prefix>.log 2>&1 
# 
# python3 Python/compute_burden_scores.py 