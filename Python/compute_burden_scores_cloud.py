import subprocess
import os
import sys
# subprocess.check_call([sys.executable, "-m", "pip", "install", "polars"])
# subprocess.check_call([sys.executable, "-m", "pip", "install", "tqdm"])
# subprocess.check_call([sys.executable, "-m", "pip", "install", "google"])
# subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])
# subprocess.check_call([sys.executable, "-m", "pip", "install", "scipy"])
from pathlib import Path
import numpy as np
import polars as pl
from scipy import sparse
from typing import Optional, Dict, List, Tuple
import glob # For finding files
from tqdm import tqdm # For progress bar
import matplotlib.pyplot as plt
import argparse
from google.cloud import storage
from io import BytesIO, StringIO

script_dir = Path(__file__).parent.resolve()
if str(script_dir) not in sys.path:
    sys.path.append(str(script_dir))

# --- Import custom functions/classes --- 
from linkage_disequilibrium import get_burden_score

# --- Define Paths --- 
dataname = 'aou_amr'

TOP_DIR = "/Users/wlu/Dropbox (Partners HealthCare)"

# Directory containing the within-gene LD matrices (e.g., .npz files)
LD_MATRIX_DIR = Path(TOP_DIR) / "within_gene_ld_ukbb"

# Directory containing variant-level data
# VARIANT_DATA_FILE = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/utility/test_aou_afr_exome_vep2.txt.bgz"
VARIANT_DATA_FILE = Path(TOP_DIR) / f"burdenEM/burdenEM_results/data/utility/test_{dataname}_exome_vep2.txt.bgz"

# Path to the gnomAD LoF metrics file (for gene ID mapping)
GNOMAD_LOF_METRICS_FILE = Path(TOP_DIR) / "burdenEM/burdenEM_results/data/utility/gnomad.v4.constraint.by_gene.tsv"

# Annotations to get from the gnomAD file
GENE_ANNOT_NAMES = ['lof.oe', 'mis_pphen.oe', 'lof.oe_ci.upper', 'cds_length', 'lof.mu', 'mis.mu']

# Output directory for results (keep for later)
OUTPUT_DIR = Path("./burden_scores_output")
OUTPUT_DIR.mkdir(exist_ok=True) # Create output directory if it doesn't exist

# >>> Define Output File for R <<<
R_OUTPUT_FILE = Path(TOP_DIR) / f"burdenEM/burdenEM_results/data/utility/{dataname}_ld_corrected_burden_scores"

# --- Define Constants ---
TRAIT = "height" # Hard-coded trait name based on user request

FUNCTIONAL_CATEGORIES = [
    'missense_benign',
    'missense_damaging',
    'pLoF', 
    'synonymous'
]
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


# --- Main Execution --- 

def main():
    print("Starting burden score computation script...")
    
    # --- Argument Parsing --- 
    parser = argparse.ArgumentParser(description="Compute LD-corrected burden scores from variant summary statistics.")
    parser.add_argument("--variant_data", default=VARIANT_DATA_FILE, help="Path to the *directory* containing summary statistics files (e.g., where genebass_TRAIT_CATEGORY_*.txt.bgz files are located)")
    parser.add_argument("--gene_map", default=GNOMAD_LOF_METRICS_FILE, help="Path to the gene symbol to Ensembl ID mapping file.")
    parser.add_argument("--ld_matrices", 
    # default=LD_MATRIX_DIR, 
    help="Path to the directory containing LD matrices.")
    parser.add_argument("--output_prefix", "-o", default=R_OUTPUT_FILE, help="Prefix for the output files (e.g., Results/ukbb_ld_corrected_burden_scores)")
    parser.add_argument("--num_genes", "-n", type=int, default=None, help="Optional: Limit processing to the first N genes.")
    parser.add_argument("--no_save", default=False, action="store_true")
    parser.add_argument('--gcs_bucket', type=str, default=None,
                    help='Name of the GCS bucket containing LD .npz and .snplist files')
    parser.add_argument('--gcs_prefix', type=str, default='',
                        help='Prefix path inside the GCS bucket (e.g., "within_gene_ld_ukbb/")')
    # Optional: Add other parameters like AF limits if needed as arguments
    # parser.add_argument("--min_af", type=float, default=0.0, help="Minimum allele frequency threshold.")
    # parser.add_argument("--max_af", type=float, default=0.001, help="Maximum allele frequency threshold.")

    args = parser.parse_args()
    
    # Initialize GCS client if requested
    if args.gcs_bucket:
        if storage is None:
            raise ImportError("google-cloud-storage is required to read from GCS. Please install it with `pip install google-cloud-storage`.")
        gcs_client = storage.Client()
        gcs_bucket = gcs_client.bucket(args.gcs_bucket)

    # Define LD matrix directory (local or GCS)
    if args.gcs_bucket:
        LD_MATRIX_DIR = None  # we will build blob paths dynamically

    OUTPUT_DIR = Path(args.output_prefix).parent
    print(VARIANT_DATA_FILE)
    gene_to_id_map, gnomad_df = load_gnomad_mapping(GNOMAD_LOF_METRICS_FILE, GENE_ANNOT_NAMES)

    # 2. Load sumstats for each functional category
    variants_df = pl.read_csv(
        VARIANT_DATA_FILE,
        separator='\t',
        has_header=True,
        comment_prefix='#',
        null_values="NA",
        ).filter((pl.col('AF') >= MIN_AF), (pl.col('AF') <= MAX_AF)) \
        .group_by('SNP').first()

    # Determine genes to process based on --num_genes
    genes_to_process = sorted([g for g in variants_df['gene'].unique() if g is not None])
    if args.num_genes is not None and args.num_genes > 0:
        genes_to_process = genes_to_process[:args.num_genes]
        print(f"\nLimiting processing to the first {len(genes_to_process)} genes.")

    genes_without_id = 0
    gene_results_list = []
    for gene_symbol in tqdm(genes_to_process, desc="Processing genes"): # Updated desc
        gene_id = gene_to_id_map.get(gene_symbol)
        print(gene_id)
        if not gene_id:
            genes_without_id += 1
            continue

        # Read LD matrix and SNP list either from GCS or local FS
        if args.gcs_bucket:
            # Build blob names
            print(args.gcs_prefix)
            npz_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.npz")
            snp_blob = gcs_bucket.blob(f"{args.gcs_prefix}{gene_id}.snplist")
            # Check existence
            if not (npz_blob.exists() and snp_blob.exists()):
                genes_without_id += 1
                continue
            # Load SNP list into Polars DataFrame
            snp_text = snp_blob.download_as_text()
            if dataname == 'ukbb':
              snplist_df = pl.read_csv(StringIO(snp_text), separator='\t', has_header=True)\
                            .with_columns(('chr' + pl.col('SNP')).name.keep())
            else:
              snplist_df = pl.read_csv(StringIO(snp_text), separator='\t', has_header=True) \
                              .with_columns((pl.col('SNP')).alias('SNP'))
            # Load sparse LD matrix from NPZ
            npz_bytes = npz_blob.download_as_bytes()
            ld_matrix = sparse.load_npz(BytesIO(npz_bytes))

            # # Visualize sparsity
            # plt.figure(figsize=(8, 8))
            # plt.spy(ld_matrix, markersize=0.5)
            # plt.title('Sparsity Pattern of the Matrix')
            # plt.xlabel('Columns')
            # plt.ylabel('Rows')
            # plt.show()
        else:
            ld_matrix_path = LD_MATRIX_DIR / f"{gene_id}.npz"
            ld_snplist_path = LD_MATRIX_DIR / f"{gene_id}.snplist"
            if not (ld_matrix_path.exists() and ld_snplist_path.exists()):
                genes_without_id += 1
                continue
            snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True) \
                            .with_columns((pl.col('SNP')).alias('SNP'))
            ld_matrix = sparse.load_npz(ld_matrix_path)

        gene_sumstats_df = variants_df\
                        .filter(pl.col("gene") == gene_symbol)\
                        .rename({"AF": "AF_annot"}) # avoid name conflict
        # snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True)\
        #                 .with_columns(('chr' + pl.col('SNP')).name.keep())

        # ld_matrix = sparse.load_npz(ld_matrix_path)

        corrected_scores_list, uncorrected_scores_list = get_burden_score(
            matrices=[ld_matrix],
            matrix_snplists=[snplist_df],
            annot_snplists=[gene_sumstats_df],
            annot_af_name="AF_annot",
            annot_names=ANNOT_NAMES,
            merge_fields=['SNP'],
            AF_rtol=AF_RTOL,
        )

        corrected_scores_dict = corrected_scores_list[0]
        uncorrected_scores_dict = uncorrected_scores_list[0]
        
        print(corrected_scores_dict)
        print(uncorrected_scores_dict)

        for category in ANNOT_NAMES:
            gene_results_list.append({
                "gene": gene_symbol,
                "gene_id": gene_id,
                "functional_category": category,
                "burden_score": corrected_scores_dict[category],
                "burden_score_no_ld": uncorrected_scores_dict[category],
            })
    print(gene_results_list)
    print("\nConverting results to DataFrame...")
    results_df = pl.DataFrame(gene_results_list)
    
    gnomad_join_df = gnomad_df.select(['gene'] + GENE_ANNOT_NAMES)
    
    print("results_df columns:", results_df.columns)
    print("gnomad_join_df columns:", gnomad_join_df.columns)
    final_results_df = results_df \
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
        if args.no_save:
            print(f"Would write {category_df.height} rows for category '{category}' to {output_path}")
        else:
            print(f"Writing {category_df.height} rows for category '{category}' to {output_path}")
            category_df.select(output_columns).write_csv(output_path, separator='\t')

if __name__ == "__main__":
    main()
