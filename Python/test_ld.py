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
from linkage_disequilibrium import get_burden_score, XtXOpt

# --- Define Paths --- 

TOP_DIR = "/Users/lukeoconnor/Dropbox"

# Directory containing the within-gene LD matrices (e.g., .npz files)
LD_MATRIX_DIR = Path(TOP_DIR) / "within_gene_ld_ukbb"

DENSE_LD_DIR = Path("/Users/lukeoconnor/Downloads")

# Directory containing variant-level data
VARIANT_DATA_DIR = Path(TOP_DIR) / "burdenEM_results/data/genebass/var_txt" # Corrected path

def main():
    
    genes_to_process = ['ENSG00000259040']
    genes_to_process = ['ENSG00000161939']
    for gene_symbol in genes_to_process: # Updated desc
        sparse_ld_path = LD_MATRIX_DIR / f"{gene_symbol}.npz"
        ld_snplist_path = LD_MATRIX_DIR / f"{gene_symbol}.snplist"
        dense_ld_path = DENSE_LD_DIR / f"chrom17_{gene_symbol}.npy"
        snplist_df = pl.read_csv(ld_snplist_path, separator='\t', has_header=True)

        # Load LD Matrix
        ld_matrix = sparse.load_npz(sparse_ld_path)
        dense_ld_matrix = np.load(dense_ld_path)
        print(ld_matrix.shape, dense_ld_matrix.shape)
        print(dense_ld_matrix)
        print(ld_matrix)
        # Print max absolute off-diagonal element of dense_ld_matrix
        off_diag_mask = ~np.eye(dense_ld_matrix.shape[0], dtype=bool)
        max_off_diag = np.abs(dense_ld_matrix[off_diag_mask]).max()
        print(f"Max absolute off-diagonal element: {max_off_diag}")

        n = ld_matrix.shape[0]
        operator = XtXOpt(XtX_by_n=ld_matrix,
                          allele_frequency=snplist_df['AF'].to_numpy(),
                          indices=np.arange(n),
                          shape=(n, n))
        print(operator.matrix.diagonal())

        # Prepare sumstats (renaming AF)
        # gene_sumstats_df_renamed = gene_sumstats_df.rename({ANNOT_AF_NAME: RENAMED_ANNOT_AF_NAME})
        # Call get_burden_score
        # corrected_scores_dict = get_burden_score(
        #     matrices=[ld_matrix],
        #     matrix_snplists=[snplist_df],
        #     annot_snplists=[gene_sumstats_df_renamed],
        #     annot_af_name=RENAMED_ANNOT_AF_NAME,
        #     annot_names=ANNOT_NAMES
        # )[0] # Assuming it returns a list/tuple with the dict as the first element



if __name__ == "__main__":
    main()