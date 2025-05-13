import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, aslinearoperator
from scipy.sparse import diags
from dataclasses import dataclass
from functools import cached_property
from typing import List, Callable, Optional, Tuple
import polars as pl


@dataclass
class XtXOpt(LinearOperator):
    """
    Represents the correlation matrix corresponding to the sparse matrix X'X.
    
    For missing variants (indicated by `indices`), acts as the identity.
    """
    XtX_by_n: sparse.spmatrix     # X'X/num_haplotypes matrix
    allele_frequency: np.ndarray  # Column means of X
    indices: np.ndarray # indices into the vector that correspond to rows/cols of XtX_by_n
    shape: Tuple[int, int]
    
    def __post_init__(self):
        if np.any(self.XtX_by_n.data > 2):
            raise ValueError(f"All entries in XtX_by_n must be <= 2")
        
        if np.any(self.XtX_by_n.data < 0):
            raise ValueError("All entries in XtX_by_n must be non-negative")
        
        if np.any((self.allele_frequency < 0) | (self.allele_frequency > 1)):
            raise ValueError("All entries in allele_frequency must be between 0 and 1")
            
        if self.allele_frequency.ndim != 1:
            raise ValueError(f"allele_frequency must be a 1D array")
            
        if self.XtX_by_n.shape != (len(self.allele_frequency), len(self.allele_frequency)):
            raise ValueError(f"Shape of XtX_by_n does not match expected shape")
        
        if np.any(self.allele_frequency == 1):
            raise ValueError("Allele frequency cannot be 1")

        # Operator acts as the identity on missing variants, which is also
        # the desired behavior for variants with zero allele frequency
        self._discard_zero_allele_frequency()

    def _discard_zero_allele_frequency(self):
        af_is_zero = self.allele_frequency == 0

        xtx_is_zero = self.XtX_by_n.diagonal() == 0
        if np.any(af_is_zero != xtx_is_zero):
            raise ValueError("Xt1 cannot be 0 if XtX is not, and vice versa")

        self.allele_frequency = self.allele_frequency[~af_is_zero]
        self.XtX_by_n = self.XtX_by_n[~af_is_zero, :][:, ~af_is_zero]
        self.indices = self.indices[~af_is_zero]
    
    @property
    def dtype(self):
        return float
    
    @property
    def matrix(self):
        return self @ np.eye(self.shape[0])
    
    @cached_property
    def _operator(self):
        # Leverages scipy LinearOperator lazy execution
        mu = aslinearoperator(self.allele_frequency.reshape(-1, 1))  # Column vector
        
        # Account for ploidy: n is num_haplotypes, not num_individuals
        std_deviations = np.sqrt(self.XtX_by_n.diagonal() - 2 * self.allele_frequency**2)
        scaling_matrix = aslinearoperator(diags(1.0 / std_deviations))
        
        XtX_by_n = aslinearoperator(self.XtX_by_n)
        
        return scaling_matrix @ (XtX_by_n - 2 * mu @ mu.T) @ scaling_matrix

    def _matvec(self, v):
        result = v.copy()
        result[self.indices] = self._operator @ v[self.indices]
        return result


def get_burden_score(matrices: List[sparse.spmatrix], 
                matrix_snplists: List[pl.dataframe],
                annot_snplists: List[pl.dataframe],
                annot_af_name: str,
                annot_names: List[str],
                merge_fields: List[str]=['CHR', 'POS'],
                AF_rtol: float=0.5):
    """
    Computes the LD-corrected burden score from a set of LD matrices and annotation data.

    The LD-corrected burden score is w'Rw, where w_i=2*AF_i*(1-AF_i) and R is the LD correlation matrix.
    
    Args:
        matrices: List of sparse matrices X'X/n, where X is the 0-1-2 matrix of genotypes.
        matrix_snplists: List of Polars DataFrames corresponding to LD matrices.
        annot_snplists: List of Polars DataFrames containing allele frequencies,
            merge_fields, and functional categories.
        annot_af_name: Name of the column in the annotation data containing allele frequencies.
        annot_names: Names of the column in the annotation data containing functional categories.
        merge_fields: Names of fields to merge on between LD matrix and annotation data.
        AF_rtol: Variants are discarded if |AF - AF_annot| / (AF + AF_annot) > AF_rtol.
    
    Returns:
        List of dictionaries containing LD-corrected burden scores for each annotation.
    """
    if annot_af_name == 'AF':
        raise ValueError("Please name the annotation allele frequencies something else")
    
    merged_snplists = [
        matrix_snplist.with_row_count('matrix_index')\
                        .join(annot_snplist, on=merge_fields, how='right')\
                        .group_by(merge_fields).first()\
                        .with_row_count('merged_index')\
                        .with_columns(
                            ((pl.col(annot_af_name) - pl.col('AF')).abs()\
                                / (pl.col(annot_af_name) + pl.col('AF'))\
                                < AF_rtol)\
                                .alias('AF_mask'),\
                            pl.col('matrix_index').is_not_null()\
                                .alias('missing_mask')
                                    )
        for matrix_snplist, annot_snplist in zip(matrix_snplists, annot_snplists)
    ]
    
    burden_scores = []
    for ld_mat, snplist in zip(matrices, merged_snplists):
        # print(f"Fraction passing AF filter: {snplist.filter('missing_mask').select('AF_mask').mean().item()}")
        # Drop variants from the LD matrix that didn't merge with the annotation data
        idx_into_matrix = snplist.filter('missing_mask', 'AF_mask').select('matrix_index').to_numpy().flatten()
        ld_mat = ld_mat[idx_into_matrix, :][:, idx_into_matrix]

        ldmat_af = snplist.filter('missing_mask', 'AF_mask').select('AF').to_numpy().flatten()
        idx_into_merged = snplist.filter('missing_mask', 'AF_mask').select('merged_index').to_numpy().flatten()
        n = snplist.shape[0]
        operator = XtXOpt(XtX_by_n=ld_mat,
                          allele_frequency=ldmat_af,
                          indices=idx_into_merged,
                          shape=(n, n))
        
        result = {}
        for annot_name in annot_names:
            af_or_0 = snplist.with_columns(
                (pl.col(annot_af_name) * (pl.col("annotation") == annot_name)).alias(annot_name)
            ).select(annot_name).to_numpy().flatten()

            weights = np.sqrt(2 * af_or_0 * (1-af_or_0)).ravel()
            result[annot_name] = np.dot(weights, operator @ weights)
        burden_scores.append(result)

    return burden_scores