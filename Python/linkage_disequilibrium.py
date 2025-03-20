import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, aslinearoperator
from dataclasses import dataclass
from functools import cached_property
from typing import List, Callable, Optional, Tuple
import polars as pl


@dataclass
class DiagOpt(LinearOperator):
    diag: np.ndarray

    @property
    def shape(self):
        return (self.diag.shape[0], self.diag.shape[0])

    @property
    def dtype(self):
        return float

    @property
    def matrix(self):
        return np.diag(self.diag.flatten()) 

    @property
    def data(self):
        return np.diag(self.diag.flatten()) 

    def _matmat(self, v):
        return self.diag.reshape(-1, 1) * v

@dataclass
class XtXOpt(LinearOperator):
    """Represents the correlation matrix corresponding to the sparse matrix X'X"""
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
        
        std_deviations = np.sqrt(self.XtX_by_n.diagonal() - self.allele_frequency**2)
        std_deviations[std_deviations==0] = 1
        scaling_matrix = DiagOpt(diag=1.0 / std_deviations)
        
        XtX_by_n = aslinearoperator(self.XtX_by_n)
        
        return scaling_matrix @ (XtX_by_n - mu @ mu.T) @ scaling_matrix

    def _matvec(self, v):
        result = v.copy()
        result[self.indices] = self._operator @ v[self.indices]
        return result


def get_burden_score(matrices: List[sparse.spmatrix], 
                matrix_snplists: List[pl.dataframe],
                num_haplotypes: int,
                annot_snplists: List[pl.dataframe],
                annot_af_name: str,
                annot_names: List[str]):
    
    if annot_af_name == 'AF':
        raise ValueError("Please name the annotation allele frequencies something else")
    
    for i, matrix_snplist in enumerate(matrix_snplists):
        matrix_snplists[i] = matrix_snplist.with_row_count('matrix_index')
    
    merged_snplists = [
        matrix_snplist.join(annot_snplist, on='SNP', how='right')
        for matrix_snplist, annot_snplist in zip(matrix_snplists, annot_snplists)
    ]
    for i, snplist in enumerate(merged_snplists):
        original_length = len(snplist)
        merged_snplists[i] = snplist.group_by('SNP').first()
        new_length = len(merged_snplists[i])
        if new_length < original_length:
            print(f"Warning: {original_length - new_length} duplicate SNPs were discarded from matrix {i}")
    
        merged_snplists[i] = merged_snplists[i].with_row_count('merged_index')
    
    # TODO maybe check for discrepencies between allele frequencies
    
    burden_scores = []
    for ld_mat, snplist in zip(matrices, merged_snplists):
        idx_into_matrix = snplist.filter(pl.col('matrix_index').is_not_null()).select('matrix_index').to_numpy().flatten()
        ld_mat = ld_mat[idx_into_matrix, :][:, idx_into_matrix]
        print(f"Number of diagonal nonzeros: {np.sum(ld_mat.diagonal()!=0)}, Number of non-zero elements: {ld_mat.nnz}")
        
        
        idx_into_merged = snplist.filter(pl.col('matrix_index').is_not_null()).select('merged_index').to_numpy().flatten()
        af = snplist.filter(pl.col('matrix_index').is_not_null()).select('AF').to_numpy().flatten()
        
        n = snplist.shape[0]

        operator = XtXOpt(XtX_by_n=ld_mat / num_haplotypes,
                          allele_frequency=af,
                          indices=idx_into_merged,
                          shape=(n, n))
        
        result = {}
        af = snplist.select(annot_af_name).to_numpy()
        sqrt_2pq = np.sqrt(2 * af * (1-af)).flatten()
        for annot_name in annot_names:
            burden_weights = sqrt_2pq * snplist.select(annot_name).to_numpy().flatten()
            result[annot_name] = np.dot(burden_weights.flatten(), (operator @ burden_weights).flatten())
        burden_scores.append(result)

    return burden_scores