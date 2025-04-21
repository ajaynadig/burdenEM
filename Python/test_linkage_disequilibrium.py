import numpy as np
from scipy import sparse
import polars as pl
import sys
import os
from typing import List, Tuple

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Python.linkage_disequilibrium import get_burden_score, XtXOpt


def correlation_from_covariance(XtX_by_n, mean):
    covariance = XtX_by_n - np.outer(mean,mean)
    v = np.sqrt(np.diag(covariance))
    v[v==0] = 1
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    # print(correlation[:10,:10])
    # correlation[covariance == 0] = 0
    return correlation


def test_get_burden_score():
    """
    Test the get_burden_score function with sample data.
    """
    # Load the sample data
    sample_npz_path = '/Users/lukeoconnor/Dropbox/within_gene_ld_ukbb/ENSG00000000419.npz'
    sample_snplist_path = '/Users/lukeoconnor/Dropbox/within_gene_ld_ukbb/ENSG00000000419.snplist'
    
    # Load the sparse matrix
    npz_data = np.load(sample_npz_path)
    matrix_shape = tuple(npz_data['shape'])
    matrix_format = npz_data['format'].tobytes().decode('utf-8')
    
    # Create the sparse matrix
    matrix = sparse.csr_matrix(
        (npz_data['data'], npz_data['indices'], npz_data['indptr']),
        shape=matrix_shape
    )
    
    # Load the snplist
    snplist = pl.read_csv(sample_snplist_path, separator='\t')
    
    # Create a modified version of the snplist for annotation
    # Drop some variants and add some new ones
    annot_snplist = snplist.clone()
    
    # Drop 10% of the variants
    num_variants = len(annot_snplist)
    drop_indices = np.random.choice(num_variants, size=int(num_variants * 0.1), replace=False)
    
    # Create a list of indices to keep
    keep_indices = [i for i in range(num_variants) if i not in drop_indices]
    annot_snplist = annot_snplist.filter(pl.col('SNP').is_in(annot_snplist.select('SNP').filter(pl.arange(0, num_variants).is_in(keep_indices)).to_series()))
    annot_snplist = annot_snplist.rename({'AF': 'annot_AF'})

    
    # Add some new variants
    new_variants = []
    for i in range(5):
        pos = int(annot_snplist['POS'].max()) + i + 1
        new_variants.append({
            'POS': pos,
            'REF': 'A',
            'ALT': 'G',
            'SNP': f'20:{pos}:G:A',
            'MISSINGNESS': 0.001,
            'annot_AF': 0.0001
        })
    
    new_variants_df = pl.DataFrame(new_variants)
    variants_with_ld = range(len(annot_snplist))
    annot_snplist = pl.concat([annot_snplist, new_variants_df])

    # Add a new column with all 1's to annot_snplist
    annot_snplist = annot_snplist.with_columns(pl.lit(1).alias('annot_column'))
    
    # Step 1: Convert the matrix into a correlation matrix
    allele_frequencies = snplist['AF'].to_numpy()
    
    # Convert to correlation matrix using the provided function
    matrix_dense = matrix.toarray()
    
    correlation_matrix = correlation_from_covariance(matrix_dense, allele_frequencies)

    # Filter the correlation matrix
    correlation_matrix = correlation_matrix[keep_indices, :][:, keep_indices]
    
    # Call the function
    matrices = [matrix]
    matrix_snplists = [snplist]
    annot_snplists = [annot_snplist]
    
    # Test the function
    burden_scores = get_burden_score(matrices, matrix_snplists, annot_snplists, 
                                    'annot_AF', ['annot_column'])
    result = burden_scores[0]['annot_column']
    print(f"Function executed successfully wtih result {result}")
    
    v = annot_snplist.select('annot_column').to_numpy().astype('float64')
    af = annot_snplist.select('annot_AF').to_numpy()
    v *= np.sqrt(2 * af * (1-af))
    Rv_true = v.copy()
    Rv_true[variants_with_ld] = correlation_matrix @ v[variants_with_ld]
    result_true = np.dot(v.flatten(), Rv_true.flatten())
    print(correlation_matrix[:5, :5])
    print(v[:10])
    print(f"True value: {result_true}, calculated value: {result}, no LD value: {np.dot(v.flatten(), v.flatten())}")
    assert np.isclose(result, result_true)

    print("All tests passed!")
    return True



if __name__ == "__main__":
    result = test_get_burden_score()
    if result:
        print("Test completed successfully")
    else:
        print("Test failed")
