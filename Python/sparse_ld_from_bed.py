#!/usr/bin/env python3
"""
Estimate within-gene LD matrices for WES data.

Operates on:
    1. plink BED files
    2. gene-to-variant mapping file: a tab-delim file with columns gene ID and variant ID
    3. sample_file: text file with sample IDs to include
(1) and (2) should be per-chromosome

AI-written code that is definitely not optimized - in particular when processing each gene, 
the entire gene-to-variant file is scanned for lines matching that gene. Bed reading logic
might be similarly suboptimal. On UKBB exomes it took ~1hr to run a long chromosome 
with 96 cores, but you could probably save 10-100x with some optimizations.
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import sparse
from pysnptools.snpreader import Bed


def read_bim_file(chrom, plink_dir="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"):
    """
    Read the BIM file for a chromosome.
    
    Args:
        chrom: Chromosome number
        plink_dir: Directory containing PLINK files
        
    Returns:
        pandas DataFrame containing BIM data
    """
    bed_file = f"{plink_dir}/ukb23158_c{chrom}_b0_v1.bed"
    bim = pd.read_csv(f"{bed_file[:-4]}.bim", sep='\t', header=None,
                      names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    return bim


def process_gene(gene, chrom, bim=None, out_dir="DATA/WES_RELATED_SPARSE_LD", 
                plink_dir="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release",
                set_file="/mnt/project/kangcheng/WES_QC/genebass_var_gene_map.tsv.gz",
                sample_file="/mnt/project/kangcheng/WES_QC/qced_related_samples.eid"):
    """
    Process a single gene to estimate LD matrix.
    
    Args:
        gene: Gene ID
        chrom: Chromosome number
        bim: Pre-loaded BIM DataFrame (if None, will be loaded)
        out_dir: Output directory
        plink_dir: Directory containing PLINK files
        set_file: File containing gene-variant mapping
        sample_file: File containing sample IDs
    """
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    # Create a temporary file to list all genes we've processed
    # This helps track progress and ensure no files are missed
    with open(f"{out_dir}/processed_genes.txt", "a") as f:
        f.write(f"{gene}\n")
    
    try:
        # Load sample IDs
        with open(sample_file, 'r') as f:
            samples = [line.strip() for line in f]
        print(f"Loaded {len(samples)} samples")
        
        # Load variant IDs for the gene
        gene_variants = []
        with gzip.open(set_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if fields[0] == gene:
                    gene_variants.append(fields[1])
        
        if len(gene_variants) == 0:
            print(f"No variants found for gene {gene}")
            
            # Create empty output files for genes with no variants
            snp_info = pd.DataFrame(columns=['POS', 'REF', 'ALT', 'SNP', 'MISSINGNESS', 'AF'])
            snp_info.to_csv(f"{out_dir}/{gene}.snplist", sep='\t', index=False, na_rep='NA')
            
            # Create an empty LD matrix of shape (0,0)
            empty_ld = sparse.csr_matrix((0, 0))
            sparse.save_npz(f"{out_dir}/{gene}.npz", empty_ld)
            return
        
        print(f"Found {len(gene_variants)} variants for gene {gene}")
        
        # Read PLINK files
        bed_file = f"{plink_dir}/ukb23158_c{chrom}_b0_v1.bed"
        
        # Load BIM file if not provided
        if bim is None:
            bim = read_bim_file(chrom, plink_dir)
        
        # Filter BIM to variants in the gene
        bim_subset = bim[bim['SNP'].isin(gene_variants)].copy()
        
        # Read genotypes
        snp_filter = np.array(bim_subset['SNP'])
        
        # Set up the Bed reader
        bed = Bed(bed_file, count_A1=True)
        
        # Convert sample IDs into PLINK (FID, IID) format if needed
        formatted_samples = []
        for s in samples:
            if ' ' in s:
                formatted_samples.append(s.split())
            else:
                formatted_samples.append([s, s])
        
        # Get indices of samples and SNPs
        sample_indices = bed.iid_to_index(formatted_samples)
        snp_indices = bed.sid_to_index(snp_filter)
        
        # Get the actual SNP IDs that were loaded
        loaded_snp_ids = bed.sid[snp_indices]
        
        # Filter bim_subset to only include loaded SNPs in the same order
        bim_subset = bim_subset[bim_subset['SNP'].isin(loaded_snp_ids)].copy()
        bim_subset = pd.DataFrame(bim_subset.set_index('SNP').loc[loaded_snp_ids]).reset_index()
        
        # Read genotype data
        gene_geno = bed[sample_indices, snp_indices].read(order='A').val
        
        # Keep original genotype values (0, 1, 2)
        genotypes = gene_geno.astype(np.float32)
        
        # Calculate missingness rate for each variant
        nan_mask = np.isnan(genotypes)
        nan_count = nan_mask.sum()
        missingness = nan_mask.mean(axis=0)
        
        if nan_count > 0:
            print(f"Warning: {gene} has {nan_count} NaN values in genotype data")
            
        # For allele frequency calculation, replace NaNs with zeros
        genotypes_for_calc = np.nan_to_num(genotypes, nan=0.0)
        
        # Calculate allele frequencies (mean of genotype / 2)
        allele_freq = genotypes_for_calc.mean(axis=0) / 2.0
        
        # Sparse LD matrix as genotypes.T @ genotypes
        # Normalize by number of haplotypes (2 * number of samples)
        num_haplotypes = 2 * genotypes.shape[0]
        sparse_ld = sparse.csr_matrix((genotypes_for_calc.T @ genotypes_for_calc) / num_haplotypes)
        
        # Write out SNP info with allele frequencies
        snp_info = bim_subset[['POS', 'A1', 'A2']].copy()
        snp_info.columns = ['POS', 'REF', 'ALT']
        snp_info['SNP'] = loaded_snp_ids
        snp_info['MISSINGNESS'] = missingness
        snp_info['AF'] = allele_freq
        
        # Debug SNP info shape before writing
        print(f"SNP info shape: {snp_info.shape}, Allele freq length: {len(allele_freq)}")
        
        # Write to file
        snp_info.to_csv(f"{out_dir}/{gene}.snplist", sep='\t', index=False, na_rep='NA')
        
        # Verify the written file
        try:
            written_df = pd.read_csv(f"{out_dir}/{gene}.snplist", sep='\t')
            if 'AF' in written_df.columns:
                print(f"NaN count in written AF: {written_df['AF'].isna().sum()}/{len(written_df)}")
        except Exception as e:
            print(f"Error verifying written file: {str(e)}")
            
        # Save sparse matrix
        sparse.save_npz(f"{out_dir}/{gene}.npz", sparse_ld)
        
        print(f"Finished processing gene {gene}")
        
    except Exception as e:
        print(f"Error processing gene {gene}: {e}")
        import traceback
        traceback.print_exc()
        
        # Create empty output files for genes that fail
        try:
            snp_info = pd.DataFrame(columns=['POS', 'REF', 'ALT', 'SNP', 'MISSINGNESS', 'AF'])
            snp_info.to_csv(f"{out_dir}/{gene}.snplist", sep='\t', index=False, na_rep='NA')
            
            # Create an empty LD matrix of shape (0,0)
            empty_ld = sparse.csr_matrix((0, 0))
            sparse.save_npz(f"{out_dir}/{gene}.npz", empty_ld)
        except:
            pass


def process_gene_in_parallel(args):
    """
    Wrapper function for parallel processing compatibility.
    
    Args:
        args: Tuple containing (gene, chrom, out_dir, set_file, sample_ids, plink_dir)
    
    Returns:
        Result of process_gene
    """
    # Unpack arguments
    gene, chrom, out_dir, set_file, sample_ids, plink_dir = args
    
    # Call process_gene with unpacked arguments
    return process_gene(
        gene=gene,
        chrom=chrom,
        out_dir=out_dir,
        plink_dir=plink_dir,
        set_file=set_file,
        sample_file=sample_ids
    )


def run(chrom, num_genes=None, num_processes=None):
    """
    Process one chromosome at a time using parallel processing.
    
    Args:
        chrom: Chromosome to process
        num_genes: If specified, limits processing to this many genes
        num_processes: Number of processes to use for parallel execution
                      (default: number of available CPU cores)
    """
    import multiprocessing
    import shutil
    
    # Set default paths
    out_dir = "DATA/WES_RELATED_SPARSE_LD"
    plink_dir = "/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release"
    set_file = "/mnt/project/kangcheng/WES_QC/genebass_var_gene_map.tsv.gz"
    sample_file = "/mnt/project/kangcheng/WES_QC/qced_related_samples.eid"
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    # Clear the processed genes tracking file
    with open(f"{out_dir}/processed_genes.txt", "w") as f:
        f.write("")
        
    # Create a temp directory for this chromosome's files
    temp_dir = f"{out_dir}/temp_chr{chrom}"
    os.makedirs(temp_dir, exist_ok=True)
    
    # Get list of genes for this chromosome
    genes = []
    with gzip.open(set_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] not in genes:
                genes.append(fields[0])
    
    # Limit number of genes if specified
    if num_genes is not None:
        genes = genes[:num_genes]
    
    print(f"Chromosome {chrom}: {len(genes)} total genes (processing up to {num_genes if num_genes else 'all'})")
    
    # Prepare arguments for parallel processing
    arg_list = [(gene, chrom, temp_dir, set_file, sample_file, plink_dir) for gene in genes]
    
    # Set number of processes
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    print(f"Using {num_processes} processes")
    
    # Process genes in parallel
    with multiprocessing.Pool(num_processes) as pool:
        pool.map(process_gene_in_parallel, arg_list)
    
    print(f"Finished processing chromosome {chrom}")
    
    # Check which genes were actually processed
    processed_genes = []
    try:
        with open(f"{temp_dir}/processed_genes.txt", "r") as f:
            processed_genes = [line.strip() for line in f.readlines()]
    except:
        pass
        
    print(f"Successfully processed {len(processed_genes)} genes")
    
    # Create zip file of results
    zip_file = os.path.join(out_dir, f"chrom{chrom}_sparse.zip")
    print(f"Creating zip file: {zip_file}")
    os.system(f"cd {temp_dir} && zip chrom{chrom}_sparse.zip *.snplist *.npz")
    
    # Move the zip file to the output directory
    shutil.move(f"{temp_dir}/chrom{chrom}_sparse.zip", zip_file)
    
    # Copy all files to the main output directory
    for file in os.listdir(temp_dir):
        if file.endswith('.snplist') or file.endswith('.npz'):
            shutil.copy(f"{temp_dir}/{file}", f"{out_dir}/{file}")
    
    # Attempt upload to DNANexus
    try:
        import dxpy
        dest_folder = "/luke/within_gene_ld"
        print(f"Uploading {zip_file} to DNANexus folder: {dest_folder}")
        dxfile = dxpy.upload_local_file(zip_file, folder=dest_folder, wait_on_close=True)
        print(f"Upload complete. File ID: {dxfile.get_id()}")
    except ImportError:
        print("dxpy not available. Please upload the results manually:")
        print(f"  Source: {zip_file}")
        print(f"  Destination: /luke/within_gene_ld")
    except Exception as e:
        print(f"Error uploading to DNANexus: {str(e)}")
        print(f"Please upload the file manually: {zip_file}")
        
    # Clean up temp directory
    shutil.rmtree(temp_dir)
