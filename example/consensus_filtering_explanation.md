# "Consensus Across Variant Classes" Filtering

## What It Does

The "consensus" step (lines 409-418 in `set_up_kim.R`) ensures that every gene in the final analysis has **valid (non-zero, non-NA) mutation rates for ALL variant classes** (PTV, Mis2, Mis1, Mis0, Syn).

## The Code

```r
kim_genes_across <- lapply(1:length(kim_data_obj$loop_vars$subsets), function(i) {
  processed_input <- get_genetic_data(i, kim_data_obj)
  input_df <- processed_input$genetic_data
  return(rownames(input_df))
})

kim_genes_consensus <- Reduce("intersect", kim_genes_across)
```

This:
1. Calls `get_genetic_data()` for each variant class (PTV, Mis2, Mis1, Mis0, Syn)
2. `get_genetic_data()` returns only genes with valid mutation rates for that class
3. Takes the **intersection** of all 5 gene lists
4. Only keeps genes that pass for ALL variant classes

## Why Are There Different Gene Counts?

**Some genes have zero mutation rates for certain variant classes:**

| Variant Class | Genes with ZERO mutation rate |
|---------------|-------------------------------|
| PTV | 0 genes |
| **Mis2** | **3,261 genes** ⚠️ |
| Mis1 | 200 genes |
| Mis0 | 0 genes |
| Syn | 0 genes |

**Why does Mis2 have so many zeros?**
- Mis2 = missense variants with MPC ≥ 2 (high-confidence damaging)
- Many genes have **no predicted high-damage missense sites** according to the mutation rate model
- These genes might:
  - Lack MPC scores for their possible missense mutations
  - Have all missense mutations predicted as benign (MPC < 2)
  - Be small genes with few missense-mutable positions

**Example scenario:**
- Gene X has valid mutation rates for PTV, Mis1, Mis0, and Syn
- But Gene X has **zero** mutation rate for Mis2 (no high-damage missense sites)
- The consensus step **excludes** Gene X from the analysis entirely

## Impact on Gene Count

Before consensus: Each variant class might have different numbers of valid genes
- PTV: ~9,000 genes
- Mis2: ~5,700 genes (excludes 3,261 with zero rates)
- Mis1: ~8,800 genes
- Others: ~9,000 genes

After consensus: **7,518 genes** (intersection of all)

**Lost in this step**: ~1,500 genes that are valid for SOME but not ALL variant classes

## Why Is This Done?

### Advantages:
1. **Consistent gene set** - All analyses use the same 7,518 genes
2. **Fair comparisons** - Can directly compare enrichment across variant classes
3. **Shared features** - LOEUF stratification is identical across all analyses
4. **Cleaner interpretation** - No need to explain why gene sets differ

### Disadvantages:
1. **Reduced power** - Loses ~1,500 genes that could be analyzed for some classes
2. **Conservative** - Excludes genes with valid data for 4/5 classes due to 1 missing class

## Is This Too Stringent?

**Potentially yes**, especially for the Mis2 class.

### Alternative approach:
Analyze each variant class with its own optimal gene set:
- PTV: 9,000 genes
- Mis2: 5,700 genes
- Mis1: 8,800 genes
- Mis0: 9,000 genes
- Syn: 9,000 genes

This would:
- Maximize power for each individual analysis
- Still allow comparisons (just note that gene sets differ slightly)
- Be especially helpful for Mis2, which loses the most genes

### Current approach is defensible because:
- The Kim 2024 cohort is underpowered anyway (680 probands)
- Adding 1,500 genes won't fundamentally change the conclusions
- Consistency across analyses simplifies interpretation
- The 7,518-gene set still includes adequate background genes

## Recommendation

**For small cohorts like Kim 2024 (680 probands):** Keep the consensus filtering. The sample size is the limiting factor, not the gene count.

**For large cohorts (>5,000 probands):** Consider class-specific gene sets to maximize power, especially for Mis2 which has many genes with zero mutation rates.

## Technical Note

The function `get_genetic_data()` (from `R/io.R` or `R/secondary_analysis_functions.R`) likely filters genes by:
1. Removing genes with NA mutation rates for that class
2. Removing genes with zero mutation rates for that class
3. Removing genes with NA expected counts
4. Possibly other quality filters

So the "consensus" is really an intersection of genes that pass these filters for ALL variant classes simultaneously.
