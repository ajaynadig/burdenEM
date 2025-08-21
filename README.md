# burdenEM RVAS

This repository implements a method for estimating the distribution of burden effect sizes from rare variant association statistics, and calculating downstream estimators. A typical workflow is to (1) fit mixture models for some set of traits, and (2) calculate a downstream estimator. Step (1) requires (a) a `studies.tsv` file listing the traits/datsets to be analyzed, (b) a `genes.tsv` file containing gene annotations and possibly LD-corrected burden scores, and (c) one or more `sumstats.txt.bgz` files containing the variant-level summary statistics for each dataset. It computes a `.rds` file containing the fitted model for each row of `studies.tsv`. Step (2) requires (a) the same `studies.tsv` file, and (b) the `.rds` files. 

For example:

```bash
# Step (1)
Rscript luke/main_cli.R data/my.studies.tsv --annotation pLoF --genes_file data/my.genes.txt

# Step (2)
Rscript luke/cli.R heritability data/my.studies.tsv --annotation pLoF

# or:
Rscript luke/cli.R effect_replication data/my.studies.tsv --primary_dataset genebass -a pLoF
```

## CLI

 The CLI has two commands, `luke/main_cli.R` and `luke/cli.R`, for steps 1 and 2 respectively.

### Step 1 options
Positional arguments:
- `studies_file`: Path to the TSV file defining studies to run (must contain columns: `identifier`, `dataset`, `sumstats_filename_pattern`). [required]

Optional flags:
- `-a`, `--annotation`: Single functional annotation category to process (e.g., 'pLoF', 'missense', 'synonymous'). [required]
- `-g`, `--genes_file`: Path to the genes file template (use `<ANNOTATION>`, `<LOWER>`, `<UPPER>`, and optionally `<DATASET>` placeholders). [required]
- `-f`, `--feature_col_name`: Optional: Name of the column in the genes file for gene features (e.g., 'oe_lof').
- `-b`, `--num_feature_bins`: Number of bins for feature column [default 5].
- `-i`, `--num_iter`: Number of EM iterations [default 5000].
- `-c`, `--num_positive_components`: Number of positive components for BurdenEM [default 10].
- `--burdenem_grid_size`: Grid size for BurdenEM [default 4].
- `--per_allele_effects`: Calculate per-allele effect sizes instead of per-gene. [default: FALSE]
- `-m`, `--binary_trait_model_type`: Model type for binary traits (betabinom, binom, nbinom, or pois). [default: betabinom]
- `--correct_for_ld`: Apply LD correction to burden scores and gamma. [default: FALSE]
- `--frequency_range`: Comma-separated min,max for allele frequency range (e.g., '0,0.001'). [default: "0,0.001"]
- `--intercept_frequency_bin_edges`: Comma-separated string of AF bin edges for intercept calculation (e.g., '0,1e-5,1e-4,1e-3'). [default: "0,1e-5,1e-4,1e-3"]
- `--no_parallel`: Disable parallel execution across studies (runs sequentially). [default: FALSE]
- `-v`, `--verbose`: Print extra details during execution. [default: FALSE]
- `-n`, `--name`: Override the output directory for model files. When provided, fitted models are written under `<studies_dir>/<name>/` while preserving the full relative `model_filename` path from the studies TSV.
- `--skip_existing`: Skip studies whose output .rds already exists. [default: FALSE]


### Step 2 subcommands
Each step 2 subcommand requires a `studies.tsv` file as its first positional argument. It outputs a single table, which will be located within a `tables` subdirectory of the directory containing the `studies.tsv` file. This table is named with the same prefix as the `studies.tsv` file and with the name of the subcommand.

    The subcommands are:
    - `heritability`: computes heritability components for each study
    - `polygenicity`: computes several polygenicity metrics
    - `distribution`: Computes the number of genes required to explain various fractions of heritability
    - `replication`: Calculates p-value replication metrics between a primary study and other studies for the same trait. Requires specifying a `--primary_dataset`.
    - `effect_replication`: Calculates effect-size replication metrics between a primary study and other studies for the same trait. Requires specifying a `--primary_dataset`. Not supported for binary traits.

### Step 2 options
- `studies_file` (required)
- `--annotation`, `-a` the variant annotation (default: pLoF)
- `--trait_type`, `-t` either 'binary', 'continuous', or 'all' (default: all)
- `--no_parallel` use sequential instead of parallel execution across studies
- `-n`, `--name` Override the output directory and table filename prefix. When provided, tables are written to `<studies_dir>/<name>/tables/` and table filenames are prefixed by `<name>`. For loading model files, the full relative `model_filename` from the studies TSV is preserved (after `<ANNOTATION>` substitution) and prefixed with `<studies_dir>/<name>/`.
- `--primary_dataset=<pd>` For `replication` and `effect_replication` subcommands, specifies the identifier of the dataset to be used as the primary study for comparisons.
- `--pvalue_threshold=<pval>` For `replication` only: two-tailed p-value threshold for significance in the primary study.
- `--help`, `-h`
- `--verbose`, `-v`
- `--version`

## File formats
- `.studies.tsv` with one line per study-phenotype pair
    - `identifier`
    - `trait_type`
    - `dataset`
    - `description`
    - `abbreviation`
    - `sumstats_filename_pattern` containing placeholder `<ANNOTATION>`
    - `model_filename` containing placeholder `<ANNOTATION>`
- `.model.rds`
- `.sumstats.txt.bgz`
    - `gene`, either a symbol or ENSG ID
    - `POS`
    - `CHR`
    - `beta`
    - `N`
    - `variant_variance`
    - `AF`
    - `trait_type`
    - `AC_cases` (if trait_type is binary)
    - `prevalence` (if trait_type is binary)
- `.genes.txt`
    - `gene`, the symbol
    - `gene_id`, the ENSG ID
    - `burden_score_corrected`
    - `burden_score_uncorrected`
    - any others, e.g. `lof.oe`

## API

### `model` object
The `model` object contains likelihoods
for each gene evaluated at each point in a prespecified grid of effect sizes. After
model fitting, it also contains the mixture weights for each component,
as well as an information matrix for these mixturue weights.

It is a list with the following components:

- `df`: A data frame containing gene-level data, including the computed likelihood matrix as a list column `likelihood` (genes x grid points), the features for each gene as a list column `features`, and usually additional columns. Currently `features` is assumed to be a one-hot encoding of a gene partition.
- `grid_effects`: A numeric vector of the grid points at which the likelihood is evaluated.
- `components`: A matrix approximating each mixture component as a discrete distribution over the grid points. Each row is a component, and each column is a grid point.
- `component_endpoints`: A numeric vector of the component endpoints.
- `null_index`: One component is required to be 'null', and this is the index of that component.
- `delta`: A matrix of the mixture weights for each feature (row) and component (column).
- `h2_function`: A function that operates on an effect size and a row of the data frame `df`, returning the heritability explained by that gene if it has that effect size.
- `pval_function`: A function that operates on a row of the data frame `df`, returning the one-tailed p-value for that row. A p-value close to zero indicates a postive effect, and a p-value close to one indicates a negative effect.
- `get_power_function`: A function that returns a power function for a given p-value threshold. This power function operates on a vector of effect sizes and a row of the data frame `df`, returning the one-tailed power for each effect size at the chosen p-value threshold.

### Fitting the model
The model fitting procedure requires two inputs: a variant-level summary statistics file and a gene-level file containing any gene annotation data and (optionally) LD-corrected burden scores. There can be multiple variant-level files, which will be concatenated. It writes a .rds file containing the fitted model. The gene-level file contains features that are binned.

The same function (`run_burdenEM_rvas` in `luke/main.R`) runs either on continuous or binary traits. For the latter it fits a negative binomial model with overdispersion parameter estimated using method of moments. For a continuous trait it runs either a per-allele or a per-sd effect size model.

### Computing arbitrary posterior means
Most estimators of interest can be expressed as the posterior mean of some function of beta. To compute
such an estimator, specify a function `my_function(beta)` or `my_function(beta, row)` where `row` is a row of the data frame `df`. Then, call `result<-posterior_expectation_with_se(model, my_function)`, and this returns a list with the posterior expectation for each feature (`result$mean`) and its standard error (`result$se`). A limitation is that you cannot currently obtain the covariance of two posterior means, which would be useful, for example, to get the standard error of their ratio.

To compute the heritability and its standard error, call `result<-posterior_expectation_with_se(model, model$h2_function)`. The idea is that this can be completely general across different likelihoods (gaussian, poisson, etc.) and also across different definitions of the effect size (per-allele vs per-s.d.).

To compute gene-specific posterior means, instead of calling `posterior_expectation_with_se`, call `posterior_expectation2(model, my_function)`. This does not return standard errors.