# burdenEM RVAS

This repository implements a method for estimating the distribution of burden effect sizes from rare variant association statistics, and calculating downstream estimators. A typical workflow is to (1) fit mixture models for some set of traits, and (2) calculate a downstream estimator. Step (1) requires (a) a `studies.tsv` file listing the traits/datasets to be analyzed, (b) optionally a `genes.tsv` file containing gene annotations and LD-corrected burden scores, and (c) one or more `sumstats.txt.bgz` files containing the variant-level summary statistics for each dataset. It computes a `.rds` file containing the fitted model for each row of `studies.tsv`. Step (2) requires (a) the same `studies.tsv` file, and (b) the `.rds` files.

## Directory Structure

Following R package conventions:
- `R/` contains all R source code organized into subdirectories
- `src/` contains compiled code (C++/C/Fortran)

```
burdenEM/
├── R/
│   ├── cli/                    # Command-line entry points
│   │   ├── fit_models.R        # Step 1: Fit BurdenEM models
│   │   └── analyze.R           # Step 2: Compute downstream estimators
│   │
│   ├── io/                     # Data I/O
│   │   ├── io.R                # Load variant files, data processing
│   │   └── variant_to_gene.R   # Aggregate variants to gene level
│   │
│   ├── estimation/             # Core model fitting
│   │   ├── model.R             # Model initialization, posterior calculations
│   │   ├── em.R                # EM algorithm implementation
│   │   ├── likelihoods.R       # Likelihood functions
│   │   ├── fit_burdenem.R      # Grid-based model fitting
│   │   ├── run_pipeline.R      # Main workflow orchestration
│   │   ├── intercept.R         # Intercept/calibration estimation
│   │   ├── gene_features.R     # Feature binning
│   │   └── binary_traits.R     # Binary trait likelihood and overdispersion
│   │
│   ├── analysis/               # Step 2: Post-fit analysis
│   │   ├── heritability.R      # Heritability estimation
│   │   ├── heritability_core.R # Core heritability functions
│   │   ├── polygenicity.R      # Polygenicity metrics
│   │   ├── distribution.R      # Effect size distribution
│   │   ├── replication.R       # P-value replication
│   │   ├── effect_replication.R # Effect size replication
│   │   ├── power.R             # Power calculations
│   │   └── qqplot.R            # QQ plot utilities
│   │
│   ├── simulations/            # Simulation framework
│   │   ├── sim_cli.R           # Simulation analysis CLI
│   │   ├── sim_heritability.R  # True heritability from simulations
│   │   ├── sim_polygenicity.R  # True polygenicity from simulations
│   │   ├── sim_distribution.R  # True distribution from simulations
│   │   └── sim_calibration.R   # Calibration metrics
│   │
│   └── plotting/               # Visualization
│       ├── plotting_functions.R
│       └── plot_*.R            # Various plotting scripts
│
├── src/
│   └── EM.cpp                  # Compiled C++ code for EM algorithm
│
├── Python/                     # Python utilities (LD computation)
├── data/                       # Input data
├── tables/                     # Output tables
├── figures/                    # Output figures
└── fitted_models/              # Fitted model .rds files
```

## Quick Start

```bash
# Step 1: Fit models
Rscript R/cli/fit_models.R data/my.studies.tsv --annotation pLoF --genes_file data/my.genes.txt

# Step 2: Compute heritability
Rscript R/cli/analyze.R heritability data/my.studies.tsv --annotation pLoF

# Or compute effect replication
Rscript R/cli/analyze.R effect_replication data/my.studies.tsv --primary_dataset genebass -a pLoF
```

## CLI Reference

### Step 1: Model Fitting (`R/cli/fit_models.R`)

Fits BurdenEM mixture models to variant-level summary statistics.

**Positional arguments:**
- `studies_file`: Path to TSV file defining studies (must contain: `identifier`, `dataset`, `sumstats_filename_pattern`, `model_filename`)

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-a`, `--annotation` | Functional annotation (e.g., 'pLoF', 'missense'). Use 'all' for simulations. | `all` |
| `-g`, `--genes_file` | Path to genes file template. Supports `<ANNOTATION>`, `<LOWER>`, `<UPPER>`, `<DATASET>` placeholders. | None (optional) |
| `-f`, `--feature_col_name` | Column name for gene features (e.g., 'oe_lof') | None |
| `-b`, `--num_feature_bins` | Number of bins for feature column | `5` |
| `-c`, `--num_positive_components` | Number of positive mixture components | `10` |
| `--burdenem_grid_size` | Grid points per component | `4` |
| `-i`, `--num_iter` | Maximum EM iterations | `5000` |
| `--per_allele_effects` | Use per-allele instead of per-gene effects | `FALSE` |
| `-m`, `--binary_trait_model_type` | Model for binary traits: betabinom, binom, nbinom, pois | `betabinom` |
| `--correct_for_ld` | Apply LD correction | `FALSE` |
| `--frequency_range` | AF range as "min,max" | `0,0.001` |
| `--intercept_frequency_bin_edges` | AF bin edges for intercept estimation | `0,1e-5,1e-4,1e-3,1` |
| `--optimizer` | Optimization method: 'EM' or 'mixsqp' | `EM` |
| `--no_parallel` | Run sequentially | `FALSE` |
| `--skip_existing` | Skip if output exists | `FALSE` |
| `-n`, `--name` | Override output directory | None |
| `-v`, `--verbose` | Print detailed output | `FALSE` |

### Step 2: Analysis (`R/cli/analyze.R`)

Computes downstream estimators from fitted models.

**Subcommands:**
- `heritability`: Heritability components per study
- `polygenicity`: Various polygenicity metrics
- `distribution`: Genes required to explain fractions of heritability
- `replication`: P-value replication (requires `--primary_dataset`)
- `effect_replication`: Effect size replication (requires `--primary_dataset`, continuous traits only)

**Options:**
| Flag | Description | Default |
|------|-------------|---------|
| `-a`, `--annotation` | Variant annotation | `pLoF` |
| `-t`, `--trait_type` | Filter: 'binary', 'continuous', or 'all' | `all` |
| `--primary_dataset` | Primary dataset for replication analyses | None |
| `--pvalue_threshold` | P-value threshold for replication | `5e-8` |
| `--no_parallel` | Run sequentially | `FALSE` |
| `-n`, `--name` | Override output directory | None |
| `-v`, `--verbose` | Print detailed output | `FALSE` |

## Running Simulations

Simulations use the same CLI as real data analysis to ensure any bugs are caught where ground truth is known.

```bash
# Fit models on simulation data
Rscript R/cli/fit_models.R sims/my.studies.tsv --frequency_range "0,1"

# Analyze simulation results (computes true vs estimated metrics)
Rscript R/simulations/sim_cli.R heritability sims/my.studies.tsv -a all
Rscript R/simulations/sim_cli.R polygenicity sims/my.studies.tsv -a all
Rscript R/simulations/sim_cli.R distribution sims/my.studies.tsv -a all
Rscript R/simulations/sim_cli.R calibration sims/my.studies.tsv -a all
```

## File Formats

### `studies.tsv`
One row per study-phenotype pair:
| Column | Description |
|--------|-------------|
| `identifier` | Unique study identifier |
| `trait_type` | `continuous` or `binary` |
| `dataset` | Dataset name |
| `description` | Human-readable description |
| `abbreviation` | Short name for grouping replicates |
| `sumstats_filename_pattern` | Path to sumstats, may contain `<ANNOTATION>` |
| `model_filename` | Path for output model, may contain `<ANNOTATION>` |

### `sumstats.txt.bgz`
Variant-level summary statistics:
| Column | Required | Description |
|--------|----------|-------------|
| `gene` | Yes | Gene symbol or ENSG ID |
| `AF` | Yes | Allele frequency |
| `beta` | Yes | Effect size estimate |
| `variant_variance` | Yes | Variance contribution |
| `trait_type` | Yes | `continuous` or `binary` |
| `N` | Binary | Sample size |
| `AC_cases` | Binary | Allele count in cases |
| `prevalence` | Binary | Disease prevalence |

### `genes.txt`
Gene-level annotations:
| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `gene_id` | ENSG ID |
| `burden_score` | LD-corrected burden score |
| Additional columns for features (e.g., `lof.oe`) |

## Model Object

The fitted model (`.rds` file) is a list containing:

| Component | Description |
|-----------|-------------|
| `df` | Gene-level data frame with `likelihood` matrix and `features` |
| `grid_effects` | Effect size grid points |
| `components` | Component distributions over grid (matrix: components x grid) |
| `component_endpoints` | Component boundary values |
| `null_index` | Index of the null component |
| `delta` | Mixture weights (matrix: features x components) |
| `h2_function` | Function(beta, row) returning heritability contribution |
| `pval_function` | Function(row) returning one-tailed p-value |
| `get_power_function` | Function(pval_threshold) returning power function |
| `information` | Information matrices for standard errors |

## API

### Computing Posterior Means

Most estimators are posterior means of functions of beta:

```r
# Heritability with standard error
result <- posterior_expectation_with_se(model, model$h2_function)
# result$mean: posterior mean per feature
# result$se: standard error per feature

# Custom function
my_fn <- function(beta, row) beta^2 * row$burden_score
result <- posterior_expectation_with_se(model, my_fn)

# Gene-specific posterior means (no SE)
gene_means <- posterior_expectation2(model, model$h2_function)
```
