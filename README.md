# burdenEM RVAS

## `model` object
The new pipeline centers on the new `model` object. This object contains likelihoods
for each gene evaluated at each point in a prespecified grid of effect sizes. After
model fitting, it also contains the mixture weights for each component,
as well ass an information matrix for these mixturue weights.

It is a list with the following components:

- `df`: A data frame containing gene-level data, including the computed likelihood matrix as a list column `likelihood` (genes x grid points), the features for each gene as a list column `features`, and usually additional columns. Currently `features` is assumed to be a one-hot encoding of a gene partition.
- `grid_effects`: A numeric vector of the grid points at which the likelihood is evaluated.
- `components`: A matrix approximating each mixture component as a discrete distribution over the grid points. Each row is a component, and each column is a grid point.
- `component_endpoints`: A numeric vector of the component endpoints.
- `null_index`: One component is required to be 'null', and this is the index of that component.
- `delta`: A matrix of the mixture weights for each feature (row) and component (column).
- `h2_function`: A function which operates on an effect size and a row of the data frame `df`, returning the heritability explained by that gene if it has that effect size.

## Fitting the model
The model fitting procedure requires two inputs: a variant-level summary statistics file and a gene-level file containing any gene annotation data and (optionally) LD-corrected burden scores. There can be multiple variant-level files, which will be concatenated. It writes a .rds file containing the fitted model. The gene-level file contains features that are binned.

The same function (`run_burdenEM_rvas` in `luke/main.R`) runs either on continuous or binary traits. For the latter it fits a Poisson model, which seems to produce inflated estimates for synonymous variants (I think this is expected). For a continuous trait it runs either a per-allele or a per-sd effect size model.

The goal is that all downstream analyses will be agnostic to the kind of model used. The key to achieving this is the `h2_function` stored inside of the model object. Any downstream analysis of heritability components should use this function to avoid needing to know the details of the model that was fitted.  

## Computing heritability components
Call the function `estimate_heritability_components` with the model object. It returns a list with the total heritability, the heritability explained by each component, and the heritability explained by positive and negative effects, each with a point estimate and standard error.

This function itself calls `posterior_expectation_with_se` with the `h2_function` stored in the model object.

## Computing arbitrary posterior means
Most estimators of interest can be expressed as the posterior mean of some function of beta. To compute
such an estimator, specify a function `my_function(beta)` or `my_function(beta, row)` where `row` is a row of the data frame `df`. Then, call `result<-posterior_expectation_with_se(model, my_function)`, and this returns a list with the posterior expectation for each feature (`result$mean`) and its standard error (`result$se`). A limitation is that you cannot currently obtain the covariance of two posterior means, which would be useful, for example, to get the standard error of their ratio.

To compute the heritability and its standard error, call `result<-posterior_expectation_with_se(model, model$h2_function)`. The idea is that this can be completely general across different likelihoods (gaussian, poisson, etc.) and also across different definitions of the effect size (per-allele vs per-s.d.).

To compute gene-specific posterior means, instead of calling `posterior_expectation_with_se`, call `posterior_expectation2(model, my_function)`. Similarly, `my_function` can take either one argument or two.

## Partial list of to do items
- Integrate negative binomial model for binary traits
- Combine trio and rvas workflows, which I think should be much easier with these changes
- Integrate new posterior expectation code with various estimators, such as power calculations
- Run across traits; I've tried height, bmi, LDL / pLoF, synonymous
- Maybe get AoU ld-corrected burden scores; LDL looks like it has inflation in synonymous
- 