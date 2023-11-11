function [coefs, mixture_params] = burdenEM(mixture_params,varargin)
%burdenEM fits a mixture model for the distribution of
%some test statistics using EM. It can fit a mixture of uniform
%distributions, or a mixture of normal distributions. It can accomodate
%either Poisson-distributed or normally distributed data (in the Poisson
%case, it only fits a mixture of uniforms, not a mixture of normals). 
% 
%   gamma_hat: effect-size estimates
%   sigmasq_grid: grid of effect-size variance parameters
%   effect_variance: sampling variance of each effect size estimate

%  Input handling
p = inputParser;

% Parameters of the mixture model: either variance parameters if model is
% normal, or edges of a uniform distribution otherwise
addRequired(p, 'mixture_params', @isrow);

% Which type of mixture model: 'uniform' or 'normal'
addParameter(p, 'model_type', 'uniform', @ischar);

% Effect-size estimates
addParameter(p, 'effect_estimate', [], @iscolumn);

% SE of effect_estimate
addParameter(p, 'effect_se', [], @iscolumn);

% Allele count in cases
addParameter(p, 'case_count', [], @iscolumn);

% Expected allele count in cases
addParameter(p, 'case_rate', [], @iscolumn);

% P-value (is this one-tailed or two-tailed?)
addParameter(p, 'p_value', [], @iscolumn);

% Features: should specify a convex combination for each test
addParameter(p, 'features', [], @ismatrix);

% Iterations for the EM algorithm
addParameter(p, 'EM_iterations', 100, @isscalar);

parse(p, mixture_params,varargin{:});

% turns variables named p.Results.x into just x
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end
clear p

no_cpts = length(mixture_params);

if ~isempty(effect_estimate)
    assert(isempty(case_count), 'Specify effect_estimate or case_count but not both')
    if isempty(effect_se)
        assert(~isempty(p_value,'must specify either effect_se or p_value with effect_estimate'))
        assert(all(size(p_value) == size(effect_estimate)))
        z_score = sqrt(chi2inv(p_value,1,'upper')) .* sign(effect_estimate);
        effect_se = effect_estimate ./ z_score;
    else
        assert(all(size(effect_se) == size(effect_estimate)))
    end
    data_type = 'normal';
    no_tests = length(effect_estimate);
else
    assert(~isempty(case_count), 'Specify either effect_estimate or case_count')
    assert(strcmpi(model_type, 'uniform'), 'Only uniform mixture model is supported with Poisson data')
    if isempty(case_rate)
        assert(~isempty(p_value,'must specify either case_rate or p_value with case_count'))
        assert(all(size(p_value) == size(case_count)))
        case_rate = poisson_rate(case_count, p_value);
    else
        assert(all(size(case_rate) == size(case_count)))
    end
    data_type = 'poisson';
    no_tests = length(case_count);
end

% Compute likelihood for each gene-component pair
if strcmpi(model_type,'uniform')
    likelihood = zeros(no_tests, no_cpts);
    mu_grid = 0.05:.1:1;
    if strcmpi(data_type,'normal')
        sigma = repmat(effect_se,1,length(mu_grid));
        for kk = 1:no_cpts
            mu = mu_grid * mixture_params(kk);
            likelihood(:,kk) = mean( normpdf(effect_estimate - mu,...
                zeros(size(sigma)),sigma), 2);
        end
    elseif strcmpi(data_type,'poisson')
        for kk = 1:no_cpts
            rate = case_rate .* exp(mu_grid * mixture_params(kk));
            likelihood(:,kk) = mean( poisspdf(case_count .* ones(1,length(mu_grid)),rate), 2);
        end
    end
elseif strcmpi(model_type,'normal')
    assert(all(mixture_params >= 0), 'Normal model should have nonnegative mixture_params, representing variances')
    sigma = sqrt(mixture_params + effect_se.^2);
    likelihood = normpdf(repmat(beta_hat,1,length(sigmasq_grid)),...
        zeros(no_tests, no_cpts), sigma);
else
    error('Options for model_type are uniform or normal')
end

if isempty(features)
    features = ones(no_tests,1);
end
assert(all(sum(features,2) == 1) && all(features(:)>=0),'Features should specify a convex combination for each SNP')
assert(size(features,1) == no_tests)

% EM algorithm
coefs = ones(size(features,2),no_cpts);
for rep = 1:EM_iterations
    weights = features * coefs;
    posteriors = weights .* likelihood;
    posteriors = posteriors ./ sum(posteriors,2);
    coefs = features \ posteriors;

end

end




