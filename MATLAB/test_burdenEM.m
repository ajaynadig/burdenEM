addpath(genpath('.'))
nn = 10000;
normal_model = 1;
% features = repmat(eye(2),nn/2,1);
features = rand(nn,1); 
features = [features, 1-features];
% features = ones(nn,1);
model_params = [-4 1 0 1 4]/10;
coefs = [.01 .04 .9 .04 .01; 0 0 1 0 0];

% Mixture component assignments
weights = features * coefs;
weights = weights ./ sum(weights,2);
cpt_param = zeros(nn,1);
for ii = 1:nn
    cpt_param(ii) = randsample(model_params,1,true,weights(ii,:));
end

% Effect sizes
beta = rand(nn,1);
beta = beta .* cpt_param;

if normal_model
    % Gaussian noise in the effect estimates
    noise_var = .01*rand(nn,1);
    beta_hat = beta + sqrt(noise_var) .* randn(nn,1);

    % Estimation
    coefs_est = burdenEM(model_params,...
        'effect_estimate', beta_hat, 'effect_se', sqrt(noise_var),...
        'features', features, 'model_type', 'uniform');
else
    % Poisson-distributed counts
    case_rate = 400 * rand(nn,1);
    lambda = case_rate .* exp(beta);
    case_count = poissrnd(lambda);

    % Estimation
    coefs_est = burdenEM(model_params,...
        'case_count', case_count, 'case_rate', case_rate,...
        'features', features, 'model_type', 'uniform');
end
disp(coefs_est)