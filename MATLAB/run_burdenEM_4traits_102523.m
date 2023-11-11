addpath(genpath('.'))
clear
figure

traits = {'Height','BMI','HbA1c','Neuroticism'};
titles = {'A Height','B BMI','C HbA1c','D Neuroticism'};
trait_numbers = [36 13 29 10];

% Each trait has a different intercept for each MAF bin (0-1e-5, 1e-5-1e-4,
% 1e-4-1e-3). These numbers were pulled from BHR supplementary tablbes.
intercept_mafbins = [3.2, 3.3, 3.4;
    2.89 2.89 2.96;
    3 3 3;
    3.31, 3.34, 3.41] * 1e-6;

% Constraint bins
baselinefile = '~/Dropbox/GitHub/bhr/reference_files/ms_baseline_oe5.txt';
T_baseline=readtable(baselinefile);

% used for extracting phenotype codes
files_group1 = dir('/Volumes/T7/data/BHR/bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds_*.txt');

for jj=1:4
    ii = trait_numbers(jj);
    
    % phenotype code for trait ii
    filename = files_group1(ii).name;
    phenocode = sscanf(filename,'bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds_%d*');
    
    % three files for each trait, corresponding to MAF bins
    files = {['bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0_high1e-05_group1.ms.munged.Rds_',num2str(phenocode),'NA.txt'];
        ['bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low1e-05_high0.0001_group2.ms.munged.Rds_',num2str(phenocode),'NA.txt'];
        ['bhr_ms_gene_ss_400k_final_withnullburden_pLoF_low0.0001_high0.001_group3.ms.munged.Rds_',num2str(phenocode),'NA.txt']
        };
    
    % load data for each MAF bin, adding just one extra column containing
    % the BHR intercept for that trait-MAFbin pair
    T = table();
    for bin = 1:length(files)
        T_mafbin = readtable(['/Volumes/T7/data/BHR/',files{bin}]);
        T_mafbin.intercept(:) = intercept_mafbins(jj,bin);
        T=[T; T_mafbin];
    end

    % group rows of table by gene (as opposed to by MAF bin). Rows of T
    % with T.gene equal to genes{k} are those with indices==k.
    [genes, ~, indices] = unique(T.gene);

    noGenes = length(genes);
    nn = mean(T.N);
    mean_intercept_nn = mean(T.intercept) * nn;

    % Merge between baseline annotations and gene list
    [~,~,idx] = intersect(genes, T_baseline.Var2, 'stable');
    assert(noGenes == length(idx)); % assert all genes are found

    % noGenes x 5 matrix where each gene is assigned to a constraint
    % quintile
    oe_bins = table2array(T_baseline(idx,3:end));
    oe_bins = [oe_bins, 1-sum(oe_bins,2)];
    
    % X'*y/n where X is the n x m genotype matrix and y is the n x 1
    % phenotype vector
    Xty = T.variant_variance.*T.beta;

    % (X*W)'*y/n where W is a variants by genes matrix of burden weights. Ie,
    % sum up X'*y for all the variants within each gene
    XWty = accumarray(indices, Xty);
    
    % Burden scores
    burdenScore = accumarray(indices, T.variant_variance);
    
    % Convert variant-wise intercept into a gene-wise intercept. Happy to
    % explain this formula next time we meet.
    intercept_vector = accumarray(indices, T.intercept .* T.variant_variance) ./ burdenScore;
    
    % Burden effect size estimates in per-allele units
    gamma_perAllele = XWty ./ burdenScore;

    % Burden effect size estimatese in per-sd units
    gamma_perSD = gamma_perAllele .* sqrt(burdenScore);

    % Endpoints of the uniform mixture components
    effectVar = mean(gamma_perSD.^2);
    unif_endpoints_perSD = effectVar * 4.^(-2:5);
    unif_endpoints_perSD = [-unif_endpoints_perSD, 0, unif_endpoints_perSD];
    
    % Run burdenEM with per-SD effect sizes
    sampling_variance = sqrt(intercept_vector);
    weights_matrix_perSD = burdenEM(unif_endpoints_perSD,...
        'effect_estimate', gamma_perSD, 'effect_se', sampling_variance,...
        'features', oe_bins, 'model_type', 'uniform');
    weights_matrix_perSD = weights_matrix_perSD';
    
    % Average mixture weights across genes (assumes equal number of genes
    % in each bin)
    weights_perSD_agg = sum(weights_matrix_perSD/5,2);
    
    % Mean burden score in each bin
    for bin = 1:5
        meanBurdenScore(bin) = mean(burdenScore(oe_bins(:,bin)==1));
    end

    % Run burdenEM with per-allele effect sizes
    unif_endpoints_perAllele = unif_endpoints_perSD / sqrt(mean(burdenScore));
    sampling_variance = intercept_vector ./ burdenScore;
    weights_matrix_perAllele = burdenEM(unif_endpoints_perAllele,...
        'effect_estimate', gamma_perAllele, 'effect_se', sqrt(sampling_variance),...
        'features', oe_bins, 'model_type', 'uniform');
    weights_matrix_perAllele = weights_matrix_perAllele';
    weights_perAllele_agg = sum(weights_matrix_perAllele/5,2);

    % Organize results into a matrix with one row for each mixture
    % component. First column = endpoints; second = weights; third =
    % heritabilty per component. 1/3 comes from E(x.^2) where x~unif(0,1).
    resultsTablePerSD = sortrows([unif_endpoints_perSD' weights_perSD_agg ...
        noGenes/3*weights_perSD_agg.*unif_endpoints_perSD'.^2]);
    resultsTablePerAllele = sortrows([unif_endpoints_perAllele' weights_perAllele_agg ...
        noGenes/3*weights_perAllele_agg.*unif_endpoints_perAllele'.^2 * mean(burdenScore)]);

    % Print out heritability estimates for the two ways of running
    % burdenEM, which should roughly agree
    h2_perSD(jj) = sum(resultsTablePerSD(:,3));
    h2_perAllele(jj) = sum(resultsTablePerAllele(:,3));
    prop_h2_positive(jj) = sum(resultsTablePerSD(resultsTablePerSD(:,1)>0,3))/h2_perSD(jj);
    disp([h2_perSD(jj) h2_perAllele(jj) prop_h2_positive(jj)])

    %% Probably no need to replicate what follows
    
    % Make QQ plot
    % Simulate from the per-SD mixture distribution
    samples = []; noSamples = 1e6;
    for kk = 1:length(unif_endpoints_perSD)
        samples = [samples; unif_endpoints_perSD(kk) * rand(floor(noSamples*weights_perSD_agg(kk)),1)];
    end

    % Number of effects to explain 50% h2
    noSamples = length(samples);
    cumh2 = cumsum(sort(samples.^2,'descend'));
    first = find(cumh2 >= cumh2(end)/2, 1);
    polygenicity(jj) = first/noSamples*noGenes;
    x = sort(samples.^2,'descend');
    median_poly(jj) = h2_perSD(jj) / x(first);
    entropy_poly(jj) = exp(-mean(samples.^2 .* log(max(eps,samples.^2)/h2_perSD(jj)))/mean(samples.^2));
    effective_poly(jj) = 3*noGenes*mean(samples.^2)^2 / mean(samples.^4);

    % Add noise
    samples = samples + randn(noSamples,1) * sqrt(mean_intercept_nn/nn);
    yy = quantile(samples * sqrt(nn / mean_intercept_nn),((1:noGenes)/noGenes)-0.5/noGenes);

    subplot(2,5,jj+1); hold on

    Z = gamma_perSD * sqrt(nn / mean_intercept_nn);
    maxZ = max(abs(Z));

    expected=norminv((1:noGenes)/noGenes);
    xylim = 4;
    h1 = plot([-xylim,xylim],[-xylim,xylim],'color',[.5 .5 .5]);
    h2 = scatter(expected,sort(Z),8,[ 0    0.4470    0.7410], 'filled');
    h3 = plot(expected, yy,'color',[ 0.9290    0.6940    0.1250]);

    title(titles{jj})
    xlabel('Expected Z')
    ylabel('Observed Z')
    if jj == 1
        legend([h1 h3 h2], {'Null','Predicted','Observed'});legend boxoff
    end
    xlim([-xylim xylim]);
    ylim([-15 15])

    N_trait(jj) = nn;
    nn_array = [1 2 4 8 16];
    threshold_array = 2.^(-4:5);
    for nval = 1:length(nn_array)
        [total_ntpr(nval,jj),total_h2gwas(nval,jj),total_power(nval,jj),power_pos(nval,jj),power_neg(nval,jj)] = NTPR(unif_endpoints_perSD,weights_perSD_agg',nn/mean_intercept_nn*nn_array(nval),.05/noGenes);
    end
    for nval = 1:length(threshold_array)
        [threshold_ntpr(nval,jj),~,threshold_power(nval,jj)] = NTPR(unif_endpoints_perSD,weights_perSD_agg',nn/mean_intercept_nn,.05/noGenes*threshold_array(nval));
    end

    disp(noGenes * [total_power,power_pos,power_neg])

    % Simulate from the per-allele mixture distribution
    samples = []; samples_h2 = []; noSamples = 1e6;
    for bin = 1:5
        for kk = 1:length(unif_endpoints_perAllele)
            newSamples = unif_endpoints_perAllele(kk) * rand(floor(noSamples*weights_matrix_perAllele(kk,bin)),1);
            samples = [samples; newSamples];
            newSamples_h2 = newSamples.^2 .* randsample(burdenScore(oe_bins(:,bin)==1),length(newSamples),true);
            samples_h2 = [samples_h2; newSamples_h2];
        end
    end
    noLargeEffectsPositive(jj) = noGenes * mean(samples > .5);
    noLargeEffectsNegative(jj) = noGenes * mean(samples < -.5);

    % Effect size of significant genes
    % sig_genes = gamma_perSD.^2 > sig_threshold;
    
end

colors = colororder();

subplot(2,4,6)
bar(1,prop_h2_positive');
ylabel('Prop. h^2_{burden}')
ylim([0 1]);xlim([.5 1.5])
set(gca,'XTick',1 + (-1.5:1.5)*.18,'XTickLabel',traits)

title('F h^2 explained by trait-increasing genes')
box off

subplot(2,4,5)
bar(1,effective_poly');
set(gca,'XTick',1:4,'XTickLabel',traits)
ylabel('No. genes')
xlim([.5 1.5])
set(gca,'XTick',1 + (-1.5:1.5)*.18,'XTickLabel',traits)
title('E Effective polygenicity')
box off

subplot(2,4,7)

plot(N_trait.*nn_array', noGenes*total_power);
set(gca,'Xscale','log','YScale','log')
ylabel('No. genes')
xlabel('Sample size')
xlim([2e5 1e7])
title('G Exome-wide significant genes')
box off

subplot(2,4,8)
plot(N_trait.*nn_array', total_h2gwas);
set(gca,'Xscale','log')
ylabel('Prop. h^2_{burden}')
xlabel('Sample size')
ylim([0 1])
xlim([2e5 1e7])
title('H h^2_{burden} of EWS genes')
box off

