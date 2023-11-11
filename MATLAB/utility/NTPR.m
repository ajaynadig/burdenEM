function [total_ntpr,total_h2gwas,total_power,total_power_pos,total_power_neg] = NTPR(bb,pp,nn,alpha)
%NTPR computes the not-by-chance true positive rate for a mixture of uniform
%   distributions with parameters [0,b] for b in bb, or [b,0] if b<0, and
%   mixture weights p in pp, at sample sizes nn and significance thresholds
%   alpha. It also computes h2gwas, the fraction of h2 that will be
%   explained by significant SNPs/genes at that significance threshold, and
%   h2power, the expected fraction of SNPs/genes that will be significant.

assert(isrow(nn));
assert(isrow(alpha));
assert(isrow(bb));
assert(isrow(pp));

chisq_threshold = chi2inv(1-alpha,1)';


% Probability that marker with effect size x will be significant at
% threshold a with sample size n
powerfn=@(x,a,n)normcdf(sqrt(n).*x-sqrt(a)) + normcdf(-sqrt(n).*x-sqrt(a));

ntprfn = @(x,a,n)(normcdf(sqrt(n).*abs(x)-sqrt(a)) - normcdf(-sqrt(n).*abs(x)-sqrt(a))) ./ powerfn(x,a,n);

% Probability that 
samesignfn=@(x,a,n)normcdf(sqrt(n)*x-sqrt(a));

[total_power,total_ntpr,total_h2gwas,total_power_pos,total_power_neg] = deal(zeros(length(alpha),length(nn)));

for ii=1:length(pp)
    xx=(0:.001:1) * bb(ii);
    pow = powerfn(reshape(xx,1,1,length(xx)),chisq_threshold,nn);
    cpt_power = mean(pow,3);
    
    ntpr = ntprfn(reshape(xx,1,1,length(xx)),chisq_threshold,nn);
    cpt_ntpr = mean(pow.*ntpr,3);
    
    cpt_h2gwas = mean(pow.*reshape(xx.^2,1,1,length(xx)),3);

    total_power = total_power + pp(ii) * cpt_power;
    total_ntpr = total_ntpr + pp(ii) * cpt_ntpr;
    total_h2gwas = total_h2gwas + pp(ii) * cpt_h2gwas;

    if bb(ii) > 0
        total_power_pos = total_power_pos + pp(ii) * cpt_power;
    elseif bb(ii) < 0
        total_power_neg = total_power_neg + pp(ii) * cpt_power;
    end
end

total_ntpr = total_ntpr ./ total_power;
total_h2gwas = total_h2gwas / sum(pp.*bb.^2/3);

end