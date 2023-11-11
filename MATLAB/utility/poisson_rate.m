function lambda = poisson_rate(x, p, steps)
%poisson_rate computes the rate, lambda, such that:
%   poisscdf(lambda,x,'upper') == p
%   steps (optional, default 100): number of Newton steps to take

if ~exist('steps')
    steps = 100;
end

lambda = x;
factx = factorial(x);
for step = 1:steps
    cdf = poisscdf(x, lambda);
    lambda = lambda - (cdf + p - 1) ./ (-exp(-lambda) .* lambda.^x ./ factx);
end
end