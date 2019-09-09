function [priorModel, Labels] = initDirichletProcessModel(X)
%   Parameters of Normal-Inverse-Wishart
%   see: https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
[dim, n] = size(X);

kappa_0 = 1;
mu_0 = mean(X, 2);
nu_0 = dim;
S_0 = eye(dim);

% suppose all the points are in one cluster
Labels = ones(1, size(X, 2));
setIdx_prior = false(1, size(X, 2));
priorModel = DirichletProcessModel(mu_0, kappa_0, nu_0, S_0, X, setIdx_prior);

end

