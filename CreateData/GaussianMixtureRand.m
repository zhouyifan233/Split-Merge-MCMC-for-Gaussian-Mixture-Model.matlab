function [X, labels, GMmodel] = GaussianMixtureRand(c, n, alpha, dim)
% Genarate samples form a Gaussian mixture model.
% Input:
%   c: number of components
%   n: number of data points
%   alpha: Parameter of Dirichlet distribution
%   dim: dimension of data (default 2)
% Output:
%   X: d x n data points
%   z: 1 x n labels
%   model: Gaussian mixture model

if nargin == 3
    dim = 2;
end

% parameters for generating covariances using inverse Wishart distribution
Sigma_iw = eye(dim);
DF_iw = dim+1;

% parameters for generating mu for Gaussian random numbers
mu_base = zeros(dim, 1);
beta_base = 9;

% generating labels for each data point using Dirichlet distribution
weights = dirichletRand(alpha*ones(1, c), 1);
% labels = discreteRnd(weights, n);
labels = randsample(1:c, n, true, weights);

X = zeros(dim ,n);
mu = zeros(dim, c);
Sigma = zeros(dim, dim, c);
for i = 1:c
    l = (labels==i);
    Sigma(:, :, i) = iwishrnd(Sigma_iw, DF_iw);
    mu(:, i) = mvnrnd(mu_base, beta_base*Sigma(:, :, i))';
    X(:, l) = mvnrnd(mu(:, i), Sigma(:, :, i), sum(l))';
end

GMmodel.mu = mu;
GMmodel.Sigma = Sigma;
GMmodel.w = weights;

