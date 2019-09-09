clear;
close all
%% Collapse Gibbs sampling for Dirichelt process gaussian mixture model

dim = 2;
k = 4;
n = 400;
alpha = 3;

randomData = true;
if randomData
    [X, label, true_model] = GaussianMixtureRand(k, n, alpha);
else
    manaulGMtest
end
% display data points
figure(1);clf;
plotClusters(X,label);

% run split-merge MCMC process
[priorModel, labels] = initDirichletProcessModel(X);
[GauMixModels, y] = runRestrictGibSM(X, labels, priorModel, alpha);
figure;
plotClusters(X,y);

% run matlab built-in GMM with known number of components
GMModel = fitgmdist(X',k);
idx = cluster(GMModel,X');
figure;
plotClusters(X,idx);

