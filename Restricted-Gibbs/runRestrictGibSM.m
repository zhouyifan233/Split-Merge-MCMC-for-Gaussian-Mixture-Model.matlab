function [GauMixModels, y] = runRestrictGibSM(X, labels, priorModel, alpha)
%RUNGIBBSGM Summary of this function goes here
%   Detailed explanation goes here
Models = [];
rejectcount=0;
n = size(X,2);
maxIter = 1000;
maxlabelval = max(labels(:))+1;
for iter = 1:maxIter
    
    % Pick two distinct observations
    RandIdx = randperm(n);
    i_idx = RandIdx(1);
    j_idx = RandIdx(2);
    
    % Get labels of i, j
    c_i = labels(i_idx);
    c_j = labels(j_idx);
    
    % prepare launch set
    if c_i == c_j
        c_launch_i = maxlabelval;
        c_launch_j = c_j;
        c_i_idx = find(labels == c_i);
        S_with_ij_idx = c_i_idx;S_idx = S_with_ij_idx;
        S_idx(S_idx == i_idx) = [];S_idx(S_idx == j_idx) = [];
        unirandidx = randi(2, 1, size(S_idx, 2));
        S_i_with_i_idx = [S_idx(unirandidx==1), i_idx];
        S_j_with_j_idx = [S_idx(unirandidx==2), j_idx];
    else
        c_launch_i = c_i;
        c_launch_j = c_j;
        c_i_idx = find(labels == c_i);c_j_idx = find(labels == c_j);
        S_with_ij_idx = [c_i_idx, c_j_idx];S_idx = S_with_ij_idx;
        S_idx(S_idx == i_idx) = [];S_idx(S_idx == j_idx) = [];
        unirandidx = randi(2, 1, size(S_idx, 2));
        S_i_with_i_idx = [S_idx(unirandidx==1), i_idx];
        S_j_with_j_idx = [S_idx(unirandidx==2), j_idx];
    end
    
    Models{1} = priorModel.copy();
    Models{1} = Models{1}.AddMore(X, S_i_with_i_idx);
    Models{2} = priorModel.copy();
    Models{2} = Models{2}.AddMore(X, S_j_with_j_idx);

    [inner_labels, Models, logprob] = innerGibbs(X, S_idx, Models, 20);
    
    if c_i == c_j
        % do split
        c_i_split = c_launch_i;
        c_j_split = c_launch_j;
        [inner_labels, Models, logprob] = innerGibbs(X, S_idx, Models, 1);
        c_i_split_idx = find(inner_labels==1);
        c_j_split_idx = find(inner_labels==2);
        inner_labels(c_i_split_idx) = c_i_split;
        inner_labels(c_j_split_idx) = c_j_split;
        labels_split = labels;
        labels_split(inner_labels~=0) = inner_labels(inner_labels~=0);
        
        N_c_i_split = size(c_i_split_idx, 2);
        N_c_j_split = size(c_j_split_idx, 2);
        N_c_i = N_c_i_split + N_c_j_split;          % N_c_i = N_c_j as c_i == c_j
        
        % calculate Prior
        Prior = log(alpha) + logfactorial([1:N_c_i_split-1, 1:N_c_j_split-1], 1:N_c_i-1);
        
        % calculate Likelihood
        Model1 = priorModel.copy();
        term1 = 0;
        for j = 1:size(c_i_split_idx, 2)
            if j == 1
                term1 = term1 + Model1.GetPrior(X, c_i_split_idx(j));
                Model1 = Model1.AddOne(X, c_i_split_idx(j));
            else
                term1 = term1 + Model1.GetPosteriorPred(X, c_i_split_idx(j));
                Model1 = Model1.AddOne(X, c_i_split_idx(j));
            end
        end
        Model2 = priorModel.copy();
        term2 = 0;
        for j = 1:size(c_j_split_idx, 2)
            if j == 1
                term2 = term2 + Model2.GetPrior(X, c_j_split_idx(j));
                Model2 = Model2.AddOne(X, c_j_split_idx(j));
            else
                term2 = term2 + Model2.GetPosteriorPred(X, c_j_split_idx(j));
                Model2 = Model2.AddOne(X, c_j_split_idx(j));
            end
        end
        Model3 = priorModel.copy();
        term3 = 0;
        for j = 1:size(S_with_ij_idx, 2)
            if j == 1
                term3 = term3 + Model3.GetPrior(X, S_with_ij_idx(j));
                Model3 = Model3.AddOne(X, S_with_ij_idx(j));
            else
                term3 = term3 + Model3.GetPosteriorPred(X, S_with_ij_idx(j));
                Model3 = Model3.AddOne(X, S_with_ij_idx(j));
            end
        end
        Likelihood = term1+term2-term3;
        
        % calculat Proposal
        Proposal = -logprob;
        
        accept_prob = min([exp(Prior+Likelihood+Proposal), 1]);
        
        if rand < accept_prob
            labels = labels_split;
            maxlabelval = maxlabelval + 1;
            rejectcount = 0;
            disp([num2str(iter), '  split accept...']);
        else
            rejectcount = rejectcount + 1;
            disp([num2str(iter), '  split reject...']);
        end
        
    else
        % do merge
        c_i_merge = c_j;
        c_j_merge = c_j;
        labels_merge = labels;
        
        labels_merge(c_i_idx) = c_i_merge;
        labels_merge(c_j_idx) = c_j_merge;
        
        N_c_i = size(c_i_idx, 2);
        N_c_j = size(c_j_idx, 2);
        N_c_i_merge = N_c_i + N_c_j;
        
        % calculate Prior
        Prior = -log(alpha) + logfactorial(1:N_c_i_merge-1, [1:N_c_i-1, 1:N_c_j-1]);
        
        % calculate Likelihood
        Model1 = priorModel.copy();
        term1 = 0;
        for j = 1:size(S_with_ij_idx, 2)
            if j == 1
                term1 = term1 + Model1.GetPrior(X, S_with_ij_idx(j));
                Model1 = Model1.AddOne(X, S_with_ij_idx(j));
            else
                term1 = term1 + Model1.GetPosteriorPred(X, S_with_ij_idx(j));
                Model1 = Model1.AddOne(X, S_with_ij_idx(j));
            end
        end
        Model2 = priorModel.copy();
        term2 = 0;
        for j = 1:size(c_i_idx, 2)
            if j == 1
                term2 = term2 + Model2.GetPrior(X, c_i_idx(j));
                Model2 = Model2.AddOne(X, c_i_idx(j));
            else
                term2 = term2 + Model2.GetPosteriorPred(X, c_i_idx(j));
                Model2 = Model2.AddOne(X, c_i_idx(j));
            end
        end
        Model3 = priorModel.copy();
        term3 = 0;
        for j = 1:size(c_j_idx, 2)
            if j == 1
                term3 = term3 + Model3.GetPrior(X, c_j_idx(j));
                Model3 = Model3.AddOne(X, c_j_idx(j));
            else
                term3 = term3 + Model3.GetPosteriorPred(X, c_j_idx(j));
                Model3 = Model3.AddOne(X, c_j_idx(j));
            end
        end
        Likelihood = term1-term2-term3;
        
        % calculat Proposal
        target_labels_local = zeros(1, size(X, 2));
        target_labels_local(c_i_idx) = 1;
        target_labels_local(c_j_idx) = 2;
        
        logprob = innerGibbs_reverse(X, S_idx, target_labels_local, Models);
        Proposal = logprob;
        
        accept_prob = min([exp(Prior+Likelihood+Proposal), 1]);
        
        if rand < accept_prob
            labels = labels_merge;
            rejectcount = 0;
            disp([num2str(iter), '  merge accept...']);
        else
            rejectcount = rejectcount + 1;
            disp([num2str(iter), '  merge reject...']);
        end
    end
    
    if rejectcount == inf
       break; 
    end
    figure(101);clf;plotClusters(X, labels);pause(0.1);
end

uni_labels = unique(labels);
y = zeros(1, size(labels, 2));
for i = 1:size(uni_labels, 2)
    y(labels==uni_labels(i)) = i;
end

GauMixModels = [];
for i = 1:size(uni_labels, 2)
    obsidx = find(y==i);
    GauMixModels{i}.mu = mean(X(:, obsidx), 2);
    GauMixModels{i}.sigma = cov(X(:, obsidx));
    GauMixModels{i}.w = size(obsidx, 2)/size(X, 2);
end

end
