function [labels, Models, logprob] = innerGibbs(X, S_idx, Models, t)
%INNERGIB Summary of this function goes here
%   Detailed explanation goes here
n = size(S_idx, 2);
logprob = 0;
logprob_all = [];
for iter = 1:t
    randpidx = randperm(n);
    for i = 1:size(randpidx, 2)
        cidx = S_idx(randpidx(i));
        k = getCompLabel(Models, cidx);
        Models{k} = Models{k}.DelOne(X, cidx);
        N = getNumsOfSample(Models);
        N_ck = N(k);
        N_ci = N(1);
        N_cj = N(2);
        
        term1 = N_ck .* exp(Models{k}.GetPosteriorPred(X, cidx));
        term2 = N_ci .* exp(Models{1}.GetPosteriorPred(X, cidx));
        term3 = N_cj .* exp(Models{2}.GetPosteriorPred(X, cidx));
        
        p_stay = term1 / (term2+term3);
        p_change = 1 - p_stay;
%         k_sample = discreteRnd([p_stay, p_change]);
        
        if rand < p_stay
            Models{k} = Models{k}.AddOne(X, cidx);
            logprob = logprob + log(p_stay);
            logprob_all = [logprob_all, log(p_stay)];
        else
            Models{3-k} = Models{3-k}.AddOne(X, cidx);
            logprob = logprob + log(p_change);
            logprob_all = [logprob_all, log(p_change)];
        end
        
    end
end

w = getNumsOfSample(Models)/n;

labels = zeros(1, size(X, 2));
for i = 1:size(Models, 2)
    labels(Models{i}.setIdx) = i;
end

end

