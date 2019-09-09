classdef DirichletProcessModel
    %see: http://rajarshd.github.io/talks/DPGMM_tutorial.pdf
    
    properties
        mu_0
        kappa_0
        nu_0
        S_0
        mu
        kappa
        nu
        S
        x_s
        S_part
        N
        setIdx
    end
    
    methods
        function obj = DirichletProcessModel(mu_0, kappa_0, nu_0, S_0, X, setIdx)
            %COMPMODEL Construct an instance of this class
            %   Detailed explanation goes here
            d = size(X, 1);
            obj.mu_0 = mu_0;
            obj.kappa_0 = kappa_0;
            obj.nu_0 = nu_0;
            obj.S_0 = S_0;
            obj.x_s = zeros(d, 1);
            obj.S_part = zeros(d, d);
            obj.setIdx = false(1, size(X, 2));
            obj.N = 0;
            
            if sum(setIdx) > 0
                obj.setIdx = setIdx;
                obj.N = sum(setIdx);
                obj.kappa = kappa_0 + obj.N;
                obj.nu = nu_0 + obj.N;
                obj.x_s = sum(X(:, setIdx), 2);
                obj.S_part = X(:, setIdx)*X(:, setIdx)';
                obj.mu = (obj.kappa_0*mu_0 + obj.x_s)/obj.kappa;
                obj.S = obj.S_0 + obj.S_part + obj.kappa_0.*mu_0*mu_0'-obj.kappa.*obj.mu*obj.mu';
            end
        end
        
        function obj = AddMore(obj, X, xidx)
            obj.N = obj.N + size(xidx, 2);
            obj.x_s = obj.x_s + sum(X(:, xidx), 2);
            obj.S_part = obj.S_part + X(:, xidx)*X(:, xidx)';
            obj.setIdx(xidx) = true;
            
            obj.kappa = obj.kappa_0 + obj.N;
            obj.nu = obj.nu_0 + obj.N;
            obj.mu = (obj.kappa_0*obj.mu_0 + obj.x_s)/obj.kappa;
            obj.S = obj.S_0 + obj.S_part + obj.kappa_0.*obj.mu_0*obj.mu_0' - obj.kappa.*obj.mu*obj.mu';
        end
        
        function obj = AddOne(obj, X, xidx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.setIdx(xidx) == false
                obj.N = obj.N + 1;
                obj.x_s = obj.x_s + X(:, xidx);
                obj.S_part = obj.S_part + X(:, xidx)*X(:, xidx)';
                obj.setIdx(xidx) = true;
                
                obj.kappa = obj.kappa_0 + obj.N;
                obj.nu = obj.nu_0 + obj.N;
                obj.mu = (obj.kappa_0*obj.mu_0 + obj.x_s)/obj.kappa;
                obj.S = obj.S_0 + obj.S_part + obj.kappa_0.*obj.mu_0*obj.mu_0' - obj.kappa.*obj.mu*obj.mu';
            else
                error('The element to be added is in the model...')
            end
        end
        
        function obj = DelOne(obj, X, xidx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if obj.setIdx(xidx) == true
                obj.N = obj.N - 1;
                obj.x_s = obj.x_s - X(:, xidx);
                obj.S_part = obj.S_part - X(:, xidx)*X(:, xidx)';
                obj.setIdx(xidx) = false;
                
                obj.kappa = obj.kappa_0 + obj.N;
                obj.nu = obj.nu_0 + obj.N;
                obj.mu = (obj.kappa_0*obj.mu_0 + obj.x_s)/obj.kappa;
                obj.S = obj.S_0 + obj.S_part + obj.kappa_0.*obj.mu_0*obj.mu_0' - obj.kappa.*obj.mu*obj.mu';
            else
                error('The element to be deleted is not in the model...')
            end
        end
        
        function PPpdf = GetPosteriorPred(obj, X, xidx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d = size(X, 1);
            v = obj.nu - d + 1;
            dist = X(:, xidx) - obj.mu;
            SIGMA = ((obj.kappa+1) ./ (obj.kappa*v)) .* obj.S;
            
            % see https://en.wikipedia.org/wiki/Multivariate_t-distribution
            term1 = gammaln((v+d)/2) - gammaln(v/2) - d/2*log(v) - d/2*log(pi) - log(det(SIGMA))/2;
            term2 = log(1 + (1/v)*(dist'*inv(SIGMA)*dist));
            
            PPpdf = term1 + -(v+d)/2 .* term2;
        end
        
        function Prpdf = GetPrior(obj, X, xidx)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            d = size(X, 1);
            v = obj.nu_0 - d + 1;
            dist = X(:, xidx) - obj.mu_0;
            SIGMA = ((obj.kappa_0+1) ./ (obj.kappa_0*v)) .* obj.S_0;
            
            % see https://en.wikipedia.org/wiki/Multivariate_t-distribution
            term1 = gammaln((v+d)/2) - gammaln(v/2) - d/2*log(v) - d/2*log(pi) - log(det(SIGMA))/2;
            term2 = log(1 + (1/v)*(dist'*inv(SIGMA)*dist));
            
            Prpdf = term1 + -(v+d)/2 .* term2;
        end
        
        function result = SearchIdx(obj, xidx)
            result = obj.setIdx(xidx);
        end
        
        function NumOfSample = getNumSample(obj)
            NumOfSample = obj.N;
        end
        
        function obj = copy(obj)
            
        end
        
    end
end

