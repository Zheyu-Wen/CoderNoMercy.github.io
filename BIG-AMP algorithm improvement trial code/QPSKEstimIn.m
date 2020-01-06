classdef QPSKEstimIn < EstimIn
    % QPKSEstimIn:  QPSK Signal input estimation function
    % Construct: QPSKEstimIn(prob1,probm1,probi,probmi)
    
    properties
        %Probability of 1,-1,1i,-1i
        prob1 = 0.25;           
        probm1 = 0.25;      
        probi = 0.25;
        probmi = 0.25;
        mean0 = 0;         % Prior mean
        var0 = 1;          % Prior variance
    end
   
 
    methods
        % Constructor
        function obj = QPSKEstimIn(prob1,probm1,probi,probmi)
            obj = obj@EstimIn;
            obj.prob1 = 0.25;           
            obj.probm1 = 0.25;
            obj.probi = 0.25;
            obj.probmi = 0.25;
            if nargin ~= 0 % Allow nargin == 0 syntax
                sum_prob = prob1 + probi + probm1 + probmi;
                obj.prob1 = prob1 ./ sum_prob;
                obj.probm1 = probm1 ./ sum_prob;
                obj.probi = probi ./ sum_prob;
                obj.probmi = probmi ./ sum_prob;
            end
            obj.mean0 = (obj.prob1 - obj.probm1 + 1i * (obj.probi - obj.probmi))/4;
            obj.var0 = obj.prob1 .* abs( 1 - obj.mean0).^2 ...
                        + obj.probm1 .* abs( -1 - obj.mean0).^2 ...
                        + obj.probi .* abs( 1i - obj.mean0).^2 ...
                        + obj.probmi .* abs( -1i - obj.mean0).^2;
        end

        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(obj)
            mean0 = obj.mean0;
            var0  = obj.var0;
            valInit = 0;
        end
       
        % Size
        function [nx,ncol] = size(obj)
            [nx,ncol] = size(obj.mean0);
        end

        % QPSK estimation function
        % Provides the mean and variance of a variable x = CN(uhat0,uvar0)
        % from an observation rhat = x + w, w = CN(0,rvar)
        function [xhat, xvar, val] = estim(obj, rhat, rvar)
            % Get prior
            prior_prob1 = obj.prob1;           
            prior_probm1 = obj.probm1;      
            prior_probi = obj.probi;
            prior_probmi = obj.probmi;
            % Compute the conditional distribution
            con_prob1 = exp( -abs( 1 - rhat).^2 ./ rvar);
            con_probm1 = exp( -abs( -1 - rhat).^2 ./ rvar);
            con_probi = exp( -abs( 1i - rhat).^2 ./ rvar);
            con_probmi = exp( -abs( -1i - rhat).^2 ./ rvar);
            sum_con = con_prob1 + con_probm1 + con_probi + con_probmi;
            sum_con((sum_con < 10^(-50))) = 10^(-50);
            con_prob1 = con_prob1 ./ sum_con;
            con_probm1 = con_probm1 ./ sum_con;
            con_probi = con_probi ./ sum_con;
            con_probmi = con_probmi ./ sum_con;
            
            
            % Compute posterior mean and variance
            post_prob1 = prior_prob1 .* con_prob1;
            post_probm1 = prior_probm1 .* con_probm1;
            post_probi = prior_probi .* con_probi;
            post_probmi = prior_probmi .* con_probmi;
            
            sum_post = post_prob1 + post_probm1 + post_probi + post_probmi;
            sum_post((sum_post < 10^(-50))) = 10^(-50);
            post_prob1 = post_prob1 ./ sum_post;
            post_probm1 = post_probm1 ./ sum_post;
            post_probi = post_probi ./ sum_post;
            post_probmi = post_probmi ./ sum_post;
            
            xhat = post_prob1 - post_probm1 + 1i * (post_probi - post_probmi);
            xvar = post_prob1 .* abs( 1 - xhat).^2 ...
                        + post_probm1 .* abs( -1 - xhat).^2 ...
                        + post_probi .* abs( 1i - xhat).^2 ...
                        + post_probmi .* abs( -1i - xhat).^2;
            
%             if (nargout >= 3)
%                % Compute the negative KL divergence
%                %   klDivNeg = \sum_i \int p(x|r)*\log( p(x) / p(x|r) )dx
%                xvar_over_uvar0 = rvar./(uvar0+rvar);
%                val =  (log(xvar_over_uvar0) + (1-xvar_over_uvar0) ...
%                         - abs(xhat-uhat0).^2./uvar0 );
%             end
            if (nargout >= 3)
                % Compute the negative KL divergence
                %   klDivNeg = \int p(x|r)*\log( p(x) / p(x|r) )dx
                klDivNeg1 = max(min(con_prob1 .* log( prior_prob1 ./ con_prob1),50),-50);
                klDivNegm1 = max(min(con_probm1 .* log( prior_probm1 ./ con_probm1),50),-50);
                klDivNegi = max(min(con_probi .* log( prior_probi ./ con_probi),50),-50);
                klDivNegmi = max(min(con_probmi .* log( prior_probmi ./ con_probmi),50),-50);
                val = klDivNeg1 + klDivNegm1 +  klDivNegi +  klDivNegmi;
            end
          end
         
        % Generate random samples
        function x = genRand(obj, outSize)
            if isscalar(outSize)
                randNum = rand();
            else
                randNum = rand(outSize);
            end
            isValue1 = (randNum < obj.prob1);
            isValuem1 = (~isValue1) & (randNum - obj.prob1 < obj.probm1);
            isValuei = (~isValue1) & (~isValuem1) & (randNum - obj.prob1 - obj.probmi < obj.probi);
            isValuemi = (~isValue1) & (~isValuem1) & (~isValuei);
            x = isValue1 - isValuem1 + 1i * (isValuei - isValuemi);
        end
        
        % Computes the likelihood p(rhat) for rhat = x + v, v = CN(0,rvar)
        % this function not used
        function py = plikey(obj,rhat,rvar)
            py = exp(-1./((obj.var0+rvar)).*abs(rhat-obj.mean0).^2);
            py = py./ (pi*(obj.var0+rvar));
        end
        
        % Computes the log-likelihood, log p(rhat), for rhat = x + v, where 
        % x = CN(obj.mean0, obj.var0) and v = CN(0,rvar)
        function logpy = loglikey(obj, rhat, rvar)
          py1 = obj.prob1 .* exp( -abs(rhat - 1).^2 ./ rvar);
          pym1 = obj.probm1 .* exp( -abs(rhat + 1).^2 ./ rvar);
          pyi = obj.probi .* exp( -abs(rhat - 1i).^2 ./ rvar);
          pymi = obj.probmi .* exp( -abs(rhat + 1i).^2 ./ rvar);
          logpy =  - log(pi) - (log(rvar)/2) + log(py1 + pym1 + pymi + pyi);
          logpy = max(min(logpy,50),-50);
%             logpy = -( log(pi) + log(obj.var0 + rvar) + ...
%                 (abs(rhat - obj.mean0).^2) ./ (obj.var0 + rvar) );
        end    
    end
    
end

