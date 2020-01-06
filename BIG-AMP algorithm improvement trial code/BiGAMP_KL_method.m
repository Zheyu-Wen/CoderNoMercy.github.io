function [estFin, optFin, estHist, state] = ...
    BiGAMP_KL_method(gX, gA, gOut, problem, opt, state,~,stop_num)
%%
% BiGAMP: Bilinear Generalized Approximate Message Passing
%
% The BiG-AMP algorithm is intended for the estimation of
% random matrices A and X observed through the Markov Chain
%
%   X,A -> Z = A*X -> Y,
%
% where the components of X and A are independent and the mapping Z -> Y is
% separable. X is NxL, A is MxN, and Z,Y are consequently MxL.
%
% INPUTS:
% -------
% gX:  An input estimator derived from the EstimIn class
%    based on the input distribution p_x_{nl}(x_nl).
% gA:  An input estimator derived from the EstimIn class
%    based on the input distribution p_a_{mn}(a_mn).
% gOut:  An output estimator derived from the EstimOut class
%    based on the output distribution p_{Y|Z}(y_ml|z_ml).
% problem: An objet of the class BiGAMPProblem specifying the problem
%   setup, including the matrix dimensions and observation locations
% opt (optional):  A set of options of the class BiGAMPOpt.
% state (optional): A structure containing all the values needed to warm
%   start BiG-AMP
%
% OUTPUTS:
% --------
% estFin: Structure containing final BiG-AMP outputs
% optFin: The BiGAMPOpt object used
% estHist: Structure containing per iteration metrics about the run
% state: The values of all parameters required to warm start the algorithm


%% Setup why define these setting
Y = gOut.y;
nuw = gOut.wvar;
nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;           % minimum number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
adaptStep = opt.adaptStep;      % adaptive step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
stepWindow = opt.stepWindow;    % step size check window size
diagnostics = opt.diagnostics;  % Save diagnostic information
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
pvarStep = opt.pvarStep;        % incldue step size in pvar/zvar
uniformVariance =...
    opt.uniformVariance;        % use scalar variances
compVal = adaptStep;            % only compute cost function for adaptive
zvarToPvarMax = opt.zvarToPvarMax;  % maximum zvar/pvar ratio

%Determine requested outputs
saveEM = opt.saveEM;
saveHist = (nargout >= 3);

it_hist = [];
Aerror_hist = [];
Xerror_hist = [];
Avar_hist = [];
Xvar_hist = [];
Zerror_hist = [];
%Get problem dimensions

sparseMode = opt.sparseMode;
xvarMin = opt.xvarMin;
Avarmin = opt.AvarMin;

%Preallocate storage for estHist if user requires it
if (saveHist)
    estHist.errZ = zeros(nit,1);
    estHist.errX = zeros(nit,1);
    estHist.errA = zeros(nit,1);
    estHist.val = zeros(nit,1);
    estHist.step = zeros(nit,1);
    estHist.pass = false(nit,1);
    estHist.timing = zeros(nit,1);
    if diagnostics %Only create these if diagnostics are requested
        estHist.pvarMin = zeros(nit,1);
        estHist.pvarMax = zeros(nit,1);
        estHist.pvarMean = zeros(nit,1);
        estHist.zvarMin = zeros(nit,1);
        estHist.zvarMax = zeros(nit,1);
        estHist.zvarMean = zeros(nit,1);
        estHist.AvarMin = zeros(nit,1);
        estHist.AvarMax = zeros(nit,1);
        estHist.AvarMean = zeros(nit,1);
        estHist.xvarMin = zeros(nit,1);
        estHist.xvarMax = zeros(nit,1);
        estHist.xvarMean = zeros(nit,1);
        estHist.svarMin = zeros(nit,1);
        estHist.svarMax = zeros(nit,1);
        estHist.svarMean = zeros(nit,1);
        estHist.qvarMin = zeros(nit,1);
        estHist.qvarMax = zeros(nit,1);
        estHist.qvarMean = zeros(nit,1);
        estHist.rvarMin = zeros(nit,1);
        estHist.rvarMax = zeros(nit,1);
        estHist.rvarMean = zeros(nit,1);
        estHist.normA = zeros(nit,1);
        estHist.normX = zeros(nit,1);
        estHist.normZ = zeros(nit,1);
        estHist.zError = zeros(nit,1);  
        estHist.testxvar = zeros(nit,1);
        estHist.testAvar = zeros(nit,1);
    end
end


%% Initialization

%Check for provided state
if nargin < 6
    state = [];
end
state = [];
if isempty(state) %if no state is provided

    xhat = opt.xhat0;
    valIn = -inf;
    xvar = opt.xvar0;
    Ahat = opt.Ahat0;
    Avar = opt.Avar0;    
    
    %Placeholder initializations- values are not used
    xhatBar = 0;
    AhatBar = 0;
    pvarOpt = 0;
    zvarOpt = 0;
    CoefficientOpt = 0;
    Coefficient = 0;
    CoefficientNew = 0;
    %Address warm starting of shat0 what mean for warm starting?
    
    %Init valOpt empty
    valOpt = [];
    
    %Set pvar_mean to unity initially
    pvar_mean = 1;
    temp = 0;
    tempOpt = 0;
else %Use the provided state information
    
    %A variables
    Ahat = state.Ahat;
    Avar = state.Avar;
    AhatBar = state.AhatBar;
    AhatOpt = state.AhatOpt;
    AhatBarOpt = state.AhatBarOpt;
    
    %X Variables
    xhat = state.xhat;
    xvar = state.xvar;
    xhatBar = state.xhatBar;
    xhatOpt = state.xhatOpt;
    xhatBarOpt = state.xhatBarOpt;
    
    
    %Cost stuff
    valIn = state.valIn;
    
    %Variance momentum terms
    pvarOpt = state.pvarOpt;
    zvarOpt = state.zvarOpt;
    
    %Step
    step = state.step;
    
    %Old cost values
    valOpt = state.valOpt;
    
    %Set pvar_mean
    pvar_mean = state.pvar_mean;
    
%     Coefficient = state.Coefficient;
%     CoefficientNew = state.Coefficient;
%     CoefficientOpt = state.CoefficientOpt;
%     temp = state.temp;
%     tempOpt = state.tempOpt;
%     tempNew = state.tempNew;
end

%Specify minimum variances
pvarMin = opt.pvarMin;

%Placeholder initializations
rhat = 0;
rvar = 0;
qhat = 0;
qvar = 0;

%Cost init
val = zeros(nit,1);
zhatOpt = 0;


%% Iterations

%Start timing first iteration
tstart = tic;

%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;


%Handle first step
step1 = 1;
% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    
    if ~uniformVariance
        %Precompute squares quantities for use, these change on every
        %iteration
        Ahat2 = abs(Ahat).^2;
        xhat2 = abs(xhat).^2;
        
        %Compute zvar   %%% corresponding to the paper    (A5) How to
        %observe the difference?
        if ~sparseMode
            zvar = Avar*xhat2 + Ahat2*xvar;
        end
        
        %Compute pvar
        if ~sparseMode
            pvar = zvar + Avar*xvar;
        end
    
    end
    
    %Include pvar step   ?
    if pvarStep
        pvar = step1*pvar + (1-step1)*pvarOpt;
        zvar = step1*zvar + (1-step1)*zvarOpt;
    end
    
    %Update zhat
    zhat = Ahat * xhat;
    % Continued output step
    phat = zhat- Coefficient.*(zvar);%...
    %+ abs(shat).^2 .* ( (Ahat .* Avar) * (xhat .* xvar)  ); %%% ???????

    % Compute log likelihood at the output and add it the total negative   
    % K-L distance at the input.
    if (compVal)
        valOut = sum(sum(gOut.logLike(zhat,pvar))); %%% gOut function ?           
        val(it) = valOut + valIn;
    end
    
    % Determine if candidate passed
    if ~isempty(valOpt)
        
        %Check against worst value in last stepWindow good steps
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        
        %Check the step
        pass = (val(it) > min(valOpt(startInd:stopInd))) ||...
            ~adaptStep || (step <= stepMin);
        
    else
        pass = true;
    end
    
    %Save the step size and pass result if history requested
    if saveHist
        estHist.step(it) = step;
        estHist.pass(it) = pass;
    end
    

    % If pass, set the optimal values and compute a new target shat and
    % snew.
    if (pass)
        
        %Slightly inrease step size after pass if using adaptive steps

        step = stepIncr*step;
        % Set new optimal values
        xhatBarOpt = xhatBar;
        xhatOpt = xhat;
        AhatBarOpt = AhatBar;
        AhatOpt = Ahat;
        pvarOpt = pvar;
        zvarOpt = zvar;
        tempOpt = temp;
        CoefficientOpt = Coefficient;

        %Bound pvar
        pvar = max(pvar, pvarMin);
        
        %We keep a record of only the succesful step valOpt values
        valOpt = [valOpt val(it)]; %#ok<AGROW>

        %Compute mean of pvar
        
        % Output nonlinear step
        [~,zvar0] = gOut.estim(phat,pvar);
    
        
        %Update the shat quantities
        
        CoefficientNew = 1 ./ (pvar+1*nuw).*(Y-phat);
        tempNew = 1 ./ (pvar+1*nuw);%.*(1-min(zvar0./pvar,zvarToPvarMax));
        %Enforce step size bounds
        step = min([max([step stepMin]) stepMax]);
        
    else
        
        %Check on failure count
        failCount = failCount + 1;
        % Decrease step size
        step = max(stepMin, stepDecr*step);      
        %Check for minimum step size
        if step < stepTol
            stop = true;
        end
    end
    
    
    
    % Save results
    if (saveHist)
        
        %Record timing information
        if it > 1
            estHist.timing(it) = estHist.timing(it-1) + toc(tstart);
        else
            estHist.timing(it) = toc(tstart);
        end
        
        %Compute the Z error only if needed
        if ~sparseMode
            estHist.errZ(it) = opt.error_function(zhat);
        end
        %estHist.errZ(it) = sum(sum((problem.Z - zhat).^2))/sum(sum((problem.Z).^2));
        estHist.errA(it) = opt.error_functionA(Ahat);
        estHist.errX(it) = opt.error_functionX(xhat);
        estHist.val(it) = val(it);
    end

    % Check for convergence if step was succesful
    if pass
        if any(isnan(zhat(:))) || any(isinf(zhat(:)))
            stop = true;
        else
            testVal = norm(zhat(:) - zhatOpt(:)) / norm(zhat(:));
            if (it > 1) && ...
                    (testVal < tol)
                stop = true;
            end
        end
        
        %Set other optimal values- not actually used by iterations
        AvarOpt = Avar;
        xvarOpt = xvar;
        zhatOpt = zhat;
        
        %Save EM variables if requested
        if saveEM
            rhatFinal = rhat + gX.update_mean;
            rvarFinal = pvar_mean*rvar;
            qhatFinal = qhat + gA.update_mean;
            qvarFinal = pvar_mean*qvar;
            zvarFinal = zvar0;
            pvarFinal = pvar;
        end
    end
    
    %Start timing next iteration
    tstart = tic;
    
    % Create new candidate shat
    if it > 1 || ~isempty(state)
        step1 = step;
    end
    Coefficient = (1-step1)*CoefficientOpt + step1*CoefficientNew;
    temp = (1 - step1)*tempOpt + step1*tempNew;

    xhatBar = (1-step1)*xhatBarOpt + step1*xhatOpt;
    AhatBar = (1-step1)*AhatBarOpt + step1*AhatOpt;
    gA.update_mean = AhatBar;
    gX.update_mean = xhatBar;
    
    rvar = 1./((abs(AhatBar).^2)'*(temp));
    rvar(rvar > opt.varThresh) = opt.varThresh;
    rhat = rvar.*((AhatBar'-(mean(mean(Avar))*nuw/(mean(mean((Y-phat).^2))*mean(mean(Avar))+mean(mean(AhatBar.^2))*mean(mean(pvar)))))*Coefficient);
    rvar = max(rvar, xvarMin);
    
    % Input linear step for A
    qvar = 1./((temp)*(abs(xhatBar).^2)');
    qvar(qvar > opt.varThresh) = opt.varThresh;
    qhat = qvar.*(Coefficient*(xhatBar'-(mean(mean(xvar))*mean(mean(nuw))/(mean(mean((Y-phat).^2))*mean(mean(xvar))+mean(mean(xhatBar.^2))*mean(mean(pvar)))))); 
    qvar = max(qvar,Avarmin);
    
    % Input nonlinear step
    [xhat,xvar,valInX] = gX.estim(rhat, rvar*pvar_mean);
    [Ahat,Avar,valInA] = gA.estim(qhat, qvar*pvar_mean);
    %xhat = step1*xhat + (1 - step1)*xhatOpt;
    %Ahat = step1*Ahat + (1 - step1)*AhatOpt;
    [Aerror,Xerror] = AXMSE_Compute(Ahat,xhat,problem);
    %% state evolution
    c = 1*abs(Xerror/Aerror);
    cA = problem.M*problem.N^2*problem.L*c;
    cB = problem.M*problem.N*norm(xhat,'fro')^2+problem.N*problem.L*c^2*norm(Ahat,'fro')^2;
    cC = problem.M*problem.L*nuw + norm(Ahat*xhat,'fro')^2 - norm(Y,'fro')^2;
    sigmaA = 1*(-cB+sqrt(cB^2-4*cA*cC))/(2*cA);
    sigmaX = sigmaA*c;
    Aerror_hist = [Aerror_hist,Aerror];
    Xerror_hist = [Xerror_hist,Xerror];
    %Avar = Avar + sigmaA;
    %xvar = xvar + sigmaX;
    Avar_hist = [Avar_hist,mean(mean(Avar))];
    Xvar_hist = [Xvar_hist,mean(mean(xvar))];
    %disp(['Aerror: ' num2str(Aerror) ', Xerror: ' num2str(Xerror) ', Avar: ',num2str(mean(mean(Avar))) ', xvar: ',num2str(mean(mean(xvar)))...
    %    ', sigmaA: ' num2str(sigmaA) ', sigmaX: ' num2str(sigmaX)]);
    %% DeSVD computation
%     [Ahat1, xhat1] = AX_DeSVD_Compute(Ahat,xhat,problem);
%     Zerror = norm(Ahat1*xhat1-problem.Z,'fro')^2/numel(problem.Z);
%     Zerror_hist = [Zerror_hist,Zerror];
%     %Update valIn
     valIn = sum( valInX(:) ) + sum ( valInA(:) );
%     it_hist = [it_hist, it];
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end
    
end
if stop_num == 2
    opt.Avar_estimate = [opt.Avar_estimate,Avar_hist];
    opt.Avar_true = [opt.Avar_true, Aerror_hist];
    opt.Xvar_estimate = [opt.Xvar_estimate,Xvar_hist];
    opt.Xvar_true = [opt.Xvar_true, Xerror_hist];
%     save('C:\Users\DELL\Desktop\Paper use code\_result\Aerror_hist','Aerror_hist');
%     save('C:\Users\DELL\Desktop\Paper use code\_result\Xerror_hist','Xerror_hist');
%     save('C:\Users\DELL\Desktop\Paper use code\_result\it_hist','it_hist');
%     save('C:\Users\DELL\Desktop\Paper use code\_result\Avar_hist','Avar_hist');
%     save('C:\Users\DELL\Desktop\Paper use code\_result\Xvar_hist','Xvar_hist');
end
%save('C:\Users\DELL\Desktop\Paper use code\_result\KL_Zerror_hist','Zerror_hist');
%% Save the final values

%Save the options object that was used
optFin = opt;

%Estimates of the two matrix factors
estFin.xhat = xhatOpt;
estFin.xvar = xvarOpt;
estFin.Ahat = AhatOpt;
estFin.Avar = AvarOpt;

%Save values useful for EM learning
if saveEM
    estFin.rhat = rhatFinal;
    estFin.rvar = rvarFinal;
    estFin.qhat = qhatFinal;
    estFin.qvar = qvarFinal;
    estFin.zvar = zvarFinal;
    estFin.phat = phat;
    estFin.pvar = pvarFinal;
end

%% Cleanup estHist

%Trim the outputs if early termination occurred
if saveHist && (it < nit)
    estHist.errZ = estHist.errZ(1:it);
    estHist.errA = estHist.errA(1:it);
    estHist.errX = estHist.errX(1:it);
    estHist.val = estHist.val(1:it);
    estHist.step = estHist.step(1:it);
    estHist.pass = estHist.pass(1:it);
    estHist.timing = estHist.timing(1:it);
    if diagnostics
        estHist.pvarMin = estHist.pvarMin(1:it);
        estHist.pvarMax = estHist.pvarMax(1:it);
        estHist.pvarMean = estHist.pvarMean(1:it);
        estHist.zvarMin = estHist.zvarMin(1:it);
        estHist.zvarMax = estHist.zvarMax(1:it);
        estHist.zvarMean = estHist.zvarMean(1:it);
        estHist.AvarMin = estHist.AvarMin(1:it);
        estHist.AvarMax = estHist.AvarMax(1:it);
        estHist.AvarMean = estHist.AvarMean(1:it);
        estHist.xvarMin = estHist.xvarMin(1:it);
        estHist.xvarMax = estHist.xvarMax(1:it);
        estHist.xvarMean = estHist.xvarMean(1:it);
        estHist.svarMin = estHist.svarMin(1:it);
        estHist.svarMax = estHist.svarMax(1:it);
        estHist.svarMean = estHist.svarMean(1:it);
        estHist.qvarMin = estHist.qvarMin(1:it);
        estHist.qvarMax = estHist.qvarMax(1:it);
        estHist.qvarMean = estHist.qvarMean(1:it);
        estHist.rvarMin = estHist.rvarMin(1:it);
        estHist.rvarMax = estHist.rvarMax(1:it);
        estHist.rvarMean = estHist.rvarMean(1:it);
        estHist.normA = estHist.normA(1:it);
        estHist.normX = estHist.normX(1:it);
        estHist.normZ = estHist.normZ(1:it);
        estHist.zError = estHist.zError(1:it);
        
        estHist.testxvar = estHist.testxvar(1:it);
        estHist.testAvar = estHist.testAvar(1:it);
    end
end

saveState = 0;
if saveState
    
    %A variables
    state.Ahat = Ahat;
    state.Avar = Avar;
    state.AhatBar = AhatBar;
    state.AhatOpt = AhatOpt;
    state.AhatBarOpt = AhatBarOpt;
    
    %X Variables
    state.xhat = xhat;
    state.xvar = xvar;
    state.xhatBar = xhatBar;
    state.xhatBarOpt = xhatBarOpt;
    state.xhatOpt = xhatOpt;

    
    %Cost stuff
    state.valIn = valIn;
    
    %Variance momentum terms
    state.pvarOpt = pvarOpt;
    state.zvarOpt = zvarOpt;
    
    %Step
    state.step = step;
    
    %Old cost values
    state.valOpt = valOpt;
    
    %Problem dimensions
    %state.N = N;
    
    %pvar_mean
    state.pvar_mean = pvar_mean;
%     state.Coefficient = Coefficient;
%     state.CoefficientNew = Coefficient;
%     state.CoefficientOpt = CoefficientOpt;
%     
%     state.temp = temp;
%     state.tempOpt = tempOpt;
%     state.tempNew = tempNew;
end






