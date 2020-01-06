function results = DLCommun_LMMSE_ref_SPCA_QPSK(optIn)

% subspace estimation before DL, i.e., Y = UDV', Z = U1D1V1'. Then D1V1' is regarded as the input of DL algorithm. 

%% =====================================================================
%
%trial_DL: This function runs several algorithms on a sample
%instance of the dictionary learning problem. The function can be run with no
%arguments. See the nargin==0 block below for the required fields to call
%the function with optional parameters. The function returns a structure of
%algorithm results.
%This code can be used to produce the results for the noise free phase
%plane plots for dictionary learning in the BiG-AMP arXiv submission. 
%(This code conducts a single trial at 1 point in the phase plane.)

%Add paths
%setup_DL

%Test case
if nargin == 0
    clc
    
    %Handle random seed
    defaultStream = RandStream.getGlobalStream;
    if 1
        savedState = defaultStream.State;
        save random_state.mat savedState;
    else
        load random_state.mat %#ok<UNRCH>
    end
    defaultStream.State = savedState;
    
    %Flags to test BiG-AMP variants
    optIn.tryBigampEM = 1;
    
    %Flags to turn on comparison algorithms
    optIn.tryKsvd = 0;
    optIn.tryErspud = 0;
    optIn.trySpams = 0;

    
    %Problem dimensions
    optIn.M = 20; %coherence time
    optIn.N = optIn.M; %transmit antenna
    optIn.L = ceil(5*optIn.N*log(optIn.N)); %received antenna why set in this way?
    
    
    %Specify dictionary usage
    optIn.K = 5;
    
    %Specify maximum allowed trials
    optIn.maxTrials = 1;
    
    %SNR
    optIn.SNR = inf;
    
    %Specify coding mode (0 for OMP, 1 for TST)
    %The specified function is used to compute a coding error for the
    %resulting dictionary for each algorithm. TST may be significantly
    %slower.
    optIn.useTST = 0;
    
    %Precondition option (enable to use a preconditioner for EM-BiG-AMP)
    optIn.precondition = 0;
    
    
end

%% Problem Setup

%Turn algs on/off
tryBigampEM = optIn.tryBigampEM;
tryKsvd = optIn.tryKsvd;
tryErspud = optIn.tryErspud;
trySpams = optIn.trySpams;

%SNR
SNR = optIn.SNR;

%Coding
useTST = optIn.useTST;

%Precondition
precondition = optIn.precondition;

%Max trials
maxTrials = optIn.maxTrials;

%Define problem dimensions
M = optIn.M;
L = optIn.L;
N = optIn.N;
K = optIn.K;
K2 = optIn.K2;

% NoiseVar = optIn.NoiseVar;
% NoiseMean = optIn.NoiseMean;
% ScalarVar = optIn.ScalarVar;
% ScalarMean = optIn.ScalarMean;

%Set options
opt = BiGAMPOpt; %initialize the options object
opt.nit = optIn.nit; %limit iterations

opt.uniformVariance = optIn.uniformVariance;

%Set sizes
problem = BiGAMPProblem();
problem.M = N; %%%%% only for subspace estimation
problem.N = N;
problem.L = L;



%% Compute coefficient vectors

%Compute true vectors with exactly K non-zeroes in each %? does X is a
%sparse matrix, which is different from the paper where the H is sparse.
optGen.N = N; % size of dictionary
optGen.Nbar = N; %initial (estimated) size of dictionary
optGen.L = L; %From Spielman et al.
optGen.spar = K/N;  % non-zero rate of each packet
optGen.SNR = SNR;  %transmit SNR = 1*spar/noise_var in dB
optGen.alphabet = [1,-1,1i,-1i]'; %QPSK alphabet
optGen.verbose = true; %set to true to disp results
optGen.maxTrials = 5;
X = generateX(optGen);
% Kcol = ceil(K*L/N);
% for ll = 1:N
%     yada = randperm(L);
%     yada2 = zeros(1,L);
%     yada2(yada(1:Kcol)) = 1;
%     X(ll,:) = X(ll,:) .* yada2;
% end

%Draw randomly
%% 1-bit pilot % A is the transmitted signal 
A = zeros(M,N);
A1 = (randn(M-1,N)+randn(M-1,N)*1i)/sqrt(2); % singal
%A1 = sqrt(M-1)*A1*diag(1 ./ sqrt(abs(diag(A1'*A1))));
for ll = 1:M
    yada = randperm(N);
    yada2 = zeros(N,1);
    yada2(yada(1:K)) = 1;
    A(ll,:) = A(ll,:) .* yada2';
end
A(1,:) = 1;
A(2:M,:) = A1;
%scaledFactor = sqrt(M*N/real(trace(A*A')));
%A = scaledFactor*A; % normalized the signal power to M*N, each element 1.

%Normalize the columns
A = sqrt(M)*A*diag(1 ./ sqrt(abs(diag(A'*A))));

%Dictionary error function ? what's mean for this function 
dictionary_error_function =...
    @(q) 20*log10(norm(A -...
    q*find_permutation(A,q),'fro')/norm(A,'fro'));


%scaledFactor = sqrt(N*L/real(trace(X*X')));
%X = scaledFactor*X;

%     X_temp = randn(N,L-N);
%     for ll = 1:L-N
%         yada = randperm(N);
%         yada2 = zeros(N,1);
%         yada2(yada(1:K)) = 1;
%         X_temp(:,ll) = X_temp(:,ll) .* yada2;
%     end
%     X = [eye(N), X_temp];

%% Form the output channel
%Compute noise free output
Z = A*X;
problem.A = A;
problem.X = X;
problem.Z = Z;
%Define the error function
error_function = @(qval) 20*log10(norm(qval - Z,'fro') / norm(Z,'fro'));

opt.error_function = error_function;


%Determine nuw
nuw = 10^(-SNR/10)*N;
nuw = max(nuw,1e-20);
%nuw = norm(reshape(Z,[],1))^2/M/L*10^(-SNR/10);

%Initialize results as empty
results = [];
PrioriIn.noiseVar = nuw;
PrioriIn.lambda = K/N;
PrioriIn.nuX = 1;
PrioriIn.spar = K;
PrioriIn.A = A;
PrioriIn.X = X;

%% EM BiG AMP

if tryBigampEM
    
    %Silence
    opt.verbose = false;
    disp('Starting EM-BiG-AMP')
    
    %Coding
    bestError = inf;
    bestSparsity = inf;
    
    %Compute preconditioner
    if precondition
        Q = chol(inv(Y*Y'));  %? what's mean for Q?
    else
        Q = 1;
    end
    tstart = tic; % start the programme 
    VarNoise = zeros(N,N);
    Aerror = 0;
    Xerror = 0;
    
    scalar = zeros(N,maxTrials);
    EqNoise2 = zeros(M,N,maxTrials);
    ErrorCount = 0;
    Noise = (randn(M,L) + randn(M,L)*1i)/sqrt(2) * sqrt(10^(-SNR/10));
    for Nsample = 1:maxTrials
        %Noisy output channel
        Y = Z + Noise;%sqrt(nuw)*randn(size(Z));
        
        % subspace estimation
        [Uy, Dy, Vy] = svd(Y);
        Un = Uy(:,1:N);
        Dn = Dy(1:N,1:N);
        Vn = Vy(:,1:N);
        problem.Un = Un;
        eqvZ = Dn*Vn';%% why write it in this way? 
        
        %eqvZ = Y;
        %Coding error
        coding_error_function = @(q) 20*log10(coding_error(eqvZ,q,K,useTST));
        
        % Define the error function
        QZ = Q*eqvZ; %Using (possibly) noisy data for error function
        %QZ = Q*Y;
        error_function2 = @(qval) 20*log10(norm(qval - QZ,'fro') / norm(QZ,'fro'));
        error_function2 = @(qval)norm(qval-QZ,'fro')/(norm(QZ,'fro'));
        opt.error_function = error_function2;
        
        for trial = 1:1
            %Run EM-BiGAMP
            %[estFinTemp,~,~,estHistEMtemp] = EMBiGAMP_DL(Q*Y,problem,opt);
            
            %[estFinTemp,~,~,estHistEMtemp] = BiGAMP_DL_fixEM(Q*Y,problem,opt,PrioriIn);
            [estFinTemp,optIn,~,estHistEMtemp] = BiGAMP_DL_fixEM_QPSK(Q*eqvZ,problem,opt,PrioriIn);
            %Correct the dictionary
            estFinTemp.Ahat = Q \ estFinTemp.Ahat;
            
            %Reset error function
            opt.error_function = error_function;
            
            %If the sparisty is better and the error is very small or better
%             if (sum(sum(estHistEMtemp.p1)) < bestSparsity) && ...
%                     ( (estHistEMtemp.errZ(end) < bestError)  ||...
%                     (estHistEMtemp.errZ(end) < -100)  )
                
                %Update
                bestSparsity = sum(sum(estHistEMtemp.p1));
                bestError = estHistEMtemp.errZ(end);
                AhatOptEM = estFinTemp.Ahat;
                XhatOptEM = estFinTemp.xhat;
                estHistEM = estHistEMtemp;
                p1EM = estHistEMtemp.p1;
                
                %Notify user
%                 disp(['Accepting new result. Error: ' num2str(bestError)...
%                     '  Average sparsity: ' num2str(bestSparsity/L)...
%                     '  Max Sparsity: ' num2str(max(sum(p1EM)))])
%            end
            %keyboard
        end
        
        
        %Save results
        tempAhatOptEM = Un*AhatOptEM;
        %tempAhatOptEM = AhatOptEM;
        
        %% LMMSE estimator for scaling ambiguity
        pMatrix2 = find_permutation(A,tempAhatOptEM);
        permMatrix = zeros(size(pMatrix2));
        permMatrix(abs(pMatrix2)>0) = 1; % determine the permutation matrix
        
        
        % 1: full MMSE to obtain the variance
        tempAhat1 = tempAhatOptEM*pMatrix2;
        tempEqNoise= tempAhat1 - A;
        tempNoiseVar = var(tempEqNoise);

        tempAhat = tempAhatOptEM*permMatrix;
        
        temp_scalar = A(1,:)./(abs(A(1,:)).^2+tempNoiseVar);
        %temp_scalar = (A(1,:)+NoiseMean)./(abs(A(1,:)+NoiseMean).^2+NoiseVar );
        scalar1 = temp_scalar.*tempAhat(1,:);
        scalar_MMSE = 1./scalar1;
        % refined LMMSE
        tempScalar = abs(scalar_MMSE);
        meanSca = mean(abs(tempScalar));
        count = 0;
        tempIndex = zeros(1,length(tempScalar));
        for si = 1:length(tempScalar)
           if tempScalar(si) > 2*meanSca
               tempIndex(count+1)= si;
               count = count+1;
               if si==1
                   meanSca = (sum(tempScalar)-tempScalar(si))/length(tempScalar);
                   tempScalar(si) = meanSca;
               else
                   tempScalar(si) = meanSca;
               end
           else
               meanSca = sum(tempScalar(1:si))/si;
           end       
        end
        meanSca= (sum(tempScalar)-sum(tempScalar(tempIndex(1:count))))/(length(tempScalar)-count);
        tempScalar(tempIndex(1:count)) = meanSca;
        scalar_MMSE = sign(scalar_MMSE).*tempScalar;
        %========================
        pMatrix = permMatrix*diag(scalar_MMSE);
        
%         Estimator = (abs(A(1,:)).*ScalarVar)./(abs(A(1,:)).^2*ScalarVar+NoiseVar);
%         EstMean = ScalarMean - Estimator.*(abs(A(1,:)).*ScalarMean);
%         scalar_LMMSE = Estimator.*abs(AhatOptEM(1,:))+EstMean;
%         scalar_LMMSE = scalar_LMMSE.*sign(AhatOptEM(1,:));
%         pMatrix = permMatrix*diag(1./scalar_LMMSE);
        
%         %% LS estimator for scaling ambiguity
%         pMatrix2 = find_permutation(A,tempAhatOptEM);
%         permMatrix = zeros(size(pMatrix2));
%         permMatrix(abs(pMatrix2)>0) = 1; % determine the permutation matrix
%         
%         tempAhat = tempAhatOptEM*permMatrix;
%         scalar_LS = A(1,:)./tempAhat(1,:);
%         pMatrix = permMatrix*diag(scalar_LS);
        
        %%
        Ahat = tempAhatOptEM*pMatrix;
        Xhat = pinv(pMatrix)*XhatOptEM;
        
        Yhat = Ahat*Xhat;
        Zerror = 20*log10(norm(Z-Yhat,'fro')^2/norm(Z,'fro')^2);
        Xhat2 = zeros(size(X));
        Xhat2(abs(Xhat)>1e-6) = Xhat(abs(Xhat)>1e-6);
        
        EqNoise = Ahat - A;
        tempvarNoise = EqNoise'*EqNoise;
        VarNoise = VarNoise.*(Nsample-1)/Nsample + tempvarNoise/Nsample; % average variance
        tempAerror = (norm(A - Ahat,'fro')^2/norm(A,'fro')^2);
        fprintf('%3.4e\n',tempAerror);
        if tempAerror>1e-3
            ErrorCount = ErrorCount + 1;
        end
        
        Aerror = Aerror.*(Nsample-1)/Nsample + tempAerror/Nsample;
        if Aerror > tempAerror
            Aerror = tempAerror;
        end
        
        
        tempXerror = (norm(X - Xhat,'fro')^2/norm(X,'fro')^2);
        Xerror = Xerror.*(Nsample-1)/Nsample + tempXerror/Nsample;
        if Xerror > tempXerror
            Xerror = tempXerror;
        end
       
%         if tempAerror>0.0001
%             kk=0;
%         end
        %
        scalar(:,Nsample) = 1./sum(pMatrix,2);
        EqNoise2(:,:,Nsample) = tempAhatOptEM - A*pinv(pMatrix);
        
    end
    tEMGAMP = toc(tstart);
    
    loc = length(results) + 1;
    results{loc}.name = 'EM-BiG-AMP'; %#ok<*AGROW>
    results{loc}.err = estHistEM.errZ(end);
    results{loc}.time = tEMGAMP;
    results{loc}.ErrorCount = ErrorCount;
    %results{loc}.errHist = estHistEM.errZ;
    %results{loc}.timeHist = estHistEM.timing;
    %results{loc}.dict = AhatOptEM;
    %results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    %results{loc}.codingError =...
    %    coding_error_function(results{loc}.dict);
    
    results{loc}.Aerror = Aerror;
    results{loc}.Xerror = Xerror;
    
    %----------- for test -------------
    %scaleAmbiguity = sum(pMatrix,2);
    %results{loc}.Scaled = scaleAmbiguity; % scalar factor
    results{loc}.VarNoise = VarNoise; % equivalent noise
    results{loc}.scalar = scalar;
    results{loc}.EqNoise2 = EqNoise2;
end

%% SPAMS


if trySpams
    
    %Build params
    spams_param = [];
    spams_param.K = N;
    spams_param.mode = 2;
    spams_param.lambda = 0.1/sqrt(N);
    spams_param.iter = 1000;
    
    
    %Trials
    bestSPAMSerror = inf;
    
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it   
        A_spamsTemp = mexTrainDL(Y,spams_param);
        
        
        SPAMSerrorTemp = dictionary_error_function(A_spamsTemp);
        
        %Update the estimate if this result was superior
        if SPAMSerrorTemp < bestSPAMSerror
            SPAMSerror = SPAMSerrorTemp;
            bestSPAMSerror = SPAMSerror;
            A_spams = A_spamsTemp;
            disp(['Updating Solution. Error was: '...
                num2str(SPAMSerror) ' dB'])
        end
        
    end
    tspams = toc(tstart);
    
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'SPAMS'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(A_spams);
    results{loc}.time = tspams;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = A_spams;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end



%% ER-SpUD


if tryErspud
    
    
    %Call it
    tic
    %[Aerspud,Xerspud] = ER_SpUD_SC(Y);
    [Aerspud,Xerspud]=dl_spud(Y);   %#ok<NASGU>
    tErspud = toc;
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'ER-SpUD (proj)'; %#ok<*AGROW>
    results{loc}.err = dictionary_error_function(Aerspud);
    results{loc}.time = tErspud;
    results{loc}.errHist = results{loc}.err;
    results{loc}.timeHist = zeros(size(results{loc}.errHist));
    results{loc}.dict = Aerspud;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
    
end


%% K-SVD

if tryKsvd
    
    %Build params
    ksvd_params = [];
    ksvd_params.data = Y;
    ksvd_params.Tdata = K;
    ksvd_params.dictsize = N;
    ksvd_params.iternum = 100;
    ksvd_params.exacty = 1;
    
    %Trials
    bestKSVDerror = inf;
    tstart = tic;
    for trial = 1:maxTrials
        
        %Do it
        
        [A_ksvdtemp,~,err_ksvdtemp] = ksvd(ksvd_params);
        
        
        %Fix error
        err_ksvdtemp = err_ksvdtemp.^2 * numel(Y);
        
        %Update the estimate if this result was superior
        if err_ksvdtemp(end) < bestKSVDerror
            bestKSVDerror = err_ksvdtemp(end);
            err_ksvd = err_ksvdtemp;
            A_ksvd = A_ksvdtemp;
            disp(['Updating Solution. Error was: '...
                num2str(10*log10(err_ksvd(end)/norm(Z,'fro')^2))])
        end
        
    end
    tksvd = toc(tstart);
    
    %Save results
    loc = length(results) + 1;
    results{loc}.name = 'K-SVD'; %#ok<*AGROW>
    results{loc}.err = 10*log10(err_ksvd(end)/norm(Z,'fro')^2);
    results{loc}.time = tksvd;
    results{loc}.errHist = 10*log10(err_ksvd/norm(Z,'fro')^2);
    results{loc}.timeHist = zeros(size(err_ksvd));
    results{loc}.dict = A_ksvd;
    results{loc}.dictError = dictionary_error_function(results{loc}.dict);
    results{loc}.codingError =...
        coding_error_function(results{loc}.dict);
    
    
end




%% Store the options structures in results
results{1}.optIn = optIn;
results{1}.trueDict = A;
results{1}.trueEncoding = X;


%% Show Results

if nargin == 0
    results{:} %#ok<NOPRT>
    disp('Note that dictError is the normalized error in dB for recovering')
    disp('the dictionary for each algorithm')
end

