function [Ahat, Xhat] = AX_DeSVD_Compute(Ahat,xhat,problem)
Q = 1;
Ahat = Q \ Ahat;
Aerror = 0;
Xerror = 0;
%Reset error function
A = problem.A;
X = problem.X;
Un = problem.Un;
Nsample = 1;
%Update

AhatOptEM = Ahat;
XhatOptEM = xhat;



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
pMatrix = permMatrix*diag(scalar_MMSE);


%%
Ahat = tempAhatOptEM*pMatrix;
try
    Xhat = pinv(pMatrix)*XhatOptEM;
catch
    Xhat = randn(problem.N,problem.L);
end

tempAerror = (norm(A - Ahat,'fro')^2/norm(A,'fro')^2);
Aerror = Aerror.*(Nsample-1)/Nsample + tempAerror/Nsample;

tempXerror = (norm(X - Xhat,'fro')^2/norm(X,'fro')^2);
Xerror = Xerror.*(Nsample-1)/Nsample + tempXerror/Nsample;


end