function X  = generateX(optIn)
    %generate X with iid locations of non-zero elements
    
    N = optIn.N;
    L = optIn.L;
    spar = optIn.spar;
    alphabet = optIn.alphabet;  %constellation

    %determin locations of non-zero elements
    X = rand(N,L) < spar;
    alLen = length(alphabet);
    
    %generate modulated signals i.i.d.
    randMod = ceil(rand(N,L).*alLen);
    X = X.* alphabet(randMod);
    for i = 1:N
        temp = find(abs(X(i,:)) > 0);
        %set the first non-zero symbol to 1
        X(i,:) = X(i,:) / X(i,temp(1)); 
    end
end
