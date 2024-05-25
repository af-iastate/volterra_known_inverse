function result = tauToFlatIndex(M, tau, p)
    
    function result = phi(M, tau, m)
        result = 0;
        result = result + sum(arrayfun(@(K) ...
            nchoosek(K+M-1, K), ...
            m ...
        ));
        result = result + tauToOffsetWithinOrder(M, tau); 
    end

    [nRows, K] = size(tau);
    if nargin<3
        p = 1:K;
    end
    m = p(p<K);
    
    phiWrapped = @(tau) phi(M, tau, m);
    
    if nRows == 1
        result = phiWrapped(tau);
        return
    end

    result = arrayfun(@(ii) phiWrapped(tau(ii, :)), 1:nRows);
    result = result(:);
end