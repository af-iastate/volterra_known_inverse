function [h, y, y_noise] = rlsOnlineSimulation(varargin)
    outCells = rlsOnlineSimulation_inner(varargin{:});
    h = outCells{1};
    y = outCells{2};
    y_noise = outCells{3};
end


function outCells = rlsOnlineSimulation_inner(x, M, p, f_filter, f_addNoise, hIn, M_true)

    if nargin < 7
        M_true = M;
    end
    
    N = length(x);
    y = zeros(N+M-1, 1);
    y_noiseFree = zeros(N+M-1, 1);
    
    % rls
    P = [];
    if nargin<6
        % Calculate filter size
        vvsize = sum(arrayfun(@(P) ...
            nchoosek(M+P-1, P), p ...
        ));
        h = [1; zeros(vvsize-1, 1)];
    else
        h = hIn;
    end
    
    % forward params
    yni = zeros(M_true-1, 1);
    x = [zeros(M-1, 1); x];
    xmem = zeros(M, 1);
    
    for n=1:N-1
        idx = n:-1:n-(M-1);
        xn_volt = voltVecGen(x(idx+M-1), p);
        yn_volt = voltVecGen(y(idx+M-1), p);
        y_inv_hat = xn_volt.' * h;
        [hNew, P] = rlsVolterraOnline(yn_volt, y_inv_hat, M, p, h, P);
        
        xn_volt = voltVecGen(x(idx+M), p);
        y_inv_hat = xn_volt.' * h;
        xmem = [y_inv_hat; xmem(2:end)];
        
        % run through forward system
        [yn, yni] = doForwardIter(y_inv_hat, f_filter, yni);
        y_noiseFree(n+M) = yn;
        y(n+M) = yn + f_addNoise(n);
        h = hNew;
        if mod(n,100)==0
            stem(h)
            pause(0.01)
        end
        
    end
    y_noiseFree = y_noiseFree(M+1:end);
    y = y(M+1:end);
    
    outCells = {h, y_noiseFree, y};
end


function [y, yni_out] = doForwardIter(x, f_filter, yni_in) 
    [y, yni_out] = f_filter(x, yni_in);
end

function [y, ynOut] = linSim(A, x, ynIn)
    % A and ynIn should be column vectors 
    A_scaled = (A / A(1));
    y = x/A(1) - A_scaled(2:end).'*ynIn;
    ynOut = [y; ynIn(1:end-1)];
end