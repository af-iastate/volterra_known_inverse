function [hOut, POut] = rlsVolterraOnline(xn, y, M, p, hIn, PIn, lambda)
    
    % Calculate filter size
    vvsize = sum(arrayfun(@(P) ...
        nchoosek(M+P-1, P), p ...
    ));

    % forgetting factor
    if nargin<7
        lambda = [];
    end
    if isempty(lambda)
        lambda = 1;
    end
    % initialize P matrix to identity
    if nargin<6
        PIn = [];
    end
    if isempty(PIn)
        PIn = eye(vvsize);
    end
    % initialize filter to all-pass linear filter (if not given)
    if nargin<5
        hIn = [];
    end
    if isempty(hIn)
        hIn = [1; zeros(vvsize-1,1)];
    end
    
    % ---------------------------------------------

    alpha = y-hIn.'*xn;

    % update gain vector
    Px = PIn*xn;
    k = Px/(lambda+xn'*Px);

    % matrix update
    POut = (1/lambda)*(PIn-k*xn'*PIn);

    % update kernel vector (check if conj(alpha) needed)
    hOut = hIn + k*alpha;
    
end

