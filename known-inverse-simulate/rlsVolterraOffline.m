function [h, P] = rlsVolterraOffline(X, y, range, M, p, hIn, P, lambda)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    
    % Calculate filter size
    vvsize = sum(arrayfun(@(P) ...
        nchoosek(M+P-1, P), p ...
    ));

    %vvsize = length(voltVecGen(zeros(M,1),p));

    % forgetting factor
    if nargin<8
        lambda = [];
    end
    if isempty(lambda)
        lambda = 1;
    end
    % initialize P matrix to identity
    if nargin<7
        P = [];
    end
    if isempty(P)
        P = eye(vvsize);
    end
    % initialize filter to all-pass linear filter (if not given)
    if nargin<6
        hIn = [];
    end
    if isempty(hIn)
        hIn = [1; zeros(vvsize-1,1)];
    end

    % Process range data
    N = length(y);
    if isempty(range)
        range = M:N;
    end

    h = hIn;
    for n=range
        yn = y(n);
        xn = X(n,:).';
        
        % (hopefully) faster alpha calc
        alpha = yn-h.'*xn;

        % update gain vector
        Px = P*xn;
        k = Px/(lambda+xn'*Px);

        % matrix update
        P = (1/lambda)*(P-k*xn'*P);

        % update kernel vector (check if conj(alpha) needed)
        h = h + k*alpha;
    end
end

