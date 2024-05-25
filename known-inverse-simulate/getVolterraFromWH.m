function h = getVolterraFromWH(varargin)
    % Process Inputs
    [A, B, M, p, alpha, beta, f] = attemptParseArgs(varargin{:});
    if length(M) > 1
        M = sum(M);
    end
    
    % Determine h
    hCells = cell(numel(p), 1);
    ii = 1;
    for m=p(:).'
        tau = genIdx(M, m);
        [N,~] = size(tau);
        htmp = arrayfun(@(jj) getHFromTau(tau(jj, 1:m), A, B, f(ii), alpha, beta), 1:N).';
        hCells{ii} = htmp;
        ii = ii + 1;
    end
    h = vertcat(hCells{:});
end


function [A, B, M, p, alpha, beta, f] = attemptParseArgs(varargin)
    
    function f = parsePolyInv(x, p)
        if isa(x, 'function_handle')
            f = x(p);
            return
        end 
        if isnumeric(x)
            [nRows, nCols] = size(x);
            if nRows == 1 || nCols == 1
                xCols = x(:);
                N = numel(xCols);
                if N == numel(p)
                    f = xCols;
                    return
                elseif N > max(p) && all(p > 1 & p <= N)
                    f = xCols(p);
                    return
                end
            end
        end
        f = 1./p;
    end
    
    % Attempt Parse System Struct
    failedParse = false;
    parser = inputParser();
    parser.addRequired('whSystem', @isstruct);
    parser.addOptional('M', [], @isnumeric)
    parser.addOptional('p', [], @isnumeric)
    parser.addOptional('alpha', [], @isnumeric);
    parser.addOptional('beta', [], @isnumeric);
    parser.addOptional('f', [], @(x)validateattributes(...
        x, ...
        {'numeric', 'function_handle'}, ...
        {'nonempty'} ...
    ));
    try
        parser.parse(varargin{:});
    catch
        failedParse = true;
    end
    if ~failedParse
        % Parser System Struct
        mSystem = parser.Results.whSystem;
        A = mSystem.Af;
        B = mSystem.Ai;
        M = parser.Results.M;
        if isempty(M)
            M = fliplr(mSystem.M);
        end
        p = parser.Results.p;
        if isempty(p)
            p = mSystem.p;
        end
        fRaw = parser.Results.f;
        if isempty(fRaw)
            fRaw = mSystem.f_polyInv;
        end
        f = parsePolyInv(fRaw, p);
        
        alpha = parser.Results.alpha;
        if isempty(alpha)
            alpha = mSystem.alpha;
        end
        beta = parser.Results.beta;
        if isempty(beta)
            beta = mSystem.beta;
        end
        return
    end
    % Attempt Parse Individual Params
    parser = inputParser();
    parser.addRequired('A');
    parser.addRequired('B');
    parser.addRequired('M', @isnumeric)
    parser.addRequired('p', @isnumeric)
    parser.addOptional('alpha', 1, @isnumeric);
    parser.addOptional('beta', 1, @isnumeric);
    parser.addOptional('f', @(x) 1/x, @(x)validateattributes(...
        x, ...
        {'numeric', 'function_handle'}, ...
        {'nonempty'} ...
    ));
    parser.parse(varargin{:});

    A = parser.Results.A(:).';
    B = parser.Results.B(:).';
    M = parser.Results.M;
    p = parser.Results.p;
    alpha = parser.Results.alpha;
    beta = parser.Results.beta;
    f = parsePolyInv(parser.Results.f, p);
end


function q = genIdx(M, P)
    function q = inner(s, M, p)
        if length(p)==1
            q = (s:M-1).';
            return
        end
        qCells = cell(M-s, 1);
        for ii=s:M-1
            res = inner(ii, M, p(2:end));
            [N, ~] = size(res);
            qtmp = [ii*ones(N, 1), res];
            qCells{ii-s+1} = qtmp; 
%             if ~isempty(q)
%                 q = [q; qtmp];
%             else
%                 q = qtmp;
%             end
        end
        q = vertcat(qCells{:});
    end
    q = inner(0, M, 1:P);
end


function h = getHFromTau(tau, A, B, f, alpha, beta)
    function m = inner(tau, M, f, A, beta)
        k = countIndices(tau, M(1));
        P = length(tau);
        if sum(k) == P
            m = (1/beta)^P * f * multinomialCoefficient(k) * prod(A.^k);
        else
            m = 0;
        end
    end
    M = cellfun(@length,{A, B});
    L = min([tau, M(2)-1]);
    g = arrayfun(@(ii) B(ii+1)*inner(tau - ii, M, f, A, beta), 0:L);
    h = (1/alpha)*sum(g);
end


function m = multinomialCoefficient(k, n)
    if nargin < 2
        n = sum(k);
    end
    num = factorial(n);
    kFact = arrayfun(@factorial, k);
    den = prod(kFact);
    m = num / den;
end


function k = countIndices(tau, M)
    k = arrayfun(@(x) sum(tau==x), 0:M-1);
end