function [HResult, MResult] = applyFIR(varargin)
    [h, M, p, a, Q] = parseLongForm(varargin{:});
    
    function result = hWrapper(tau)
        if any(tau>M-1)
            result = 0;
            return
        end
        result = h(tauToFlatIndex(M, tau, p) + 1);
    end
    
    MResult = M + Q - 1;
    hResultCells = cell(length(p), 1);
    
    baseOffset = 1;
    for ii=1:length(p)
        k = p(ii);
        tauArr = generateTauUpperTri(MResult, k);
        [nRows,~] = size(tauArr);
        hArrK = zeros(nRows, 1);
        for jj=0:nRows-1
            tau = tauArr(jj+1, 1:k);
            Lk = min([tau, Q-1]);
            hVal = sum(arrayfun(@(n) a(n+1)*hWrapper(tau - n), 0:Lk));
            hArrK(jj+1) = hVal;
        end
        hResultCells{ii} = hArrK;
        baseOffset = baseOffset + numel(nRows);
    end
    
    
    HResult = vertcat(hResultCells{:});
end


function [h, M, p, a, Q] = parseLongForm(varargin)
    parser = inputParser();
    parser.addRequired('h', @isnumeric);
    parser.addRequired('M', @isnumeric)
    parser.addRequired('p', @isnumeric)
    parser.addRequired('a', @(x)validateattributes(...
        x, ...
        {'numeric'}, ...
        {'nonempty'} ...
    ));
    parser.parse(varargin{:});
    h = parser.Results.h;
    M = parser.Results.M;
    p = parser.Results.p;
    a = parser.Results.a;
    Q = length(a);
end