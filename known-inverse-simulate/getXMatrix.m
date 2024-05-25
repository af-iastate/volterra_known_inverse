function X = getXMatrix(x, M, p, varargin)
% X = GETXMATRIX(x, M, p, dataRange)  Get volterra matrix from input x.
% Inputs:
%   x:         Input vector x.
%   M:         System memory M.
%   p:         Vector of term orders for which to compute the matrix.
%   dataRange: Range of x indices from which to construct X (default=M:end).
% Outputs:
%   X: Volterra matrix constructed from x
    
    function X = f_inner(x, M, p, varargin)
        % Parse arguments
        parser = inputParser;
        parser.addRequired('x');
        parser.addRequired('M');
        parser.addRequired('p');
        parser.addOptional('dataRange', []);
        parser.parse(x, M, p, varargin{:});
        dataRange = parser.Results.dataRange;

        % Calculate filter size
        vvsize = length(voltVecGen(zeros(M, 1), p)); 

        % Process dataRange data
        N = length(x);
        if isempty(dataRange)
            dataRange = 1:N;
        elseif length(dataRange) == 1
            dataRange = dataRange:N;
        end

        % construct X array
        x = [zeros(M-1, 1); x(:)]; 
        X = zeros(length(dataRange), vvsize);
        ii = 1;
        for n=dataRange
            % get new data chunks
            idx = n-1:-1:n-M;
            xn = x(idx + M);
            xn = voltVecGen(xn, p);
            X(ii,:) = xn(:).';
            ii = ii + 1;
        end  
    end
    X = f_inner( x, M, p, varargin);
%     X = rememberPastCalls(@f_inner, 'getXMatrix', x, M, p, varargin);
end