function func_inv = getInverseFromPoly(polynomial, varargin)
    p = inputParser;
    p.addRequired('poly', @(x) isnumeric(x) && isvector(x) );
    p.addOptional('bounds', [], @(x) isnumeric(x) && length(x) == 2 && all(diff(x) > 0));
    p.addOptional('numSamples', 1001, @(x) isnumeric(x) && isscalar(x));
    p.parse(polynomial, varargin{:});
    bounds = p.Results.bounds;
    if isempty(bounds)
        bounds = [-1, 1];
    end
    numSamples = p.Results.numSamples;
    
    xs = linspace(bounds(1), bounds(2), numSamples);
    ys = polyval(polynomial, xs);
    
    func_inv = @(y) interp1(ys, xs, y, 'linear', 'extrap');
    
end