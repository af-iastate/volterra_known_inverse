function [h] = pseudoInverseSimulation(x, Y, gamma)
    if nargin < 3
        gamma = 0;
    end
    if gamma > 0
        h = ridge(x, Y, gamma, 0);
        h = h(2:end);
    else
        h = pinv(Y)*x(:);
    end
end