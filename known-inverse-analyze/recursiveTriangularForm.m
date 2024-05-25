function [y] = recursiveTriangularForm(x,p,M,xp,mp)
%recursiveTriangularForm Produces input vector for Volterra filter
%   recursiveTriangularForm(x,P) produces a vector y containing all the
%   polynomial combinations of the vector x necessary for implementing a
%   Pth order (homogeneous) Volterra filter.  It is assumed that x is in
%   the form
%           x = [ x(n) x(n-1) ... x(n-M) ]'
%   where M is the memory depth of the filter.  For example:
%           recursiveTriangularForm((1:3)',2)
%   returns
%           [ x^2(n) = 1,
%             x(n)x(n-1) = 2,
%             x(n)x(n-2) = 3,
%             x^2(n-1) = 4,
%             x(n-1)x(n-2) = 6,
%             x^2(n-2) = 9 ]
%
%   As the name implies, the coefficients are in triangular form and are
%   calculated recursively.

% Andrew Bolstad, Oct. 15, 2019

% to start with, p = 1:P.  take the first one off each time

% initialize
y = [];
if nargin<3
    xp = 1;
    mp = 1;
    M = length(x);
    
    % assume user enters, e.g., p = 3 for third order
    p = 1:max(p);
end

if length(p)==1
    % innermost loop (really mp:M-1 with 1-index)
    y = xp*x(mp:M);
    y = y(:);
else
    for m = mp:M
        % multiply by relevant x
        xcum = xp*x(m);
        
        % recurse
        ytmp = recursiveTriangularForm(x,p(2:end),M,xcum,m);
        
        y = [y; ytmp];
    end
end


end

