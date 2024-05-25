clc; clearvars; close all;
LpSp_fileName = getMostRecentSimulation('lpspn2_overmem');

% FORWARD SYSTEM (NEEDED BECAUSE FUNCTION HANDLES NOT SAVED IN MAT)
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
 % initial linear sys (inverse)
g = [0.98, 0.1, -0.3, 0.2, 0, 0, 0, 0]; % initial linear sys (inverse)
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
% G = tf([1 zeros(1, length(g)-1)], g(:).', 1/fs); 
% f_filter = @(t, x) lsim(G, (inv_f(x)), t);
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
f_filter = @(t, x) filter(GB, GA, inv_f(x));

analyzeRound(fileName, f_filter)

%% AUX FUNCTIONS
function u = thirdOrderActivationFunc(t, a1, a3)
    if a1 ~= 0 && a3 ~= 0
        tmpf = nthroot(-9*a3^2*t + sqrt(12*a1^3*a3^3 + 81*a3^4*t.^2), 3);
        u = (2/3)^(1/3)*a1./tmpf - tmpf./(2^(1/3)*3^(2/3)*a3);
    elseif a3 ~= 0
        u = nthroot(t./a3, 3);
    elseif a1 ~= 0
        u = (1/a1).*t;
    else
        u = 0.*t;
    end
end
