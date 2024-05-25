close all; clc; clearvars;

if ~isfile('data.mat')
    make_test_data()
%     clc; close all; clearvars;
end
load('data.mat')

len = @length;
fs = 125e6;

% SYSTEM PARAMETERS
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
 % initial linear sys (inverse)
% g = 0.999*gen_lin_ir(17, fs, 20e6, 40e6);  % bandstop character
g = [0.98, 0.1, -0.3, 0.2]; % initial linear sys (inverse)

M = length(g); % initial memory
K = 1; % initial order

% TRUE INVERSE VOLTERRA
[h_true, K] = applyPolymap(g, M, K, f);  % apply polymap
[h_true, K] = applyPolymap(h_true, M, K, f);  % apply polymap again

% FORWARD SYSTEM
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
% G = tf([1 zeros(1, length(g)-1)], g(:).', 1/fs); 
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
f_filter = @(t, x) filter(GB, GA, inv_f(inv_f(x)));

simulateRidgeRegressionRound(...
        M, K, 'lpdp_g0', f_filter,  ...
        xTwoTone, noiseLvls, fs, 0);
    
    
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