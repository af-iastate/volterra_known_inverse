clc; clearvars; close all;
fileName = getMostRecentSimulation('bpspn1');
% FORWARD SYSTEM (NEEDED BECAUSE FUNCTION HANDLES NOT SAVED IN MAT)
load(fileName, 'fs')
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
g = 0.999*gen_lin_ir(17, fs, 20e6, 40e6);  % bandstop character
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
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

function h = gen_lin_ir(sz, Fs, f0, f1, tp)
    if nargin<5
        tp = 'stop';
    end
    [B, A] = butter(1, [f0,f1]*2*pi, tp, 's');
    H = tf(B, A);
    G = c2d(H, 1/Fs, 'tustin');

    t = (0:sz-1)/Fs;
    [h, ~] = impulse(G, t);
    h = h/Fs;
end
