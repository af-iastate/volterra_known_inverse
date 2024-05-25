clearvars; close all; clc;

len = @length;
fs = 125e6;

% SYSTEM PARAMETERS
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
 % initial linear sys (inverse)
g = 0.999*gen_lin_ir(17, fs, 20e6, 40e6);  % bandstop character
M_true = length(g); % initial memory
K = 1; % initial order

% TRUE INVERSE VOLTERRA
[h_true, K] = applyPolymap(g, M_true, K, f);  % apply polymap
% [h_true, K] = applyPolymap(h_true, M, K, f);  % apply polymap again

% GENERATE INITIAL DATA
if ~isfile('data.mat')
    make_test_data()
end
load('data.mat')

% FORWARD SYSTEM
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
%G = tf([1 zeros(1, length(g)-1)], g(:).', 1/fs); 
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
% f_filter = @(t, x) lsim(G, inv_f(inv_f(x)), t);
f_filter = @(t, x) filter(GB, GA, inv_f(x));
f_filter2 = @(t, x, xi) filter(GB, GA, inv_f(x), xi);

theNoise = noise_sig{2};
M = 9;

simulationRound(h_true, M, K, 'bpspn2_undermem', f_filter, f_filter2, ...
    xGauss, xMLN43, xTwoTone, theNoise, fs, M_true)


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


function result = addNoise(n, N, noise_sig)
    result = noise_sig(n:n+N-1);
end