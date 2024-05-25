clc; clearvars; close all;
% H = 0;
% M = 3;
% K = 3;
% a = [1 1 1];
% HOut = applyFIR(H, M, K, a);


% Get H without FIR
% A, B, M, p, alpha, beta, f

Ma = 6;
Mb = 4;
A = getFIR(Ma, [25e6, 40e6], 'stop');
B = getFIR(Mb, [10e6], 'high');

f = [1, -0.5, .25, 0.3];

p = 1:4;

h1 = getVolterraFromWH(A, [1].', Ma, p, 1, 1, f);
h2 = getVolterraFromWH(A, B, Ma+Mb-1, p, 1, 1, f);

% stem(h1)
% hold on;
% stem(h2);

h2Hat = applyFIR(h1,Ma,p,B);

% stem(h2)
% hold on;
% stem(h2Hat);

stem(h2Hat-h2);

function h = getFIR(N, fz, type)
    Fs = 125e6;

    [B, A] = butter(4, fz, type, 's');
    H = tf(B, A);


    G = c2d(H, 1/Fs, 'zero');

    t = (0:N-1)/Fs;
    [h, ~] = impulse(G, t);
    h = h/Fs;
end
