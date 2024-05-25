clc; close all; clearvars;
% For Pseudo alg *only*, 
% plot SFDRs (good) for each sys
% and norm err for each sys

LpSp_fileName = getMostRecentSimulation('lpspn5');
BpSp_fileName = getMostRecentSimulation('bpspn5');
LpDp_fileName = getMostRecentSimulation('lpdpn5');


% load common params
if ~isfile('data.mat')
    make_test_data()
end
result = load('data.mat', ...
    'xTwoTone', 'xTwoToneAlt', ...
    'freqsTwoTone', 'freqsTwoToneAlt' ...
);
xTwoTone = result.xTwoTone;
xTwoToneAlt = result.xTwoToneAlt;
freqsTwoTone = result.freqsTwoTone;
freqsTwoToneAlt = result.freqsTwoToneAlt;

result = load(LpSp_fileName, ...
    't', 'fs' ...
);
t = result.t;
fs = result.fs;

%% -- Load Specifics ---
result = load(LpSp_fileName, ...
    'h_pseudo_two_tone', 'M', 'K', 'h_true' ...
);
h_LpSp = result.h_pseudo_two_tone;
M_LpSp = result.M;
K_LpSp = result.K;
h_true_LpSp = result.h_true;

result = load(BpSp_fileName, ...
    'h_pseudo_two_tone', 'M', 'K', 'h_true' ...
);
h_BpSp = result.h_pseudo_two_tone;
M_BpSp = result.M;
K_BpSp = result.K;
h_true_BpSp = result.h_true;

result = load(LpDp_fileName, ...
    'h_pseudo_two_tone', 'M', 'K', 'h_true' ...
);
h_LpDp = result.h_pseudo_two_tone;
M_LpDp = result.M;
K_LpDp = result.K;
h_true_LpDp = result.h_true;

% -- sys functions ---
% LpSp
f = [1, 0, 0.1];
g = [0.98, 0.1, -0.3, 0.2]; 
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
f_LpSp_filter = @(t, x) filter(GB, GA, inv_f(x));

% BpSp
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
g = 0.999*gen_lin_ir(17, fs, 20e6, 40e6);  % bandstop character
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
f_BpSp_filter = @(t, x) filter(GB, GA, inv_f(x));

% LpDp
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
g = [0.98, 0.1, -0.3, 0.2]; % initial linear sys (inverse)
inv_f = @(x) thirdOrderActivationFunc(x, f(1), f(3));
GB = [1 zeros(1, length(g)-1)];
GA = g(:).';
f_LpDp_filter = @(t, x) filter(GB, GA, inv_f(inv_f(x)));

%% -- Calc SFDRs ---
sfdr_LpSp = [];
sfdr_BpSp = [];
sfdr_LpDp = [];

curFig = figure(7);
% -- trained --
x = xTwoTone;
for n=1:6
    % LpSp
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp{n}, M_LpSp, K_LpSp, freqsTwoTone, fs);
    sfdr_LpSp = [sfdr_LpSp; a];
    hold off
    
    % BpSp
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp{n}, M_BpSp, K_BpSp, freqsTwoTone, fs);
    sfdr_BpSp = [sfdr_BpSp; a];
    hold off
    
    % LpDp
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp{n}, M_LpDp, K_LpDp, freqsTwoTone, fs);
    sfdr_LpDp = [sfdr_LpDp; a];
    hold off   
end

sfdr_LpSp_alt = [];
sfdr_BpSp_alt = [];
sfdr_LpDp_alt = [];
x = xTwoToneAlt;
for n=1:6
    % LpSp
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp{n}, M_LpSp, K_LpSp, freqsTwoToneAlt, fs);
    sfdr_LpSp_alt = [sfdr_LpSp_alt; a];
    hold off
    
    % BpSp
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp{n}, M_BpSp, K_BpSp, freqsTwoToneAlt, fs);
    sfdr_BpSp_alt = [sfdr_BpSp_alt; a];
    hold off
    
    % LpDp
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp{n}, M_LpDp, K_LpDp, freqsTwoToneAlt, fs);
    sfdr_LpDp_alt = [sfdr_LpDp_alt; a];
    hold off   
end
close all


%% SFDR Plot
figure(1)
resize_typeset(3, 1.7*2)

tix = [50,100,150,200,250];
axx = [1.9, 6.1, 50, 250];

subplot(311)
hold off
plotSfdr(sfdr_LpSp, 'SFDR System 1')
hold on
plotSfdr(sfdr_LpSp_alt, 'SFDR System 1', 'x')
legend({'Training', 'Test'}, 'location', 'southeast', 'Orientation', 'horizontal')
axis(axx)
yticks(tix)
    
subplot(312)
hold off
plotSfdr(sfdr_BpSp, 'SFDR System 1')
hold on
plotSfdr(sfdr_BpSp_alt, 'SFDR System 2', 'x')
axis(axx)
yticks(tix)

subplot(313)
hold off
plotSfdr(sfdr_LpDp, 'SFDR System 1')
hold on
plotSfdr(sfdr_LpDp_alt, 'SFDR System 3', 'x')
axis(axx)
yticks(tix)
print('noiseFreeSfdr', outName, '-r900', '-dpdf');

%% H Norm Err Plot
nNf = @(f,g,n) (arrayfun(@(ii) norm(f{ii}-g), n));
db = @(x) 20*log10(x);

figure(2)
resize_typeset(3, 1.7*0.75)
hold off
stem(db(nNf(h_LpSp, h_true_LpSp, 2:6)))
hold on
stem(db(nNf(h_BpSp, h_true_BpSp, 2:6)))
stem(db(nNf(h_LpDp, h_true_LpDp, 2:6)))
legend({'S. 1', 'S. 2', 'S. 3'}, 'location', 'southwest', 'Orientation', 'horizontal')
xlabel('$n$th iteration')
ylabel('dB')
axis([xlim,-350, 0])
print('noiseFreeErr', outName, '-r900', '-dpdf');

%% -- AUX FUNCTIONS --
function result = analyzeSfdr(x, t, f_filter, n, h, M, K, freqs, fs)
    x_preinv = applyVolterra(x, h, M, 1:K);
    if ~nanappear(x_preinv)
        y = f_filter(t, x_preinv);
    else
        y = zeros(size(x_preinv));
    end

    [~,sfdr] = calcSFDR(y, [1, 1250-0, 1667], freqs, [0,fs/2], fs);
    sfdr = sfdr(:);
    a = n*ones(size(sfdr));
    result = [a, sfdr];
end

function plotSfdr(X, ttl, ss)
    if nargin < 3
        ss = 'o';
    end
    idx = X(:,1) > 1;
    scatter(X(idx,1), X(idx,2), ss)
    title(ttl)
    xlabel('$n$th iteration')
    ylabel('dB')
    xticks(unique(2:6))  % todo generalize...
    xticklabels({'1', '2', '3', '4', '5'})
%     grid on
end

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

function r=nanappear(n)
  r=(any(isnan(n(:))));
end