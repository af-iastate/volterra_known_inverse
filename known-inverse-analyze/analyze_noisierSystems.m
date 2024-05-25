clc; close all; clearvars;
% For Pseudo alg *only*, 
% plot SFDRs (bad) for each sys
% and norm err for each sys


LpSp_fileName = getMostRecentSimulation('lpspn2');
BpSp_fileName = getMostRecentSimulation('bpspn2');
LpDp_fileName = getMostRecentSimulation('lpdpn2');


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
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone', ...
    'M', 'K', 'h_true' ...
);
h_LpSp_pseudo_two_tone      = result.h_pseudo_two_tone;
h_LpSp_rls_offline_two_tone	= result.h_rls_offline_two_tone;
h_LpSp_rls_online_two_tone  = result.h_rls_online_two_tone;
M_LpSp = result.M;
K_LpSp = result.K;
h_true_LpSp = result.h_true;

result = load(BpSp_fileName, ...
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone', ...
    'M', 'K', 'h_true' ...
);
h_BpSp_pseudo_two_tone = result.h_pseudo_two_tone;
h_BpSp_rls_offline_two_tone	= result.h_rls_offline_two_tone;
h_BpSp_rls_online_two_tone  = result.h_rls_online_two_tone;
M_BpSp = result.M;
K_BpSp = result.K;
h_true_BpSp = result.h_true;

result = load(LpDp_fileName, ...
    'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone', ...
    'M', 'K', 'h_true' ...
);
h_LpDp_pseudo_two_tone = result.h_pseudo_two_tone;
h_LpDp_rls_offline_two_tone	= result.h_rls_offline_two_tone;
h_LpDp_rls_online_two_tone  = result.h_rls_online_two_tone;
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
% -- trained --
sfdr_LpSp_pseudo = [];
sfdr_LpSp_rls_offline = [];
sfdr_LpSp_rls_online = [];
sfdr_true_LpSp = [];

sfdr_BpSp_pseudo = [];
sfdr_BpSp_rls_offline = [];
sfdr_BpSp_rls_online = [];
sfdr_true_BpSp = [];

sfdr_LpDp_pseudo = [];
sfdr_LpDp_rls_offline = [];
sfdr_LpDp_rls_online = [];
sfdr_true_LpDp = [];

curFig = figure(7);
x = xTwoTone;
for n=1:6
    % LpSp
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_pseudo_two_tone{n}, M_LpSp, K_LpSp, freqsTwoTone, fs);
    sfdr_LpSp_pseudo = [sfdr_LpSp_pseudo; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_rls_offline_two_tone{n}, M_LpSp, K_LpSp, freqsTwoTone, fs);
    sfdr_LpSp_rls_offline = [sfdr_LpSp_rls_offline; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_rls_online_two_tone{n}, M_LpSp, K_LpSp, freqsTwoTone, fs);
    sfdr_LpSp_rls_online = [sfdr_LpSp_rls_online; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_true_LpSp, M_LpSp, K_LpSp, freqsTwoTone, fs);
    sfdr_true_LpSp = [sfdr_true_LpSp; a];
    hold off
    
    % BpSp
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_pseudo_two_tone{n}, M_BpSp, K_BpSp, freqsTwoTone, fs);
    sfdr_BpSp_pseudo = [sfdr_BpSp_pseudo; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_rls_offline_two_tone{n}, M_BpSp, K_BpSp, freqsTwoTone, fs);
    sfdr_BpSp_rls_offline = [sfdr_BpSp_rls_offline; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_rls_online_two_tone{n}, M_BpSp, K_BpSp, freqsTwoTone, fs);
    sfdr_BpSp_rls_online = [sfdr_BpSp_rls_online; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_true_BpSp, M_BpSp, K_BpSp, freqsTwoTone, fs);
    sfdr_true_BpSp = [sfdr_true_BpSp; a];
    hold off
    
    % LpDp
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_pseudo_two_tone{n}, M_LpDp, K_LpDp, freqsTwoTone, fs);
    sfdr_LpDp_pseudo = [sfdr_LpDp_pseudo; a];
    hold off
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_rls_offline_two_tone{n}, M_LpDp, K_LpDp, freqsTwoTone, fs);
    sfdr_LpDp_rls_offline = [sfdr_LpDp_rls_offline; a];
    hold off
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_rls_online_two_tone{n}, M_LpDp, K_LpDp, freqsTwoTone, fs);
    sfdr_LpDp_rls_online = [sfdr_LpDp_rls_online; a];
    hold off
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_true_LpDp, M_LpDp, K_LpDp, freqsTwoTone, fs);
    sfdr_true_LpDp = [sfdr_true_LpDp; a];
    hold off   
end


% -- trained --
sfdr_LpSp_pseudo_alt = [];
sfdr_LpSp_rls_offline_alt = [];
sfdr_LpSp_rls_online_alt = [];
sfdr_true_LpSp_alt = [];

sfdr_BpSp_pseudo_alt = [];
sfdr_BpSp_rls_offline_alt = [];
sfdr_BpSp_rls_online_alt = [];
sfdr_true_BpSp_alt = [];

sfdr_LpDp_pseudo_alt = [];
sfdr_LpDp_rls_offline_alt = [];
sfdr_LpDp_rls_online_alt = [];
sfdr_true_LpDp_alt = [];

x = xTwoToneAlt;
for n=1:6
    % LpSp
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_pseudo_two_tone{n}, M_LpSp, K_LpSp, freqsTwoToneAlt, fs);
    sfdr_LpSp_pseudo_alt = [sfdr_LpSp_pseudo_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_rls_offline_two_tone{n}, M_LpSp, K_LpSp, freqsTwoToneAlt, fs);
    sfdr_LpSp_rls_offline_alt = [sfdr_LpSp_rls_offline_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_LpSp_rls_online_two_tone{n}, M_LpSp, K_LpSp, freqsTwoToneAlt, fs);
    sfdr_LpSp_rls_online_alt = [sfdr_LpSp_rls_online_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_LpSp_filter, n, h_true_LpSp, M_LpSp, K_LpSp, freqsTwoToneAlt, fs);
    sfdr_true_LpSp_alt = [sfdr_true_LpSp_alt; a];
    hold off
    
    % BpSp
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_pseudo_two_tone{n}, M_BpSp, K_BpSp, freqsTwoToneAlt, fs);
    sfdr_BpSp_pseudo_alt = [sfdr_BpSp_pseudo_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_rls_offline_two_tone{n}, M_BpSp, K_BpSp, freqsTwoToneAlt, fs);
    sfdr_BpSp_rls_offline_alt = [sfdr_BpSp_rls_offline_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_BpSp_rls_online_two_tone{n}, M_BpSp, K_BpSp, freqsTwoToneAlt, fs);
    sfdr_BpSp_rls_online_alt = [sfdr_BpSp_rls_online_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_BpSp_filter, n, h_true_BpSp, M_BpSp, K_BpSp, freqsTwoToneAlt, fs);
    sfdr_true_BpSp_alt = [sfdr_true_BpSp_alt; a];
    hold off
    
    % LpDp
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_pseudo_two_tone{n}, M_LpDp, K_LpDp, freqsTwoToneAlt, fs);
    sfdr_LpDp_pseudo_alt = [sfdr_LpDp_pseudo_alt; a];
    hold off
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_rls_offline_two_tone{n}, M_LpDp, K_LpDp, freqsTwoToneAlt, fs);
    sfdr_LpDp_rls_offline_alt = [sfdr_LpDp_rls_offline_alt; a];
    hold off 
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_LpDp_rls_online_two_tone{n}, M_LpDp, K_LpDp, freqsTwoToneAlt, fs);
    sfdr_LpDp_rls_online_alt = [sfdr_LpDp_rls_online_alt; a];
    hold off 
    a = analyzeSfdr(x, t, f_LpDp_filter, n, h_true_LpDp, M_LpDp, K_LpDp, freqsTwoToneAlt, fs);
    sfdr_true_LpDp_alt = [sfdr_true_LpDp_alt; a];
    hold off   
end
close all


%% SFDR Plots
% -- Ideal Graphs ---
figure(1)
resize_typeset(3, 1.7*1.9)

tix = [-25,0,25,50,75];
axx = [-.5, .5, 260, 320];

subplot(311)
plotSfdr([sfdr_true_LpSp; sfdr_true_LpSp_alt], 'System 1', 'x');
axis(axx)
xticks([])
xlabel('')
    
subplot(312)
plotSfdr([sfdr_true_BpSp; sfdr_true_BpSp_alt], 'System 2', 'x');
axis(axx)
xticks([])
xlabel('')

subplot(313)
plotSfdr([sfdr_true_LpDp; sfdr_true_LpDp_alt], 'System 3', 'x');

axis(axx)
xticks([])
xlabel('')
print('idealSfdr', outName, '-r900', '-dpdf');


% -- System 1 Graphs ---
figure(2)
resize_typeset(3, 1.7*2.4)
xxlim = [-.1 5.1];

subplot(311)
hold off
plotSfdr(sfdr_LpSp_pseudo)
hold on
plotSfdr(sfdr_LpSp_pseudo_alt, 'Pseudo', 'x')
axis([xxlim, 15, 85]);
yticks([20 40 60 80])

subplot(312)
hold off
plotSfdr(sfdr_LpSp_rls_offline)
hold on
plotSfdr(sfdr_LpSp_rls_offline_alt, 'RLS Offline', 'x')
axis([xxlim, 15, 85]);
yticks([20 40 60 80])

subplot(313)
hold off
plotSfdr(sfdr_LpSp_rls_online)
hold on
plotSfdr(sfdr_LpSp_rls_online_alt, 'RLS Online', 'x')
axis([xxlim, 15, 85]);
yticks([20 40 60 80])

legend({'Training', 'Test'}, 'orientation', 'horizontal');
p = gcf;
p = p.Children(1);
p.Position = [0.355710663646233,0.009826398289926,0.552083327434957,0.04591836628257];
p.Color = 'None';
p.EdgeColor = 'None';

%print('sysOneNSfdr', outName, '-r900', '-dpdf');


% -- System 2 Graphs ---
figure(3)
resize_typeset(3, 1.7*2.4)
xxlim = [-.1 5.1];

subplot(311)
hold off
plotSfdr(sfdr_BpSp_pseudo)
hold on
plotSfdr(sfdr_BpSp_pseudo_alt, 'Pseudo', 'x')
axis([xxlim, -40, 80]);
yticks([-40 0 40 80])

subplot(312)
hold off
plotSfdr(sfdr_BpSp_rls_offline)
hold on
plotSfdr(sfdr_BpSp_rls_offline_alt, 'RLS Offline', 'x')
axis([xxlim, -40, 80]);
yticks([-40 0 40 80])

subplot(313)
hold off
plotSfdr(sfdr_BpSp_rls_online)
hold on
plotSfdr(sfdr_BpSp_rls_online_alt, 'RLS Online', 'x')
axis([xxlim, -40, 80]);
yticks([-40 0 40 80])

legend({'Training', 'Test'}, 'orientation', 'horizontal');
p = gcf;
p = p.Children(1);
p.Position = [0.355710663646233,0.009826398289926,0.552083327434957,0.04591836628257];
p.Color = 'None';
p.EdgeColor = 'None';

% print('sysTwoNSfdr', outName, '-r900', '-dpdf');


% -- System 3 Graphs ---
figure(4)
resize_typeset(3, 1.7*2.4)
xxlim = [-.1 5.1];

subplot(311)
hold off
plotSfdr(sfdr_LpDp_pseudo)
hold on
plotSfdr(sfdr_LpDp_pseudo_alt, 'Pseudo', 'x')
axis([xxlim, 0, 80]);
yticks([0 20 40 60 80])

subplot(312)
hold off
plotSfdr(sfdr_LpDp_rls_offline)
hold on
plotSfdr(sfdr_LpDp_rls_offline_alt, 'RLS Offline', 'x')
axis([xxlim, 0, 80]);
yticks([0 20 40 60 80])

subplot(313)
hold off
plotSfdr(sfdr_LpDp_rls_online)
hold on
plotSfdr(sfdr_LpDp_rls_online_alt, 'RLS Online', 'x')
axis([xxlim, 0, 85]);
yticks([0 20 40 60 80])

legend({'Training', 'Test'}, 'orientation', 'horizontal');
p = gcf;
p = p.Children(1);
p.Position = [0.355710663646233,0.009826398289926,0.552083327434957,0.04591836628257];
p.Color = 'None';
p.EdgeColor = 'None';

% print('sysThreeNSfdr', outName, '-r900', '-dpdf');

%% H Norm Err Plot
% nNf = @(f,g,n) (arrayfun(@(ii) norm(f{ii}-g), n));
% 
% figure(2)
% resize_typeset(3, 1.7*0.75)
% semilogy([1, 2])
% hold off
% db = @(x) 20*log10(x);
% stem(db(nNf(h_LpSp_pseudo_two_tone, h_true_LpSp, 2:6)))
% hold on
% stem(db(nNf(h_BpSp_pseudo_two_tone, h_true_BpSp, 2:6)))
% stem(db(nNf(h_LpDp_pseudo_two_tone, h_true_LpDp, 2:6)))
% legend({'S. 1', 'S. 2', 'S. 3'}, 'location', 'northwest', 'Orientation', 'horizontal')
% xlabel('$n$th iteration')
% ylabel('dB')
% axis([xlim,0,220])
% % print('noisierErr', outName, '-r900', '-dpdf');

%% -- AUX FUNCTIONS --
function result = analyzeSfdr(x, t, f_filter, n, h, M, K, freqs, fs)
    x_preinv = applyVolterra(x, h, M, 1:K);
    if ~nanappear(x_preinv)
        y = f_filter(t, x_preinv);
    else
        y = zeros(size(x_preinv));
    end

    [~,sfdr] = calcSFDR(y, [1, 1250-0, 2500], freqs, [0,fs/2], fs);
    sfdr = sfdr(:);
    a = n*ones(size(sfdr));
    result = [a, sfdr];
end

function p = plotSfdr(X, ttl, ss)
    if nargin < 3
        ss = 'o';
    end
    if nargin < 2
        ttl = '';
    end
    p = scatter(X(:,1)-1, X(:,2), ss);
    title(ttl)
    xlabel('$n$th iteration')
    ylabel('dB')
    xticks(unique(0:5))  % todo generalize...
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