clearvars; close all; clc;

len = @length;
fs = 125e6;

% SYSTEM PARAMETERS
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
 % initial linear sys (inverse)
% g = 0.999*gen_lin_ir(17, fs, 20e6, 40e6);  % bandstop character
g = [0.98, 0.1, -0.3, 0.2, 0, 0, 0, 0]; % initial linear sys (inverse)

M = length(g); % initial memory
K = 1; % initial order

% TRUE INVERSE VOLTERRA
[h_true, K] = applyPolymap(g, M, K, f);  % apply polymap
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

simulationRound(h_true, M, K, 'lpspn2_overmem', f_filter, f_filter2, ...
    xGauss, xMLN43, xTwoTone, theNoise, fs)

% h_passthrough = zeros(size(h_true));
% h_passthrough(1) = 1;
% 
% %% MEASURE PERFORMANCE 
% tic
% 
% %% --- GAUSSIAN ----
% N = 5;  % Repeat 5 times
% 
% h_pseudo_g = cell(N+1,1);
% h_rls_offline_g = cell(N+1,1);
% h_rls_online_g = cell(N+1,1);
% h_rla_g = cell(N+1,1);
% 
% snr_pseudo_g = cell(N,1);
% snr_rls_offline_g = cell(N,1);
% snr_rls_online_g = cell(N,1);
% snr_rla_g = cell(N,1);
% 
% % The initial assumption is that h is passthrough
% h_pseudo_g{1} = h_passthrough;
% h_rls_offline_g{1} = h_passthrough;
% h_rls_online_g{1} = h_passthrough;
% h_rla_g{1} = h_passthrough;
% 
% x = xGauss;
% t = (0:length(x)-1)/fs.';
% 
% for ii=2:N+1 
%     % -------- CALC H ----------------
%     % Pseudo
%     h_cur = h_pseudo_g{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_pseudo_g{ii} = pseudoInverseSimulation(x_predist, Y);
%     snr_pseudo_g{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Offline
%     h_cur = h_rls_offline_g{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_rls_offline_g{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
%     snr_rls_offline_g{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Online
%     h_cur = h_rls_online_g{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     [h_rls_online_g{ii}, y] = rlsOnlineSimulation(x_predist, M, 1:K, ...
%         @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, theNoise), h_cur);
%     snr_rls_online_g{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLA
%     h_cur = h_rla_g{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     h_rla_g{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
% end
% 
% 
% %% --- MLN ----
% N = 5;  % Repeat 5 times
% 
% h_pseudo_mln = cell(N+1,1);
% h_rls_offline_mln = cell(N+1,1);
% h_rls_online_mln = cell(N+1,1);
% h_rla_mln = cell(N+1,1);
% 
% snr_pseudo_mln = cell(N,1);
% snr_rls_offline_mln = cell(N,1);
% snr_rls_online_mln = cell(N,1);
% snr_rla_mln = cell(N,1);
% 
% % The initial assumption is that h is passthrough
% h_pseudo_mln{1} = h_passthrough;
% h_rls_offline_mln{1} = h_passthrough;
% h_rls_online_mln{1} = h_passthrough;
% h_rla_mln{1} = h_passthrough;
% 
% x = xMLN43;
% t = (0:length(x)-1)/fs.';
% 
% for ii=2:N+1
%     % -------- CALC H ----------------
%     % Pseudo
%     h_cur = h_pseudo_mln{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_pseudo_mln{ii} = pseudoInverseSimulation(x_predist, Y);
%     snr_pseudo_mln{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Offline
%     h_cur = h_rls_offline_mln{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_rls_offline_mln{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
%     snr_rls_offline_mln{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Online
%     h_cur = h_rls_online_mln{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     [h_rls_online_mln{ii}, y] = rlsOnlineSimulation(x_predist, M, 1:K, ...
%         @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, theNoise), h_cur);
%     snr_rls_online_mln{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLA
%     h_cur = h_rla_mln{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     h_rla_mln{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
% end
% 
% %% --- TWO TONE -------------
% N = 5;  % Repeat 5 times
% 
% h_pseudo_two_tone = cell(N+1,1);
% h_rls_offline_two_tone = cell(N+1,1);
% h_rls_online_two_tone = cell(N+1,1);
% h_rla_two_tone = cell(N+1,1);
% 
% snr_pseudo_two_tone = cell(N,1);
% snr_rls_offline_two_tone = cell(N,1);
% snr_rls_online_two_tone = cell(N,1);
% snr_rla_two_tone = cell(N,1);
% 
% % The initial assumption is that h is passthrough
% h_pseudo_two_tone{1} = h_passthrough;
% h_rls_offline_two_tone{1} = h_passthrough;
% h_rls_online_two_tone{1} = h_passthrough;
% h_rla_two_tone{1} = h_passthrough;
% 
% x = xTwoTone;
% t = (0:length(x)-1)/fs.';
% 
% for ii=2:N+1
% % -------- CALC H ----------------
%     % Pseudo
%     h_cur = h_pseudo_two_tone{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_pseudo_two_tone{ii} = pseudoInverseSimulation(x_predist, Y);
%     snr_pseudo_two_tone{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Offline
%     h_cur = h_rls_offline_two_tone{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     y = f_filter(t, x_predist);
%     y_noise = y + addNoise(1, numel(y), theNoise);
%     Y = getXMatrix(y_noise, M, 1:K);
%     h_rls_offline_two_tone{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
%     snr_rls_offline_two_tone{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLS Online
%     h_cur = h_rls_online_two_tone{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     [h_rls_online_two_tone{ii}, y] = rlsOnlineSimulation(x_predist, M, 1:K, ...
%         @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, theNoise), h_cur);
%     snr_rls_online_two_tone{ii-1} = norm(y)./norm(theNoise);
%     
%     % RLA
%     h_cur = h_rla_two_tone{ii-1};
%     x_predist = applyVolterra(x, h_cur, M, 1:K);
%     h_rla_two_tone{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
% end
% %% MEASURE PERFORMANCE
% total_sim_time = toc;
% disp(total_sim_time)
% 
% %% SAVE RESULTS
% file_name = ['sim_v301_', datestr(now, 'yyyymmddTHHMMSS')];
% field_names = {
%     'h_true'
%     'h_pseudo_g'
%     'h_rls_offline_g'
%     'h_rls_online_g'
%     'h_rla_g'
%     'h_pseudo_mln'
%     'h_rls_offline_mln'
%     'h_rls_online_mln'
%     'h_rla_mln'
%     'h_pseudo_two_tone'
%     'h_rls_offline_two_tone'
%     'h_rls_online_two_tone'
%     'h_rla_two_tone'
%     'total_sim_time'
% }.';
% save([file_name, '_min.mat'], field_names{:})
% save([file_name, '_full.mat'])

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