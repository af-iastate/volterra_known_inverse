function simulationRound(h_true, M, K, theName, f_filter, f_filter2, xGauss, xMLN, xTwoTone, noiseLvl, fs, M_true, gamma)
    if nargin < 13
        gamma = 0;
    end
    if nargin < 12
        M_true = M;
    end
    
    h_under = zeros(size(voltVecGen(1:M, 1:K)));
    if numel(h_under) == numel(h_true)
        h_under = h_true;
    end
    
    h_passthrough = zeros(size(h_under));
    h_passthrough(1) = 1;
    
    nX = length(xTwoTone);
    noise_signal = cell(5,1);
    rng(23)
    for i=1:5
        noise_signal{i} = makeNoise(nX, noiseLvl);
    end
    
    %% MEASURE PERFORMANCE 
    tic

    %% --- GAUSSIAN ----
    N = 5;  % Repeat 5 times

    h_pseudo_g = cell(N+1,1);
    h_rls_offline_g = cell(N+1,1);
    h_rls_online_g = cell(N+1,1);
    h_rla_g = cell(N+1,1);

    x_predist_pseudo_g = cell(N,1);
    x_predist_rls_offline_g = cell(N,1);
    x_predist_rls_online_g = cell(N,1);
    x_predist_rla_g = cell(N,1);   

    y_pseudo_g = cell(N,1);
    y_rls_offline_g = cell(N,1);
    y_rls_online_g = cell(N,1);
    y_rla_g = cell(N,1);   
    
    y_noise_pseudo_g = cell(N,1);
    y_noise_rls_offline_g = cell(N,1);
    y_noise_rls_online_g = cell(N,1);
    y_noise_rla_g = cell(N,1);   
    
    snr_pseudo_g = cell(N,1);
    snr_rls_offline_g = cell(N,1);
    snr_rls_online_g = cell(N,1);
    snr_rla_g = cell(N,1);

    % The initial assumption is that h is passthrough
    h_pseudo_g{1} = h_passthrough;
    h_rls_offline_g{1} = h_passthrough;
    h_rls_online_g{1} = h_passthrough;
    h_rla_g{1} = h_passthrough;

    x = xGauss;
    t = (0:length(x)-1)/fs.';

    for ii=2:N+1 
        % -------- CALC H ----------------
        % Pseudo
        h_cur = h_pseudo_g{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_pseudo_g{ii} = pseudoInverseSimulation(x_predist, Y, gamma);
        
        x_predist_pseudo_g{ii-1} = x_predist;
        y_pseudo_g{ii-1} = y;
        y_noise_pseudo_g{ii-1} = y_noise;
        snr_pseudo_g{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Offline
        h_cur = h_rls_offline_g{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_rls_offline_g{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
        
        x_predist_rls_offline_g{ii-1} = x_predist;
        y_rls_offline_g{ii-1} = y;
        y_noise_rls_offline_g{ii-1} = y_noise;
        snr_rls_offline_g{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Online
        h_cur = h_rls_online_g{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        [h_rls_online_g{ii}, y, y_noise] = rlsOnlineSimulation(x_predist, M, 1:K, ...
            @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, noise_signal{ii-1}), h_cur, M_true);
        
        x_predist_rls_online_g{ii-1} = x_predist;
        y_rls_online_g{ii-1} = y;
        y_noise_rls_online_g{ii-1} = y_noise;
        snr_rls_online_g{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLA
%         h_cur = h_rla_g{ii-1};
%         x_predist = applyVolterra(x, h_cur, M, 1:K);
%         h_rla_g{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
    end


    %% --- MLN ----
    N = 5;  % Repeat 5 times

    h_pseudo_mln = cell(N+1,1);
    h_rls_offline_mln = cell(N+1,1);
    h_rls_online_mln = cell(N+1,1);
    h_rla_mln = cell(N+1,1);
    
    x_predist_pseudo_mln = cell(N,1);
    x_predist_rls_offline_mln = cell(N,1);
    x_predist_rls_online_mln = cell(N,1);
    x_predist_rla_mln = cell(N,1);  
    
    y_pseudo_mln = cell(N,1);
    y_rls_offline_mln = cell(N,1);
    y_rls_online_mln = cell(N,1);
    y_rla_mln = cell(N,1);   
    
    y_noise_pseudo_mln = cell(N,1);
    y_noise_rls_offline_mln = cell(N,1);
    y_noise_rls_online_mln = cell(N,1);
    y_noise_rla_mln = cell(N,1);   

    snr_pseudo_mln = cell(N,1);
    snr_rls_offline_mln = cell(N,1);
    snr_rls_online_mln = cell(N,1);
    snr_rla_mln = cell(N,1);

    % The initial assumption is that h is passthrough
    h_pseudo_mln{1} = h_passthrough;
    h_rls_offline_mln{1} = h_passthrough;
    h_rls_online_mln{1} = h_passthrough;
    h_rla_mln{1} = h_passthrough;

    x = xMLN;
    t = (0:length(x)-1)/fs.';

    for ii=2:N+1
        % -------- CALC H ----------------
        % Pseudo
        h_cur = h_pseudo_mln{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_pseudo_mln{ii} = pseudoInverseSimulation(x_predist, Y, gamma);
        
        x_predist_pseudo_mln{ii-1} = x_predist;
        y_pseudo_mln{ii-1} = y;
        y_noise_pseudo_mln{ii-1} = y_noise;
        snr_pseudo_mln{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Offline
        h_cur = h_rls_offline_mln{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_rls_offline_mln{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
        
        x_predist_rls_offline_mln{ii-1} = x_predist;
        y_rls_offline_mln{ii-1} = y;
        y_noise_rls_offline_mln{ii-1} = y_noise;
        snr_rls_offline_mln{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Online
        h_cur = h_rls_online_mln{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        [h_rls_online_mln{ii}, y, y_noise] = rlsOnlineSimulation(x_predist, M, 1:K, ...
            @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, noise_signal{ii-1}), h_cur, M_true);
        
        x_predist_rls_online_mln{ii-1} = x_predist;
        y_rls_online_mln{ii-1} = y;
        y_noise_rls_online_mln{ii-1} = y_noise;
        snr_rls_online_mln{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLA
%         h_cur = h_rla_mln{ii-1};
%         x_predist = applyVolterra(x, h_cur, M, 1:K);
%         h_rla_mln{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
    end

    %% --- TWO TONE -------------
    N = 5;  % Repeat 5 times

    h_pseudo_two_tone = cell(N+1,1);
    h_rls_offline_two_tone = cell(N+1,1);
    h_rls_online_two_tone = cell(N+1,1);
    h_rla_two_tone = cell(N+1,1);
    
    x_predist_pseudo_two_tone = cell(N,1);
    x_predist_rls_offline_two_tone = cell(N,1);
    x_predist_rls_online_two_tone = cell(N,1);
    x_predist_rla_two_tone = cell(N,1);  
    
    y_pseudo_two_tone = cell(N,1);
    y_rls_offline_two_tone = cell(N,1);
    y_rls_online_two_tone = cell(N,1);
    y_rla_two_tone = cell(N,1);   
    
    y_noise_pseudo_two_tone = cell(N,1);
    y_noise_rls_offline_two_tone = cell(N,1);
    y_noise_rls_online_two_tone = cell(N,1);
    y_noise_rla_two_tone = cell(N,1);   

    snr_pseudo_two_tone = cell(N,1);
    snr_rls_offline_two_tone = cell(N,1);
    snr_rls_online_two_tone = cell(N,1);
    snr_rla_two_tone = cell(N,1);

    % The initial assumption is that h is passthrough
    h_pseudo_two_tone{1} = h_passthrough;
    h_rls_offline_two_tone{1} = h_passthrough;
    h_rls_online_two_tone{1} = h_passthrough;
    h_rla_two_tone{1} = h_passthrough;

    x = xTwoTone;
    t = (0:length(x)-1)/fs.';

    for ii=2:N+1
    % -------- CALC H ----------------
        % Pseudo
        h_cur = h_pseudo_two_tone{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_pseudo_two_tone{ii} = pseudoInverseSimulation(x_predist, Y, gamma);
        
        x_predist_pseudo_two_tone{ii-1} = x_predist;
        y_pseudo_two_tone{ii-1} = y;
        y_noise_pseudo_two_tone{ii-1} = y_noise;
        snr_pseudo_two_tone{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Offline
        h_cur = h_rls_offline_two_tone{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        y = f_filter(t, x_predist);
        y_noise = y + addNoise(1, numel(y), noise_signal{ii-1});
        Y = getXMatrix(y_noise, M, 1:K);
        h_rls_offline_two_tone{ii} = rlsOfflineSimulation(x_predist, Y, M, 1:K, h_cur);
        
        x_predist_rls_offline_two_tone{ii-1} = x_predist;
        y_rls_offline_two_tone{ii-1} = y;
        y_noise_rls_offline_two_tone{ii-1} = y_noise;
        snr_rls_offline_two_tone{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLS Online
        h_cur = h_rls_online_two_tone{ii-1};
        x_predist = applyVolterra(x, h_cur, M, 1:K);
        [h_rls_online_two_tone{ii}, y] = rlsOnlineSimulation(x_predist, M, 1:K, ...
            @(x, xi) f_filter2(t, x, xi), @(n) addNoise(n, 1, noise_signal{ii-1}), h_cur, M_true);
        
        x_predist_rls_online_two_tone{ii-1} = x_predist;
        y_rls_online_two_tone{ii-1} = y;
        y_noise_rls_online_two_tone{ii-1} = y_noise;
        snr_rls_online_two_tone{ii-1} = norm(y)./norm(noise_signal{ii-1});

        % RLA
%         h_cur = h_rla_two_tone{ii-1};
%         x_predist = applyVolterra(x, h_cur, M, 1:K);
%         h_rla_two_tone{ii} = rlaRefineSimulation(x_predist, M, 1:K, @(x) f_filter(t, x));
    end
    %% MEASURE PERFORMANCE
    total_sim_time = toc;
    disp(total_sim_time)

    %% SAVE RESULTS
    file_name = ['sim_v302_', theName, '_', datestr(now, 'yyyymmddTHHMMSS')];
    save([file_name, '_full.mat'])
end

%% AUX FUNCTIONS
function result = addNoise(n, N, noise_sig)
    result = noise_sig(n:n+N-1);
end

function result = makeNoise(N, sigma, mu)
        if nargin < 3
            mu = 0;
        end
        if nargin < 2
            sigma = 0.005/3;
        end
        result = sigma*randn([N, 1]) + mu;
end