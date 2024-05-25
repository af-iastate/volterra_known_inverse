function analyzeRound(dataFileName, f_filter)
    
    if ~isfile('data.mat')
        make_test_data()
    end
    load('data.mat', 'freqsTwoTone', 'xTwoToneAlt', 'freqsTwoToneAlt')
        
    load( ...
        dataFileName, ...
        'h_true', 'M', 'K', 'N', 'Y', 'fs',  ...
        'h_cur', 'h_passthrough', ...
        'h_pseudo_g', 'h_pseudo_mln', 'h_pseudo_two_tone', ...
        'h_rla_g', 'h_rla_mln', 'h_rla_two_tone', ...
        'h_rls_offline_g', 'h_rls_offline_mln', 'h_rls_offline_two_tone', ...
        'h_rls_online_g', 'h_rls_online_mln', 'h_rls_online_two_tone', ...
        't', 'noise_signal',...
        'xGauss', 'xMLN', 'xTwoTone', ...
        'h_pseudo_g', 'h_rls_offline_g', 'h_rls_online_g', 'h_rla_g', ...
        'h_pseudo_mln', 'h_rls_offline_mln', 'h_rls_online_mln', 'h_rla_mln', ...
        'h_pseudo_two_tone', 'h_rls_offline_two_tone', 'h_rls_online_two_tone', 'h_rla_two_tone', ...
        'x_predist_pseudo_g', 'x_predist_rls_offline_g', 'x_predist_rls_online_g', 'x_predist_rla_g', ...
        'x_predist_pseudo_mln', 'x_predist_rls_offline_mln', 'x_predist_rls_online_mln', 'x_predist_rla_mln', ...
        'x_predist_pseudo_two_tone', 'x_predist_rls_offline_two_tone', 'x_predist_rls_online_two_tone', 'x_predist_rla_two_tone', ...
        'y_pseudo_g', 'y_rls_offline_g', 'y_rls_online_g', 'y_rla_g', ...
        'y_pseudo_mln', 'y_rls_offline_mln', 'y_rls_online_mln', 'y_rla_mln', ...
        'y_pseudo_two_tone', 'y_rls_offline_two_tone', 'y_rls_online_two_tone', 'y_rla_two_tone', ...
        'snr_pseudo_g', 'snr_rls_offline_g', 'snr_rls_online_g', 'snr_rla_g', ...
        'snr_pseudo_mln', 'snr_rls_offline_mln', 'snr_rls_online_mln', 'snr_rla_mln', ...
        'snr_pseudo_two_tone', 'snr_rls_offline_two_tone', 'snr_rls_online_two_tone', 'snr_rla_two_tone', ...
        'y_noise_pseudo_g', 'y_noise_rls_offline_g', 'y_noise_rls_online_g', 'y_noise_rla_g', ...
        'y_noise_pseudo_mln', 'y_noise_rls_offline_mln', 'y_noise_rls_online_mln', 'y_noise_rla_mln', ...
        'y_noise_pseudo_two_tone', 'y_noise_rls_offline_two_tone', 'y_noise_rls_online_two_tone', 'y_noise_rla_two_tone' ...
    )

    x_predist_two_tone_true = applyVolterra(xTwoTone, h_true, M, 1:K);
    
    if numel(h_pseudo_g{1}) == numel(h_true)
        h_under = h_true;
    else
        h_under = zeros(size(h_pseudo_g{1}));
    end

    for n=6
        figure(1)
        subplot(311)
        analyze_result(h_under, h_pseudo_g{n}, ['Gaussian Pseudo Iteration ', num2str(n-1)])
        subplot(312)
        analyze_result(h_under, h_rls_offline_g{n}, ['Gaussian RLS Off Iteration ', num2str(n-1)])
        subplot(313)
        analyze_result(h_under, h_rls_online_g{n}, ['Gaussian RLS On Iteration ', num2str(n-1)])
    end

    for n=6
        figure(2)
        subplot(311)
        analyze_result(h_under, h_pseudo_mln{n}, ['MLN Pseudo Iteration ', num2str(n-1)])
        subplot(312)
        analyze_result(h_under, h_rls_offline_mln{n}, ['MLN RLS Off Iteration ', num2str(n-1)])
        subplot(313)
        analyze_result(h_under, h_rls_online_mln{n}, ['MLN RLS On Iteration ', num2str(n-1)])
    end

    for n=6
        figure(3)
        subplot(311)
        analyze_result(h_under, h_pseudo_two_tone{n}, ['5T Pseudo Iteration ', num2str(n-1)])
        subplot(312)
        analyze_result(h_under, h_rls_offline_two_tone{n}, ['5T RLS Off Iteration ', num2str(n-1)])
        subplot(313)
        analyze_result(h_under, h_rls_online_two_tone{n}, ['5T RLS On Iteration ', num2str(n-1)])
    end

    %%
    
    if (norm(h_under) > 0 && numel(h_pseudo_g{1}) == numel(h_true))
        figure(4)
        nNf = @(f,n) (arrayfun(@(ii) norm(f{ii}-h_true), n));

        % 
        subplot(311)
        stem(nNf(h_pseudo_g, 2:6))
        title('$||\hat{h}-h||_2$ Gaussian', 'interpreter', 'latex')
        hold on
        stem(nNf(h_rls_offline_g, 2:6))
        stem(nNf(h_rls_online_g, 2:6))
        legend({'Gaussian Pseudo', 'Gaussian RLS Off', 'Gaussian RLS On'})
        xlabel('nth iteration')

        subplot(312)
        stem(nNf(h_pseudo_mln, 2:6))
        title('$||\hat{h}-h||_2$ MLN', 'interpreter', 'latex')
        hold on
        stem(nNf(h_rls_offline_mln, 2:6))
        stem(nNf(h_rls_online_mln, 2:6))
        legend({'MLN Pseudo', 'MLN RLS Off', 'MLN RLS On'})
        xlabel('nth iteration')

        subplot(313)
        stem(nNf(h_pseudo_two_tone, 2:6))
        title('$||\hat{h}-h||_2$ Two Tone', 'interpreter', 'latex')
        hold on
        stem(nNf(h_rls_offline_two_tone, 2:6))
        stem(nNf(h_rls_online_two_tone, 2:6))
        legend({'5T Pseudo', '5T RLS Off', '5T RLS On'})
        xlabel('nth iteration')
    end

    %% ------------------------------------ SFDR ------------------------------

    % ------Gaussian----
    sfdr_pseudo_g = [];
    sfdr_rls_offline_g = [];
    sfdr_rls_online_g = [];

    x = xTwoTone;
    
    figure(5)
    for n=1:6
        % PSEUDO
        a = analyzeSfdr(x, t, f_filter, n, h_pseudo_g{n}, M, K, freqsTwoTone, fs);
        sfdr_pseudo_g = [sfdr_pseudo_g; a];
        hold off

        % RLS OFFLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_offline_g{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_offline_g = [sfdr_rls_offline_g; a];
        hold off

        % RLS ONLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_online_g{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_online_g = [sfdr_rls_online_g; a];
        hold off
    end

    % SFDR Plot
    subplot(311)
    plotSfdr(sfdr_pseudo_g, 'Gaussian Pseudo SFDR')
    subplot(312)
    plotSfdr(sfdr_rls_offline_g, 'Gaussian RLS Offline SFDR')
    subplot(313)
    plotSfdr(sfdr_rls_online_g, 'Gaussian RLS Online SFDR')

    % ------MLN----
    sfdr_pseudo_mln = [];
    sfdr_rls_offline_mln = [];
    sfdr_rls_online_mln = [];

    x = xTwoTone;
    figure(6)
    for n=1:6
        % PSEUDO
        a = analyzeSfdr(x, t, f_filter, n, h_pseudo_mln{n}, M, K, freqsTwoTone, fs);
        sfdr_pseudo_mln = [sfdr_pseudo_g; a];
        hold off

        % RLS OFFLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_offline_mln{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_offline_mln = [sfdr_rls_offline_mln; a];
        hold off

        % RLS ONLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_online_mln{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_online_mln = [sfdr_rls_online_mln; a];
        hold off
    end

    % SFDR Plot
    subplot(311)
    plotSfdr(sfdr_pseudo_mln, 'MLN Pseudo SFDR')
    subplot(312)
    plotSfdr(sfdr_rls_offline_mln, 'MLN RLS Offline SFDR')
    subplot(313)
    plotSfdr(sfdr_rls_online_mln, 'MLN RLS Online SFDR')


    % ------Two Tone----
    sfdr_pseudo_two_tone = [];
    sfdr_rls_offline_two_tone = [];
    sfdr_rls_online_two_tone = [];

    x = xTwoTone;
    figure(7)
    for n=1:6
        % PSEUDO
        a = analyzeSfdr(x, t, f_filter, n, h_pseudo_two_tone{n}, M, K, freqsTwoTone, fs);
        sfdr_pseudo_two_tone = [sfdr_pseudo_two_tone; a];
        hold off

        % RLS OFFLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_offline_two_tone{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_offline_two_tone = [sfdr_rls_offline_two_tone; a];
        hold off

        % RLS ONLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_online_two_tone{n}, M, K, freqsTwoTone, fs);
        sfdr_rls_online_two_tone = [sfdr_rls_online_two_tone; a];
        hold off
    end
    
    sfdr_pseudo_two_tone_alt = [];
    sfdr_rls_offline_two_tone_alt = [];
    sfdr_rls_online_two_tone_alt = [];
    x = xTwoToneAlt;
    for n=1:6
        % PSEUDO
        a = analyzeSfdr(x, t, f_filter, n, h_pseudo_two_tone{n}, M, K, freqsTwoToneAlt, fs);
        sfdr_pseudo_two_tone_alt = [sfdr_pseudo_two_tone_alt; a];
        hold off

        % RLS OFFLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_offline_two_tone{n}, M, K, freqsTwoToneAlt, fs);
        sfdr_rls_offline_two_tone_alt = [sfdr_rls_offline_two_tone_alt; a];
        hold off

        % RLS ONLINE
        a = analyzeSfdr(x, t, f_filter, n, h_rls_online_two_tone{n}, M, K, freqsTwoToneAlt, fs);
        sfdr_rls_online_two_tone_alt = [sfdr_rls_online_two_tone_alt; a];
        hold off
    end

    % SFDR Plot
    tix = [0,50,100,150];
    
    subplot(311)
    plotSfdr(sfdr_pseudo_two_tone, '5T Pseudo SFDR')
    hold on
    plotSfdr(sfdr_pseudo_two_tone_alt, '5T Pseudo SFDR', 'x')
%     axis([xlim, -10, 155])
%     yticks(tix)
    
    subplot(312)
    plotSfdr(sfdr_rls_offline_two_tone, '5T RLS Offline SFDR')
    hold on
    plotSfdr(sfdr_rls_offline_two_tone_alt, '5T RLS Offline SFDR', 'x')
%     axis([xlim, -10, 155])
%     yticks(tix)
    
    subplot(313)
    plotSfdr(sfdr_rls_online_two_tone, '5T RLS Online SFDR')
    hold on
    plotSfdr(sfdr_rls_online_two_tone_alt, '5T RLS Online SFDR', 'x')
%     axis([xlim, -10, 155])
%     yticks(tix)

    %% looking at Two Tone x_predist
    figure(8)
    hold off
    plot(x_predist_two_tone_true)
    hold on
    plot(x_predist_pseudo_two_tone{5})
    plot(x_predist_rls_offline_two_tone{5})
    plot(x_predist_rls_online_two_tone{5})
    hold off
    legend({   ...
        'Truth', ...
        'Pseudo', ...
        'RLS Offline', ...
        'RLS Offline' ...
    }, 'interpreter', 'latex')
    title('$\hat{x}$ for 5 Tone', 'interpreter', 'latex')

    %% looking at Two Tone x_predist with Gauss
    x = xTwoToneAlt;% + 10*noise_sig{1};

    figure(9)
    hold off
    plot(applyVolterra(x, h_true, M, 1:K))
    hold on
    plot(applyVolterra(x, h_pseudo_two_tone{5}, M, 1:K))
    plot(applyVolterra(x, h_rls_offline_two_tone{5}, M, 1:K))
    plot(applyVolterra(x, h_rls_online_two_tone{5}, M, 1:K))
    hold off
    legend({   ...
        'Truth', ...
        'Pseudo' ...
        'RLS Offline', ...
        'RLS Offline' ...
    }, 'interpreter', 'latex')
    title('$\hat{x}$ for 5 Tone for resp to Gauss', 'interpreter', 'latex')

    %% SNR
    f_db = @(x) 20*log10([x{:}]);

    figure(10)

    % gauss
    subplot(311)
    hold off;
    stem(f_db(snr_pseudo_g))
    hold on;
    stem(f_db(snr_rls_offline_g))
    stem(f_db(snr_rls_online_g))
    ylabel('dB')
    title('SNR Guassian')
    xlabel('nth iteration')
    legend({'pseudo', 'rls offline', 'rls online'})

    % mln
    subplot(312)
    hold off;
    stem(f_db(snr_pseudo_mln))
    hold on;
    stem(f_db(snr_rls_offline_mln))
    stem(f_db(snr_rls_online_mln))
    ylabel('dB')
    title('SNR MLN')
    xlabel('nth iteration')
    legend({'pseudo', 'rls offline', 'rls online'})

    % two tone
    subplot(313)
    hold off;
    stem(f_db(snr_pseudo_two_tone))
    hold on;
    stem(f_db(snr_rls_offline_two_tone))
    stem(f_db(snr_rls_online_two_tone))
    ylabel('dB')
    title('SNR 5T')
    xlabel('nth iteration')
    legend({'pseudo', 'rls offline', 'rls online'})



end

%% AUX FUNCTIONS
function analyze_result(h, h_hat, ttl)
    if nargin<3
        ttl = [];
    end
    if isempty(ttl)
        ttl = 'none';
    end
    % figures of merit
    h_err = norm(h - h_hat);
    
    stem(h)
    hold on;
    stem(h_hat)
    hold off;
    
    title_str = sprintf('%s $||\\hat{h}-h||_2 = %g$', ttl, h_err);
    title(title_str, 'interpreter', 'latex')
    
end

function result = analyzeSfdr(x, t, f_filter, n, h, M, K, freqs, fs)
    x_preinv = applyVolterra(x, h, M, 1:K);
    if ~nanappear(x_preinv)
        y = f_filter(t, x_preinv);
    else
        y = zeros(size(x_preinv));
    end
    
    % add splits analyzer
    [~,sfdr] = calcSFDR(y, [1, 1250-0, 2500], freqs, [0,fs/2], fs);
    sfdr = sfdr(:);
    a = n*ones(size(sfdr));
    result = [a, sfdr];
end

function plotSfdr(X, ttl, ss)
    if nargin < 3
        ss = 'o';
    end
    scatter(X(:,1)-1, X(:,2), ss)
    title(ttl)
    xlabel('nth iteration')
    ylabel('dB')
    xticks(unique(0:5))  % todo generalize...
    grid on
end

function r=nanappear(n)
  r=(any(isnan(n(:))));
end