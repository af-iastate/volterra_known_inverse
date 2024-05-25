function analyze_increasingNoiseConditionNumRound(fileNames)

    snrs = zeros(5,1);
    
    yPseudos = cell(5,1);
    yRlsOfflines = cell(5,1);
    yRlsOnlines = cell(5,1);

    kappaPseudos = zeros(5,1);
    kappaRlsOfflines = zeros(5,1);
    kappaRlsOnlines = zeros(5,1);

    x_lbls = cell(5,1);
    f_db = @(x) 20*log10(x);

    for i=1:5
        fName = fileNames{i};
        result = load(fName, 'h_true', 'y_noise_pseudo_two_tone', ...
            'y_noise_rls_offline_two_tone', 'y_noise_rls_online_two_tone', ...
            'x_predist_rls_offline_two_tone', 'x_predist_rls_online_two_tone', ...
            'snr_pseudo_two_tone', 'snr_rls_offline_two_tone', ...
            'snr_rls_online_two_tone', 'M', 'K');
        M = result.M;
        K = result.K;
        
%         x_predist_offline_two_tone = result.x_predist_rls_offline_two_tone{1};
%         x_predist_online_two_tone = result.x_predist_rls_online_two_tone{1};
        
%         figure(2)
%         plot(x_predist_offline_two_tone)
%         hold on
%         plot(x_predist_online_two_tone)
%         hold off
        
        
        yPseudos{i} = result.y_noise_pseudo_two_tone{1};
%         yRlsOfflines{i} = result.y_noise_rls_offline_two_tone{1};
%         yRlsOnlines{i} = result.y_noise_rls_online_two_tone{1};
        
%         L = numel(yRlsOfflines{i}) - numel(yRlsOnlines{i});
%         if L > 0
%             yRlsOnlines{i} = [yRlsOnlines{i}; zeros(L, 1)];
%         end
        
        YPseudo = getXMatrix(yPseudos{i}, M, 1:K);
%         YRLSOffline = getXMatrix(yRlsOfflines{i}, M, 1:K);
%         YRLSOnline = getXMatrix(yRlsOnlines{i}, M, 1:K);

        kappaPseudos(i) = rememberPastCalls(@cond, 'cond', (YPseudo));
%         kappaRlsOfflines(i) = rememberPastCalls(@cond, 'cond', (YRLSOffline));
%         kappaRlsOnlines(i) = rememberPastCalls(@cond, 'cond', (YRLSOnline));

        running_sum = 0;
        for jj=1:5
            running_sum =  running_sum+(result.snr_pseudo_two_tone{jj} + ...
            result.snr_rls_offline_two_tone{jj} + ...
            result.snr_rls_online_two_tone{jj})/3;
        end
        snrs(i) = running_sum/5;
        x_lbls{i} = sprintf('%0.3g', f_db(snrs(i))) ;
    end

    plot((kappaPseudos), '--o')
%     hold on;
%     plot((kappaRlsOfflines), '--x')
%     plot((kappaRlsOnlines), '--*')

    xticks(1:5);
    xticklabels(x_lbls)
%     legend({'Pseudo Inverse', 'RLS Offline', 'RLS Online'},'location','northwest')
    xlabel('SNR in dB')
end




