function analyze_increasingNoiseRound(fileNames, n)
    if nargin<2
        n = 6;
    end
    
    snrs = zeros(5,1);
    hTrues = cell(5,1); % should all be the same
    hPseudos = cell(5,1);
    hRlsOfflines = cell(5,1);
    hRlsOnlines = cell(5,1);

    hPseudoErr = zeros(5,1);
    hRlsOfflineErr = zeros(5,1);
    hRlsOnlineErr = zeros(5,1);

    x_lbls = cell(5,1);
    f_db = @(x) 20*log10(x);

    for i=1:5
        fName = fileNames{i};
        result = load(fName, 'h_true', 'h_pseudo_two_tone', ...
            'h_rls_offline_two_tone', 'h_rls_online_two_tone', ...
            'snr_pseudo_two_tone', 'snr_rls_offline_two_tone', ...
            'snr_rls_online_two_tone');
        hTrues{i} = result.h_true;
        hPseudos{i} = result.h_pseudo_two_tone{n};
        hRlsOfflines{i} = result.h_rls_offline_two_tone{n};
        hRlsOnlines{i} = result.h_rls_online_two_tone{n};

        hPseudoErr(i) = norm(hPseudos{i} - hTrues{i})./norm(hTrues{i});
        hRlsOfflineErr(i) = norm(hRlsOfflines{i} - hTrues{i})./norm(hTrues{i});
        hRlsOnlineErr(i) = norm(hRlsOnlines{i} - hTrues{i})./norm(hTrues{i});

        running_sum = 0;
        for jj=1:5
            running_sum =  running_sum+(result.snr_pseudo_two_tone{jj} + ...
            result.snr_rls_offline_two_tone{jj} + ...
            result.snr_rls_online_two_tone{jj})/3;
        end
        snrs(i) = running_sum/5;
        x_lbls{i} = sprintf('%0.3g', f_db(snrs(i))) ;
    end
    
    fprintf('Plot Data\n')
    A = [(hPseudoErr).';(hRlsOfflineErr).';(hRlsOnlineErr).'];
    disp(A)

    plot((hPseudoErr), '--o')
    hold on;
    plot((hRlsOfflineErr), '--x')
    plot((hRlsOnlineErr), '--*')
    
    xticks(1:5);
    xticklabels(x_lbls)
    legend({'Pseudo Inverse', 'RLS Offline', 'RLS Online'},'location','northwest')
    xlabel('SNR in dB')
    axis([.75, 5.25, min(ylim), max(ylim)+diff(ylim)*0.333])
end




