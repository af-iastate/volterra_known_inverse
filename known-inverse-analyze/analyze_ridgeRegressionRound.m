function analyze_ridgeRegressionRound(fileNames, h_true, byIteration, offset)
    if nargin < 4
        offset = 1;
    end
    if nargin < 3
        byIteration = true;
    end

    err = zeros(numel(fileNames),5);
    for ii=1:numel(fileNames)
        result = load(fileNames{ii}, 'h_pseudo_two_tone', 'kappas', 'gamma');
        h_hats = result.h_pseudo_two_tone; kappas = result.kappas; gamma = result.gamma;
        if ~byIteration
            tmp = [cell2mat(h_hats(:,1 + offset).') - h_true*ones(1,5)];
        else
            tmp = [cell2mat(h_hats(offset,2:end)) - h_true*ones(1,5)];
        end
        err(ii,:) = arrayfun(@(jj) norm(tmp(:, jj))./norm(h_true), 1:5);
%         disp(gamma)
    end
    
    if ~byIteration
        err = fliplr(err);
    end
    f_db = @(x) 20.*log10(x);
    
    fprintf('Plot Data\n')
    err = err./norm(h_true);
    disp(err)

    plot(err(1,:), 'o--')
    hold on
    plot(err(2,:), 'x--')
    plot(err(3,:), '*--')

    xticks([1, 2, 3, 4, 5])

    % xlabel('$n$th Iteration')
    snr_labels = {'\infty', '100', '75', '45', '35'};
    
    if ~byIteration
        fprintf('At iteration=%d\n', offset)
        xticklabels(snr_labels)
        xlabel('SNR in dB')
    else
        sprintf('At SNR=%s', snr_labels{6-offset})
        xlabel('$n$th iteration', 'interpreter', 'latex')
    end
%     ylabel('dB')
    axis([.75, 5.25, min(ylim), max(ylim)+diff(ylim)*0.333])
    l = legend({'$\gamma=10^{-4}$', '$\gamma=2$', '$\gamma=75$'}, 'location', 'northeast', 'orientation', 'horizontal', 'interpreter', 'latex');
    l.EdgeColor='none';
end