function simulateRidgeRegressionRound(...
        M, K, theName, f_filter,  ...
        xTwoTone, noiseLvls, fs, gamma)
    
    h_under = zeros(size(voltVecGen(1:M, 1:K)));
    h_passthrough = zeros(size(h_under));
    h_passthrough(1) = 1;
    
    N = 5;
    nX = length(xTwoTone);
    noise_signal = cell(numel(noiseLvls), N);
    rng(23)
    for ii=1:N
        noiseLvl = noiseLvls{ii};
        for jj=1:numel(noiseLvls)
            noise_signal{ii, jj} = makeNoise(nX, noiseLvl);
        end
    end
    
    h_pseudo_two_tone = cell(numel(noiseLvls), N+1);
    kappas = zeros(numel(noiseLvls), N);
    
    
    x = xTwoTone;
    t = (0:length(x)-1)/fs.';
%     gamma = 2;
    for ii=1:numel(noiseLvls)
        h_pseudo_two_tone{ii, 1} = h_passthrough;
        for jj=2:N+1
            % -------- CALC H ----------------
                % Pseudo
                h_cur = h_pseudo_two_tone{ii, jj-1};
                x_predist = applyVolterra(x, h_cur, M, 1:K);
                y = f_filter(t, x_predist);
                y_noise = y + addNoise(1, numel(y), noise_signal{ii, jj-1});
                Y = getXMatrix(y_noise, M, 1:K);
                h = rememberPastCalls(@pseudoInverseSimulation, 'pseudoInverseSimulation', x_predist, Y, gamma);
                h_pseudo_two_tone{ii, jj} = h;
                kappas(ii,jj-1) = rememberPastCalls(@cond, 'cond', Y);
                disp(kappas(ii,jj-1))
%                 disp(kappas(ii,jj-1))
                stem(h)
                pause(0.01)
        end
    end
    file_name = ['sim_v302_', theName, '_', datestr(now, 'yyyymmddTHHMMSS')];
    save([file_name, '.mat'])
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