% Cmp 3rd order inverse
clearvars; clc; close all;
% Reset Rng
rng(23)

% System Parameters
mDelays = 3;
volterraOrders = [1,3];
fs = 125e6;
durationTwoToneAdditional = [0.01 0.01]*1e-3;

% TT
bucketSep = 1/durationTwoToneAdditional(1);
numBuckets = floor(fs/2/bucketSep);
poolSize = floor(16*5*2.5);
freqsPool = round(bucketSep)*floor(linspace(1,numBuckets-3, poolSize));
[freqsTwoTone, xTwoTone, freqsTwoToneAlt, xTwoToneAlt, splitIdx] = ...
    makeTonePulseDataSets(freqsPool, durationTwoToneAdditional, fs, 5, 16);

% TODO rename everything TT
% % ------------------ Get X -----------------------
% freqsTwoTone = [7 7.2; 7 7.4; 7 7.6; 7 7.8; 7 8; ...
%     7.2 7.4; 7.2 7.6; 7.2 7.8; 7.2 8; ...
%     7.4 7.6; 7.4 7.8; 7.4 8; ...
%     7.6 7.8; 7.6 8; 7.8 8]*1e5;
% 
% 
% freqsTwoToneAlt = ([...
%     5 7.2; 5 7.4; 5 7.6; 5 7.8; 5 8; ...
%     5.2 7.4; 5.2 7.6; 5.2 5.8; 5.2 8; ...
%     5.4 5.6; 5.4 5.8; 5.4 8; ...
%     5.6 5.8; 5.6 8; 5.8 8]*5)*1e6;
% 
% 
% durationTwoToneAdditional = [0.01 0.01]*1e-3;
% 
% [xTwoTone, splitsTwoTone] = getInputs( ...
%     freqsTwoTone, ...
%     durationTwoToneAdditional, ...
%     fs ...
% );
% 
% [xTwoToneAlt, splitsTwoToneAlt] = getInputs( ...
%     freqsTwoToneAlt, ...
%     durationTwoToneAdditional, ...
%     fs ...
% );

nX = length(xTwoTone);
qLevels = linspace(-1, 1, 3);
xMLN43 = makeMLN(4, 3, nX);
xMLN63 = makeMLN(6, 3, nX);
xMLN49 = makeMLN(4, 9, nX);
xGauss = randn(nX, 1);
xMLN = xGauss;

noiseLvls = cell([5, 1]);
noiseLvls{1} = 1e-2/3;
noiseLvls{2} = 5e-3/3;
noiseLvls{3} = 1e-4/3;
noiseLvls{4} = 1e-5/3;
noiseLvls{5} = 0;


save('data.mat', 'xTwoTone', 'xMLN43', 'xMLN63', 'xMLN49', 'xMLN', 'xGauss', ...
    'noiseLvls', 'freqsTwoTone', 'xTwoToneAlt', 'freqsTwoToneAlt' ...
)

%%
function [x, splitIdx] = getInputs(freqsTwoTone, durationTwoToneAdditional, fs)
    [x, ndivs] = genPulsedToneSeries( ...
        freqsTwoTone, ...
        durationTwoToneAdditional, ...
        'dutyCycle', 0.75, ...
        'fs', fs ...
    );
    N = floor(length(x)/ndivs);
    splitIdx = (1:N) * ndivs;
end

function x = makeMLN(M, K, n)
    x_seq = permn(linspace(-1, 1, K+1), M);
    x_seq = reshape(x_seq.', [numel(x_seq), 1]);
    nRepeats = floor(n/numel(x_seq));
    bound = (nRepeats*numel(x_seq));
    x = zeros(n, 1);
    x(1:bound) = repmat(x_seq, nRepeats, 1);
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
