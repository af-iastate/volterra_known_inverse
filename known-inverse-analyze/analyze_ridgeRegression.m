clc; clearvars; close all;
% SYSTEM PARAMETERS
f = [1, 0, 0.1]; % memoriless non-linear (inverse)
 % initial linear sys (inverse)
g = [0.98, 0.1, -0.3, 0.2]; % initial linear sys (inverse)

M = length(g); % initial memory
K = 1; % initial order

% TRUE INVERSE VOLTERRA
[h_true, K] = applyPolymap(g, M, K, f);  % apply polymap
[h_true, K] = applyPolymap(h_true, M, K, f);  % apply polymap again

fileNames = cell(3,1);
% ADD Ridge Simulation for gamma = 0
fileNames{1} = getMostRecentSimulation('lpdp_gen4');
fileNames{2} = getMostRecentSimulation('lpdp_g2');
fileNames{3} = getMostRecentSimulation('lpdp_g75');


figure(1)
resize_typeset(3.5, 1.7)
analyze_ridgeRegressionRound(fileNames, h_true, true, 5)

figure(2)
resize_typeset(3.5, 1.7)
analyze_ridgeRegressionRound(fileNames, h_true, true, 2)

figure(3)
resize_typeset(3.5, 1.7)
analyze_ridgeRegressionRound(fileNames, h_true, true, 1)

figure(4)
resize_typeset(3.5, 1.7)
analyze_ridgeRegressionRound(fileNames, h_true, false, 5)

figure(5)
resize_typeset(3.5, 1.7)
analyze_ridgeRegressionRound(fileNames, h_true, false, 1)

% print('ridgeAnalysis', outName, '-r900', '-dpdf');