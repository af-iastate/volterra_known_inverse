clc; close all; clearvars;

% Sys 1
figure(1)
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('lpspn5');
fileNames{2} = getMostRecentSimulation('lpspn4');
fileNames{3} = getMostRecentSimulation('lpspn3');
fileNames{4} = getMostRecentSimulation('lpspn2');
fileNames{5} = getMostRecentSimulation('lpspn1');
analyze_increasingNoiseConditionNumRound(fileNames);
outName ='sysOneSnrIncKappa';
print(gcf, outName, '-r900', '-dpdf');
% 
% % Sys 2
figure(2)
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('bpspn5');
fileNames{2} = getMostRecentSimulation('bpspn4');
fileNames{3} = getMostRecentSimulation('bpspn3');
fileNames{4} = getMostRecentSimulation('bpspn2');
fileNames{5} = getMostRecentSimulation('bpspn1');
analyze_increasingNoiseConditionNumRound(fileNames);
% % axis([xlim, 0, 1250])
q=gcf;
% q.Children(1).Location = 'southwest';
outName ='sysTwoSnrIncKappa';
print(gcf, outName, '-r900', '-dpdf');
% 
% Sys 3
figure(3)
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('lpdpn5');
fileNames{2} = getMostRecentSimulation('lpdpn4');
fileNames{3} = getMostRecentSimulation('lpdpn3');
fileNames{4} = getMostRecentSimulation('lpdpn2');
fileNames{5} = getMostRecentSimulation('lpdpn1');

analyze_increasingNoiseConditionNumRound(fileNames);
outName ='sysThreeSnrIncKappa';
print(gcf, outName, '-r900', '-dpdf');