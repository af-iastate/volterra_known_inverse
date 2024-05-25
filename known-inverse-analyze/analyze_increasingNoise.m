clc; close all; clearvars;


%% ------------------- N Iteration = 5 ------------------------
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
analyze_increasingNoiseRound(fileNames);
outName ='sysOneSnrInc';
print(gcf, outName, '-r900', '-dpdf');

% Sys 2
figure(2)
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('bpspn5');
fileNames{2} = getMostRecentSimulation('bpspn4');
fileNames{3} = getMostRecentSimulation('bpspn3');
fileNames{4} = getMostRecentSimulation('bpspn2');
fileNames{5} = getMostRecentSimulation('bpspn1');
analyze_increasingNoiseRound(fileNames);
% axis([xlim, 0, 1250])
q=gcf;
q.Children(1).Location = 'northwest';
outName ='sysTwoSnrInc';
print(gcf, outName, '-r900', '-dpdf');

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

analyze_increasingNoiseRound(fileNames);
outName ='sysThreeSnrInc';
print(gcf, outName, '-r900', '-dpdf');


%% ------------------- N Iteration = 1 ------------------------
% Sys 1
figure(4)
cla
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('lpspn5');
fileNames{2} = getMostRecentSimulation('lpspn4');
fileNames{3} = getMostRecentSimulation('lpspn3');
fileNames{4} = getMostRecentSimulation('lpspn2');
fileNames{5} = getMostRecentSimulation('lpspn1');
analyze_increasingNoiseRound(fileNames, 2);
axis([.75, 5.25, min(ylim), max(ylim)+diff(ylim)*0.66555])
outName ='sysOneSnrIncFirst';
print(gcf, outName, '-r900', '-dpdf');

% Sys 2
figure(5)
cla
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('bpspn5');
fileNames{2} = getMostRecentSimulation('bpspn4');
fileNames{3} = getMostRecentSimulation('bpspn3');
fileNames{4} = getMostRecentSimulation('bpspn2');
fileNames{5} = getMostRecentSimulation('bpspn1');
analyze_increasingNoiseRound(fileNames, 2);
axis([.75, 5.25, min(ylim), max(ylim)+diff(ylim)*0.333])
% axis([xlim, 0, 1250])
q=gcf;
q.Children(1).Location = 'northwest';
outName ='sysTwoSnrIncFirst';
print(gcf, outName, '-r900', '-dpdf');

% Sys 3
figure(6)
cla
resize_typeset(3, 1.7)
fileNames = cell(5,1);
% higher filenames = less noise
fileNames{1} = getMostRecentSimulation('lpdpn5');
fileNames{2} = getMostRecentSimulation('lpdpn4');
fileNames{3} = getMostRecentSimulation('lpdpn3');
fileNames{4} = getMostRecentSimulation('lpdpn2');
fileNames{5} = getMostRecentSimulation('lpdpn1');

analyze_increasingNoiseRound(fileNames, 2);

outName ='sysThreeSnrIncFirst';
print(gcf, outName, '-r900', '-dpdf');