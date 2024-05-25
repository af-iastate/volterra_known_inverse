clc; clearvars; close all;
fileName = getMostRecentSimulation('lpdpn2_g75');
result = load(fileName, 'h_true', 'h_pseudo_two_tone', ...
    'h_rls_offline_two_tone', 'h_rls_online_two_tone' ...
    );

figure(1)
resize_typeset(3, 1.7*1.9)
% 
% tix = [-25,0,25,50,75];
% axx = [-.5, .5, 260, 320];

subplot(311)
stem(result.h_pseudo_two_tone{6} - result.h_true)
axis([1 numel(result.h_pseudo_two_tone{6}) ylim])
% xticks([])
yticks(linspace(-10, 10, 5))

subplot(312)
stem(result.h_rls_offline_two_tone{6} - result.h_true)
axis([1 numel(result.h_pseudo_two_tone{6}) -.01 .02])
% xticks([])

subplot(313)
stem(result.h_rls_online_two_tone{6} - result.h_true)
axis([1 numel(result.h_pseudo_two_tone{6}) ylim])


print('overFitH', outName, '-r900', '-dpdf');