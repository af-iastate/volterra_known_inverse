clc; clearvars; close all;
a1 = 0.2; 
a3 = 0.3;
a5 = 0.51;
h = [a5, 0, a3, 0, a1, 0];
M = 6; 
K = 1;

qthPower(h, M, K, 3)