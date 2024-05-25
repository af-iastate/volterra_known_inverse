clc; close all; clearvars;


tau = generateTauUpperTri(4,4)
offset = tauToOffsetWithinOrder(4,tau)
idx = tauToFlatIndex(4,tau)
