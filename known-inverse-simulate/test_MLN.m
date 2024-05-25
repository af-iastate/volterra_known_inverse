clc; close all; clearvars;

x = makeMLN(4, 9, 50000);

function x = makeMLN(M, K, n)
    x_seq = permn(linspace(-1, 1, K+1), M);
    x_seq = reshape(x_seq.', [numel(x_seq), 1]);
    nRepeats = floor(n/numel(x_seq));
    bound = (nRepeats*numel(x_seq));
    x = zeros(n, 1);
    x(1:bound) = repmat(x_seq, nRepeats, 1);
end