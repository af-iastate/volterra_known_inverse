clearvars; close all; clc;
len = @length;
fs = 125e6;

h = [0.98, 0.1, -0.3, 0.1, -0.2, 0.01];
h2 = [1, -0.3, 0.3];
M = length(h); 
K = 1;

global f
f = [1, 0, 0.1];

hHat1 = applyPolymap(h, M, 1, f);
M = M;
K = 3;

% hHat2 = hHat1;
hHat2 = applyPolymap(hHat1, M, K, f);
K = 9;

% hTrue = getVolterraFromWH(h, h2, M, 1:K, 1, 1, f);

figure(1)
% stem(hTrue)
% hold on
stem(hHat2)

% forward
G = tf([1 zeros(1, length(h)-1)], h(:).', 1/fs); 
G2 = tf([1 zeros(1, length(h2)-1)], h2(:).', 1/fs); 
f_activation = @thirdOrderActivationFunc;
f_filter = @(t, x) lsim(G, f_activation(f_activation(x)), t);

% input

nX = 10001;
rng(200);
xGauss = randn(nX, 1);
t = (0:length(xGauss)-1)/fs;
x = 0.5*sin(0.2*fs*t/4) + 0.05*xGauss(:)' + 1.2;
% x = 0.05*xGauss;
% x = 0.5*sin(0.2*fs*t/4);
x = x(:);

% y
y = f_filter(t, x);
Y = getXMatrix(y, M, 1:K);

x_approx = Y*hHat2;
x_approx = x_approx((0)+1:end);
x = x(1:end-(0));

nffts = @(x) 20*log10((1/sqrt(len(x)))*abs(fft(x)));
nfft = @(x) ((1/length(x))*abs(fft(x)));
figure(2)

plot(x)
hold on;
plot(x_approx)

% plot(x-x_approx)
% plot(nfft(x-x_approx))
% hold on
% plot(nfft(x_approx))




%%




% nX = 1001;
% rng(200);
% xGauss = randn(nX, 1);
% t = (0:length(xGauss)-1)/fs;
% y = f_filter(t, xGauss);
% 
% % Y = getXMatrix(y, 6, 1:9);
% % hHat3 = pinv(Y)*xGauss(:);
% 
% figure(3)
% stem(hHat3)

function u = thirdOrderActivationFunc(t)
    global f
    a1 = f(1); a3 = f(3);
    %a1 = abs(a1) * (sign(a3) + (a3 == 0));
    if a1 ~= 0 && a3 ~= 0
        tmpf = nthroot(-9*a3^2*t + sqrt(12*a1^3*a3^3 + 81*a3^4*t.^2), 3);
        u = (2/3)^(1/3)*a1./tmpf - tmpf./(2^(1/3)*3^(2/3)*a3);
    elseif a3 ~= 0
        u = nthroot(t./a3, 3);
    elseif a1 ~= 0
        u = (1/a1).*t;
    else
        u = 0.*t;
    end
end