clc; clearvars; close all;
a1 = 0.2; 
a3 = 0.3;
a5 = 0.51;
a = [a5, 0, a3, 0, a1, 0];
t = linspace(-1, 1, 1001);
YoT = @(t) polyval(a, t);
XoT = getInverseFromPoly(a);

figure
hold on
plot(t, XoT(t))

plot(YoT(t), t)
axis([-1 1 -2 2])