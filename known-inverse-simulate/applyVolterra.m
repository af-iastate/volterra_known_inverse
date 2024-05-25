function y = applyVolterra(x, h, M, p)
    X = getXMatrix(x, M, p);
    y = X * h;
end