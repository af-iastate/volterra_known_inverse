function [h, y] = rlaRefineSimulation(x, M, p, f_filter, numRefine)
    
    if nargin < 5
        numRefine = 100;
    end
    
    y = f_filter(x);
    ek = y - x;
    d_prev = x - ek;
    d = d_prev;

    for ii = 1:numRefine

        
        y = (f_filter(d));
        
%         subplot(211)        
%         plot(d)
%         title('x hat')
%         subplot(212)
%         plot(ek)
%         title('error')
% %         plot(y0)
% %         hold on;
% %         plot(y1)
%         hold off;

        
        ek = (y - x);
        d = d - 0.3*ek;
        
%         norm(ek)/norm(x)
        
%         h = pinv(X)*d;
%         stem(h)
%         pause(0.05)
        
        
    end

    X = getXMatrix(x, M, p);
    h = pinv(X)*d;

end