clc; close all; clearvars;

% L(4, 3, 0) == L3(4, 3, 0) 
% L(4, 3, 1) == L3(4, 3, 1) 
% L(4, 3, 2) == L3(4, 3, 2) 
% L(4, 3, 3) == L3(4, 3, 3) 

phi(5,[3,3,3])
reverseIdx(5, [3,3,3])

function p = L(M, K, jj)
    if K == 1
        p = M - jj;
        return;
    end
    p = 0;
    for ii=jj:M-1
        p = p + L(M, K-1, ii);
    end
end



% function p = L2(M, K, z)
%     if nargin < 3
%         z = M - 1;
%     end
%     if K == 0 || z==0
%         p = 1;
%         return
%     end
%     p = L2(M,K-1,z) + L2(M,K,z-1);
% end

function p = L2(K, z)
%     if nargin < 3
%         z = M - 1;
%     end
    if K == 0 || z==0
        p = 1;
        return
    end
    p = L2(K-1,z) + L2(K,z-1);
end




function p = L3(M, K, jj)
    z = M-jj;
    p = nchoosek(z+K-1, K);
end

function result = phi(M, tau, carryTerm)
    if nargin<3
        carryTerm = 0;
    end
    
    K = length(tau);
    if K == 1
        result = tau(1) - carryTerm;
        return;
    end
    
    subsetSize = @(M, K, tau0) nchoosek(M-tau0+K-1, K);
    result = subsetSize(M,K, carryTerm) - subsetSize(M,K, tau(1));
    
    result = result + phi(M, tau(2:end), tau(1));
end

function result = reverseIdx(M, tau)
    K = numel(tau);
    result = 0;
    result = result + sum(arrayfun(@(K) ...
        nchoosek(K+M-1, K), ...
        1:K-1 ...
    ));
    result = result + phi(M, tau);
end