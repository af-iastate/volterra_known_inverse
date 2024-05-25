function [h] = rlsOfflineSimulation(x, Y, M, p, hIn)
    if nargin<5
        hIn = [];
    end
    [h, ~] = rlsVolterraOffline(Y, x, [], M, p, hIn);
end