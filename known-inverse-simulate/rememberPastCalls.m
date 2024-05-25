function result = rememberPastCalls(f_func, idStr, varargin)
    persistent allCallMemory;
    if isempty(allCallMemory)
        allCallMemory = containers.Map;
    end
    
    if ~isKey(allCallMemory, idStr)
        allCallMemory(idStr) = containers.Map;
    end
    callMemory = allCallMemory(idStr);
    
    argsHash = DataHash(varargin);
    if isKey(callMemory, argsHash)
        result = callMemory(argsHash);
        return;
    end
    
    % otherwise, actually perform the call
    result = f_func(varargin{:});
    % and store the result
    callMemory(argsHash) = result;
    allCallMemory(idStr) = callMemory;
    
end