function result = getMostRecentSimulation(typeStr)
    matNames = {dir('*.mat').name}.';
    simRegexStr = 'sim_v(?<version>\d+)_(?<simtype>\w+)_(?<datestr>\w+)_full\.mat';
    matchesFull = regexp(matNames, simRegexStr, 'match');
    matchesFull = matchesFull(cellfun(@(x) ~isempty(x), matchesFull));
    matchesFull = cellfun(@(x) x{1}, matchesFull, 'UniformOutput', false);
    matchesStruct = regexp(matchesFull, simRegexStr, 'names');
    % parse strings to nums/dates/etc
    simCandidates = cell(numel(matchesFull), 1);
    for k=1:numel(matchesFull)
        matchFull = matchesFull{k};
        matchStruct = matchesStruct{k};
        curSim = struct;
        curSim.fileName = matchFull;
        curSim.version = str2double(matchStruct.version);
        curSim.simType = matchStruct.simtype;
        curSim.date = datetime({matchStruct.datestr}, 'InputFormat', 'yyyyMMdd''T''HHmmSS');
        
        simCandidates{k} = curSim;
    end
    
    % filtering
    simCandidates = simCandidates(cellfun(@(x) strcmp(x.simType, typeStr) == 1, simCandidates));
    
    versions = cellfun(@(x) x.version, simCandidates);
    versionMax = max(versions);
    simCandidates = simCandidates(cellfun(@(x) x.version >= versionMax, simCandidates));
    dates = cellfun(@(x) x.date, simCandidates);
    dateMax = max(dates);
    simCandidates = simCandidates(cellfun(@(x) x.date >= dateMax, simCandidates));
    resultSim = simCandidates{1};
    result = resultSim.fileName;
end