function [regions indexWindow relativeTime] = getPeriEventRegions( traces, eventInds, window, varargin)
% traces is ntime x ntraces, eventInds is cell{ntraces}, window is [-before after] or cell array of [-before after] for each trace
% specify sis as an argument to get relativeTime in real time units

par.onlyIndexWindow = 0;
par.sis = 1;
assignargs(par, varargin);

nTraces = size(traces,2);
regions = cell(nTraces,1);

if(iscell(window))
    windowList = window;
else
    for iTrace = 1:nTraces
        windowList{iTrace} = window;
    end
end

for iTrace = 1:nTraces
    if(iscell(eventInds))
        events = eventInds{iTrace};
    else
        events = eventInds;
    end
    
    window = windowList{iTrace};
    
    if(min(events) + window(1) < 1 || max(events)+window(2) > size(traces,1))
        error('Specified event index plus surrounding window is out of bounds.');
    end
    
    if(~onlyIndexWindow)
        grabRegion = @(iTrace, iEvent) traces(iEvent+window(1):iEvent+window(2), iTrace);
        regions{iTrace} = cell2mat(arrayfun(@(iEvent) grabRegion(iTrace, iEvent), events, 'UniformOutput', false)');
    end
    
    if(size(events,2) > size(events, 1)); events = events'; end;
    indexWindow{iTrace} = [events+window(1) events+window(2)];
end

relativeTime = [window(1) : window(2)]' * sis;

