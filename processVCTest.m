function out = processVCTest( fnameOrData, varargin )
% computes accessResistance and membraneResistance (in megaohms) from a voltage step voltage clamp protocol

%% Parameter initialization
par.showPlots = 1;
par.savePlots = 0;
par.chIName = 'IN [01]';
par.chVCName = 'IN [23]';
par.postStepSteadyStateBegin = 0.3;
par.plotTrace = 1;
assignargs(par, {});

%% Data load if necessary

if(~exist('fnameOrData', 'var'))
    fnameOrData = '';
end

chMatchRegexes = {chIName, chVCName};
[d sis h chMatches time fname fnamelast] = abfDataLoad(fnameOrData, chMatchRegexes);
nTraces = size(d,3);
chI = chMatches{1};
chVC = chMatches{2};

out = [];

% Find voltage step times
[ stepInds stepLevels preStepLevel] = stepDetect(chVC);
stepStartInds = cellfun(@(inds) inds(1), stepInds, 'UniformOutput', false);
minStartInd = min(cellfun(@(inds) inds(1), stepInds));
stepWidthMin = min(cellfun(@(inds) inds(2) - inds(1) - 1, stepInds));

% Find average pre-pulse baseline region
regions = getPeriEventRegions(chI, 1, [0 minStartInd-2]);
meanBaseline = mean(cell2mat(regions));

% Average the region from each trace immediately following the command step
[regions, ~, relativeTime] = getPeriEventRegions(chI, stepStartInds, [0 stepWidthMin], 'sis', sis);
meanTrace = mean(cell2mat(regions'), 2);

% compute the mean command step height
vDelta = mean(cellfun(@(levels, pre) levels(1) - pre, stepLevels, preStepLevel));

[val iPeak] = min(meanTrace);
currentPeak = val -  meanBaseline;
steadyStateSkip = round(postStepSteadyStateBegin / sis);
currentSteadyState = mean(meanTrace(steadyStateSkip:end)) - meanBaseline;

% Isteady = Vdelta / (Rmembrane + Raccess);
% volts*10^-3 / amps * 10^-12 = ohms / 10^9 = Mohms / 10^3
totalResistance = vDelta / currentSteadyState * 1000;
out.accessResistance = vDelta / currentPeak * 1000;
out.membraneResistance = totalResistance - out.accessResistance;

if(plotTrace)
    figure(1), clf;
    plot(relativeTime*1000, meanTrace-meanBaseline,'k-')
    hold on
    plot(relativeTime(iPeak)*1000, currentPeak, 'r.');
    plot([0 max(relativeTime)*1000], [currentSteadyState currentSteadyState], 'r-');
    box off;
    xlim([0 10]);
    xlabel('Time Relative to Step (ms)');
    ylabel('Current - Baseline (pA)');
end

end

