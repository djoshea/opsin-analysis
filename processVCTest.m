function out = processVCTest( fnameOrData, varargin )
% computes accessResistance and membraneResistance (in megaohms) from a voltage step voltage clamp protocol

%% Parameter initialization
par.chIName = '';
par.chVCName = '';
par.chINum = 1;
par.chVCNum = 2;
par.postStepSteadyStateBegin = 0.3;
par.voltageDelta = []; % used to override the extracted value
par.showPlots = 1;
par.savePlots = 0;
par.holding = -70;
par.yl = [];
par.protocolParams = []; % for overriding params detected from channel 2
assignargs(par, varargin);

%% Data load if necessary

if(~exist('fnameOrData', 'var'))
    fnameOrData = '';
end

if(isempty(chIName) || isempty(chVCName))
    chMatchRegexes = {}; % request them numerically
else
    chMatchRegexes = {chIName, chVCName};
end
    
abfData = abfDataLoad(fnameOrData, chMatchRegexes);
sis = abfData.sis;
time = abfData.time; 
nTraces = abfData.nTraces;
chI = abfData.channels{1};
chVC = abfData.channels{2};

out = [];

% Find voltage step times
if(isempty(protocolParams))
    % detect the protocol params from ch 2
    [ stepInds stepLevels preStepLevel] = stepDetect(chVC);
    stepStartInds = cellfun(@(inds) inds(1), stepInds, 'UniformOutput', false);
    minStartInd = min(cellfun(@(inds) inds(1), stepInds));
    stepWidthMin = min(cellfun(@(inds) inds(2) - inds(1) - 1, stepInds));
else
    stepStartInds = mat2cell(protocolParams.stepStartInd*ones(nTraces,1), ones(nTraces,1));
    minStartInd = protocolParams.stepStartInd;
    stepWidthMin = protocolParams.stepWidth;
end

% pass out detected protocol parameters
out.stepStartInd = minStartInd;
out.stepWidth = stepWidthMin;

% Find average pre-pulse baseline region
regions = getPeriEventRegions(chI, 1, [0 minStartInd-2]);
meanBaseline = mean(cell2mat(regions));

% Average the region from each trace immediately following the command step
[regions, ~, relativeTime] = getPeriEventRegions(chI, stepStartInds, [0 stepWidthMin], 'sis', sis);
meanTrace = mean(cell2mat(regions'), 2);

% compute the mean command step height
if(isempty(voltageDelta))
    % use the measured values from the traces
    vDelta = mean(cellfun(@(levels, pre) levels(1) - pre, stepLevels, preStepLevel));
else
    % override the measured values (in case something is wrong with the gain)
    vDelta = voltageDelta;
end
% vDelta = 10;
[val iPeak] = min(meanTrace);
currentPeak = val -  meanBaseline;
steadyStateSkip = round(postStepSteadyStateBegin / sis);
currentSteadyState = mean(meanTrace(steadyStateSkip:end)) - meanBaseline;

% Isteady = Vdelta / (Rmembrane + Raccess);
% volts*10^-3 / amps * 10^-12 = ohms / 10^9 = Mohms / 10^3
totalResistance = vDelta / currentSteadyState * 1000;
out.accessResistance = vDelta / currentPeak * 1000;
out.membraneResistance = totalResistance - out.accessResistance;
out.totalResistance = totalResistance;
out.meanBaseline = meanBaseline;
out.voltageDelta = vDelta;
% calculate resting potential (I = 0): MOhms * pA / 1000 = mV
out.resting = holding - 10^-3*meanBaseline*totalResistance; 

if(showPlots)
    figure(1), clf;
    plot(relativeTime*1000, meanTrace-meanBaseline,'k-')
    hold on
    plot(relativeTime(iPeak)*1000, currentPeak, 'r.');
    plot([0 max(relativeTime)*1000], [currentSteadyState currentSteadyState], 'r-');
    box off;
    xlim([0 600]);
    xlabel('Time Relative to Step (ms)');
    ylabel('Current - Baseline (pA)');
    if(~isempty(yl))
        ylim(yl);
    end
end

end

