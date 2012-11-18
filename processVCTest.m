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
chI = abfData.channels{chINum};
chVC = abfData.channels{chVCNum};

out = [];

% Find voltage step times
if(isempty(protocolParams))
    % detect the protocol params from ch 2
    [ stepInds stepLevels preStepLevel] = stepDetect(chVC);
    stepStartInds = cellfun(@(inds) inds(1), stepInds, 'UniformOutput', false);
    minStartInd = min(cellfun(@(inds) inds(1), stepInds));
    stepWidthMin = min(cellfun(@(inds) inds(2) - inds(1) - 1, stepInds));
else
    if ~isempty(voltageDelta)
    elseif isfield(protocolParams, 'voltageDelta')
        voltageDelta = protocolParams.voltageDelta;
    else
        error('Please provide voltage delta when manually specifying protocol parameters');
    end
        
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
%meanTrace = smooth(meanTrace, 3);

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
out.steadyStateFromBaseline = currentSteadyState;
out.peakFromBaseline = currentPeak;
out.voltageDelta = vDelta;
% calculate resting potential (I = 0): MOhms * pA / 1000 = mV
out.resting = holding - 10^-3*meanBaseline*totalResistance; 

if(showPlots)
    
    windowMs = [-10 50];
    windowInds = round(windowMs / 1000 /sis);
    % Average the region from each trace immediately following the command step
    [regionsExpanded, ~, relativeTimeExpanded] = getPeriEventRegions(chI, stepStartInds, windowInds, 'sis', sis);
    meanTraceExpanded = mean(cell2mat(regionsExpanded'), 2);
    %meanTraceExpanded = smooth(meanTraceExpanded, 5);
    
    t = sprintf('PeakDelta: %.0f pA, SteadyDelta: %0.f pA (from Baseline: %.0f pA)', ...
        out.peakFromBaseline, out.steadyStateFromBaseline, out.meanBaseline); 
    fprintf('%s\n', t);
    
    figure(1), clf;
    plot(relativeTimeExpanded*1000, meanTraceExpanded-meanBaseline,'k-')
    hold on
    plot(relativeTime(iPeak)*1000, currentPeak, 'r.');
    plot(windowMs, [currentSteadyState currentSteadyState], 'r-');
    plot(windowMs, [0 0], 'b-');
   
    box off;
    xlim(windowMs);
    xlabel('Time Relative to Step (ms)');
    ylabel('Current - Baseline (pA)');
    title(t);
    if(~isempty(yl))
        ylim(yl);
    end
end

end


