function out = processLightTrains(fnameOrData, varargin)
% out = processLightTrains(data,samplingInMicrosec,varargin) or pLT(fname, varargin)

%% Parameter initialization
par.showPlots = 1;
par.savePlots = 0;
par.chVmName = 'IN [01]';
par.chLaserName = 'IN [23]';
par.elicitedSpikeWindow = [1e-3 8e-3]; % where valid spikes occur
par.snippetWindow = [-1e-3 13e-3]; % to display on the snippets plot
par.laserSignalType = 'toggle'; % 'toggle' or 'state', toggle means each rising edge toggle on/off
assignargs(par, varargin);

% spike detection params
spikeDetectParams.threshHigh = -10;
spikeDetectParams.threshLow = -30;
spikeDetectParams.refractoryMs = 2;
spikeDetectParams.maxSpikeWidthMs = 6;
spikeDetectParams = structargs(spikeDetectParams, varargin);

%% Data load if necessary

if(~exist('fnameOrData', 'var'))
    fnameOrData = '';
elseif(isempty(fnameOrData))
    fnameOrData = '';
end

chMatchRegexes = {chVmName, chLaserName};
[d sis h chMatches time fname fnamelast] = abfDataLoad(fnameOrData, chMatchRegexes);
nTraces = size(d,3);
chVm = chMatches{1};
chLaser = chMatches{2};

if(isempty(fnamelast)) % don't know what filename to save to without prefix
    savePlots = 0;
end

% convert times into index units
spikeDetectParams.refractory = round(spikeDetectParams.refractoryMs / 1000 / sis);
spikeDetectParams.maxSpikeWidthMS = round(spikeDetectParams.maxSpikeWidthMs / 1000 / sis);
snippetWindowInds = round(snippetWindow / sis);
elicitedSpikeWindowInds = elicitedSpikeWindow / sis;

% Find laser pulse times as cell{iTrace}
[laserPulseInds laserPulseOffInds] = getLaserSignal(chLaser, laserSignalType);
    
% Calculate laser pulse duration and frequency for each trace
laserPulseDurations = arrayfun(@(i) mean(laserPulseOffInds{i}(1) - laserPulseInds{i}(1)), 1:nTraces) * sis;
meanLaserPulseDuration = mean(laserPulseDurations);
laserPulseIntervals = arrayfun(@(i) mean(laserPulseInds{i}(2:end) - laserPulseInds{i}(1:end-1)), 1:nTraces) * sis;
laserPulseFreqHz = 1./laserPulseIntervals;

% Find all spike times
spikeInds = spikeDetect(chVm, spikeDetectParams);

% Grab acceptable windows (indices) for spiking around each pulse time
% cell{iTrace} = nspikes x 2  == [min index, max index];
[~, pulseElicitedSpikeIndWindows] = getPeriEventRegions(chVm, laserPulseInds, elicitedSpikeWindowInds, 'onlyIndexWindow', 1);

% Determine which pulseElicitedWindows contain detected spikes
% Calculate spike-success rate for each trace separately
spikesFoundInWindow = cell(nTraces, 1);
pulsesWithSpike = zeros(nTraces,1);
for iTrace = 1:nTraces
    pulseWindows = pulseElicitedSpikeIndWindows{iTrace};
    % number of spikes found in each window {iTrace}(pulseNum)
    spikesFoundInWindow{iTrace} = arrayfun(@(iMin, iMax) nnz(spikeInds{iTrace} >= iMin & spikeInds{iTrace} <= iMax), ...
        pulseWindows(:,1), pulseWindows(:,2));
    
    proportionPulsesWithSpike(iTrace) = mean(spikesFoundInWindow{iTrace} > 0);
end

% Grab larger regions near laser pulse times for display only
[pulseSnippets ~] = getPeriEventRegions(chVm, laserPulseInds, snippetWindowInds);

% Collect successful and failure snippets for visual display
successSnippets = []; failureSnippets = [];
for iTrace = 1:nTraces
    successSnippets = [successSnippets, pulseSnippets{iTrace}(:, spikesFoundInWindow{iTrace} >0)];
    failureSnippets = [failureSnippets, pulseSnippets{iTrace}(:, spikesFoundInWindow{iTrace}==0)];
end

%% Print calculated results
fprintf('Mean Laser Pulse Duration(ms) \t%.1f\n\n', meanLaserPulseDuration*1000);
fprintf('Pulse Frequency (Hz)\tProportion of Successful Pulses\n');
for iTrace = 1:nTraces
    fprintf('%.0f\t%.2f\n', laserPulseFreqHz(iTrace), proportionPulsesWithSpike(iTrace));
end

out.params = par;
out.fname = fname;
out.meanLaserPulseDuration = meanLaserPulseDuration;
out.laserPulseDurations = laserPulseDurations;
out.laserPulseFreqHz = laserPulseFreqHz;
out.pulseSnippets = pulseSnippets;
out.spikesFoundInWindow = spikesFoundInWindow;
out.proportionPulsesWithSpike = proportionPulsesWithSpike;

%% Show extracted spike snippets
if(showPlots)

    tSnippetWindowMs = 1000*sis*(snippetWindowInds(1):snippetWindowInds(2));
    hfigSnippets = figure(5); clf;
    set(hfigSnippets, 'Name', 'Pulse-Triggered Waveforms', 'NumberTitle', 'off', 'Color', [1 1 1]);
    axh = gca;
    ymin = -80;
    ymax = 30;
    hSpikeRegionOverlay = fill(1000*[elicitedSpikeWindow(1) elicitedSpikeWindow(2) elicitedSpikeWindow(2) elicitedSpikeWindow(1)], ...
     [ymin+0.01 ymin+0.01 ymax ymax], [0.8 1 1], 'EdgeColor', 'none');
    hold on
    if(~isempty(successSnippets))
        plot(axh, tSnippetWindowMs, successSnippets, 'Color',[0.39 0.67 0.82], 'LineWidth', 1);
    end
    if(~isempty(failureSnippets))
        plot(axh, tSnippetWindowMs, failureSnippets, 'Color',[0.5 0.5 0.5], 'LineWidth', 1);
    end
    plot([tSnippetWindowMs(1) tSnippetWindowMs(end)], [spikeDetectParams.threshHigh spikeDetectParams.threshHigh], '-', 'Color', [0.7 0.7 0.7]);
    plot([tSnippetWindowMs(1) tSnippetWindowMs(end)], [spikeDetectParams.threshLow spikeDetectParams.threshLow], '--', 'Color', [0.7 0.7 0.7]);
    
    xlabel('Time After Pulse Onset (ms)');
    ylabel('Voltage (mV)');
    xlim([tSnippetWindowMs(1) tSnippetWindowMs(end)]);
    ylim([ymin ymax]);
    box off
    title('Pulse-triggered Voltage Waveforms');
    
    if(savePlots)
        saveas(hfigSnippets, [fnamelast '.snippets.png'], 'png');
        saveas(hfigSnippets, [fnamelast '.snippets.fig'], 'fig');
    end
    
    %% Show spike efficacy
    hfigEfficacy = figure(6); clf;
    set(hfigEfficacy, 'Name', 'Pulse Efficacy', 'NumberTitle', 'off', 'Color', [1 1 1]);
    plot(laserPulseFreqHz, proportionPulsesWithSpike, 'ks', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    xlabel('Pulse Frequency (Hz)');
    ylabel('Probability of Elicing Spike');
    title(sprintf('Pulse Efficacy, Duration = %d ms', meanLaserPulseDuration * 1000));
    ylim([0.9*min(proportionPulsesWithSpike) 1.05]);
    box off
    
    if(savePlots)
        saveas(hfigEfficacy, [fnamelast '.snippets.png'], 'png');
        saveas(hfigEfficacy, [fnamelast '.snippets.png'], 'fig');
    end

end
