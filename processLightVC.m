function out = processLightVC( fnameOrData, varargin )

% throughout the code regions of time are associated with four variables:
%     somethingWindow is the lower and upper bounds on time (s) relative to some triggering event
%     somethingWindowInds are the indexing offsets relative to the index of the event
%     tvSomething is the time vector relative to trace start at 0 of this window
%     tvSomethingOffset is the time vector relative to the event at 0 of this window

%% Parameter initialization
par.showPlots = 1;
par.savePlots = 0;
par.chVmName = 'IN [01]';
par.chLaserName = 'IN [23]';
par.minLaserPulseInactivation = 0.02; % minimum laser pulse duration to compute inactivation tau
par.forcePeakDeactivation = 1;
par.shadeLaserRegion = 1;
assignargs(par, {});

%% Data load if necessary

if(~exist('fnameOrData', 'var'))
    fnameOrData = '';
end

chMatchRegexes = {chVmName, chLaserName};
[d sis h chMatches time fname fnamelast] = abfDataLoad(fnameOrData, chMatchRegexes);
nTraces = size(d,3);
chVm = chMatches{1};
chLaser = chMatches{2};

out = [];

% Find laser pulse times
laserThreshs = linterp([0 1], [min(chLaser(:)) max(chLaser(:))], [0.3 0.7]);
laserTriggerPulses = spikeDetect(chLaser, 'threshLow', laserThreshs(1), 'threshHigh', laserThreshs(2));
  
for iTrace = 1:nTraces

    if(isempty(fnamelast)) % don't know what filename to save to without prefix
        savePlots = 0;
    end
    
    % retrieve laser pulse times
    laserOnInd = laserTriggerPulses{iTrace}(1);
    tLaserOn = time(laserOnInd);
    laserOffInd = laserTriggerPulses{iTrace}(2);
    tLaserOff = time(laserOffInd);
    laserPulseDuration = tLaserOff - tLaserOn;
    out(iTrace).laserPulseDuration = laserPulseDuration;
    
    % find baseline holding from pre-pulse region
    prePulseWindowInds = [0 laserOnInd-1];
    regions = getPeriEventRegions(chVm(:,iTrace), 1, prePulseWindowInds);
    baselineCurrent = -mean(regions{1});
    out(iTrace).baselineCurrent = baselineCurrent;
    
    % extract region in which to search for maximum
    pulseWindow = [0 laserPulseDuration + 0.5];
    pulseWindowInds = floor(pulseWindow / sis);
    tvPulseWindowOffset = (pulseWindowInds(1)*sis:sis:pulseWindowInds(2)*sis)';
    tvPulseWindow = tvPulseWindowOffset + tLaserOn;
    regions = getPeriEventRegions( chVm(:, iTrace), laserOnInd, pulseWindowInds);
    duringPulse = regions{1};
    
    % find peak current and time at which it occurs
    [iMax indMax] = min(duringPulse);
    iMax = -iMax;
    tPeak = tvPulseWindow(indMax);
    tPeakInd = indMax + laserOnInd;
    out(iTrace).tFlashToPeak = tPeak - tLaserOn;
    out(iTrace).peakCurrent = iMax - baselineCurrent;
    
    computeInactivation = laserPulseDuration >= minLaserPulseInactivation;
    if(computeInactivation)
        % extract region during pulse on
        pulseWindow = [0 laserPulseDuration];
        pulseWindowInds = floor(pulseWindow / sis);
        tvPulseWindowOffset = (pulseWindowInds(1)*sis:sis:pulseWindowInds(2)*sis)';
        tvPulseWindow = tvPulseWindowOffset + tLaserOn;
        regions = getPeriEventRegions( chVm(:, iTrace), laserOnInd, pulseWindowInds);
        duringPulse = regions{1};
        
        % find inactivation kinetics
        % fit a monoexponential with known peak value to the inactivation (peak to pulse off)
        postPeak = duringPulse(indMax:end);
        tvPostPeak = tvPulseWindow(indMax:end);
        tvPostPeakOffset = tvPostPeak - tPeak;
        optFit = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0], ...
            'Upper', [laserPulseDuration, iMax], ...
            'Startpoint', [laserPulseDuration/5, mean(postPeak(:))]);
        monoexpFun = fittype('-(peak-steadyState)*exp(-t / tau) - steadyState', 'options', optFit, ...
            'independent', 't', 'coefficients', {'tau', 'steadyState'}, 'problem', 'peak');
        fitInactivation = fit(tvPostPeakOffset, postPeak, monoexpFun, 'problem', iMax);
     
        out(iTrace).steadyState = fitInactivation.steadyState - baselineCurrent;
        out(iTrace).tauInactivation = fitInactivation.tau;
    end
    
    % find deactivation kinetics: we take the region following the light pulse off
    % or the peak time, whichever is later (short laser pulses have a alpha fn shape)
    deactivationWindow = [0 0.5]; % determines how far beyond the pulse off/peak to consider
    deactivationWindowInds = round(deactivationWindow / sis);
    tvDeactivationOffset = (deactivationWindowInds(1)*sis:sis:deactivationWindowInds(2)*sis)';
    
    if(tLaserOff > tPeak)
        % use the period following laser pulse off    
        tvDeactivation = tvDeactivationOffset + tLaserOff;
        regions = getPeriEventRegions( chVm(:,iTrace), laserOffInd, deactivationWindowInds);
        deactivation = regions{1};
    else
        % use a period following the peak value
        tvDeactivation = tvDeactivationOffset + tPeak;
        regions = getPeriEventRegions( chVm(:,iTrace), tPeakInd, deactivationWindowInds);
        deactivation = regions{1};
    end

    % should we fit the peak parameter in the exponential or just the tau?
    if(forcePeakDeactivation)
        optFit = fitoptions('Method', 'NonlinearLeastSquares', ...
            'Lower', [0 0], ...
            'Upper', [deactivationWindow(2), iMax], ...
            'Startpoint', [laserPulseDuration/5, baselineCurrent]);
        monoexpFun = fittype('-(peak-steadyState)*exp(-t / tau) - steadyState', 'options', optFit, ...
            'independent', 't', 'coefficients', {'tau', 'steadyState'}, 'problem', 'peak');
        fitDeactivation = fit(tvDeactivationOffset, deactivation, monoexpFun, 'problem', -deactivation(1));
    else
        optFit = fitoptions('Method', 'NonlinearLeastSquares', ...
            'Lower', [0 0 0], ...
            'Upper', [iMax deactivationWindow(2), iMax], ...
            'Startpoint', [baselineCurrent laserPulseDuration/5, baselineCurrent]);
        monoexpFun = fittype('-(peak-steadyState)*exp(-t / tau) - steadyState', 'options', optFit, ...
            'independent', 't', 'coefficients', {'peak', 'tau', 'steadyState'});
        fitDeactivation = fit(tvDeactivationOffset, deactivation, monoexpFun);
    end
    out(iTrace).tauDeactivation = fitDeactivation.tau;
    
    % extract surrounding region
    if(showPlots)
        
        % grab the region of the trace to display
        showBeforePulseOn = max(min(0.1,2*fitDeactivation.tau), 0.2*laserPulseDuration);
        showAfterPulseOff = max(min(0.4,10*fitDeactivation.tau), 0.2*laserPulseDuration);
        displayWindow = [-showBeforePulseOn laserPulseDuration+showAfterPulseOff];
        tDisplayOffsetMS = 1000*tLaserOn; % subtract this from real-time to plot on this figure
        
        displayWindowInds = round(displayWindow / sis);
        [regions] = getPeriEventRegions( chVm(:,iTrace), laserOnInd, displayWindowInds);
        response = regions{1};
        tvDisplayWindowOffset = displayWindowInds(1)*sis:sis:displayWindowInds(2)*sis;
       
        figTrace = figure(14+iTrace); clf; set(figTrace, 'Color', [1 1 1], 'Name', 'Photocurrent Trace', 'NumberTitle', 'off');
        trLim = [min(response) max(response)];
        range = trLim(2) - trLim(1);
        yl = [trLim(1)-0.1*range trLim(2)+0.2*range];
        
        % shade the laser on portion
        if(shadeLaserRegion)
            fill([0 0 1000*laserPulseDuration 1000*laserPulseDuration], ...
                [yl(1)+1 yl(2)-1 yl(2)-1 yl(1)+1], ...
                [0.9 0.9 1], 'EdgeColor', 'none');
        end
 
        hold on
        
        % plot region of trace
        plot(tvDisplayWindowOffset*1000, response, 'k-');
        
        % plot the laser signal
        blueLightLevel = -baselineCurrent + 0.1*range;
        plot([0 1000*laserPulseDuration], [blueLightLevel blueLightLevel], '-', 'Color', [0.3 0.5 1], 'LineWidth', 3);
        
        % plot the inactivation fit
        if(laserPulseDuration >= minLaserPulseInactivation)
            plot(tvPostPeak*1000 - tDisplayOffsetMS, fitInactivation(tvPostPeakOffset), 'r-','LineWidth',2);
        end
        
        % plot the deactivation fit
        plot(tvDeactivation*1000 - tDisplayOffsetMS, fitDeactivation(tvDeactivationOffset), 'g-', 'LineWidth', 2);
        
        % plot the peak location
        plot(tPeak*1000 - tDisplayOffsetMS, -iMax, 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        
        box off
        ylim(yl);
        xlim(1000*displayWindow);
        xlabel('Time Relative to Light Onset (ms)');
        ylabel('Current (pA)');
        title(sprintf('Light Pulse duration: %.0f ms', laserPulseDuration*1000));
    end

end

