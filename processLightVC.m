function out = processLightVC( fnameOrData, varargin )

% throughout the code regions of time are associated with four variables:
%     somethingWindow is the lower and upper bounds on time (s) relative to some triggering event
%     somethingWindowInds are the indexing offsets relative to the index of the event
%     tvSomething is the time vector relative to trace start at 0 of this window
%     tvSomethingOffset is the time vector relative to the event at 0 of this window

%% Parameter initialization
par.showPlots = 1;
par.savePlots = 0;
par.chVmName = '';
par.chLaserName = '';
par.minLaserPulseInactivation = 0.02; % minimum laser pulse duration to compute inactivation tau
par.forcePeakDeactivation = 1;
par.shadeLaserRegion = 1;
par.maxPeakLatency = 0.2;
par.laserSignalType = 'toggle';
par.processVCTestOnly = 0;
par.currentPolarity = 'depolarizing'; % or hyperpolarizing

% some of these traces have a quick voltage delta at the beginning
% must feed in the timing of the peak and steady state regions
% since the secondary channel is typically laser, not command voltage
par.searchWindowVC = [];
par.searchWindowVCBaselinePrePeak = [.05 0.01];
par.searchWindowVCSteadyPostPeak = 0.01; % these are search windows that depend on the protocol being used
par.searchWindowVCSteadyPreMax = 0.01;
par.searchWindowVCMaxPostPeak = 0.15;

par.smoothPeakRegionTime = [];
par.smoothInactivationRegionTime = [];
par.smoothDeactivationRegionTime = [];

par.baselineWindowWidth = 0.5;
par.postPeakFitInactivationOffset = 0.002;
par.fitDeactivationOffset = 0.002;

par.baselinePostWindowWidth = 0.5;
par.fitDeactivationToPostBaseline = 0;
par.endPulseCurrentWindowWidth = 0.1;

par.computeLevelCrossingTimesOnset = [0.15 0.85];
par.computeLevelCrossingTimesOffset = [0.15 0.85];
par.imagesPath = '.';

assignargs(par, varargin);

%% Data load if necessary

if(~exist('fnameOrData', 'var'))
    fnameOrData = '';
end

if(isempty(chVmName) || isempty(chLaserName))
    chMatchRegexes = {}; % request them numerically
else
    chMatchRegexes = {chVmName, chLaserName};
end
    
abfData = abfDataLoad(fnameOrData, chMatchRegexes);
sis = abfData.sis;
time = abfData.time; 
nTraces = abfData.nTraces;
chVm = abfData.channels{1};
chLaser = abfData.channels{end};
fnamelast = abfData.fnamePart;
out = [];

% Find laser pulse times

if(strcmp(laserSignalType,'auto'))
    laserSignalType = 'toggle';
    [laserPulseInds laserPulseOffInds] = getLaserSignal(chLaser, 'toggle');
    if(isempty(laserPulseOffInds{1}))
        % if only one pulse, must be toggle
        [laserPulseInds laserPulseOffInds] = getLaserSignal(chLaser, 'state');
        laserSignalType = 'state';
    end
else
   [laserPulseInds laserPulseOffInds] = getLaserSignal(chLaser, laserSignalType); 
end

if(strcmp(currentPolarity, 'depolarizing'))
    currentPolaritySign = -1;
else
    currentPolaritySign = 1;
end

for iTrace = 1:nTraces

    if(isempty(fnamelast)) % don't know what filename to save to without prefix
        savePlots = 0;
    end
    
    % process initial VC probe if asked
    if(~isempty(searchWindowVC))
        % extract the region around the search window
        [regionVC indRegionVC tvRegionVC] = getPeriEventRegions(chVm(:,iTrace), 1, round(searchWindowVC/sis), 'sis', sis);
        
        % search for min and max values to figure out when the step happens
        [peakVC indPeak] = min(regionVC{1});
        tPeakVC = tvRegionVC(indPeak);
        
        regionSearchMax = regionVC{1};
        regionSearchMax(tvRegionVC <= tPeakVC) = -Inf;
        regionSearchMax(tvRegionVC >= tPeakVC+searchWindowVCMaxPostPeak) = -Inf;
        [peakMaxVC indPeakMax] = max(regionSearchMax);
        tPeakMaxVC = tvRegionVC(indPeakMax);
        
        searchWindowVCBaseline = [max([searchWindowVC(1), tPeakVC-searchWindowVCBaselinePrePeak(1)]) tPeakVC-searchWindowVCBaselinePrePeak(2)];
        regionBaseline = getPeriEventRegions(chVm(:,iTrace), 1, round(searchWindowVCBaseline/sis));
        
        searchWindowVCSteady = [tPeakVC+searchWindowVCSteadyPostPeak tPeakMaxVC-searchWindowVCSteadyPreMax];
        regionSteady = getPeriEventRegions(chVm(:,iTrace), 1, round(searchWindowVCSteady/sis));
    
        out(iTrace).vcBaseline = mean(regionBaseline{1});
        out(iTrace).vcPeakDelta = peakVC - mean(regionBaseline{1});
        out(iTrace).vcSteadyDelta = mean(regionSteady{1}) - mean(regionBaseline{1});
    end
    
    if(~processVCTestOnly)
        % retrieve laser pulse times
        laserOnInd = laserPulseInds{iTrace}(1);
        tLaserOn = time(laserOnInd);
        laserOffInd = laserPulseOffInds{iTrace}(1);
        tLaserOff = time(laserOffInd);
        laserPulseDuration = tLaserOff - tLaserOn;
        out(iTrace).laserPulseDuration = laserPulseDuration;

        % find baseline holding from pre-pulse region
        prePulseWindow = [max(0, tLaserOn-baselineWindowWidth) tLaserOn-sis];
        prePulseWindowInds = floor(prePulseWindow / sis);
        regions = getPeriEventRegions(chVm(:,iTrace), 1, prePulseWindowInds);
        baselineCurrent = mean(regions{1});
        out(iTrace).baselineCurrent = baselineCurrent;

        % extract region in which to search for maximum
        periPulseWindow = [0 maxPeakLatency];
        periPulseInds = floor(periPulseWindow / sis);
        tvPeriPulseOffset = (periPulseInds(1)*sis:sis:periPulseInds(2)*sis)';
        tvPeriPulseWindow = tvPeriPulseOffset + tLaserOn;
        regions = getPeriEventRegions( chVm(:, iTrace), laserOnInd, periPulseInds);
        periPeakRegion = regions{1};

        % smooth peri-peak region to obtain more robust estimates of the peak?
        if(~isempty(smoothPeakRegionTime))
            smoothLength = round(smoothPeakRegionTime / sis);
            periPeakRegion = smooth(periPeakRegion, smoothLength, 'rlowess');
        end

        % find peak current and time at which peak occurs in this region
        % currentPolaritySign = -1 for down currents, +1 for up currents, make them upwards
        [iMax indMax] = max(currentPolaritySign*periPeakRegion(1:round(maxPeakLatency/sis)));
        iMax = currentPolaritySign*iMax; % use its original sign
        tPeak = tvPeriPulseWindow(indMax);
        tPeakInd = indMax + laserOnInd;
        out(iTrace).tFlashToPeak = tPeak - tLaserOn;
        out(iTrace).peakCurrent = iMax - baselineCurrent;
        
        % compute percentage threshold level crossing times for rise to peak
        normalizedRise = (periPeakRegion-baselineCurrent)/(iMax-baselineCurrent);
        for iLevel = 1:length(computeLevelCrossingTimesOnset)
            threshLevel = computeLevelCrossingTimesOnset(iLevel);
            crossInd = find(normalizedRise > threshLevel, 1, 'first');
            tCross = tvPeriPulseWindow(crossInd);
            fldName = sprintf('tOnsetCrossLevel%02dp', round(threshLevel*100));
            out(iTrace).(fldName) = tCross - tLaserOn;
            levelCrossingTimesOnset(iLevel) = tCross;
            valueCrossingTimesOnset(iLevel) = periPeakRegion(crossInd);
        end

        computeInactivation = laserPulseDuration >= minLaserPulseInactivation;
        if(computeInactivation)
            % extract region during pulse on
            pulseWindow = [0 laserPulseDuration];
            pulseWindowInds = floor(pulseWindow / sis);
            tvPulseWindowOffset = (pulseWindowInds(1)*sis:sis:pulseWindowInds(2)*sis)';
            tvPulseWindow = tvPulseWindowOffset + tLaserOn;
            regions = getPeriEventRegions( chVm(:, iTrace), laserOnInd, pulseWindowInds);
            duringPulse = currentPolaritySign*regions{1};

            % smooth if requested
            if(~isempty(smoothInactivationRegionTime))
                duringPulse = smooth(duringPulse, smoothInactivationRegionTime, 'rlowess');
            end

            % find inactivation kinetics
            % fit a monoexponential with known peak value to the inactivation (peak to pulse off)

            % begin fitting with first value coming from smoothed periPulseRegion trace at appropriate offset
            postPeakFitOffsetInds = round(postPeakFitInactivationOffset/sis);
            indStartFit = indMax+postPeakFitOffsetInds;
            if(indStartFit < length(periPeakRegion))
                valStartFit = currentPolaritySign*periPeakRegion(indStartFit);
            else
                valStartFit = duringPulse(indStartFit);
            end
            tStartFit = tPeak + postPeakFitInactivationOffset;

            postPeak = duringPulse(indStartFit:end);
            tvPostPeak = tvPulseWindow(indStartFit:end);
            tvPostPeakOffset = tvPostPeak - tStartFit;

            % bound steady state to baseline and peak
            steadyStateBounds = [currentPolaritySign*baselineCurrent currentPolaritySign*iMax];
            if(steadyStateBounds(2) < steadyStateBounds(1)) % flip if nonmonotonic
                steadyStateBounds = fliplr(steadyStateBounds);
            end
            optFit = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 steadyStateBounds(1)], ...
                'Upper', [laserPulseDuration, steadyStateBounds(2)], ...
                'Startpoint', [laserPulseDuration/5, mean(postPeak(:))]);
            monoexpFun = fittype('(peak-steadyState)*exp(-t / tau) + steadyState', 'options', optFit, ...
                'independent', 't', 'coefficients', {'tau', 'steadyState'}, 'problem', 'peak');
            fitInactivation = fit(tvPostPeakOffset, postPeak, monoexpFun, 'problem', valStartFit);

            out(iTrace).steadyState = currentPolaritySign*fitInactivation.steadyState - baselineCurrent;
            out(iTrace).tauInactivation = fitInactivation.tau;
            
            % find ending value (could differ from steady state for really long tau)
            % extract region during pulse on
            endPulseWindowInds = length(duringPulse)-floor(endPulseCurrentWindowWidth/sis):length(duringPulse);
            endPulse = duringPulse(endPulseWindowInds);
            tvEndPulse = tvPulseWindow(endPulseWindowInds);
            endPulseCurrent = mean(endPulse);
            out(iTrace).endPulseCurrent = endPulseCurrent;
            
        else
            out(iTrace).steadyState = NaN;
            out(iTrace).tauInactivation = NaN;
            out(iTrace).endPulseCurrent = NaN;
        end

        % find deactivation kinetics: we take the region following the light pulse off
        % or the peak time, whichever is later (short laser pulses have a alpha fn shape)
        preOffset = 0.005; % grab a bit before offset for smoothing purposes
        deactivationWindow = [-preOffset fitDeactivationOffset+baselinePostWindowWidth]; % determines how far beyond the pulse off/peak to consider
        deactivationWindowInds = round(deactivationWindow / sis);
        tvDeactivationOffset = (0:sis:(deactivationWindowInds(2)-deactivationWindowInds(1))*sis)' - preOffset;
        
        if(tLaserOff > tPeak)
            % use the period following laser pulse off    
            tvDeactivation = tvDeactivationOffset + tLaserOff + fitDeactivationOffset + preOffset;
            regions = getPeriEventRegions( chVm(:,iTrace), laserOffInd, deactivationWindowInds);
            deactivation = regions{1};
        else
            % use a period following the peak value
            tvDeactivation = tvDeactivationOffset + tPeak + fitDeactivationOffset + preOffset;
            regions = getPeriEventRegions( chVm(:,iTrace), tPeakInd, deactivationWindowInds);
            deactivation = regions{1};
        end
        
        % make upwards
        deactivation = currentPolaritySign*deactivation;

        % smooth if requested
        if(~isempty(smoothDeactivationRegionTime))
            deactivation = smooth(deactivation, smoothDeactivationRegionTime, 'rlowess');
        end

        % isolate region fitDeactivationOffset after the laser offset for fitting purposes (leads to better fit)
        deactivationForFitting = deactivation(tvDeactivationOffset > fitDeactivationOffset);
        tvDeactivationForFitting = tvDeactivationOffset(tvDeactivationOffset > fitDeactivationOffset) - fitDeactivationOffset;
        
        % should we fit the peak parameter in the exponential or just the tau?
        if(forcePeakDeactivation)
            optFit = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Lower', [0 -Inf], ...
                'Upper', [deactivationWindow(2), currentPolaritySign*iMax], ...
                'Startpoint', [laserPulseDuration/5, baselineCurrent]);
            monoexpFun = fittype('(peak-steadyState)*exp(-t / tau) + steadyState', 'options', optFit, ...
                'independent', 't', 'coefficients', {'tau', 'steadyState'}, 'problem', 'peak');
            fitDeactivation = fit(tvDeactivationForFitting, deactivationForFitting, monoexpFun, 'problem', deactivationForFitting(1));
        else
            optFit = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Lower', [-Inf 0 -Inf], ...
                'Upper', [Inf deactivationWindow(2), iMax], ...
                'Startpoint', [baselineCurrent laserPulseDuration/5, baselineCurrent]);
            monoexpFun = fittype('(peak-steadyState)*exp(-t / tau) + steadyState', 'options', optFit, ...
                'independent', 't', 'coefficients', {'peak', 'tau', 'steadyState'});
            fitDeactivation = fit(tvDeactivationForFitting, deactivationForFitting, monoexpFun);
        end
        postBaselineCurrent = currentPolaritySign*fitDeactivation.steadyState;
        
        out(iTrace).tauDeactivation = fitDeactivation.tau;
        out(iTrace).postBaselineCurrent = postBaselineCurrent;
        
        % compute level crossing times for offset
        normalizedRise = 1-(deactivation-postBaselineCurrent)/(endPulseCurrent-postBaselineCurrent);
        for iLevel = 1:length(computeLevelCrossingTimesOffset)
            threshLevel = computeLevelCrossingTimesOffset(iLevel);
            crossInd = find(normalizedRise > threshLevel, 1, 'first');
            tCross = tvDeactivation(crossInd);
            fldName = sprintf('tOffsetCrossLevel%02dp', round(threshLevel*100));
            out(iTrace).(fldName) = tCross - tLaserOff;
            levelCrossingTimesOffset(iLevel) = tCross;
            valueCrossingTimesOffset(iLevel) = deactivation(crossInd);
        end

        
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
            plot(tvDisplayWindowOffset*1000, response, '-', 'Color', [0.5 0.5 0.5]);

            % plot smoothed peak region
            if(~isempty(smoothPeakRegionTime))
                plot(tvPeriPulseWindow*1000-tDisplayOffsetMS, periPeakRegion, 'k-', 'LineWidth', 2);
            end
            
            % plot the smoothed deactivation region
            if(~isempty(smoothDeactivationRegionTime))
                plot(tvDeactivation*1000-tDisplayOffsetMS, deactivation, 'k-', 'LineWidth', 2);
            end
            
            % plot the laser signal
            blueLightLevel = yl(2)-0.08*(yl(2)-yl(1));
            plot([0 1000*laserPulseDuration], [blueLightLevel blueLightLevel], '-', 'Color', [0.3 0.5 1], 'LineWidth', 3);

            % plot the pre baseline level
            plot(prePulseWindow*1000-tDisplayOffsetMS, currentPolaritySign*[baselineCurrent baselineCurrent], 'b-', 'LineWidth', 2);
            % and post baseline level
            plot([tvDeactivation(1) tvDeactivation(end)]*1000-tDisplayOffsetMS, currentPolaritySign*[postBaselineCurrent postBaselineCurrent], 'b-', 'LineWidth', 2);
            
            % plot the inactivation fit
            if(computeInactivation)
                plot(tvPostPeak*1000 - tDisplayOffsetMS, fitInactivation(tvPostPeakOffset), 'r-','LineWidth',2);
            end

            % plot the deactivation fit
            plot(tvDeactivation*1000 - tDisplayOffsetMS, fitDeactivation(tvDeactivationOffset), 'g-', 'LineWidth', 2);
    %         plot(tvDeactivation*1000 - tDisplayOffsetMS, deactivation, 'c-');

            % plot the onset level crossing times
            for iLevel = 1:length(computeLevelCrossingTimesOnset)
               plot(levelCrossingTimesOnset(iLevel)*1000-tDisplayOffsetMS, valueCrossingTimesOnset(iLevel), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
            end
            % plot the offset level crossing times
            for iLevel = 1:length(computeLevelCrossingTimesOffset)
               plot(levelCrossingTimesOffset(iLevel)*1000-tDisplayOffsetMS, valueCrossingTimesOffset(iLevel), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
            end
    
            % plot the peak location
            plot(tPeak*1000 - tDisplayOffsetMS, iMax, 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

            % plot the endPulseCurrent level
            if(computeInactivation)
                plot([tvEndPulse(1) tvEndPulse(end)]*1000 - tDisplayOffsetMS, [endPulseCurrent endPulseCurrent], 'b-', 'LineWidth', 2);
            end
            
            box off
            if(yl(2)>yl(1))
                ylim(yl);
            end
            xlim(1000*displayWindow);
            xlabel('Time Relative to Light Onset (ms)');
            ylabel('Current (pA)');
            title(sprintf('Light Pulse duration: %.0f ms', laserPulseDuration*1000));

            if(savePlots)
                [path name] = fileparts(fnameOrData);
                figName = sprintf('%s/%s.trace%d.png', imagesPath, name, iTrace);
                saveas(gcf, figName, 'png');
            end
        end
    end
    
    if(showPlots && ~isempty(searchWindowVC))
        figTrace = figure(34+iTrace); clf; set(figTrace, 'Color', [1 1 1], 'Name', 'VC Test Probe', 'NumberTitle', 'off');
        plot(tvRegionVC, regionVC{1}, '-', 'Color', [0.4 0.4 0.4]);
        hold on
        plot(searchWindowVCSteady, [1 1]*mean(regionSteady{1}), 'g-', 'LineWidth', 2);
        plot(searchWindowVCBaseline, [1 1]*mean(regionBaseline{1}), 'b-', 'LineWidth', 2);
        plot(searchWindowVCBaseline, [1 1]*mean(regionSteady{1}), 'g-', 'LineWidth', 2);
        plot(tPeakVC, peakVC, 'r+', 'MarkerSize', 8);
        box off
        title(sprintf('VC Probe: Baseline = %.0f, dPeak = %.0f, dSteady = %.0f', ...
            out(iTrace).vcBaseline, out(iTrace).vcPeakDelta, out(iTrace).vcSteadyDelta));

        xl = [searchWindowVCBaseline(1) searchWindowVCSteady(2)+2*searchWindowVCSteadyPreMax];
        if(xl(2) > xl(1))
            xlim(xl);
        end
        xlabel('Time (s)');
        ylabel('Current (pA)');

        if(savePlots)
            [path name] = fileparts(fnameOrData);
            figName = sprintf('%s/%s.trace%d.vcTest.png', imagesPath, name, iTrace);
            saveas(figTrace, figName, 'png');
        end
    end

end

