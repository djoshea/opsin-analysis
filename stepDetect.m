function [ stepStartInds stepLevels preStepLevel] = stepDetect(traces, varargin)
% [spikeInds spikeFallingInds] = spikeDetect(traces, varagin)
% params: stdThresh
%
% Criteria: finds points where the time derivative exceeds some multiple of
%     the std of the derivative signal. Measures mean value following each step change.

par.stdThresh = 10;
par.asCell = true;
par.plotTrace = 0;
assignargs(par, varargin);

nsweeps = size(traces,2);
stepStartInds = cell(nsweeps,1);
stepLevels = cell(nsweeps, 1);
preStepLevel = cell(nsweeps, 1);

for iSweep = 1:nsweeps
  data = traces(:,iSweep);
  df = diff(data);
 
  % find inds where derivative > stdThresh * std of the derivative
  evStartInds = find(abs(df) > stdThresh*std(df));
  % find the end of each step interval for mean calculation
  evEndInds = [evStartInds(2:end)-1; length(data)];
  % find the mean signal value within each step
  meanLevels = arrayfun(@(istart, iend) mean(data(istart:iend)), evStartInds, evEndInds);
 
  stepStartInds{iSweep} = evStartInds;
  stepLevels{iSweep} = meanLevels;
  preStepLevel{iSweep} = mean(data(1:evStartInds-1));
end

if(plotTrace)
    % plot one trace (plotTrace) along with markers indicating inds and levels at each step
    figure(), clf;
    plot(traces(:,plotTrace),'k-'), hold on;
    plot([1; stepStartInds{plotTrace}], [preStepLevel{plotTrace}; stepLevels{plotTrace}], ...
        'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    box off;
end

if(nsweeps == 1 && ~asCell)
    stepStartInds = stepStartInds{1};
    stepLevels = stepLevels{1};
    preStepLevel = preStepLevel{1};
end


   