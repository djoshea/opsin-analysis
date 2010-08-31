function [ spikeInds spikeLowInds ] = spikeDetect(traces, varargin)
% [spikeInds spikeFallingInds] = spikeDetect(traces, varagin)
% params: threshHigh, threshLow, refractory, maxWidth
%
% Criteria: a spike must rise above threshHigh (at index spikeInds) and then fall 
%   below threshLow (at index spikeLowInds) within maxWidth inds, with no 
%   subsequent spikes allowed before refractory inds have elapsed since the
%   threshHigh crossing of the previous spike

par.threshHigh = 0;
par.refractory = []; % in units of sampling
par.threshLow = -30;
par.maxSpikeWidth = Inf; % in units of sampling
assignargs(par, varargin);

nsweeps = size(traces,2);
spikeInds = {};
spikeLowInds = {};

for sweepno = 1:nsweeps
  data = traces(:,sweepno);
  
  % times when trace crosses the  or repol thresholds in the correct direction
  threshUpCross =   [0; diff(data >= threshHigh)] == 1;
  threshDownCross = [0; diff(data > threshLow)] == -1;
  
  % find times of upward threshold crossings separated by at least a refractory period
  upTimes = find(threshUpCross);
  nextDownTimes = [];
  if(~isempty(upTimes))
      if(refractory)
          upTimes = upTimes([true; diff(upTimes) > refractory]);
      end
      
      % reject spikes that don't cross repolarization threshold within maxspikewidth
      downTimes = find(threshDownCross);
      nextDownTime = @(ut) min(downTimes(downTimes > ut) - ut) + ut;
      nextDownTimes = arrayfun(nextDownTime, upTimes);
      nextUpTimes = [upTimes(2:end); Inf];

      % select spikes that cross lowThresh within maxSpikeWidth and before the threshHighCrossing
      selectedUpTimes = nextDownTimes < (upTimes + maxSpikeWidth) & nextDownTimes < nextUpTimes;
      upTimes = upTimes(selectedUpTimes);
      nextDownTimes = nextDownTimes(selectedUpTimes);
  end
  
  spikeInds{sweepno} = upTimes;
  spikeLowInds{sweepno} = nextDownTimes;
end

% if(nsweeps == 1)
%     spikeInds = spikeInds{1};
% end