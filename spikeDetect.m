function spikeInds = spikeDetect(traces, varargin)
% spiketimes = spikeDetect(traces, varagin)
% params: threshHigh, threshLow, refractory, maxWidth

par.threshHigh = 0;
par.refractory = []; % in units of sampling
par.threshLow = -30;
par.maxSpikeWidth = Inf; % in units of sampling
par.keepEveryOther = 0; % e.g. for on/off toggling signal
par.keepEveryOtherEven = 0; % e.g. for on/off toggling signal, overrides above
assignargs(par, varargin);

nsweeps = size(traces,2);
spiketimes = {};

for sweepno = 1:nsweeps
  data = traces(:,sweepno);
  
  % times when trace crosses the  or repol thresholds in the correct direction
  threshUpCross =   [0; diff(data >= threshHigh)] == 1;
  threshDownCross = [0; diff(data > threshLow)] == -1;
  
  % find times of upward threshold crossings separated by at least a refractory period
  upTimes = find(threshUpCross);
  if(~isempty(upTimes))
      if(refractory)
          upTimes = upTimes([true; diff(upTimes) > refractory]);
      end
      
      % reject spikes that don't cross repolarization threshold within maxspikewidth
      downTimes = find(threshDownCross);
      nextDownTime = @(ut) min(downTimes(downTimes > ut) - ut) + ut;
      nextDownTimes = arrayfun(nextDownTime, upTimes);
      nextUpTimes = [upTimes(2:end); Inf];

      upTimes = upTimes(nextDownTimes < upTimes + maxSpikeWidth & ...
          nextDownTimes < nextUpTimes);
  end
  
  if(keepEveryOtherEven)
      upTimes = upTimes(2:2:end);
  elseif(keepEveryOther)
      upTimes = upTimes(1:2:end);
  end
  
  spikeInds{sweepno} = upTimes;
end

% if(nsweeps == 1)
%     spikeInds = spikeInds{1};
% end