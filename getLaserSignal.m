function [ laserOnInds laserOffInds ] = getLaserSignal( chLaser, laserSignalType )
% [ laserOnInds laserOffInds ] = getLaserSignal( chLaser, laserSignalType )
% laserSignalType is either 'toggle' or 'state'

% Find laser pulse times
laserThreshs = linterp([0 1], [min(chLaser(:)) max(chLaser(:))], [0.3 0.7]);
[laserSignalRisingInds laserSignalFallingInds] = spikeDetect(chLaser, 'threshHigh', laserThreshs(2), 'threshLow', laserThreshs(1));

if(strcmp(laserSignalType, 'toggle'))
    % Using rising edges as toggle signal
    laserOnInds = cellfun(@(inds) inds(1:2:end), laserSignalRisingInds, 'UniformOutput', false);
    laserOffInds = cellfun(@(inds) inds(2:2:end), laserSignalRisingInds, 'UniformOutput', false);
else
    % Use rising edge as turning on and falling edge as turning off
    laserOnInds = laserSignalRisingInds;
    laserOffInds = laserSignalFallingInds;
end

end

