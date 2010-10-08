includeList = zeros(size(cellList));
includeList([3 17]) = 1;
def.includeList = includeList;

def.holdings = [-60];
def.ndFilters = 0;
def.pulseWidths = [1 2 5];
def.printIncludedOnly = 1;
def.showPlots = 1;
def.metrics = {'proportionPulsesWithSpike', 'proportionPulsesWithMultipleSpikes', 'spikesElicited', 'spontaneousSpikes'};
def.metricNames = {'Pulse Efficacy', 'Multiple Spike Likelihood', 'Total Spikes Elicited', 'Spontaneous Spikes'};
def.metricMultiplier = [40, 40 1 1];
def.metricLimits = [0 40; 0 40; 0 200; 0 10];
def.metricPlot = [];
def.metricPrintDigits = [];
def.printCellDataFields = 0;
def.printData = 1;
def.savePlots = 0;
def.showPlots = 0;

[usableCellList results] = analysisLightTrains(cellList, cellData, def);
