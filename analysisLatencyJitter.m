%% Gather the data by pulse width
pulseWidths = [5 2 1];
holdings = [-60 -70];
def.ndFilters = 0;

def.pulseWidths = pulseWidths;
def.holdings = holdings;
def.includeAll = 1;
def.printData = 1;
def.savePlots = 0;
def.printCellDataFields = 0;

def.metrics = {'meanSpikeLatency', 'stdSpikeLatency'};
def.metricNames = {'Spike Latency', 'Spike Jitter'};
def.metricMultiplier = [1000, 1000];
def.metricLimits = [0 8; 0 3];
def.metricPrintDigits = [3 3];

[usableCellList, results] = analysisLightTrains(cellList, cellData, def);


% %% Print the results: one line for each cell
% 
% % first header row, print pulse widths and holding
% fprintf('\t'); % skip column for name
% for iPulseWidth = 1:length(pulseWidths)
%     pulseWidth = pulseWidths(iPulseWidth);
%     freqHz = results{iPulseWidth}(1).laserPulseFreqHz;
%     for iHolding = 1:length(holdings)
%         holding = holdings(iHolding);
%         fprintf('%dms, %dmV', pulseWidth, holding);
%         
%         for f = 1:length(freqHz)
%             for s= 1:length(subPeriods)
%                 fprintf('\t');
%             end
%         end
%     end
% end
% fprintf('\n');
% 
% % second header row
% fprintf('Name\t');
% for iHolding = 1:length(holdings)
%     holding = holdings(iHolding);
%      for iPulseWidth = 1:length(pulseWidths)
%         pulseWidth = pulseWidths(iPulseWidth);
%         freqHz = results{iPulseWidth}(1).laserPulseFreqHz;
%         for f = 1:length(freqHz)
%             for s= 1:length(subPeriods)
%                 fprintf('%.0fHz P%d\t', freqHz(f), s);
%             end
%         end
%     end
% end
% fprintf('\n');
%     
% % data rows, one for each cell, protocols concatenated horizontally
% for iCell = 1:length(cellData)
%     fprintf('%s\t',cellData(iCell).Cell_Name);
%    
%     for iHolding = 1:length(holdings)
%         holding = holdings(iHolding);
%          for iPulseWidth = 1:length(pulseWidths)
%             pulseWidth = pulseWidths(iPulseWidth);
%             freqHz = results{iPulseWidth}(1).laserPulseFreqHz;
%             inList = find(strcmp(cellData(iCell).Cell_Name, {usableCellList{iPulseWidth,1,iHolding}.name}));
%             if(~isempty(inList))
%                 result = results{iPulseWidth,1,iHolding}(inList(1));
%                 for f = 1:length(freqHz)
%                     for subperiod = 1:length(subPeriods)
%                         fprintf('%2d\t', result.proportionPulsesWithSpikeSubperiod{f}(subperiod)*10);
%                     end
%                 end
%             end
%         end
%     end
%     fprintf('\n');
% end
% 
% %% Plot the results: one row of plots for each freq, one column for each 
% 
% for iHolding = 1:length(holdings)
%     hFig = 27+iHolding;
%     figure(hFig), clf, set(hFig,'Color', [1 1 1]);
%     
%     holding = holdings(iHolding);
% 
%     for iPulseWidth = 1:length(pulseWidths)
%         pulseWidth = pulseWidths(iPulseWidth);
%         freqHz = results{iPulseWidth}(1).laserPulseFreqHz;
% 
%         res = {results{1,1,1}.proportionPulsesWithSpikeSubperiod}';
%         for f = 1:length(freqHz)            
%             subplot(length(freqHz), length(pulseWidths), (f-1)*length(pulseWidths)+iPulseWidth);
%             
%             for s = 1:length(subPeriods)
%                 propValues = arrayfun(@(c) res{c}{f}(s), 1:length(res));
%                 meanEfficacy(iPulseWidth, iHolding, f) = mean(propValues);
%                 semEfficacy(iPulseWidth, iHolding, f) = std(propValues) / sqrt(length(propValues));
%             end
%             
%         end
%         
%      
%     end
% end
% % figure(hFig), clf, set(hFig,'Color', [1 1 1]);
