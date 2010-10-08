function [usableCellList results] = analysisLightTrains(cellList, cellData, varargin)
% analyze all light train data matching a template

def.pulseWidths = [5 2 1];
% ndFilters = [0 0.2 0.4];
def.ndFilters = [0];
def.holdings = [-60 -70];
def.includeList = [];
def.printData = 0;
def.savePlots = 0;
def.showPlots = 1;
def.printCellDataFields = 0;
def.printIncludedOnly = 0;
            
def.subPeriods = {1:10, 11:20, 21:30, 31:40};

def.metrics = {'proportionPulsesWithSpike', 'proportionPulsesWithMultipleSpikes'};
def.metricNames = {'Pulse Efficacy', 'Multiple Spike Likelihood'};
def.metricMultiplier = [40, 40];
def.metricLimits = [0 40; 0 40];
def.metricPlot = [];
def.metricPrintDigits = [];

assignargs(def, varargin);

if(isempty(includeList))
    includeList = ones(size(cellList));
end
    

for iPulseWidth = 1:length(pulseWidths)
    pulseWidth = pulseWidths(iPulseWidth);
    for iNd = 1:length(ndFilters)
        nd = ndFilters(iNd);
        for iHolding = 1:length(holdings)
            holding = holdings(iHolding);
            
            par.pulseWidth = pulseWidth;
            par.ndFilter = nd;
            par.holding = holding;
            par.cellData = cellData;
            par.subPeriods = subPeriods;
            
            par.printData = printData;
            par.printCellDataFields = printCellDataFields;
            par.includeList = includeList;
            par.printIncludedOnly = printIncludedOnly;
            
            par.savePlots = savePlots;
            par.showPlots = showPlots;
                       
            par.metrics = metrics;
            par.metricNames = metricNames;
            par.metricMultiplier = metricMultiplier;
            par.metricLimits = metricLimits;
            par.metricPrintDigits = metricPrintDigits;
            par.metricPlot = metricPlot;
            
%             par.metrics = {'proportionPulsesWithSpike', 'proportionPulsesWithMultipleSpikes'};
%             par.metricNames = {'Pulse Efficacy', 'Multiple Spike Likelihood'};
%             par.metricMultiplier = [40, 40];
%             par.metricLimits = [0 40; 0 40];
%             
%             par.metrics = {'spikesElicited'};
%             par.metricNames = {'Total Spikes Elicited'};
%             par.metricMultiplier = [1];
%             par.metricLimits = [0 200];
%             
%             par.metrics = {'meanSpikeLatency', 'stdSpikeLatency'};
%             par.metricNames = {'Spike Latency', 'Spike Jitter'};
%             par.metricMultiplier = [1000, 1000];
%             par.metricLimits = [0 8; 0 3];
%             par.metricPrintDigits = [3 3];
            
%             par.metrics = {'spontaneousSpikes'};
%             par.metricNames = {'Spontaneous Spikes'};
%             par.metricMultiplier = 1;
%             par.metricLimits = [0 10];

%             par.metrics = {'spikesInWindow'};
%             par.metricNames = {'Spikes In Window'};
            
            [usableCellList{iPulseWidth, iNd, iHolding} results{iPulseWidth, iNd, iHolding}] = ...
                doAnalysisLightTrains(cellList, par);
        end
    end
end

end

function [usableCellList results] = doAnalysisLightTrains(cellList, varargin)

def.pulseWidth = 2;
def.ndFilter = 0;
def.holding = -60;
def.plotPath = '/Users/djoshea/Documents/Research/dLab/Opsin Comparison/Plots/Light Trains/';
def.printData = 0;
def.includeList = [cellList.include];
def.cellData = [];
def.savePlots = 1;
def.showPlots = 1;
def.printCellDataFields = 1;
def.printIncludedOnly = 0;

def.displayMultipleMatchWarning = 0;
def.displayPValueWarning = 0;

def.metrics = {'proportionPulsesWithSpike', 'proportionPulsesWithMultipleSpikes', 'spikesElicited'};
def.metricNames = {'Pulse Efficacy', 'Multiple Spike Likelihood', 'Total Spike Count'};
def.metricMultiplier = [40, 40, 1];
def.metricLimits = [0 1; 0 1; 0 40*5];
def.metricPrintDigits = [];
def.metricPlot = [];
assignargs(def, varargin);

if(isempty(metricPrintDigits))
    metricPrintDigits = zeros(length(metrics), 1);
end

if(isempty(metricPlot))
    metricPlot = ones(size(metrics));
end

if(pulseWidth == 5)
    maxLaserFreq = 100;
else
    maxLaserFreq = 200;
end

elicitedSpikeWindow = [1.5e-3 50e-3];
laserSignalTypes = {'toggle', 'state'};

ttl = sprintf('%d ms Width, %.1f ND, Holding %d', pulseWidth, ndFilter, holding);
protoName = sprintf('protocol%dms_train', pulseWidth);

% fprintf('%s\n', ttl);
usableCellList = [];
results = [];
constructNumList = [];

constructList = unique({cellList.construct});
nConstructs = length(constructList);
constructCmap = [0.4 1 0.4; 0.6 0.47 1];
[~, constructNum] = ismember({cellList.construct}, constructList);

% fprintf('Cell ');
for iCell = 1:length(cellList)
    if(~includeList(iCell))
        continue;
    end
    
    protoList = cellList(iCell).(protoName);
    % find the last protocol matching this cell
    foundMatch = 0;
    for iProto = length(protoList):-1:1
        if(protoList(iProto).ndFilter == ndFilter && ...
           protoList(iProto).holding == holding)
%            fprintf('%d ', iCell);
            % use this cell, try both laserTypes until reasonable pulse lengths are obtained
            if(foundMatch && displayMultipleMatchWarning)
                fprintf('Warning: Multiple matches for cell %d: %s\n', iCell, cellList(iCell).name);
            end
            for iLaserType = 1:length(laserSignalTypes);
                out = processLightTrains(protoList(iProto).fname, ...
                    'showPlots', 0, 'printOutput', 0, 'laserSignalType', laserSignalTypes{iLaserType}, 'elicitedSpikeWindow', elicitedSpikeWindow);
                
                if(out.meanLaserPulseDuration < 0.01 && out.meanLaserPulseDuration > 9e-4)
                    % laser pulse acceptable, keep this one                
                    break;
                end
            end
            
            if(length(out.laserPulseFreqHz) < 7)
                fprintf('\nWarning: cell %d has only %d traces, skipping.\n', iCell, length(out.laserPulseFreqHz));
                continue;
            end
            if(round(max(out.laserPulseFreqHz)) < maxLaserFreq || round(max(out.laserPulseFreqHz)) > maxLaserFreq)
                fprintf('\nWarning: cell %d laser pulse frequency estimates still messed up!\n', iCell);
                continue
            end
            
            if(~foundMatch)
                if(isempty(results))
                    results = out;
                    constructNumList = constructNum(iCell);
                    usableCellList = cellList(iCell);
                else
                    results(end+1) = out;
                    constructNumList(end+1) = constructNum(iCell);
                    usableCellList(end+1) = cellList(iCell);
                end
            end
            
            foundMatch = 1;
            laserPulseFreqHz = out.laserPulseFreqHz;
%             break; % only process one protocol
        end
    end
end

fprintf('\n');

%% Print data by cell along with cellData fields if requested

if(printData)
       infoFields = fieldnames(cellData);
       infoFields = infoFields(~strcmp(infoFields, 'Cell_Name'));
    
    for metricNum = 1:length(metrics)
       fprintf('%s: %s\n', metricNames{metricNum}, ttl);
       
       if(~isempty(results))
           fprintf('Name\t');
       
           if(printCellDataFields)
               for iField = 1:length(infoFields);
                   fprintf('%s\t',infoFields{iField});
               end
           end
       
           for j = 1:length(laserPulseFreqHz)
               fprintf('%.0f\t', laserPulseFreqHz(j));
           end
       end
       fprintf('\n');
       
       if(~isempty(results))
           for iCell = 1:length(cellData)
               inList = find(strcmp({cellData(iCell).Cell_Name}, {usableCellList.name}));
               if(~printIncludedOnly || includeList(iCell))
                   fprintf('%s\t', cellData(iCell).Cell_Name);
                   if(printCellDataFields)
                       for iField = 1:length(infoFields)
                           fprintf('%s', num2str(cellData(iCell).(infoFields{iField})));
                           fprintf('\t');
                       end
                   end
                   
                   if(~isempty(inList))
                       for j = 1:length(laserPulseFreqHz)
                           fprintf(['%.' num2str(metricPrintDigits(metricNum)) 'f\t'], ...
                               metricMultiplier(metricNum)*results(inList).(metrics{metricNum})(j));
                       end
                   end
                   fprintf('\n');
               end
           end
       end
       fprintf('\n');
    end
end

%% Plot individual and summary plots for each metric

if(showPlots && ~isempty(results))
    % build legend for construct with n= value
    for cn = 1:nConstructs
        constructLegend{cn} = sprintf('%s (n=%d)', constructList{cn}, nnz(constructNumList == cn));
    end
    
    for metricNum = 1:length(metrics)
        if(~metricPlot(metricNum))
            continue;
        end
        
        % draw individual cells doublet plot
        hfigEfficacy = figure(12+metricNum); clf;
        set(hfigEfficacy, 'Name', [metricNames{metricNum} ' Individuals'], 'NumberTitle', 'off', 'Color', [1 1 1]);
        %     set(gcf, 'Position', [ 2031         445         825         439]);
        ax = gca;
        topOffset = 0;
        botOffset = 0;
        
        mMax = metricLimits(metricNum, 2);
        mMin = metricLimits(metricNum, 1);
        mSpan = mMax-mMin;
        
        for cn = 1:nConstructs
            h(cn) = plot(ax, -1, -1, 's-',  'Color', constructCmap(cn, :), 'MarkerFaceColor', constructCmap(cn, :));
            hold on
        end
        legend(constructLegend, 'Location', 'EastOutside');
        legendboxoff
        box off
        xlim([0 10+maxLaserFreq]);
        ylim([mMin+botOffset mMax+topOffset]);
        
        for iCell = 1:length(results)
            prop = metricMultiplier(metricNum)*results(iCell).(metrics{metricNum});
            
            if(nnz(prop > 0.995*mSpan + mMin))
                topOffset = topOffset + 0.004*mSpan;
                prop(prop > 0.995*mSpan + mMin) = prop(prop > 0.995*mSpan + mMin) + topOffset;
                ylim([mMin+botOffset mMax+topOffset]);
            end
            if(nnz(prop < 0.005*mSpan + mMin))
                botOffset = botOffset - 0.004*mSpan;
                prop(prop < 0.005*mSpan + mMin) = prop(prop < 0.005*mSpan + mMin) + botOffset;
                ylim([mMin+botOffset mMax+topOffset]);
            end
            
            % add this cell to the plot
            plot(ax, out.laserPulseFreqHz, prop,...
                's-', 'Color', constructCmap(constructNumList(iCell), :), 'MarkerFaceColor', constructCmap(constructNumList(iCell), :));
            hold on
            drawnow
        end
        
        xlim([0 10+maxLaserFreq]);
        xlabel('Pulse Frequency (Hz)');
        ylabel(metricNames{metricNum});
        box off
        plotTitle = [metricNames{metricNum} ' Individual - ' ttl];
        title(plotTitle);
        if(savePlots)
            saveas(gcf, [plotPath plotTitle '.png'], 'png');
        end
        
        % draw summary plot
        prop = [];
        for cn = 1:nConstructs
            r = results(constructNumList == cn);
            prop{cn} = metricMultiplier(metricNum)*reshape([r.(metrics{metricNum})], [], nnz(constructNumList==cn))';
            
            if(nnz(isnan(prop{cn})))
                for f = 1:size(prop{cn},2)
                    meanProp{cn}(f) = mean(prop{cn}(~isnan(prop{cn}(:, f)),f));
                    semProp{cn}(f) = std(prop{cn}(~isnan(prop{cn}(:, f)),f)) / sqrt(nnz(~isnan(prop{cn}(:,f))));
                end
            else
                meanProp{cn} = mean(prop{cn},1);
                semProp{cn} = std(prop{cn}, 1) / sqrt(nnz(constructNumList == cn));
            end
        end
        
        if(length(unique(constructNumList)) >= 2)
            freqList = out.laserPulseFreqHz;
            for f = 1:length(freqList)
                [~, pValues(f)] = ttest2(prop{1}(:,f), prop{2}(:,f));
            end
            if(min(pValues) < 0.05 && displayPValueWarning)
                fprintf('Significant p-value encountered!\n');
            end
            
            hfigMeans = figure(22+metricNum); clf;
            set(hfigMeans, 'Name', [metricNames{metricNum} ' Summary'], 'NumberTitle', 'off', 'Color', [1 1 1]);
            
            title(sprintf('Pulse Efficacy, %d ms Width, %d ND, Holding %d', pulseWidth, ndFilter, holding));
            h = [];
            for cn = 1:nConstructs
                h(cn) = errorbar(freqList, meanProp{cn}, semProp{cn}, 's-', 'Color', constructCmap(cn, :), 'MarkerFaceColor', constructCmap(cn,:));
                hold on
            end
            ylim([mMin mMax]);
            xlim([0 10+maxLaserFreq]);
            legend(h, constructLegend, 'Location', 'EastOutside');
            legendboxoff
            box off
            ax = gca;
            title(ttl);
            %         set(gcf, 'Position', [2530         191         825         439]);
            
            for iP = 1:length(pValues)
                if(pValues(iP) < 0.005)
                    text(freqList(iP)-5, meanProp{cn}(iP)+semProp{cn}(iP)+0.05, 'p<.005');
                elseif(pValues(iP) < 0.01)
                    text(freqList(iP)-5, meanProp{cn}(iP)+semProp{cn}(iP)+0.05, 'p<.01');
                elseif(pValues(iP) < 0.05)
                    text(freqList(iP)-5, meanProp{cn}(iP)+semProp{cn}(iP)+0.05, 'p<.05');
                end
            end
            
            xlabel('Pulse Frequency (Hz)');
            ylabel(metricNames{metricNum});
            
            plotTitle = [metricNames{metricNum} ' Summary - ' ttl];
            title(plotTitle);
            if(savePlots)
                saveas(gcf, [plotPath plotTitle '.png'], 'png');
            end
        end
        
    end
end

end
