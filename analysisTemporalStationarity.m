%% Gather the data by pulse width
pulseWidths = [5 2 1];
holdings = [-60 -70];
def.ndFilters = 0;

def.pulseWidths = pulseWidths;
def.holdings = holdings;
def.includeList = [cellList.include];
def.printData = 0;
def.savePlots = 0;
def.showPlots = 0;
def.printCellDataFields = 0;
def.metrics = {};

subPeriods = {1:10, 11:20, 21:30, 31:40};
def.subPeriods = subPeriods;

[usableCellList, results] = analysisLightTrains(cellList, cellData, def);

%% Print the results: one line for each cell

% first header row, print pulse widths and holding
fprintf('\t'); % skip column for name
for iHolding = 1:length(holdings)
    holding = holdings(iHolding);
    for iPulseWidth = 1:length(pulseWidths)
        pulseWidth = pulseWidths(iPulseWidth);
        freqHz = results{iPulseWidth,1,iHolding}(1).laserPulseFreqHz;
        fprintf('%dms, %dmV', pulseWidth, holding);
        
        for f = 1:length(freqHz)
            for s= 1:length(subPeriods)
                fprintf('\t');
            end
        end
    end
end
fprintf('\n');

% second header row
fprintf('Name\t');
for iHolding = 1:length(holdings)
    holding = holdings(iHolding);
     for iPulseWidth = 1:length(pulseWidths)
        pulseWidth = pulseWidths(iPulseWidth);
        freqHz = results{iPulseWidth,1,iHolding}(1).laserPulseFreqHz;
        for f = 1:length(freqHz)
            for s= 1:length(subPeriods)
                fprintf('%.0fHz P%d\t', freqHz(f), s);
            end
        end
    end
end
fprintf('\n');
    
% data rows, one for each cell, protocols concatenated horizontally
for iCell = 1:length(cellData)
    fprintf('%s\t',cellData(iCell).Cell_Name);
   
    for iHolding = 1:length(holdings)
        holding = holdings(iHolding);
         for iPulseWidth = 1:length(pulseWidths)
            pulseWidth = pulseWidths(iPulseWidth);
            freqHz = results{iPulseWidth,1,iHolding}(1).laserPulseFreqHz;
            inList = find(strcmp(cellData(iCell).Cell_Name, {usableCellList{iPulseWidth,1,iHolding}.name}));
            if(~isempty(inList))
                result = results{iPulseWidth,1,iHolding}(inList(1));
                for f = 1:length(freqHz)
                    for subperiod = 1:length(subPeriods)
                        fprintf('%2d\t', result.proportionPulsesWithSpikeSubperiod{f}(subperiod)*10);
                    end
                end
            end
        end
    end
    fprintf('\n');
end

%% Plot the results: one row of plots for each freq, one column for each 

meanEfficacy = cell(nConstructs,1);
semEfficacy = cell(nConstructs,1);
propValues = cell(nConstructs,1);

nFreqAll = length(results{3,1,1}(1).laserPulseFreqHz);

for iHolding = 1:length(holdings)
    hFig = 27+iHolding;
    figure(hFig), clf, set(hFig,'Color', [1 1 1]);
    set(hFig,'Position',[680 44 916 1034]);
    
    holding = holdings(iHolding);

    for iPulseWidth = 1:length(pulseWidths)
        pulseWidth = pulseWidths(iPulseWidth);
        freqHz = results{iPulseWidth}(1).laserPulseFreqHz;

        resultsFiltered = {results{iPulseWidth,1,iHolding}.proportionPulsesWithSpikeSubperiod}';
        
        for cn = 1:nConstructs
            cellsThisConstruct = strcmp(constructList{cn}, {usableCellList{iPulseWidth,1,iHolding}.construct});
            resultsFilteredByConstruct{cn} = resultsFiltered(cellsThisConstruct);
        end
        
        for freq = 1:length(freqHz)            
            subplot(nFreqAll, length(pulseWidths), (freq-1)*length(pulseWidths)+iPulseWidth);
            
            for sub = 1:length(subPeriods)
                for cn = 1:nConstructs
                    propValues{cn} = arrayfun(@(iCell) resultsFilteredByConstruct{cn}{iCell}{freq}(sub), 1:length(resultsFilteredByConstruct{cn}));
                    meanEfficacy{cn,iPulseWidth, iHolding, freq}(sub) = mean(propValues{cn});
                    semEfficacy{cn, iPulseWidth, iHolding, freq}(sub) = std(propValues{cn}) / sqrt(length(propValues{cn}));
                    nval(cn) = length(propValues{cn});
                end
                
                [~, pValues{iPulseWidth,iHolding,freq}(sub)] = ttest2(propValues{1},propValues{2});
            end
            
            for cn = 1:nConstructs
                errorbar(1:length(subPeriods), meanEfficacy{cn,iPulseWidth, iHolding, freq}, ...
                    semEfficacy{cn,iPulseWidth, iHolding, freq}, '-', 'Color', constructCmap(cn,:), 'LineWidth', 2);
                hold on
            end
            box off
            ylim([0 1]);
            xlim([0.5 4.5]);
            
            nlabel = sprintf('%s=%d,%s=%d',constructList{1},nval(1),constructList{2},nval(2));
            text(3,0.15,nlabel);
            
            if(iPulseWidth==1)
               ylabel(sprintf('%d Hz', freq)); 
            end
            if(freq==1)
                title(sprintf('%d ms Pulse', pulseWidth));
            end
            
            for sub = 1:length(subPeriods)
                pv = pValues{iPulseWidth,iHolding,freq}(sub);
                if(pv < 0.005)
                    text(sub-0.3, 1, 'p<.005');
                elseif(pv < 0.01)
                    text(sub-0.3, 1, 'p<.01');
                elseif(pv < 0.05)
                    text(sub-0.3, 1, 'p<.05');
                end
            end
        end
      
    end
    h = uicontrol('Style','text','String',sprintf('Temporal Stationarity, Holding = %d mV',holding));
    set(h,'Units','normalized','Position',[0.1 0.95 0.8 0.02],'BackgroundColor', [1 1 1],'FontSize',12);
end
% figure(hFig), clf, set(hFig,'Color', [1 1 1]);
