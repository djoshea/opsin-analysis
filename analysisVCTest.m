%% Process VC tests

ncells = length(cellList);
vcTestShort = generateName('vc test', 'protocolShort');
data = [];

for c = 1:ncells
    vctest = cellList(c).(vcTestShort)(1);
    if(isempty(vctest))
        fprintf('Warning: Cell %s has no vc test protocol\n', cellList(c).name);
        continue;
    end
    
    % override the vc test voltage delta based on the known value used on each rig
    if(strcmp(cellList(c).rig, 'Leica'))
        voltageDelta = -5;
    elseif(strcmp(cellList(c).rig, 'Dynamic'))
        voltageDelta = -10;
    end
    
    abfname = vctest.fname;
    
    out = processVCTest(abfname, 'voltageDelta', voltageDelta, 'holding', vctest.holding, 'showPlots', 0, 'yl', [-500 10]);
    if(isempty(data))
        data = out;
    else
        data(end+1) = out;
    end
end

%% Determine cell health

threshRa = 30;
threshLeak = 200;
threshResting = -50;

Ra = [data.accessResistance]';
Rm = [data.membraneResistance]';
leak = -[data.meanBaseline]';
resting = [data.resting]';
construct = {cellList.construct}';
constructList = unique(construct);
nConstructs = length(constructList);
constructCmap = [0.4 1 0.4; 0.6 0.47 1];
[~, constructNum] = ismember(construct, constructList);

filterAccessResistance = 1;
filterLeak = 0;
filterResting = 1;
scatterAxes = {Ra, resting};
scatterAxesNames = {'Access Resistance (MOhms)', 'Resting Potential (mV)'};
scatterAxesThreshDirection = {-1, -1}; % negative means below, postive means above is healthy
scatterAxesThresh = {threshRa, threshResting};

% filter for cell health
healthy = ones(size(Ra));
if(filterAccessResistance)
    healthy = healthy & Ra <= threshRa;
end
if(filterLeak)
    healthy = healthy & leak <= threshLeak;
end
if(filterResting)
    healthy = healthy & resting <= threshResting;
end

%% Plot summary scatter plot

figure(2), clf, hList = [];
for c = 1:ncells
    color = constructCmap(constructNum(c), :);
    if(healthy(c))
        facecolor = color;
    else
        facecolor = 'w';
    end
    
    hList(constructNum(c)) = plot(scatterAxes{1}(c), scatterAxes{2}(c), ...
        's', 'Color', color, 'MarkerFaceColor', facecolor, 'MarkerSize', 5);
    hold on
end

% draw health thresholds
axis tight
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
if(scatterAxesThreshDirection{2} < 0)
    plot(scatterAxesThresh{1}*ones(1,2), [yl(1) scatterAxesThresh{2}], '--', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
else
    plot(scatterAxesThresh{1}*ones(1,2), [scatterAxesThresh{2} yl(2)], '--', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
end
if(scatterAxesThreshDirection{1} < 0)
    plot([xl(1) scatterAxesThresh{1}], scatterAxesThresh{2}*ones(1,2), '--', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
else
    plot([scatterAxesThresh{1} xl(2)], scatterAxesThresh{2}*ones(1,2), '--', 'LineWidth', 1, 'Color', [0.4 0.4 0.4]);
end

xlabel(scatterAxesNames{1})
ylabel(scatterAxesNames{2});
box off

legend(hList, constructList, 'Location', 'Best');
legendboxoff

title('Cell Health Summary');

%% Dump full data table

fprintf('%28s\t%8s\t%8s\t%8s\t%8s\t%8s\n', 'Cell Name', 'Ra', 'Rm', 'Resting', 'Leak', 'Healthy');
for c = 1:ncells
    fprintf('%28s\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8d\n', ...
        cellList(c).name, Ra(c), Rm(c), resting(c), leak(c), healthy(c));
end
fprintf('\n');

%% Dump summary

fprintf('Cell Health Summary:\n');
nHealthy = zeros(nConstructs, 1);
nTotal = zeros(nConstructs, 1);
for cn = 1:nConstructs
    nHealthy(cn) = nnz(constructNum == cn & healthy);
    nTotal(cn) = nnz(constructNum == cn);
    fprintf('Construct %s: %2d / %2d healthy cells\n', constructList{cn}, nHealthy(cn), nTotal(cn));
end

fprintf('\nPatcher Summary:\n');
patcherList = unique({cellList.patcher});
nPatchers = length(patcherList);
[~, patcherNum] = ismember({cellList.patcher}', patcherList);
nHealthy = zeros(nPatchers, 1);
nTotal = zeros(nPatchers, 1);
for cn = 1:nPatchers
    nHealthy(cn) = nnz(patcherNum == cn & healthy);
    nTotal(cn) = nnz(patcherNum == cn);
    fprintf('Patcher %s: %2d / %2d healthy cells\n', patcherList{cn}, nHealthy(cn), nTotal(cn));
end


