ncells = length(cellList);
vcTestShort = generateName('1000ms vc', 'protocolShort');
data = [];
cellListUsable = [];

healthyList = [];

for c = 1:ncells
    prs = cellList(c).(vcTestShort);
    
    proto = [];
    for p = 1:length(prs)
        if(prs(p).ndFilter == 0 & prs(p).holding==-70)
            proto = prs(p);
        end
    end
    if(isempty(proto))
        fprintf('Warning: Cell %s has no 1000ms vc protocol\n', cellList(c).name);
        continue;
    end
    
    if(isempty(cellListUsable))
        cellListUsable = cellList(c);
    else
        cellListUsable(end+1) = cellList(c);
    end
        
    abfname = proto.fname;
    
    if(strcmp(cellList(c).rig, 'Leica'))
        laserSignalType = 'state';
    else
        laserSignalType = 'toggle';
    end
    out = processLightVC(abfname, 'laserSignalType', laserSignalType, 'showPlots', 1);
    if(isempty(data))
        data = out;
    else
        data(end+1) = out;
    end
    
    if(isempty(healthyList))
        healthyList = healthy(c);
    else
        healthyList(end+1) = healthy(c);
    end
end

return

%% Plot Summary scatter plot
ncells = length(data);
tauDeactivation = [data.tauInactivation]';
steadyState = [data.steadyState]';
construct = {cellListUsable.construct}';
constructList = unique(construct);
nConstructs = length(constructList);
constructCmap = copper(nConstructs);
[~, constructNum] = ismember(construct, constructList);

figure(2), clf, hList = [];
for c = 1:ncells
    color = constructCmap(constructNum(c), :);
    if(~healthyList(c))
        continue;
    end
%     if(healthy(c))
%         facecolor = color;
%     else
%         facecolor = 'w';
%     end
    hList(constructNum(c)) = plot(tauDeactivation(c)*1000, steadyState(c), ...
        's', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 5);
    hold on
end

xlabel('Tau Deactivation (ms)')
ylabel('Peak Current (pA)');
box off

legend(hList, constructList, 'Location', 'Best');
legendboxoff

title('1000ms VC Summary');

%%

fprintf('Cell Health Summary:\n');
meanSteadyState = zeros(nConstructs, 1);
stdSteadyState = zeros(nConstructs, 1);
semSteadyState = zeros(nConstructs, 1);
for cn = 1:nConstructs
    meanSteadyState(cn) = mean(steadyState(constructNum == cn));
    stdSteadyState(cn) = std(steadyState(constructNum == cn));
    semSteadyState(cn) = stdSteadyState(cn) / sqrt(nnz(constructNum == cn));
    
    fprintf('Construct %s: %.2f +- %.2f pA SteadyState Photocurrent\n', constructList{cn}, meanSteadyState(cn), stdSteadyState(cn));
end

