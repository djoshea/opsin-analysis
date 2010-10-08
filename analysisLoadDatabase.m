protocolList = {'vc test', 'cciv', '1000ms vc', '1000ms cc', '1ms', '2+3ms', ...
                '1ms train', '2ms train', '5ms train', 'white light'};

% specify acceptable values for some fields
valid.rig = {'Dynamic', 'Leica'};
valid.construct = {'T/R', 'A', 'Chief'};
valid.ndFilter = { '0', '0.2', '0.4' };
valid.protocol = protocolList;
valid.date.dateFormat = 'yyyy.mm.dd';
valid.holding = {'-60', '-70'};

% assign meaning to the fields in the database
inclusionField = 'include'; % boolean inclusion flag
% certain fields are specific and consistent to the cell, other so the file
fieldsWithCell = {'date', 'patcher', 'rig', 'mouseID', 'construct', 'cellNum', 'sliceNum'};
fieldsWithFile = {'fileNum', 'ndFilter', 'clamp', 'notes', 'holding'};

par.valid = valid;
par.inclusionField = inclusionField;
par.fieldsWithCell = fieldsWithCell;
par.fieldsWithFile = fieldsWithFile;
par.protocolList = protocolList;

% load the databases from each rig
fnames = {};
fnames{1} = '/Users/djoshea/Documents/Research/dLab/Opsin Comparison/Data/opsin project leica.csv';
fnames{2} = '/Users/djoshea/Documents/Research/dLab/Opsin Comparison/Data/opsin project dynamic.csv';
[fileList cellList] = loadABFDatabase(fnames, par);

cellDataFname = '/Users/djoshea/Documents/Research/dLab/Opsin Comparison/Data/opsin project cell list.csv';
cellData = loadGridDatabase(cellDataFname, 'skipHeaderLines', 1);
cellInclude = strcmp({cellData.Include},'y');

fprintf('\n');
nCells = length(cellList);
if(nCells ~= length(cellData))
    fprintf('Cell counts do not match!\n');
else
    for i = 1:nCells
        % ensure names match
        if(~strcmp(cellList(i).name, cellData(i).Cell_Name))
            fprintf('Error: Cell %d name "%s" differs from "%s"\n', i, cellList(i).name, cellData(i).Cell_Name);
        end
        cellList(i).include = cellInclude(i);
    end
end

% mark cells for inclusion using cell health database, print health summary by construct 
fprintf('Cell Inclusion Summary: Including %d / %d cells\n', nnz(cellInclude), nCells);
constructList = unique({cellList.construct});
nConstructs = length(constructList);
constructCmap = [0.4 1 0.4; 0.6 0.47 1];
[~, constructNum] = ismember({cellList.construct}, constructList);
for cn = 1:nConstructs
    nHealthy(cn) = nnz(constructNum == cn & [cellList.include]);
    nTotal(cn) = nnz(constructNum == cn);
    fprintf('\tConstruct %s: \t%2d / %2d healthy cells\n', constructList{cn}, nHealthy(cn), nTotal(cn));
end


