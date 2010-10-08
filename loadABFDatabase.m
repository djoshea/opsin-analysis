function [fileList cellList] = loadABFDatabase(fnames, varargin)
% takes in the ATF list CSV file with columns as below and converts it into a struct array
% colList = {'date', 'patcher', 'rig', 'mouseID', 'construct', 'cellNum', 'fileNum',...
%            'sliceNum', 'ndFilter', 'clampType', 'protocol', 'include', 'notes'};

if(~exist('fnames', 'var'))
    fnames = [];
end

par.protocolList = {};
par.valid = struct([]);
par.inclusionField = 'include';
par.fieldsWithCell = {};
par.fieldsWithFile = {};
par.filterList = {'*.csv', 'CSV File'};
assignargs(par, varargin);

if(isempty(protocolList) && isfield(valid, 'protocol'))
    protocolList = valid.protocol;
end

if(isempty(fnames))
    % interactively request for the file(s) if not passed
   [fnamesPart path] = uigetfile(filterList, ...
        'Choose the Database Files to Load','MultiSelect', 'on');
    if isequal(fnamesPart,0) || isequal(path,0)
        fprintf('Warning: No CSV file chosen. Aborting.\n');
        fileList = [];
        cellList = [];
        return;
    end
    
    if(~iscell(fnamesPart))
        fnamesPart = {fnamesPart};
    end
    
    for i = 1:length(fnamesPart)
        fnames{i} = strcat(path,fnamesPart{i});
        
    end
end

if(~iscell(fnames))
    fnames = {fnames};
end

fileList = struct([]);
cellList = struct([]);
rowsInvalid = 0;
filesMissing = 0;

for iFile = 1:length(fnames) 
    % load in the csv database
    [fileListPart rowsInvalidPart] = loadGridDatabase(fnames{iFile}, 'delimiter', ',', 'valid', valid);
    
    rowsInvalid = rowsInvalid + rowsInvalidPart;

    % build file names for each row in this list and check for existence
    filesMissingPart = 0;
    for i = 1:length(fileListPart)
        try % in case there's an error generating the file name (from datenum usage)
            fileRowName = generateName(fileListPart(i), 'dataFileName');
            fileListPart(i).fname = fileRowName;
            fileListPart(i).fileExists = exist(fileRowName, 'file');
        catch
            fileListPart(i).fname = '';
            fileListPart(i).fileExists = 0;
        end

        if(~fileListPart(i).fileExists)
            filesMissingPart = filesMissingPart+1;
            fprintf('Warning: Line %3d: Could not locate file %s\n', fileListPart(i).lineNum, fileRowName);
            fileListPart(i).fileExists = 0;
        end
    end
    
    if(filesMissingPart)
        fprintf('\nWarning: Including %d rows with missing ABF files. Check database accuracy and correct path generation.\n\n', filesMissingPart);
    end
    
    filesMissing = filesMissing + filesMissingPart;
    
    % now we build the list of cells, attaching protocols to each cell
    % or creating a new cell as we go along 
    for i = 1:length(fileListPart)
        row = fileListPart(i);
        
        if(~strcmp(num2str(row.(inclusionField)),'1') || ~row.isValid)
            continue;
        end
        
        name = generateName(row, 'cellName');
        fileListPart(i).cellName = name;
        
        if(isempty(cellList) || ~ismember(name, {cellList.name}))
            % add this cell's base info to the list
            newcell = struct();
            newcell.name = name;
            newcell.protocolList = {};
            
            % copy over fields that are cell-specific
            for f = 1:length(fieldsWithCell)
                newcell.(fieldsWithCell{f}) = row.(fieldsWithCell{f});
            end
            
            % copy empty fields for each protocol in list (so that the array is self-consistent)
            for p = 1:length(protocolList)
                pShort = generateName(protocolList{p}, 'protocolShort');
                newcell.(pShort) = [];
            end
            
            if(isempty(cellList))
                cellList = newcell;
            else
                cellList(end+1) = newcell;
            end
            ic = length(cellList);
        else
            ic = find(strcmp({cellList.name}, name));
        end
        
        % add this to the list of found protocols
        protocolShort = generateName(row.protocol, 'protocolShort');
        cellList(ic).protocolList = unique(cat(1,cellList(ic).protocolList,{row.protocol}));
        
        % copy over fields that are file-specific
        thisFile = [];
        for f = 1:length(fieldsWithFile)
            thisFile.(fieldsWithFile{f}) = row.(fieldsWithFile{f});
        end
        % and include the ABF name as well
        thisFile.fname = generateName(row, 'dataFileName');
        
        % now add this row's data to the list as an example of the particular protocol
        if(isempty(cellList(ic).(protocolShort)))
            cellList(ic).(protocolShort) = thisFile;
        else
            cellList(ic).(protocolShort)(end+1) = thisFile;
        end
        
    end
    
    % concatenate this onto the full fileList
    if(isempty(fileList))
        fileList = fileListPart;
    else
        fileList = [fileList fileListPart];
    end
    
end

fprintf('\nSummary: Loaded %d records from %d database files comprising %d cells\n', length(fileList), length(fnames), length(cellList));

if(length(fnames) > 1)
    % print multiple file error count summaries
    if(rowsInvalid)
        fprintf('Warning: Excluding %d invalid rows found within ABF Database files. Check database accuracy and validation rules.\n', rowsInvalid);
    end
    if(filesMissing)
        fprintf('Warning: Including %d rows with missing ABF files. Check database accuracy and correct path generation.\n', filesMissing);
    end
end



end

