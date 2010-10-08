function [list rowsInvalid] = loadGridDatabase(fname, varargin)
% loads in a generic grid databse (like a CSV file) and does basic validation

% valid is a struct where each field's name is a field name in the file's column list
% and the value is a cell array of acceptable strings for this field to take
par.valid = struct(); 
par.filterList = {'*.csv', 'CSV File'};
par.delimiter = ',';
par.skipHeaderLines = 0;
par.colList = []; % specify this to use this list and process the first line as data
par.inclusionField = 'include';
assignargs(par, varargin);

% interactively request for the file if not passed
if(~exist('fname','var') || isempty(fname))
   [file path] = uigetfile(filterList, ...
        'Choose the Database File to Load');
    if isequal(file,0) || isequal(path,0)
        fprintf('Warning: No CSV file chosen. Aborting.');
        list = [];
        allValid = 0;
        return;
    end

    fname = strcat(path,file);
end

fprintf('Loading database from %s...\n', fname);
fid = fopen(fname,'r');

% skip the first lines if requested
for skip = 1:skipHeaderLines
    fgetl(fid);
end

% if not specified already, get list of column names / fields from first line
if(isempty(colList))
    header = fgetl(fid);
    colList = {};
    while(~isempty(strtrim(header)))
        [fldName header] = strtok(strtrim(header), delimiter);
        colList{end+1} = safeFieldName(fldName);
    end
end

list = [];
rowsInvalid = 0;

% loop through each line and place each value into the field for that column
invalidFieldList = {};
lineNum = skipHeaderLines+1;
listIndex = 0;
while true
    lineNum = lineNum + 1;
    listIndex = listIndex+1;
    origln = fgetl(fid);
    ln = origln;
    if ~ischar(ln), break, end

    list(listIndex).lineNum = lineNum;
    if(~isempty(valid))
        list(listIndex).isValid = 1;
    end
    
    % split the line by commas
    for fld = 1:length(colList)
        if(length(ln) >= 2 && strcmp(ln(1:2), ',,'))
            % double delimiters present a problem for this approach to strtok parse usage
            part = '';
            ln = ln(2:end);
        else
            [part ln] = strtok(strtrim(ln), delimiter);
        end
        if(fld == 1 && isempty(part)) % empty lines with just commas?
            break;
        end

        % convert to number if possible
        num = str2double(part);
        if(~isnan(num)), part = num; end

        fldName = colList{fld};
        list(listIndex).(fldName) = part;
    end

    % check the inlcusion field to see if this should get included
    if(isfield(list(listIndex), inclusionField) && ~strcmp('1', num2str(list(listIndex).(inclusionField))))
        listIndex = listIndex - 1;
        continue;
    end
       
    % perform data validation?
    if(~isempty(valid))
        validateFields = fieldnames(valid);
        for v = 1:length(validateFields)
            fldName = validateFields{v};
            validValues = valid.(fldName);
            
            fldValue = num2str(list(listIndex).(fldName));
            if(iscell(validValues))
                % see if the field value is one of the prescribed list of valid values
                isValid = ismember(fldValue, validValues);
            elseif(isstruct(validValues) && isfield(validValues, 'dateFormat'))
                % test to see if the field value is acceptable given the dateFormat
                try
                    datenum(fldValue, validValues.dateFormat);
                    isValid = 1;
                catch
                    isValid = 0;
                end
            else
                error('Unrecognized validation rule on field %s', fldName);
            end
            
            % check the validity of this field's value against valid struct
            if(~isValid)
                % invalid value: display warning
                fprintf('Warning: Line %4d: Invalid "%s" value "%s"\n', ...
                    lineNum, fldName, num2str(fldValue));
                
                % add to list of invalid fields
                if(~ismember(fldName, invalidFieldList))
                    invalidFieldList{end+1} = fldName;
                end
                
                % mark row as invalid
                list(listIndex).isValid = 0;
                
                % increment the invalid rows counter
                rowsInvalid = rowsInvalid + 1;
            end
        end
    
        if(~isempty(invalidFieldList))
            fprintf('\n%d invalid rows. Expected values for invalid fields:\n', rowsInvalid);
            invalidFieldList = sort(invalidFieldList);
            for i = 1:length(invalidFieldList)
                validValues = valid.(invalidFieldList{i});
                fprintf('%s:\n', invalidFieldList{i})
                if(isstruct(validValues) && isfield(validValues, 'dateFormat'))
                    fprintf('\tDate matching format %s\n', validValues.dateFormat);
                elseif(iscell(validValues))
                    cellfun(@(str) fprintf('\t%s\n', str), validValues);
                end
                fprintf('\n');
            end
        end
    end
end

fclose(fid);


end

