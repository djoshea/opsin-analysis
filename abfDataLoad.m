function [d sis h chMatches time fname fnamelast] = abfDataLoad(fnameOrData, chMatchRegexes)
% [d sis h fname fnamelast] = abfDataLoad(fnameOrData)
% either loads data from fname or extracts it from a struct
% sis is si but in seconds rather than microseconds
% chMatchRegexes is a cell array of string regexes that can retrieve channels whose names match, e.g. "IN 0" or "IN [02]"

% load the file if necessary or extract the data form the struct
if(exist('fnameOrData','var') && ischar(fnameOrData) && ~isempty(fnameOrData))
    % load data directly
    fname = fnameOrData;
    [d si h] = abfload(fname);
    
elseif(exist('fnameOrData','var') && ~isempty(fnameOrData))
    % data and si included as first 2 args
    d = fnameOrData.d;
    si = fnameOrData.si;
    h = fnameOrData.h;
    fname = [];
    
else
   % prompt user for file
   [file path] = uigetfile({'*.abf', 'Axon Binary File (*.abf)'}, ...
        'Choose the ABF to Load');
    if isequal(file,0) || isequal(path,0)
        error('No ABF file chosen. Aborting.');
        return % canceled
    end
    fname = strcat(path,file);
    fprintf('File: %s\n', fname);
    [d si h] = abfload(fname);
    
end

% some utility variables
sis = si * 10^-6;
fnamelast = fname(max([strfind(fname, '/') strfind(fname, '\')])+1:end);
time = (0:size(d,1)-1) * si*10^-6;

% find channels by names or regexes 
findChannel = @(nameRegex) find(cellfun(@(chName) ~isempty(regexp(chName, nameRegex)), h.recChNames));
chMatches = cell(length(chMatchRegexes),1);
for i = 1:length(chMatchRegexes)
    chMatchInd = findChannel(chMatchRegexes{i});
    if(isempty(chMatchInd))
        error('Could not find channel matching regex "%s". Check the channel names specified.', chMatchRegexes{i});
    end
    chMatches{i} = squeeze(d(:,chMatchInd,:));
end


