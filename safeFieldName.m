function safe = safeFieldName( str )

remove = {'(', ')', '?'};
replace = {' ', '_'; '=', 'eq'; '.', '_'; '%', 'percent'; '-', 'to'};

safe = strtrim(str);

for i = 1:length(remove)
    safe = strrep(safe, remove{i}, '');
end

for i = 1:size(replace,1)
    safe = strrep(safe, replace{i,1}, replace{i,2});
end

% prefix with v if it starts with a number
if(~isnan(str2double(safe(1))))
    safe = strcat('v', safe);
end

end

