function safe = safeFieldName( str )

remove = {'(', ')', '?'};
replace = {' ', '_'; '=', 'eq'};

safe = strtrim(str);

for i = 1:length(remove)
    safe = strrep(safe, remove{i}, '');
end

for i = 1:size(replace,1)
    safe = strrep(safe, replace{i,1}, replace{i,2});
end

end

