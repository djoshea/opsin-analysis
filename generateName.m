function str = generateName(data, type, varargin)
% generates a name from the field values in data for a specific purpose
% the point of this function is to put name generation all into one file

basePath = '/Users/djoshea/Documents/Research/dLab/Opsin Comparison';
dbDateFormat = 'yyyy.mm.dd';

switch(type)
    case 'cellName'
        str = sprintf('%s %s %s %02d', data.construct, data.date, data.rig, data.cellNum);
    case 'dataFileName'
        [y m d] = datevec(datenum(data.date, dbDateFormat));
        path = sprintf('%s/Data/%s %04d_%02d_%02d', basePath, data.rig, y, m, d);
        if(strcmp(data.rig, 'Dynamic'))
            str = sprintf('%s/%04d_%02d_%02d_%04d.abf', path, y, m, d, data.fileNum);
        elseif(strcmp(data.rig, 'Leica'))
            str = sprintf('%s/%02d%d%02d%03d.abf', path, mod(y,100), m, d, data.fileNum);
        else
            error('Unknown rig "%s"', data.rig);
        end
    case 'protocolShort'
        str = ['protocol' strrep(strrep(data, ' ', '_'),'+','plus')];
    otherwise
        error('Unknown type');
end

end

