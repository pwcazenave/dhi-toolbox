function dfs02nc(filename, varargin)
% dfs02c(filename) converts a dfs0 into a netCDF file.
%
% Optionally supply argument pair:
%   'overwrite' - set to 'yes' to overwrite existing files
%
% Pierre Cazenave 2017/03/14 v1.0 - Based on dfs02nc.
%

if nargin ~= 1 && nargin ~= 3
    error('Incorrect number of arguments specified. Please supply the complete path to a dfs0 file.')
end

clobber = false;
for aa = 1:2:nargin - 1
    arg = varargin{aa};
    switch lower(arg)
        case 'overwrite'
            if strcmpi(varargin{aa + 1}, 'yes') || strcmpi(varargin{aa + 1}, 'y')
                clobber = true;
            end
    end
end

% Create output file name.
[fpath, fname, ~] = fileparts(filename);
fout = fullfile(fpath, [fname, '.nc']);
clear fpath fname

if exist(fout, 'file') == 2 && ~clobber
    warning('File %s exists. Not overwriting.', fout)
    return
end

fprintf('Saving to %s\n', fout)

% Load the MATLAB interface for dfs files.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfs0.*;

% Open the file for the times.
dfs0 = dfsTSO(filename);
all_time = readTimes(dfs0);
time.Times = datestr(all_time, 'yyyy-mm-ddTHH:MM:SS.FFF');
time.time = datenum(all_time);
clear all_time
dfs0.close()

% Reopen for the rest of the data.
dfs0 = DfsFileFactory.DfsGenericOpen(filename);
% Grab all the data. Since this is time series data, the first column is
% the time (in units of?), the remaining columns are the time series.
dd = double(Dfs0Util.ReadDfs0DataDouble(dfs0));
all_data = dd(:, 2:end);
clear dd

% We need the metadata to figure out how to partition everything up for the
% netCDF variables.
for i = 0:dfs0.ItemInfo.Count - 1
   item = dfs0.ItemInfo.Item(i);
   items.name{i+1} = char(item.Name);
   items.unit{i+1} = char(item.Quantity.Unit);
   items.abbreviation{i+1} = char(item.Quantity.UnitAbbreviation);
end
% Find the unique names for the metadata.
vars.name = {};
vars.units = {};
vars.abbreviation = {};
for i = 1:length(items.name)
    name = strsplit(items.name{i}, ': ');
    name = name(end);
    % Replace illegal characters for the netCDF variable.
    name = strrep(name{1}, ' ', '_');
    name = strrep(name, ':', '-');
    name = strrep(name, '(', '');
    name = strrep(name, ')', '');
    name = strrep(name, '\.', '');
    name = strrep(name, ',', '');
    % Remove repeated underscores.
    name = regexprep(name, '([\s._])\1+', '_');
    if ~any(ismember(vars.name, name))
        fprintf('Found %s\n', name)
        current_index = length(vars.name) + 1;
        vars.name{current_index} = name;
        vars.units{current_index} = items.unit{i};
        vars.abbreviation{current_index} = items.abbreviation{i};
    end
end

% Make the dimensions struct.
dims.site = size(all_data, 2) / length(vars.name);
dims.time = length(time.time);

% Now preallocate all the variable arrays and extract the relevant data
% from the all_data array. We have to do this after the loop identifying
% variable names as we don't know how many sites we have until that point.
for i = 1:length(vars.name)
    data.(vars.name{i}) = nan(dims.time, dims.site);
end
% Extract the relevant data from the all_data array. This makes writing to
% netCDF much more straightforward.
for name_index = 0:length(vars.name) - 1
    for site_index = 1:dims.site
        data_index = (name_index * dims.site) + site_index;
        data.(vars.name{name_index + 1})(:, site_index) = all_data(:, data_index);
    end
end
clear all_data

% Clean up after ourselves
dfs0.Close();

% Now we've read in all the data we need, create the netCDF output file.
nc = netcdf.create(fout, 'clobber');

% Add global attributes.
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'title', char(dfs0.FileInfo.FileTitle))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'type', char(dfs0.FileInfo.FileType))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'projection', char(dfs0.FileInfo.Projection.WKTString))
netcdf.putAtt(nc, netcdf.getConstant('NC_GLOBAL'), 'history', 'File created using dfs02nc (Pierre Cazenave pica@pml.ac.uk)')

% Define the dimensions.
site_dimid = netcdf.defDim(nc, 'site', dims.site); %#ok<NASGU>
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
Times_dimid = netcdf.defDim(nc, 'DateStrLen', 26);

% Create the data variables.
for k = 1:length(vars.name)
    eval([vars.name{k}, '_varid = netcdf.defVar(nc, ''', vars.name{k}, ''', ''NC_FLOAT'', [site_dimid, time_dimid]);']);
    eval(['netcdf.putAtt(nc, ', vars.name{k}, '_varid, ''units'', ''', char(vars.abbreviation{k}), ''')']);
    eval(['netcdf.putAtt(nc, ', vars.name{k}, '_varid, ''standard_units'', ''', char(vars.units{k}), ''')']);
    eval(['netcdf.putAtt(nc, ', vars.name{k}, '_varid, ''short_name'', ''', char(vars.name{k}), ''')']);
end

% Add the fixed data we know we'll have.
time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'standard_name', 'Time in days')
netcdf.putAtt(nc, time_varid, 'long_name', 'days since 0000-01-01 00:00:00')
netcdf.putAtt(nc, time_varid, 'units', 'days')
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC')
netcdf.putAtt(nc, time_varid, 'format', 'MATLAB time')

Times_varid = netcdf.defVar(nc, 'Times', 'NC_CHAR', [Times_dimid, time_dimid]);
netcdf.putAtt(nc, Times_varid, 'time_zone', 'UTC')

netcdf.endDef(nc);

% Add the data.
netcdf.putVar(nc, time_varid, 0, dims.time, time.time)
nStringOut = char();
for i = 1:dims.time
    nStringOut = [nStringOut, sprintf('%-26s', time.Times(i, :))];
end
netcdf.putVar(nc, Times_varid, nStringOut)

for k = 1:length(vars.name)
    eval(['netcdf.putVar(nc, ',vars.name{k}, '_varid, [0, 0], [', ...
        num2str(dims.site), ', ', num2str(dims.time),'], data.', ...
        vars.name{k}, ''');'])
end

netcdf.close(nc)
