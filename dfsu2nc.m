function dfsu2nc(filename, varargin)
% DFSU2NC(FILENAME) converts a 2D DFSU into a NetCDF file.
%
% Optionally supply argument pair:
%   'overwrite' - set to 'yes' to overwrite existing files
%
% Pierre Cazenave 2012/12/10 v1.0
%                 2014/04/14 v1.1 - Add option to overwrite existing files.
%

if nargin ~= 1 && nargin ~= 3
    error('Incorrect number of arguments specified. Please supply the complete path to a dfsu file.')
end

clobber = false;
for aa = 2:2:nargin
    arg = varargin{aa};
    switch aa
        case {'overwrite', 'Overwrite', 'OVERWRITE'}
            if strcmpi(arg, 'yes') || strcmpi(arg, 'y')
                clobber = true;
            end
    end
end

% Create output file name.
[fpath, fname, ~] = fileparts(filename);
fout = fullfile(fpath, [fname, '.nc']);
clear fpath fname

if exist(fout, 'file') == 2 && ~clobber
    fprintf('\n')
    warning('File %s exists. Not overwriting.', fout)
    return
end

% Load the MATLAB interface for dfs files.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

% Open the file.
dfsu = DfsFileFactory.DfsuFileOpen(filename);

% Extract the dimensions of the outputs in the DFSU file.
dims.node = dfsu.NumberOfNodes;
dims.elem = dfsu.NumberOfElements;
dims.lay = dfsu.NumberOfLayers;
dims.siglay = dfsu.NumberOfSigmaLayers;
dims.time = dfsu.NumberOfTimeSteps;

% Spatial data coordinates.
xn = double(dfsu.X);
yn = double(dfsu.Y);
zn = double(dfsu.Z);

tn = mzNetFromElmtArray(dfsu.ElementTable);
[xe,ye,ze] = mzCalcElmtCenterCoords(tn,xn,yn,zn);

% Time data.
time.timestep = double(dfsu.TimeStepInSeconds);
time.start = double([dfsu.StartDateTime.Year,dfsu.StartDateTime.Month,...
    dfsu.StartDateTime.Day,dfsu.StartDateTime.Hour,...
    dfsu.StartDateTime.Minute,dfsu.StartDateTime.Second]);
time.end = datevec(datenum(time.start) + ...
    (double(dims.time-1) * (time.timestep/(24*60*60))));
time.time = datenum(time.start):time.timestep/86400:datenum(time.end);
time.Times = datestr(time.time, 'yyyy-mm-ddTHH:MM:SS.FFF');

% Variables to be exported.
items=cell(dfsu.ItemInfo.Count,3);
for i=0:dfsu.ItemInfo.Count - 1
   item = dfsu.ItemInfo.Item(i);
   % Strip out invalid characters.
   items{i+1,1} = regexprep(regexprep(regexprep(regexprep(regexprep(char(item.Name), ':', '-'), '(', ''), ')', ''), '\.', ''), ',', '');
   items{i+1,2} = char(item.Quantity.Unit);
   items{i+1,3} = char(item.Quantity.UnitAbbreviation);
end

% Iterate through the items in the DFSU file and save to a struct ready to
% be exported to NetCDF.
data = struct;
for k = 1:dfsu.ItemInfo.Count
    name = regexprep(regexprep(regexprep(items(k, 1), ' ', '_'), '-', '_'), '___', '_');
    % This assumes we're doing 2D DFSU files. If we were doing 3D DFSU,
    % then there'd be a dims.lay or dims.siglay in there too.
    Z = nan(dims.elem, dims.time);
    for i = 0:dims.time - 1
        Z(:, i + 1) = dfsu.ReadItemTimeStep(k, i).Data.single';
    end
    data.(name{1}).data = Z;
    data.(name{1}).shortname = items(k, 2);
    data.(name{1}).units = items(k, 3);
end
clear Z

% Clean up after ourselves
dfsu.Close();



% Now we've read in all the data we need, create the NetCDF output file.
nc = netcdf.create(fout, 'clobber');

% Add global attributes.
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'title', char(dfsu.FileTitle))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'type', char(dfsu.DfsuFileType))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'projection', char(dfsu.Projection.WKTString))
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'), 'history', 'FILE CREATED using dfsu2nc (Pierre Cazenave pica@pml.ac.uk)')

% Define the dimensions.
elem_dimid = netcdf.defDim(nc, 'nele', dims.elem);
node_dimid = netcdf.defDim(nc, 'node', dims.node);
lay_dimid = netcdf.defDim(nc, 'lay', dims.lay + 1);
siglay_dimid = netcdf.defDim(nc, 'siglay', dims.siglay + 1);
time_dimid = netcdf.defDim(nc, 'time', netcdf.getConstant('NC_UNLIMITED'));
Times_dimid = netcdf.defDim(nc, 'DateStrLen', 26);
three_dimid = netcdf.defDim(nc, 'Three', 3);

% Create the data variables.
fnames = fieldnames(data);
for k = 1:size(fnames, 1)
    eval([fnames{k}, '_varid = netcdf.defVar(nc, ''', fnames{k}, ''', ''NC_FLOAT'', [elem_dimid, time_dimid]);']);
    eval(['netcdf.putAtt(nc, ', fnames{k}, '_varid, ''units'', ''', char(data.(fnames{k}).units), ''')']);
    eval(['netcdf.putAtt(nc, ', fnames{k}, '_varid, ''short_name'', ''', char(data.(fnames{k}).shortname), ''')']);
end

% Add the fixed data we know we'll have.
lonc_varid = netcdf.defVar(nc, 'lonc', 'NC_FLOAT', elem_dimid);
netcdf.putAtt(nc, lonc_varid, 'standard_name', 'longitude')
netcdf.putAtt(nc, lonc_varid, 'long_name', 'zonal longitude')
netcdf.putAtt(nc, lonc_varid, 'units', 'degrees_east')

latc_varid = netcdf.defVar(nc, 'latc', 'NC_FLOAT', elem_dimid);
netcdf.putAtt(nc, latc_varid, 'standard_name', 'latitude')
netcdf.putAtt(nc, latc_varid, 'long_name', 'zonal latitude')
netcdf.putAtt(nc, latc_varid, 'units', 'degrees_north')

lon_varid = netcdf.defVar(nc, 'lon', 'NC_FLOAT', node_dimid);
netcdf.putAtt(nc, lon_varid, 'standard_name', 'longitude')
netcdf.putAtt(nc, lon_varid, 'long_name', 'nodal longitude')
netcdf.putAtt(nc, lon_varid, 'units', 'degrees_east')

lat_varid = netcdf.defVar(nc, 'lat', 'NC_FLOAT', node_dimid);
netcdf.putAtt(nc, lat_varid, 'standard_name', 'latitude')
netcdf.putAtt(nc, lat_varid, 'long_name', 'nodal latitude')
netcdf.putAtt(nc, lat_varid, 'units', 'degrees_north')

h_varid = netcdf.defVar(nc, 'h', 'NC_FLOAT', node_dimid);
netcdf.putAtt(nc, h_varid, 'standard_name', 'nodal_sea_floor_depth')
netcdf.putAtt(nc, h_varid, 'long_name', 'Nodal Bathymetry')
netcdf.putAtt(nc, h_varid, 'positive', 'up')
netcdf.putAtt(nc, h_varid, 'units', 'metres') % this might be wrong in some weird cases

hc_varid = netcdf.defVar(nc, 'hc', 'NC_FLOAT', elem_dimid);
netcdf.putAtt(nc, hc_varid, 'standard_name', 'zonal_sea_floor_depth')
netcdf.putAtt(nc, hc_varid, 'long_name', 'Zonal Bathymetry')
netcdf.putAtt(nc, hc_varid, 'positive', 'up')
netcdf.putAtt(nc, hc_varid, 'units', 'metres') % this might be wrong in some weird cases

time_varid = netcdf.defVar(nc, 'time', 'NC_FLOAT', time_dimid);
netcdf.putAtt(nc, time_varid, 'standard_name', 'Time in days')
netcdf.putAtt(nc, time_varid, 'long_name', 'days since 0000-01-01 00:00:00')
netcdf.putAtt(nc, time_varid, 'units', 'days')
netcdf.putAtt(nc, time_varid, 'time_zone', 'UTC')
netcdf.putAtt(nc, time_varid, 'format', 'MATLAB time')

Times_varid = netcdf.defVar(nc, 'Times', 'NC_CHAR', [Times_dimid, time_dimid]);
netcdf.putAtt(nc, Times_varid, 'time_zone', 'UTC')

% We need to store the triangulation table in the NetCDF file because the
% mesh files don't store it in the same order. Thus, if we want to plot the
% data in anything other than MIKE, we need this.
tri_varid = netcdf.defVar(nc, 'tri', 'NC_INT', [three_dimid, elem_dimid]);
netcdf.putAtt(nc, tri_varid, 'standard_name', 'Triangulation table')

netcdf.endDef(nc);

% Add the data.
netcdf.putVar(nc, lon_varid, xn);
netcdf.putVar(nc, lat_varid, yn);
netcdf.putVar(nc, lonc_varid, xe);
netcdf.putVar(nc, latc_varid, ye);
netcdf.putVar(nc, h_varid, zn);
netcdf.putVar(nc, hc_varid, ze);
netcdf.putVar(nc, time_varid, 0, dims.time, time.time)
netcdf.putVar(nc, tri_varid, tn')
nStringOut = char();
for i = 1:dims.time
    nStringOut = [nStringOut, sprintf('%-26s', time.Times(i, :))];
end
netcdf.putVar(nc, Times_varid, nStringOut)

for k = 1:size(fnames, 1)
    eval(['netcdf.putVar(nc,',fnames{k}, '_varid, [0, 0], [', ...
        num2str(dims.elem), ',', num2str(dims.time),'], data.', ...
        fnames{k}, '.data);'])
end

netcdf.close(nc)