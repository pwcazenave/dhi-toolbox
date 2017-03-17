% Convert the MIKE dfs0 files into netCDFs.

base = 'C:\Users\IEUser\Nextcloud\Papers\roughness\results';
files = fullfile(base, 'round_8_palaeo', 'rms_calibration', '*', '*.dfs0');

cd(base)
filelist = glob(files);
for file = 1:length(filelist)
    dfs02nc(filelist{file}, 'overwrite', 'yes');
    % Spatial data coordinates. This is a pain as the dfs0 doesn't contain
    % location data! Use the corresponding .txt file from the model
    % configuration to add the spatial data.
    [fpath, fname, ~] = fileparts(filelist{file});
    fout = fullfile(fpath, [fname, '.nc']);
    [~, runname] = fileparts(fpath);
    
    model_filename = fullfile('..', 'model', 'round_8_palaeo', ...
        'rms_calibration', runname, [fname, '.txt']);
    if exist(model_filename, 'file') ~= 2
        continue
    end
    fid = fopen(model_filename, 'r');
    assert(fid > 0, 'Error opening position data file.')
    positions = textscan(fid, ...
        '%s%f%f%f%[^\n\r]', ...
        'Delimiter', ',', ...
        'TextType', 'string', ...
        'ReturnOnError', false);
    fclose(fid);
    
    % Re-open the netCDF and add the coordinate data.
    nc = netcdf.open(fout, 'NC_WRITE');
    site_dimid = netcdf.inqDimID(nc, 'site');
    
    netcdf.reDef(nc)
    
    lon_varid = netcdf.defVar(nc, 'lon', 'NC_FLOAT', site_dimid);
    netcdf.putAtt(nc, lon_varid, 'standard_name', 'longitude')
    netcdf.putAtt(nc, lon_varid, 'long_name', 'longitude')
    netcdf.putAtt(nc, lon_varid, 'units', 'degrees_east')
    
    lat_varid = netcdf.defVar(nc, 'lat', 'NC_FLOAT', site_dimid);
    netcdf.putAtt(nc, lat_varid, 'standard_name', 'latitude')
    netcdf.putAtt(nc, lat_varid, 'long_name', 'latitude')
    netcdf.putAtt(nc, lat_varid, 'units', 'degrees_north')
    
    h_varid = netcdf.defVar(nc, 'h', 'NC_FLOAT', site_dimid);
    netcdf.putAtt(nc, h_varid, 'standard_name', 'sea_floor_depth')
    netcdf.putAtt(nc, h_varid, 'long_name', 'Bathymetry')
    netcdf.putAtt(nc, h_varid, 'positive', 'up')
    netcdf.putAtt(nc, h_varid, 'units', 'metres') % this might be wrong in some weird cases
    
    netcdf.endDef(nc);

    netcdf.putVar(nc, lon_varid, positions{2});
    netcdf.putVar(nc, lat_varid, positions{3});
    netcdf.putVar(nc, h_varid, -positions{4}); % depth is negative down
    netcdf.close(nc)
end