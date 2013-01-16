function [lon,lat,tide]=readTPXO(filename)
% READTPXO(FILENAME) reads in a TPXO netcdf grid file and
% extracts the eight tidal constituents.
%
% Pierre Cazenave 02/09/2010 v1.0

ncid=netcdf.open(filename,'NC_NOWRITE');

for i=1:100 % just go way beyond the number of parameters...
    try
        % Never know how many variables there will be, so can't
        % preallocate. Should only be small anyway...
        [varname{i},xtype(i),dimid(i,1:2),numatts(i)]=netcdf.inqVar(ncid,i-1); %#ok<NASGU,AGROW>
    catch %#ok<CTCH>
        break
    end
end

lon=netcdf.getVar(ncid,1);
lat=netcdf.getVar(ncid,2);

tide=netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{4}));

% Let's put the Greenwich meridian in the middle
lon=lon-180;
halftide=floor(size(tide,2)/2);
% Checked this puts the right z values in the correct place in the array
% (0deg=1081)
tide=[tide(:,halftide+1:end,:),tide(:,1:halftide,:)];

% Close the open grid file
netcdf.close(ncid)
