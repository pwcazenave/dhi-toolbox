function [x,y,z]=readPeltier(filename,attribute)
% READPELTIER(FILENAME,ATTRIBUTE) reads in a Peltier netcdf grid file and
% extracts ATTRIBUTE, where ATTRIBUTE is a string, either:
%   'altitude' - extracts the heights (including the ice sheet).
%   'thickness' - extracts the ice sheet thickness.
%   'mask' - extracts a binary ice/no ice map.
%   'bathy' - subtracts the thickness from the altitude to leave
%   palaeo-bathymetry.
%
% Pierre Cazenave   17/06/2010 v1.0
%                   26/10/2010 v1.1 Changed to support new version of the
%                   grids with 500 year intervals.

ncid=netcdf.open(filename,'NC_NOWRITE');

% A bit of a hack...
for i=1:100 % just go (hopefully) way beyond the number of parameters.
    try
        % Never know how many variables there will be, so can't
        % preallocate. Should only be small anyway.
        [varname{i},xtype(i),dimid(i,1:2),numatts(i)]=netcdf.inqVar(ncid,i-1); %#ok<NASGU,AGROW>
    catch %#ok<CTCH>
        break
    end
end

splitname=regexpi(filename,'/','split');
splitname=regexpi(splitname,'_','split');

% Need to differentiate between Peltier's grid files and those from
% Rosemarie. Peltier's end-member files have 7 variables in the grid;
% Rosemarie's only have 5 (they're missing bounds_lat and bounds_long).

% Amended because I now have new 500 year interval grids from 0.5 ka BP to
% 14.5 ka BP, and then 1ka BP interval grids from there to 20 ka BP. So,
% the check need only find the 00 and 21 ka BP grids, which are still
% Peltier's original ones, rather than Rosemarie's ones.
if strcmpi(splitname{end}(3),'00.0KBP')==1 || strcmpi(splitname{end}(3),'21.0KBP')==1
    % Grids from Peltier's site
    x=netcdf.getVar(ncid,2);
    y=netcdf.getVar(ncid,0);
    correction=0;
else
    % Rosemarie's grids
    % Lat and long are stored in the opposite order, other values are
    % offset by two positions.
    x=netcdf.getVar(ncid,0);
    y=netcdf.getVar(ncid,1);
    correction=2;
end

if strcmpi(attribute,'altitude')==1
    z=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{5-correction})))';
elseif strcmpi(attribute,'thickness')==1
    z=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{6-correction})))';
elseif strcmpi(attribute,'mask')==1
    z=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{7-correction})))';
elseif strcmpi(attribute,'bathy')==1
    zalt=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{5-correction})))';
    zthick=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,varname{6-correction})))';
    z=zalt-zthick;
end

% Let's put the Greenwich meridian in the middle
x=x-180;
halfz=floor(size(z,2)/2);
% Checked this puts the right z values in the correct place in the array
% (0deg=1081)
z=[z(:,halfz+1:end),z(:,1:halfz)];

% Close the open grid file
netcdf.close(ncid)
