function varargout=getRoughness(formulation,bathyFile,grainFile,wavelengthFile,heightFile)
% GETROUGHNESS Calculates Manning's roughness across a domain.
%   [M]=GETROUGHNESS(FORMULATION,BATHY,GRAINSIZE) takes water depth
%   (BATHY) and combines it with grain size data (GRAINSIZE) to calculate
%   the Manning's number across the domain. FORMULATION determines values
% 	of the constants for the drag coefficients.
%
%   'ms' - Manning-Strickler
%   'dawson' - Dawson et al. (1984)
%   'soulsby' - Soulsby (1997)
%
%   'Best' values are those of Soulsby (1997).
%
%   [M]=GETROUGHNESS(FORMULATION,BATHY,GRAINSIZE,WAVELENGTH,HEIGHT) is as
%   above, but uses van Rijn's relationship between grain size and the
%   bedform dimensions.
%
%   NOTE: All two or three meshes must be spatially identical: I perform no
%   interpolation here. The number of elements must match perfectly; their
%   values can, obviously, differ.
%
%   The Manning's number is calculated as follows:
%
%   1. Roughness length (z0) can be calculated from either the grain size
%   (z0s), or from the ripple wavelength and height (z0f):
%       a. Grain size:
%           z0s=d50/12 (d50 in mm)
%       b. Bedforms:
%           ks=(1.1*H)*(1-exp(-25*H/w))+(3*d90)
%           where
%               H=bedform height (m)
%               w=bedform wavelength (m)
%               d90=90th percentile grain size. Since we most often have
%               d50, we've assumed a Gaussian grain size distribution,
%               thus, d90=d50*1.8.
%     2. Convert the larger of the two (GRAINSIZE and WAVELENGTH/HEIGHT) to
%     a measurement of the drag coefficient in combination with the water
%     depth (BATHY).
%       Cd=a*(z0{s,f}/d)^B
%       where
%           d=water depth (m)
%           a=0.0474, 0.0190 or 0.0415
%           B=1/3, 1.208 or 2/7
%           First set are from the Manning-Strickler law, the second from
%           Dawson et al. (1983), the third from Soulsby (1997).
%     3. Given the drag coefficient above (which varies by depth), we can
%     convert it to the Manning's roughness with the following:
%       M=sqrt(g/Cd)/d^(1/6)
%       where
%           g=acceleration due to gravity 9.81ms^{-2}
%     (From Equation 2.83 in the DHI 2009 scientific documentation entitled
%     "MIKE 21 and MIKE 3 Flow Model FM").
%
%     The values of M across the entire domain are stored in M.
%
%     An alternative calculation of M can be made, but without the depth
%     dependence, which is a nice addition, I think.
%
%     Given a value of z0{s,f}, the Nikuradse roughness length is simply
%       ks=30*z0{s,f}
%     From this, the Manning's number is
%       M=25.4/ks^(1/6)
%
%     This is not, at this time, implemented as an export, however (see
%     ks_grain and altM for those values).
%
%     Pierre Cazenave 2010-09-29 v1.0 pwc101 [at] soton {dot} ac <dot> uk
%       v1.1 2011-05-05 Updated the help above slightly
%       v2.0 2011-06-10 Added the bedform dimension analysis too. All dates
%       now yyyy-mm-dd in ChangeLog.
%       v3.0 2011-07-05 Changed the calculation of dune derived roughness
%       to that of van Rijn (1984c) which combines the grain and bedform
%       derived roughnesses into a single value of z0. Also updated the
%       code to read dfsu files using the new MATLAB interface in MIKE
%       2011.
%       v3.1 2011-07-06 Added code to check wavelength and height
%       measurements don't fall outside of Flemming (1988) power law. If
%       they do, then predict the wavelength based on the height. This
%       allows widely spaced, semi-independent dunes to have more
%       appropriate wavelengths based on their sensible height. This is
%       controlled by the checkWavelength value (one or zero/on or off).

% Let's get it on...

if nargin==3;
    % Grain size only.
    both=0;
elseif nargin==5;
    % We've got both grain size and wavelength and height data, and we'll
    % use both.
    both=1;
else
    error('Something''s amiss - check the number of arguments to this function')
end

% Acceleration due to gravity ms^{-2}.
g=9.81;

% Set up MIKE MATLAB interface.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

% Read in the bathy dfsu file.
dfsuBathy=DfsFileFactory.DfsuFileOpen(bathyFile);
bathyItems={};
for i=0:dfsuBathy.ItemInfo.Count-1
   item=dfsuBathy.ItemInfo.Item(i);
   bathyItems{i+1,1}=char(item.Name);
   bathyItems{i+1,2}=char(item.Quantity.Unit);
   bathyItems{i+1,3}=char(item.Quantity.UnitAbbreviation);
end
bathyNum=find(strcmpi('Bathymetry',bathyItems));
zBathy=double(dfsuBathy.ReadItemTimeStep(bathyNum,0).Data);
zBathy=abs(zBathy); % Positive depths only, please.

if both
    % Read in the bedform wavelength data.
    dfsuBedformWavelength=DfsFileFactory.DfsuFileOpen(wavelengthFile);
    bedformWavelengthItems={};
    for i=0:dfsuBedformWavelength.ItemInfo.Count-1
       item=dfsuBedformWavelength.ItemInfo.Item(i);
       bedformWavelengthItems{i+1,1}=char(item.Name);
       bedformWavelengthItems{i+1,2}=char(item.Quantity.Unit);
       bedformWavelengthItems{i+1,3}=char(item.Quantity.UnitAbbreviation);
    end
    bedformWavelengthNum=find(strcmpi('Wavelength',bedformWavelengthItems));
    zBedformWavelength=double(dfsuBedformWavelength.ReadItemTimeStep(bedformWavelengthNum,0).Data);

    % Read in the bedform height data.
    dfsuBedformHeight=DfsFileFactory.DfsuFileOpen(heightFile);
    bedformHeightItems={};
    for i=0:dfsuBedformHeight.ItemInfo.Count-1
       item=dfsuBedformHeight.ItemInfo.Item(i);
       bedformHeightItems{i+1,1}=char(item.Name);
       bedformHeightItems{i+1,2}=char(item.Quantity.Unit);
       bedformHeightItems{i+1,3}=char(item.Quantity.UnitAbbreviation);
    end
    bedformHeightNum=find(strcmpi('Height',bedformHeightItems));
    zBedformHeight=double(dfsuBedformHeight.ReadItemTimeStep(bedformHeightNum,0).Data);

    % Close the opened dfsu files.
    dfsuBedformHeight.Close();
    dfsuBedformWavelength.Close();
end

% Work on the grain size data, which will always be given.
% WARNING: grain size data is almost certainly in mm - make sure you
% convert it to metres for the calculations.
warning('Check input grain size is in millimetres: I''m converting to metres here.') %#ok<WNTAG>
dfsuGrain=DfsFileFactory.DfsuFileOpen(grainFile);
grainItems={};
for i=0:dfsuGrain.ItemInfo.Count-1
   item=dfsuGrain.ItemInfo.Item(i);
   grainItems{i+1,1}=char(item.Name);
   grainItems{i+1,2}=char(item.Quantity.Unit);
   grainItems{i+1,3}=char(item.Quantity.UnitAbbreviation);
end
grainNum=find(strcmpi('Sediments',grainItems));
zGrain=double(dfsuGrain.ReadItemTimeStep(grainNum,0).Data);
zGrain=zGrain/1000; % Convert grain size to metres.

% Close the opened dfsu files.
dfsuBathy.Close();
dfsuGrain.Close();

% Let's convert this grain size into a measurement of roughness
z0s=zGrain./12; % in m
% Then get Nikuradse's roughness length
% ks_grain=30.*z0s; % in m
ks_grain=2.5.*(zGrain); % in m

% Finally, an initial set of Manning's numbers (non-depth dependent). This
% is currently unused.
altM=25.4./ks_grain.^(1/6); % in m^{1/3}s^{-1}

% Calculate the Nikuradse roughness based on the wavelength and heights.
% Use the drag coefficient with the grain roughness length (z0s) to get
% depth dependent drag. If we have both grain and bedforms, use bedform
% roughness (ks_beform) to get Manning's; otherwise, use the drag
% coefficient and depth to get Manning's.

% First, we need to get a decent value of the drag coefficient (Cd). We
% need to choose between the Manning-Strickler law, Dawson et al.
% (1983) and Soulsby's formulations.
if strcmpi(formulation,'ms')==1
    % Manning-Strickler wins!
    a=0.0474;
    B=1/3;
elseif strcmpi(formulation,'dawson')==1
    % Dawson et al. (1984) win!
    a=0.019;
    B=0.208;
elseif strcmpi(formulation,'soulsby')==1
    % Soulsby (1997) wins!
    a=0.0415;
    B=0.285714;
else
    error('No drag coefficient formulation chosen. Please specify either ''ms'', ''dawson'' or ''soulsby'' in the variable FORMULATION')
end

% We'll calculate Cd irrespective of whether we have bedform data as well
% as grain size data based on the grain size data only. The choice of which
% to use comes later and is based on the presence of bedform data. This is
% essentially the roughness from grains only.
Cd=a*(z0s./zBathy).^B; % dimensionless

% Do the final calculations.
if both
    ks_bedform=nan(size(zBedformHeight,2),1); % in m
    M=nan(size(zBedformHeight,2),1);
    for i=1:size(zBedformHeight,2)
        % We've set the NaN value in the dfsu files as -10, so we need to skip
        % anything with a value of < 0. Need both wavelength and height to be
        % greater than zero because the interpolation of the wavelength can
        % cause a "fringe" to develop where values are greater than zero, but
        % still nonsense. This is less apparent on the height data beceause the
        % values are lower, so the fringe is smaller.
        if zBedformHeight(i)>0 && zBedformWavelength(i)>0
            % We have bedform measurements, so use van Rijn's (1984c) approach.
            % We need d90, so convert d50 to d90 assuming log-normal grain size
            % distribution (in metres).
            d90=zGrain(i)*3.46;
            ks_bedform(i)=(1.1*zBedformHeight(i))*(1-exp(-25*(zBedformHeight(i)/zBedformWavelength(i))))+(3*d90); % in metres
            M(i)=25.4./ks_bedform(i).^(1/6); % in m^{1/3}s^{-1}
        else
            % We have only grain size data, so use the drag coefficient from
            % grain size roughness length and depth.
            M(i)=sqrt(g/Cd(i))/zBathy(i)^(1/6); % in m^{1/3}s^{-1}
        end
    end
else
    % We have only grain size data.
    M=sqrt(g./Cd)./zBathy.^(1/6); % in m^{1/3}s^{-1}
end

outputs = {M, Cd', altM', ks_grain'};
varargout = outputs(1:nargout);
