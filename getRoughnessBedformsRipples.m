function [M]=getRoughness(formulation,bathyFile,grainFile,wavelengthFile,heightFile)
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
%   above, but uses the greater of grain size or bedform data to calculate
%   the roughness (after recommendataion in Harris, C. K. and Wiberg, P. L.
%   (2001) A two-dimensional, time-dependent model of suspended sediment
%   transport and bed reworking for continental shelves, Computers
%   Geosciences, volume 27, p 675-690).
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
%           z0f=a_t*(H^2/w)
%           where
%               a_t=0.3<a_t<3, typically 1
%               H=bedform height (m)
%               w=bedform wavelength (m)
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
%       k_s=30*z0{s,f}
%     From this, the Manning's number is
%       M=25.4/k_s^(1/6)
%
%     Pierre Cazenave 2010-09-29 v1.0 pwc101 [at] soton {dot} ac <dot> uk
%       v1.1 2011-05-05 Updated the help above slightly
%       v2.0 2011-06-10 Added the bedform dimension analysis too. All dates
%       now yyyy-mm-dd in ChangeLog.

% TODO: Write M to a dfsu file using the bathy mesh positions.
%         DONE - See depth_dependent_M.m for output to dfsu.
%       Check the values in Cd - something's making complex numbers.
%         DONE - negative depth values are the culprit.

% Let's get it on...

if nargin==3;
    % Grain size only.
    both=0;
elseif nargin==5;
    % We've got both grain size and wavelength and height data, and we'll
    % compare the two.
    both=1;
else
    error('Something''s amiss - check the number of arguments to this function')
end

% Read in the bathy dfsu file (can't do with a mesh file, I don't think).
dfsuBathy=dfsManager(bathyFile);
bathyItems=get(dfsuBathy,'items');
bathyNum=find(strcmpi('Bathymetry',bathyItems));
% tnBathy=readElmtNodeConnectivity(dfsuBathy);
% [xnBathy,ynBathy,znBathy]=readNodes(dfsuBathy); % Node centre coordinates
% [xeBathy,yeBathy,zeBathy]=readElmts(dfsuBathy); % Element vertex coordinates
% xBathy(:,1:3)=xnBathy(tnBathy(:,1:3));
% yBathy(:,1:3)=ynBathy(tnBathy(:,1:3));
zBathy(:,1)=dfsuBathy(bathyNum,0);
zBathy=zBathy*-1; % Flip the depth sign.

% Read in the bedform wavelength data.
dfsuBedformWavelength=dfsManager(wavelengthFile);
bedformWavelengthItems=get(dfsuBedformWavelength,'items');
bedformWavelengthNum=find(strcmpi('Wavelength',bedformWavelengthItems));
% tnBedformWavelength=readElmtNodeConnectivity(dfsuBedformWavelength);
% [xnBedformWavelength,ynBedformWavelength,znBedformWavelength]=readNodes(dfsuBedformWavelength); % Node centre coordinates
% [xeBedformWavelength,yeBedformWavelength,zeBedformWavelength]=readElmts(dfsuBedformWavelength); % Element vertex coordinates
% xBedformWavelength(:,1:3)=xnBedformWavelength(tnBedformWavelength(:,1:3));
% yBedformWavelength(:,1:3)=ynBedformWavelength(tnBedformWavelength(:,1:3));
zBedformWavelength(:,1)=dfsuBedformWavelength(bedformWavelengthNum,0);

% Read in the bedform height data.
dfsuBedformHeight=dfsManager(heightFile);
bedformHeightItems=get(dfsuBedformHeight,'items');
bedformHeightNum=find(strcmpi('Height',bedformHeightItems));
% tnBedformHeight=readElmtNodeConnectivity(dfsuBedformHeight);
% [xnBedformHeight,ynBedformHeight,znBedformHeight]=readNodes(dfsuBedformHeight); % Node centre coordinates
% [xeBedformHeight,yeBedformHeight,zeBedformHeight]=readElmts(dfsuBedformHeight); % Element vertex coordinates
% xBedformHeight(:,1:3)=xnBedformHeight(tnBedformHeight(:,1:3));
% yBedformHeight(:,1:3)=ynBedformHeight(tnBedformHeight(:,1:3));
zBedformHeight(:,1)=dfsuBedformHeight(bedformHeightNum,0);

% Work on the grain size data, which will always be given.
% WARNING: grain size data is almost certainly in mm - make sure you
% convert it to metres for the calculations.
warning('Check input grain size is in millimetres: I''m converting to metres here.') %#ok<WNTAG>
dfsuGrain=dfsManager(grainFile);
grainItems=get(dfsuGrain,'items');
grainNum=find(strcmpi('Sediments',grainItems));
% tnGrain=readElmtNodeConnectivity(dfsuGrain);
% [xnGrain,ynGrain,znGrain]=readNodes(dfsuGrain); % Node centre coordinates
% [xeGrain,yeGrain,zeGrain]=readElmts(dfsuGrain); % Element vertex coordinates
% xGrain(:,1:3)=xnBathy(tnGrain(:,1:3));
% yGrain(:,1:3)=ynBathy(tnGrain(:,1:3));
zGrain(:,1)=dfsuGrain(grainNum,0);
zGrain=zGrain/1000; % Convert grain size to metres.

% Let's convert this grain size into a measurement of roughness
z0s=zGrain./12;
% then get Nikuradse's roughness length
k_s=30.*z0s;
% finally, an initial set of Manning's numbers (non-depth dependent). This
% is currently unused.
altM=25.4./k_s.^(1/6);

% Calculate the roughness based on the wavelength and heights supplied.
a_t=1; % See Soulsby (1997) for more information on this factor.
z0f=a_t.*((zBedformHeight.^2)./zBedformWavelength);


% Now for a more comprehensive calculation of Manning's roughness.
% First, we need to get a decent value of the drag coefficient (Cd). We
% need to choose between the Manning-Strickler law and Dawson et al.
% (1983)'s formulation.
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

if both==0;
    % Only have grain size data
    Cd=a*(z0s./zBathy).^B;
elseif both==1;
    % We have grain size and bedform data, so we can compare the two
    % calculated roughnesses (z0s and z0f), extract the larger of the two,
    % and use that instead.
    % Can't think of an elegant way to do this now, so we'll just loop over
    % every element.
    % Need to check the value of z0f is greater than zero too (negative is
    % how I represent NaN in the input dfsu file.

    % Just a bit of info: print out the percentage of times when z0f is
    % greater than z0s (i.e. bedforms win over grain size).
    totalGreater=0;
    totalCompared=0;

    Cd=nan(size(z0s));
    for i=1:size(z0f,1);
        % No need to check if z0f(i) is greater than zero because it's
        % greater than z0s(i), which can only be positive.
        if z0f(i)>z0s(i)
            Cd(i)=a*(z0f(i)/zBathy(i))^B;
            totalGreater=totalGreater+1;
            totalCompared=totalCompared+1;
        % Form drag is smaller than skin friction, but greater than zero.
        elseif z0f(i)<z0s(i) && z0f(i)>0
            Cd(i)=a*(z0s(i)./zBathy(i)).^B;
            totalCompared=totalCompared+1;
%             Cd(i)=a*(z0s(i)/zBathy(i))^B;
        % Pick z0s(i) since we've got a negative z0f(i) i.e. NaN.
        elseif z0f(i)<z0s(i) && z0f(i)<=0
            Cd(i)=a*(z0s(i)./zBathy(i)).^B;
%             Cd(i)=a*(z0s(i)/zBathy(i))^B;
        else
            error('Something''s amiss, you shouldn''t see this...')
        end
    end
    fprintf('Bedform derived roughness used %.2f%% of the time\n',(totalGreater/totalCompared)*100)
else
    error('Ooooh... something''s wrong with your input arguments.')
end

% Now we've got something approaching a sensible distribution of values, we
% can calculate out final Manning's numbers.
g=9.81; % acceleration due to gravity
M=sqrt(g./Cd)./zBathy.^(1/6);

% However, we have twice used the depth now (once in the calculation of Cd
% and once to get M)...

