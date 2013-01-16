function array2dfs1(data,filename,startTime)
% ARRAY2DFS1 exports a 2D array of tidal data to a dfs1 file for use in
%   MIKE21's modelling software.
%
%   ARRAY2DFS1(DATA,FILENAME,STARTTIME) takes DATA, which contains the time
%   vector in the first column, then the data along the boundary in the
%   remaining columns, and outputs a dfs1 starting from the time and date
%   in STARTTIME to a new dfs1 file FILENAME.
%
%   Pierre Cazenave 12/10/2010 v1.0
%

% Let's get it on...

% Set up and generate the dfs1 file.
timeStep = floor((data(2,1)-data(1,1))*86400);
t=data(:,1);
dfs1 = dfsManager(filename,1);
createTemporalDef(dfs1,'Equidistant_Calendar',startTime,size(t,1),timeStep);
setSpatialDef(dfs1,'LONG/LAT',[0,0],0,[0,1,size(data(1,2:end),2)]);
addItem(dfs1,'Predicted','Surface Elevation','m');
create(dfs1);

for i=0:size(data,1)-1
    dfs1(1,i)=double(data(i+1,2:end));
%     Alternative output function (should be the same effect as above).
%     writeItemTimestep(dfs1,1,i,double(data(i+1,2:end)));
end

% Aaaand we're done.
saveAndClose(dfs1);