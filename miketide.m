function [time, heights, info] = miketide(resultsFile)
% MIKETIDE Extract the height, metainformation and time from a MIKE dfs0
% file.
%   [TIME, HEIGHTS, INFO] = MIKETIDE(DFS0) reads in a MIKE dfs0 file with
%   surface elevation results in it and outputs TIME (days since the start
%   of the year) and the surface heights (HEIGHTS). INFO contains the
%   output types from the model as well as site names.
%
% Pierre Cazenave 13/01/2010 v1.0.

% Let's get it on...

dfs0=dfsTSO(resultsFile);
timeRaw=readTimes(dfs0);
% Grab the year from the model data (hope it doesn't start in December!).
yr=mode(timeRaw(:,1));
time=[datenum(timeRaw)-datenum([yr 00 01]),timeRaw];
clear timeRaw

% Metainformation
itemnames=get(dfs0,'itemnames');

% Preallocate arrays
info=cell(length(itemnames),2);
heights=nan(length(itemnames),size(time,1));

for getInfo=1:length(itemnames)
    info(getInfo,:)=regexp(itemnames{getInfo},': ','split');
    heights(getInfo,:)=dfs0(itemnames{getInfo});
end
clear getInfo

close(dfs0)
