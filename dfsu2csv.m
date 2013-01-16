function [X,Y,Z]=dfsu2csv(filename,outname)
% DFSU2CSV(FILENAME,OUTNAME) converts a single item single timestep dfsu
%   file to a three column CSV file with X,Y,Z values. Values are stored as
%   singles internally otherwise there are some precision issues (i.e. added
%   extraneous decimal points at very small magnitudes.
%
%   [X,Y,Z]=DFSU2CSV(FILENAME,OUTNAME) creates a CSV file and also outputs
%   the necessary information for plotting the dfsu file with PATCH.
%
%   DFSU2CSV makes use of SAVECSV to create the CSV file, so you must also
%   have my version of that in your PATH.
%
% Pierre Cazenave 2011/09/12 v1.0
%

% OK soldier, let's go...

if nargin~=2
    error('Incorrect number of arguments specified. Please supply the complete path to a dfsu file.')
end

% Load the MATLAB interface for dfs files.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

% Open the file.
dfsu=DfsFileFactory.DfsuFileOpen(filename);
% Node coordinates.
xn=double(dfsu.X);
yn=double(dfsu.Y);
zn=double(dfsu.Z);

% Create element table in Matlab format.
tn=mzNetFromElmtArray(dfsu.ElementTable);

if nargout==3
    % General position data for the patch command.
    X(:,1)=xn(tn(:,1));
    X(:,2)=xn(tn(:,2));
    X(:,3)=xn(tn(:,3));
    Y(:,1)=yn(tn(:,1));
    Y(:,2)=yn(tn(:,2));
    Y(:,3)=yn(tn(:,3));
end

% Element centres
[xe,ye,~]=mzCalcElmtCenterCoords(tn,xn,yn,zn);

% Data from the first item at the first timestep.
Z=single(dfsu.ReadItemTimeStep(1,0).Data)';

% Write out the CSV file using the element centres rather than the mesh
% positions. Output those for use with patch() instead.
savecsv(outname,[xe,ye,Z]);