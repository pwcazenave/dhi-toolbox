function dfsu2mesh(dfsu)
% DFSU2MESH Convert a dfsu file to a mesh file.
%   DFSU2MESH(DFSU,MESH) converts the DFSU file to the corresponding mesh
%   file based on the spatial information in MESH.
%
%   In instances where you need to work on a file in MATLAB, doing so with
%   dfsu files is easier than with mesh files. As such, this function will
%   convert the dfsu file to a mesh file.
%
%   Requires two inputs:
%       1. A dfsu file.
%       2. A mesh file whose spatial information matching the dfsu's.
%   The mesh file should be the bathymetry mesh file so that the boundaries
%   have the same codes.
%
%   Output mesh file is in the same directory as the input dfsu file.
%
%   This is particularly useful for me when I'm merging the variable
%   Manning's meshes together. Since the depth_dependent_M.m file creates
%   dfsu files only and the mergemesh() function only works with mesh
%   files, I need to be able to convert the variable Manning's dfsu file to
%   mesh files first.
%
%   This is mostly lifted from read_dfsu_2D.m and create_dfsu_2D.m supplied
%   in the Examples directory of the DHI MATLAB tools.
%
% Pierre Cazenave v1.0 25/11/2011

% Vamos!

% For testing.
% [~,base]=checkos;
% dfsu=[base,'/modelling/data/bedforms/round_8_palaeo/iridis/combo_v3_middle_bedforms_wavelength.dfsu'];
% mesh=[base,'/modelling/data/bedforms/round_8_palaeo/iridis/combo_v3_middle_bedforms_wavelength.mesh'];

% We're assuming we're using the 2011 DHI MATLAB tools.
NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;

dfsu2 = DfsFileFactory.DfsuFileOpen(dfsu);

nsteps = dfsu2.NumberOfTimeSteps;
if nsteps>1
    warning('Number of timesteps exceeds one. For a mesh file, a single timestep only is supported. Defaulting to the first timestep.')
end

% Node coordinates.
xn=double(dfsu2.X);
yn=double(dfsu2.Y);
zn=double(dfsu2.Z);

[~,Nodes,proj] = mzReadMesh(mesh);

% Create a new mesh file from the dfsu info
NodesNew = [xn',yn',zn',Nodes(:,4)];

% Write out the result
mzWriteMesh([dfsu(1:end-5),'.mesh'],te,NodesNew,proj);



