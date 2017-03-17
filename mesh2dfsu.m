function mesh2dfsu(mesh)
% MESH2DFSU Read in a mesh file and output data to a single timestep dfsu
%   file.
%   MESH2DFSU(MESH) will convert MESH to a dfsu file with a single item
%   (WaterDepth) in meters. From a model run perspective, this shouldn't
%   make any difference so long as you point the right file to the right
%   bit of the model (i.e. roughness file to the roughness input section
%   etc.).
%
%   Output dfsu file uses the same basename as the input .mesh file.
%
%   Pierre Cazenave v1.0 25/11/2011

% Vamos!

warning('This doesn''t work as a function for some reason... needs to be run manually.')

% Set up the MIKE tools for creating dfsu files
NET.addAssembly('DHI.Generic.MikeZero.DFS');
NET.addAssembly('DHI.Generic.MikeZero.EUM');
import DHI.Generic.MikeZero.DFS.*;
import DHI.Generic.MikeZero.DFS.dfsu.*;
import DHI.Generic.MikeZero.*

% Output dfsu filename
filename=[mesh(1:end-5),'.dfsu'];

% Read in the mesh
[Elmts,Nodes,proj] = mzReadMesh(mesh);
X = Nodes(:,1);
Y = Nodes(:,2);
Z = Nodes(:,3);
code = Nodes(:,4);

% Create a new empty dfsu 3D file object
factory = DfsFactory();
builder = DfsuBuilder.Create(DfsuFileType.Dfsu2D);

% Create a temporal definition matching input file
startDate = [2000,1,1,00,00,00]; % Is this important?
start = System.DateTime(startDate(1),...
    startDate(2),...
    startDate(3),...
    startDate(4),...
    startDate(5),...
    startDate(6));
builder.SetTimeInfo(start, 3600);
% Populate the builder handle with the mesh spatial information
builder.SetNodes(NET.convertArray(single(X)),...
    NET.convertArray(single(Y)),...
    NET.convertArray(single(Z)),...
    NET.convertArray(int32(code)));
builder.SetElements(mzNetToElmtArray(Elmts));
builder.SetProjection(factory.CreateProjection(proj))

% Add an item
builder.AddDynamicItem('Depth',eumQuantity(eumItem.eumIWaterDepth,eumUnit.eumUmeter));

% Create the file - make it ready for data
dfs = builder.CreateFile(filename);

[~,~,Ze] = mzCalcElmtCenterCoords(Elmts,X,Y,Z);

% Add the Z values from the mesh file
dfs.WriteItemTimeStepNext(0, NET.convertArray(single(Ze)));

dfs.Close();

