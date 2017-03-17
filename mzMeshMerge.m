function mzMeshMerge(corenodes,coreelements,extensionnodes,extensionelements,boundary,outputfile,surrounding)
% MZMESHMERGE merges two meshes with a common edge (either adjacent meshes,
%   or a mesh within another mesh).
%
%   Usage:
%
%   MZMESHMERGE(CORENODES,COREELEMENTS,EXTENSIONNODES,EXTENSIONELEMENTS,...
%       BOUNDARY,OUTPUT,SURROUNDING)
%
%   CORENODES is the nodes part of the file forming the base of the merge.
%   COREELEMENTS is the elements of that file.
%   EXTENSIONNODES is the nodes part of the mesh file being attached to the
%   core file.
%   EXTENSIONELEMENTS is the elements part of the extension mesh file.
%   BOUNDARY is the common points between the two mesh files, either along
%   the boundary or along the edge, if the extension is being put within
%   the core mesh.
%   OUTPUT is the name of the output file.
%   SURROUNDING either 1 or 0 for a nested merge or an adjacent merge
%   respectively.
%
%   The files for the CORENODES, COREELEMENTS, EXTENSIONNODES and
%   EXTENSIONELEMENTS can be created with the mzMeshPrep function. The
%   BOUNDARY file can be created with the mzMeshBoundary function.
%
%   PURPOSE:
%   - merges two meshes which touch each other
%   FEATURES:
%   - keeps node and element IDs from the first mesh
%   - drops double nodes in mesh file 2 by searching for them in an radius
%   HOW IT WORKS:
%   1. look for double nodes using x y z data of the touching points (it
%   could be done by using the code values included in the mesh file also
%   but then mesh 2 must surround mesh 1. This solution works for all
%   touching meshes)
%   2. makes a translation matrix of old and new node IDs including x and y
%   coordinates
%   3. adds highest node ID of file 1 to all mesh nodes of file 2, but
%   drops double nodes
%   4. gives new IDs to the nodes in the element section of file 2
%   5. writes merged mesh file
%   NOTE:
%   The merged mesh file is set to have LONG/LAT coordinate system.     !!!
%   Please change it manually if necessary.                             !!!
%
%   WHAT YOU NEED:
%   Use mzMeshPrep to prepare 5 files:
%   1) Node section of mesh file 1 including header (1 line)
%   2) Element section of mesh file 1 including header (1 line)
%   3) Node section of mesh file 2 including header (1 line)
%   4) Element section of mesh file 2 including header (1 line)
%   5) Boundary file with x y z data of the points where both meshes touch
%   each other
%   Both meshes have to based on the same projection.
%
%   TEST
%   this script has been tested for merging 2 mesh files
%             file 1       file 2
%   Nodes     164100       38389
%   Elements  327804       76029
%   boundary file: 394 xyz values
%   run time on a T60: 13 min
%
% Original code by sle.
% Turned into a function by Pierre Cazenave pwc101@soton.ac.uk
%   v1.0 02/11/2011:
%         - Added two functions to prepare data for this function
%         (mzMeshPrep and mzMeshBoundary).
%         - Added check to make sure boundaries along each input mesh match
%         in number of elements.
%         - Some cosmetic code cleanup
%   v1.1 16/01/2012:
%         - BUG: As it stands, meshes will retain their original arc
%         attribute values. I have left code in which _should_ replace the
%         attribute with user-defined one, but it doesn't work all the
%         time, so I've commented it out for the time being. Maybe one day
%         I'll find the time to figure out what's wrong.
%
%

% Right, let's do this...

if nargin~=7
    error('Incorrect number of arguments specified. See ''help MZMESHMERGE'' for more information.')
end

% INPUT SECTION

% node files to read
InfileCoreNodes = corenodes; % containing nodes of mesh 1
InfileExtensionNodes = extensionnodes; % containing nodes of mesh 2
% structure of node files
FileColumnsNodes = 5; % No. of columns of node files
FileHeaderRowsNodes = 1; % No. of header rows in node file

% boundary file to read, containing boundary nodes of mesh 1
InfileBorder = boundary;
% structure of boundary file
FileColumnsBorder = 3; % No. of data columns of border file
FileHeaderRowsBorder = 0; % No. of header rows in border file

% element files to read
InfileCoreElements = coreelements;
InfileExtensionElements = extensionelements;
% structure of element files
FileColumnsElements = 4;
FileHeaderRowsElements = 1;

% processing parameters
% minimum angle for mesh generation for search radius calculation
% (default: 26°)
angle_min = 26;
% search radius for double nodes (if 0, it is calculated automatically)
radius_constant = 0;

% files to write
mergedfile = outputfile;

% =========================================================================
% NODES SECTION

% read nodes of tile 1
fid=fopen(InfileCoreNodes);
clear filename M
if fid == -1
    msgbox(sprintf('Cannot find file ''%s''.',InfileCoreNodes),...
        'Error!','error')
end
% number of rows in the head
s = '%s';
for i = 1:FileColumnsNodes-1
    s = [s, '%s'];
end
% reads every column one by one
m = textscan(fid, s, 'headerLines', FileHeaderRowsNodes);
clear s
CORE = zeros(length(m{1}),FileColumnsNodes); % memory for extension pts
% conversion from cell array
for i = 1:FileColumnsNodes
    for j = 1:length(m{1})
        % if there is a problem of ./,
        s = sscanf(char(strrep(m{i} (j,1),',','.')),'%f');
        if isnumeric(s)
            % structure of CORE: ID x y z 0/1
            CORE(j,i) = s;
        else
            CORE(j,i) = NaN;
        end
    end
end
Core_No = length(m{1}); % store no. of core nodes
fclose(fid);
clear m

% read nodes of file 2
fid = fopen(InfileExtensionNodes);
clear filename M
if fid == -1
    msgbox(sprintf('Cannot find file ''%s''.',InfileExtensionNodes),...
        'Error!','error')
end
% number of rows in the head
s = '%s';
for i = 1:FileColumnsNodes-1
    s = [s, '%s'];
end
% reads every column one by one
m = textscan(fid, s, 'headerLines', FileHeaderRowsNodes);
clear s
EXTENSION = zeros(length(m{1}),FileColumnsNodes);
% conversion from cell array
for i = 1:FileColumnsNodes
    for j = 1:length(m{1})
        % if there is a problem of ./,
        s = sscanf(char(strrep(m{i} (j,1),',','.')),'%f');
        if isnumeric(s)
            % structure of EXTENSION: ID x y z 0/1
            EXTENSION(j,i) = s;
        else
            EXTENSION(j,i) = NaN;
        end
    end
end
Extension_No = length(m{1}); % store no. of extension nodes
fclose(fid);
clear FileHeaderRowsNodes FileColumnsNodes InfileCoreNodes

% open file
fid = fopen(InfileBorder);
clear filename M
if fid == -1
    msgbox(sprintf('Cannot find file ''%s''.',InfileBorder),...
        'Error!','error')
end
% number of rows in the head
s = '%s';
for i=1:FileColumnsBorder-1
    s = [s, '%s'];
end
% reads every column one by one
m=textscan(fid, s, 'headerLines', FileHeaderRowsBorder);
clear s
BORDER=zeros(length(m{1}),FileColumnsBorder);
% conversion from cell array
for i = 1:FileColumnsBorder
    for j = 1:length(m{1})
        % if there is a problem of ./,
        s = sscanf(char(strrep(m{i} (j,1),',','.')),'%f');
        if isnumeric(s)
            % structure of BORDER:
            BORDER(j,i) = s;
        else
            BORDER(j,i) = NaN;
        end
    end
end
Border_No = length(m{1});
fclose(fid);
clear FileHeaderRowsBorder FileColumnsBorder InfileExtensionNodes

% search radius for double nodes
if radius_constant == 0 % calculates radius for each node individually
    radius = zeros(Border_No,1)+1000000000; % set radius high
    h = waitbar(0,'Calculating search radius');
    % look for minumum distance in x and y between the neighbour nodes
    xradius = min(abs(BORDER(1,1)-BORDER(Border_No,1)),...
        abs(BORDER(2,1)-BORDER(1,1)));
    yradius = min(abs(BORDER(1,2)-BORDER(Border_No,2)),...
        abs(BORDER(2,2)-BORDER(1,2)));
    % take this distance and calculate maximum distance perpendicular to
    % border using the minimum allowed angle for mesh generation
    radius(1) = cosd(angle_min)*max(xradius, yradius)/2;
    waitbar(1/Border_No);
    for i = 2:Border_No-1 % here the same for all nodes
        xradius = min(abs(BORDER(i,1)-BORDER(i-1,1)),...
            abs(BORDER(i+1,1)-BORDER(i,1)));
        yradius = min(abs(BORDER(i,2)-BORDER(i-1,2)),...
            abs(BORDER(i+1,2)-BORDER(i,2)));
        radius(i) = cos(angle_min)*max(xradius, yradius)/2;
        waitbar(i/Border_No);
    end
    % and here the same for last node
    xradius = min(abs(BORDER(Border_No,1)-BORDER(Border_No-1,1)),...
        abs(BORDER(1,1)-BORDER(Border_No,1)));
    yradius = min(abs(BORDER(Border_No,2)-BORDER(Border_No-1,2)),...
        abs(BORDER(1,2)-BORDER(Border_No,2)));
    radius(Border_No) = cosd(angle_min)*max(xradius, yradius)/2;
    waitbar(1);
    close(h);
else % or forget it and give a constant search radius
    radius = zeros(Border_No,1)+radius_constant;
end

% LOOK FOR DOUBLE POINTS
% compares EXTENSION with BORDER and marks double points in
% Marker vector with 1

% DOUBLE structure: Core_ID Extension_ID Core_x Core_y Extension_x Extension_y
DOUBLE = zeros(Border_No,6); % for collecting IDs of double points
Marker = zeros(Extension_No,1); % for marking double points

% Counts to check the number of values checked at each border (CORE and
% EXTENSION). If the two don't match, the output mesh will probably be
% subtly invalid or completely invalid.
CORE_COUNT = 0;
EXTENSION_COUNT = 0;

h = waitbar(0,'1. Marking double nodes');
for i = 1:Border_No % loop over border nodes
    for j = 1:Core_No % compare border nodes with core nodes
        if ((BORDER(i,1)<CORE(j,2)+radius(i)) && ...
                (BORDER(i,1)>CORE(j,2)-radius(i)) && ...
                (BORDER(i,2)<CORE(j,3)+radius(i)) && ...
                (BORDER(i,2)>CORE(j,3)-radius(i)))
            CORE_COUNT = CORE_COUNT+1;
            DOUBLE(i,1) = CORE(j,1);
            DOUBLE(i,3) = CORE(j,2);
            DOUBLE(i,4) = CORE(j,3);
            break; % save time
        end
    end
    for j = 1:Extension_No % compare border nodes with extension nodes
        if ((BORDER(i,1)<EXTENSION(j,2)+radius(i)) && ...
                (BORDER(i,1)>EXTENSION(j,2)-radius(i)) && ...
                (BORDER(i,2)<EXTENSION(j,3)+radius(i)) && ...
                (BORDER(i,2)>EXTENSION(j,3)-radius(i)))
            EXTENSION_COUNT = EXTENSION_COUNT+1;
            DOUBLE(i,2) = EXTENSION(j,1);
            DOUBLE(i,5) = EXTENSION(j,2);
            DOUBLE(i,6) = EXTENSION(j,3);
            Marker(j) = 1;
            break; % save time
        end
    end
    waitbar(i/Border_No);
end
close(h);
clear BORDER

% We can't create a valid mesh if the CORE_COUNT and EXTENSION_COUNT values
% don't match, so we'll abort here. This is often the case if the mesh
% generation has had to fill a corner with a number of small elements,
% which creates new nodes along the common boundary. To find these points,
% you can check the content of DOUBLE and find the difference between
% DOUBLE(:,3) and DOUBLE(:,5) or DOUBLE(:,4) and DOUBLE(:,6) and look for
% non-zero values. Those points should be identical, but because DOUBLE is
% initialised with zeros(), you can use this to identify the points which
% have only occurred in one of the mesh borders. An alternative might be to
% initialised DOUBLE with nan(), at which point a simple
% sum(isnan(DOUBLE))~=0 && error('Quitting: unequal borders') might work.
if CORE_COUNT ~= EXTENSION_COUNT
    error('Unequal number of points on the two mesh boundaries which will produce an incorrect mesh. Check your .mdf file for elements between nodes along the boundary being merged.')
end

% TRANSLATION MATRIX
% store old IDs with corresponding new IDs in a matrix
Merged_No = Core_No+Extension_No-Border_No; % length of TRANSLATION
% TRANSLATION structure: old_ID new_ID x y z boundary_ID core=1/extension=-1
TRANSLATION = zeros(Merged_No, 7);

h = waitbar(0,'2. Write ID translation matrix');
for i = 1:Core_No
    TRANSLATION(i,1) = CORE(i,1);
    TRANSLATION(i,2) = CORE(i,1);
    TRANSLATION(i,3) = CORE(i,2);
    TRANSLATION(i,4) = CORE(i,3);
    TRANSLATION(i,5) = CORE(i,4);
    if surrounding == 1
        TRANSLATION(i,6) = 0;
    else
        % This is the original code, which does not try to replace the
        % boundary value with zero, instead leaving it as whatever it was
        % originally.
        TRANSLATION(i,6)=CORE(i,5);
        % This code tries to replace a specified arc attribute (defined as
        % boundary_swap) with zeros, effectively deleting that boundary,
        % which is what you want since the boundary will otherwise be in
        % the middle of the mesh. However, this doesn't seem to work for
        % some reason.
        % I have since written mzMeshFix(), which replaces any specified
        % arc attribute with zero to fix this problem.
        % if CORE(i,5)==boundary_swap
        %     TRANSLATION(i,6)=0;
        % else
        %     TRANSLATION(i,6)=CORE(i,5);
        % end
    end
    TRANSLATION(i,7) = 1;
    waitbar(i/(Core_No+Extension_No));
end
clear CORE
j = 0; % counter for new IDs
for i = 1:Extension_No
    if Marker(i) == 0
        j = j+1;
        TRANSLATION(Core_No+j,1) = EXTENSION(i,1);
        TRANSLATION(Core_No+j,2) = Core_No+j;
        TRANSLATION(Core_No+j,3) = EXTENSION(i,2);
        TRANSLATION(Core_No+j,4) = EXTENSION(i,3);
        TRANSLATION(Core_No+j,5) = EXTENSION(i,4);
        TRANSLATION(Core_No+j,6) = EXTENSION(i,5);
        TRANSLATION(Core_No+j,7) = -1;
    end
    waitbar((i+Core_No)/(Core_No+Extension_No));
end
close(h);
clear EXTENSION Marker


% =========================================================================
% ELEMENT SECTION

% read core element file
fid = fopen(InfileCoreElements);
clear filename M
if fid == -1
    msgbox(sprintf('Cannot find file ''%s''.',InfileCoreElements),...
        'Error!','error')
end
% number of rows in the head
s = '%s';
for i = 1:FileColumnsElements-1
    s = [s, '%s'];
end
% reads every column one by one
m = textscan(fid, s, 'headerLines', FileHeaderRowsElements);
clear s
CORE = zeros(length(m{1}),FileColumnsElements);
% conversion from cell array
for i = 1:FileColumnsElements
    for j = 1:length(m{1})
        % if there is a problem of ./,
        s = sscanf(char(strrep(m{i} (j,1),',','.')),'%f');
        if isnumeric(s)
            % structure of CORE: ID x y z 0/1
            CORE(j,i) = s;
        else
            CORE(j,i) = NaN;
        end
    end
end
Core_Element_No = length(m{1}); % store no. of core elements
fclose(fid);
clear m

% read extension element file
fid=fopen(InfileExtensionElements);
clear filename M
if fid == -1
    msgbox(sprintf('Cannot find file ''%s''.',InfileExtensionElements),'Error!','error')
end
% number of rows in the head
s = '%s';
for i = 1:FileColumnsElements-1
    s = [s, '%s'];
end
% reads every column one by one
m = textscan(fid, s, 'headerLines', FileHeaderRowsElements);
clear s
EXTENSION = zeros(length(m{1}),FileColumnsElements);
% conversion from cell array
for i = 1:FileColumnsElements
    for j = 1:length(m{1})
        % if there is a problem of ./,
        s = sscanf(char(strrep(m{i} (j,1),',','.')),'%f');
        if isnumeric(s)
            % structure of CORE: ID x y z 0/1
            EXTENSION(j,i) = s;
        else
            EXTENSION(j,i) = NaN;
        end
    end
end
Extension_Element_No = length(m{1}); % store no. of core points
fclose(fid);
clear m

% Updates the node IDs
% marks translated double nodes so that they will not be translated again
ELEMENT_MARKER = zeros(Extension_Element_No,3);
h = waitbar(0,'3. Updating node IDs for elements');
for i = 1:Extension_Element_No
    for j = 2:FileColumnsElements
        for k = 1:Border_No % translation loop for double nodes
            if EXTENSION(i,j) == DOUBLE(k,2)
                EXTENSION(i,j) = DOUBLE(k,1);
                ELEMENT_MARKER(i,j-1) = 2;
                break; % save time
            end
        end
    end
    waitbar(i/(Extension_Element_No*2));
end
for i = 1:Extension_Element_No
    for j = 2:FileColumnsElements
        for k = Core_No+1:Merged_No % translation loop for unique nodes
            if ((EXTENSION(i,j) == TRANSLATION(k,1)) && ...
                (ELEMENT_MARKER(i,j-1) == 0))
                EXTENSION(i,j) = TRANSLATION(k,2);
                ELEMENT_MARKER(i,j-1) = 1; % only for debugging
                break; % save time
            end
        end
    end
    waitbar((i+Extension_Element_No)/(Extension_Element_No*2));
end
close(h);

% WRITE MERGED POINTS INTO NEW FILE
h = waitbar(0,'4. Write mesh file');
fid = fopen(mergedfile,'w');
fprintf(fid,'%d %s\n',Merged_No,'LONG/LAT'); % header
for i = 1:Merged_No
    fprintf(fid,'%d %6.3f %7.3f %4.3f %d\n',TRANSLATION(i,2),...
        TRANSLATION(i,3),TRANSLATION(i,4),TRANSLATION(i,5),...
        TRANSLATION(i,6));
    waitbar(i/...
        (Core_Element_No+Extension_Element_No+Merged_No+Border_No));
end
fprintf(fid,'%d 3 21\n',Core_Element_No+Extension_Element_No);
for i = 1:Core_Element_No
    fprintf(fid,'%d %d %d %d\n',CORE(i,1),CORE(i,2),CORE(i,3),CORE(i,4));
    waitbar((i+Merged_No+Border_No)/...
        (Core_Element_No+Extension_Element_No+Merged_No+Border_No));
end
for i = 1:Extension_Element_No
    fprintf(fid,'%d %d %d %d\n',EXTENSION(i,1)+Core_Element_No,...
        EXTENSION(i,2),EXTENSION(i,3),EXTENSION(i,4));
    waitbar((i+Core_Element_No+Merged_No+Border_No)/...
        (Core_Element_No+Extension_Element_No+Merged_No+Border_No));
end
fclose(fid);
close(h);
