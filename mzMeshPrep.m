function mzMeshPrep(varargin)
% Takes a mesh file and prepares the files needed for MZMESHMERGE.
%
% Simple extraction of the node and element information into two separate
% files whose name is the mesh file suffixed with _nodes.txt and
% _elements.txt. By default, output files are not overwritten if they
% already exist. Override this behaviour by adding a second argument
% OVERRIDE with a value of 'yes'. Any other value and files will not be
% overwritten, likewise no second argument will not overwrite the files.
%
% Pierre Cazenave 03/11/2011 v1.0

% We're just going to read over each line in the file, check the number of
% columns and spit those out to the relevant files.
%
% Two columns is the node header.
% Five columns is the node spatial information (xyz).
% Three columns is the element header.
% Four columns is the element connectivity.

% It's time...

if size(varargin,2) == 2
    overwrite = varargin{2};
    mesh = varargin{1};
elseif size(varargin,2) == 1
    overwrite = 'no';
    mesh = varargin{1};
else
    error('Incorrect number of input arguments.')
end

fid = fopen(mesh);
if fid == -1
    error('Could not open input file %s.',mesh)
end

% Check we're not overwriting files when we shouldn't be.
if exist([mesh(1:end-5),'_nodes.txt'],'file') == 2 || exist([mesh(1:end-5),'_elements.txt'],'file') == 2
    % Output files exist, so check whether we're overwriting or not.
    if strcmpi(overwrite,'yes') ~= 1
        % Exit because we haven't been explicitly told to overwrite the
        % output files. This way, we default to safe behaviour.
        error('Aborting due to existence of output files.')
    end
end

% Create output files
fidNodes = fopen([mesh(1:end-5),'_nodes.txt'],'w');
fidElements = fopen([mesh(1:end-5),'_elements.txt'],'w');
if fidNodes == -1
    error('Could not create output file %s.',[mesh(1:end-5),'_nodes.txt'])
end
if fidElements == -1
    error('Could not create output file %s.',[mesh(1:end-5),'_elements.txt'])
end

% Count which line we're on.
counter = 1;

tline = fgets(fid);
while ischar(tline)
    if counter == 1
        % This is the header for the nodes. Since it has text in it, we'll
        % just spit it out as is. This is ignored in MERGEMESHES anyway, so
        % it's not too important that this is preserved perfectly.
        fprintf(fidNodes,'%s',tline);
    else
        test = textscan(tline,'%s');
        if numel(test{1}) == 5
            % Node spatial info
            fprintf(fidNodes,'%s',tline);
        elseif numel(test{1}) == 3
            % Element header
            fprintf(fidElements,'%s',tline);
        elseif numel(test{1}) == 4
            % Element connectivity
            fprintf(fidElements,'%s',tline);
        else
            error('Something is amiss: we have too many or too few fields. Check the mesh file is valid.')
        end
    end
    tline = fgets(fid);
    counter = counter+1;
end

fclose(fid);
fclose(fidNodes);
fclose(fidElements);