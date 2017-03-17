function mzMeshBoundary(varargin)
% MZMESHBOUNDARY finds points along mesh file boundaries for use with
% mzMeshMerge. A single or multiple boundaries can be extracted. Output is
% saved as x, y, z in a text file.
%
% MZMESHBOUNDARY(MESH,BOUNDARY) will find all boundary points in file MESH
% which occur along the specified boundary ID BOUNDARY (as specified by the
% arc attribute defined in the mesh generator). Output will be written to a
% file whose name is based on the mesh file suffixed with
% BOUNDARY_boundary.txt. If an existing file for this mesh and boundary
% exists, MZMESHBOUNDARY will exit. BOUNDARY can be an array of multiple
% values to match many arc attributes.
%
% MZMESHBOUNDARY(MESH,BOUNDARY,OVERWRITE) will find the boundary points in
% MESH and, if OVERWRITE is 'yes', will overwrite any existing output file
% (MESH name suffixed with BOUNDARY_boundary.txt). Any other value of
% OVERWRITE will not overwrite the existing output file.
%
% Pierre Cazenave 03/11/2011
% v1.0
%   - Initial version.
% v1.1 17/01/2012
%   - BOUNDARY can now contain multiple values which will all be searched for
%   and any matching nodes output accordingly.

% Let's do this

% Sort out the inputs
if ~isnumeric(varargin{2})
    warning('MIKE:Boundary:InvalidName','Specified arc attribute must be numeric (as specified when creating the mesh file).')
    boundary = varargin{2};
else
    boundary = varargin{2};
end
mesh = varargin{1};

if size(varargin,2) == 3
    overwrite=varargin(3);
elseif size(varargin,2) == 2
    overwrite = 'no';
else
    error('Incorrect number of input arguments.')
end

% Check output file existence and OVERWRITE status
if exist([mesh(1:end-5),'_',regexprep(num2str(boundary),'  ','_'),'_boundary.txt'],'file') == 2
    % Output file exists, so check whether we're overwriting or not.
    if strcmpi(overwrite,'yes') ~= 1
        % Exit because we haven't been explicitly told to overwrite the
        % output files. This way, we default to safe behaviour.
        error('Aborting due to existence of output files.')
    end
end

fid = fopen(mesh);
if fid == -1
    error('Could not open input file %s.',mesh)
end
% Create output files
fidBoundary = fopen([mesh(1:end-5),'_',regexprep(num2str(boundary),'  ','_'),'_boundary.txt'],'w');
if fidBoundary == -1
    error('Could not create output file %s.',[mesh(1:end-5),'_',regexprep(num2str(boundary),'  ','_'),'_boundary.txt'])
end
% Count which line we're on.
counter = 1;
% To count where in the output array we've got to.
i = 1;

tline = fgets(fid);
while ischar(tline)
    % Skip the header, we're not interested in that.
    if counter ~= 1
        test = textscan(tline,'%s');
        % Here, we need to check how many arc attributes we're looking for.
        if numel(boundary)==1 % easy, only one value to check.
            % Node spatial info (5 columns). Filter based on BOUNDARY value.
            if numel(test{1}) == 5 && str2double(char(test{1}(end))) == boundary
                % We can't preallocate because we don't know how many boundary
                % points we have. Also, this is fugly...
                oline(i,1) = str2double(char(test{1}(2)));
                oline(i,2) = str2double(char(test{1}(3)));
                oline(i,3) = str2double(char(test{1}(4)));
                i = i+1;
            end
        else % More complicated, but nothing a for-loop can't solve.
            for checkBoundary=1:numel(boundary)
                if numel(test{1}) == 5 && str2double(char(test{1}(end))) == boundary(checkBoundary)
                    % We can't preallocate because we don't know how many boundary
                    % points we have. Also, this is fugly...
                    oline(i,1) = str2double(char(test{1}(2)));
                    oline(i,2) = str2double(char(test{1}(3)));
                    oline(i,3) = str2double(char(test{1}(4)));
                    i = i+1;
                end
            end
        end
    end
    tline = fgets(fid);
    counter = counter+1;
end

% The reason we haven't written the file as we're going along is we need to
% sort it first: without sorting the merge doesn't produce a sensible
% output file.
oline = sortrows(oline);
for i = 1:size(oline,1)
    fprintf(fidBoundary,'%f %f %f\n',oline(i,:));
end

fclose(fid);
fclose(fidBoundary);