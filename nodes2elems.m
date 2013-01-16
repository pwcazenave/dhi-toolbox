function elems = nodes2elems(nodes, tri)
% Transfer a field from nodes to element centres.
%
% elems = nodes2elems(nodes, tri)
%
% DESCRIPTION:
%    Transfer values in nodes to the element centres described by the 
%    triangulation matrix tri.
%
% INPUT
%    nodes = 1D or 2D array of nodal values
%    tri = triangulation matrix
%
% OUTPUT:
%    elems = element centre based field
%
% EXAMPLE USAGE
%    elems = nodes2elems(nodes, tri)
%
% Author(s):
%    Pierre Cazenave (University of Southampton) based on nodes2elems from
%    fvcom-toolbox by Geoff Cowles (University of Massachusetts Dartmouth).
%
% Revision history
%
%==============================================================================

nelem = size(tri, 1);

if ndims(elems) == 1
    elems = zeros(nelem, 1);

    parfor i=1:nelem
        elems(i) = mean(nodes(tri(i, :)));
    end

elseif ndims(elems) == 2
    elems = zeros(nelem, size(elems, 2));

    parfor i=1:nelem
        elems(i) = mean(nodes(tri(i, :), :), 1);
    end

else
    error('Unsupported number of dimensions (maximum of two)')
end