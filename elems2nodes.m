function nodes = elems2nodes(elems, tri, nvert)
% Transfer a field from elements to vertices.
%
% nodes = elems2nodes(elems, tri, nvert)
%
% DESCRIPTION:
%    Transfer values in elems to the nodes in tri.
%
% INPUT
%    elems = 1D or 2D array of element centre values
%    tri = triangulation matrix
%    nvert = number of nodes in the converted values
%
% OUTPUT:
%    nodes = vertex-based field
%
% EXAMPLE USAGE
%    nodes = smoothfield(elems, tri, 1000)
%
% Author(s):
%    Pierre Cazenave (University of Southampton) based on elems2nodes from
%    fvcom-toolbox by Geoff Cowles (University of Massachusetts Dartmouth).
%
% Revision history
%
%==============================================================================

nelem = size(tri, 1);

count = zeros(nvert, 1);

if ndims(elems) == 1
    nodes = zeros(nvert, 1);

    for i=1:nelem
        n1 = tri(i, 1);
        n2 = tri(i, 2);
        n3 = tri(i, 3);
        nodes(n1) = nodes(n1) + elems(i);
        nodes(n2) = nodes(n2) + elems(i);
        nodes(n3) = nodes(n3) + elems(i);
        count(n1) = count(n1) + 1;
        count(n2) = count(n2) + 1;
        count(n3) = count(n3) + 1;
    end
    nodes = nodes./count;

elseif ndims(elems) == 2
    nodes = zeros(nvert, size(elems, 2));

    for i=1:nelem
        n1 = tri(i, 1);
        n2 = tri(i, 2);
        n3 = tri(i, 3);
        nodes(n1, :) = nodes(n1) + elems(i, :);
        nodes(n2, :) = nodes(n2) + elems(i, :);
        nodes(n3, :) = nodes(n3) + elems(i, :);
        count(n1) = count(n1) + 1;
        count(n2) = count(n2) + 1;
        count(n3) = count(n3) + 1;
    end
    nodes =  nodes ./ repmat(count, 1, size(elems, 2));

else
    error('Unsupported number of dimensions (maximum of two)')
end

