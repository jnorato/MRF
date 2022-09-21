function gendoubleLbracketmesh()
% This function generates a double L-bracket mesh of specified element size 
% and dimensions.
%
% W is the entire width of the double L-bracket, S is the width of the 
% cutout on each side. The bracket height is H.
% h is the element side.

% This function uses a different strategy than genLbracketmesh to generate
% the mesh.

global FE

W = FE.mesh_input.L_width;
H = FE.mesh_input.L_height;
S = FE.mesh_input.L_cutout;
h = FE.mesh_input.L_element_size;

% Check that W and S are multiples of h
if mod(W,h) > 0 || mod(S,h) || mod(H,h) > 0
    disp 'Element size is not an exact divisor of L-bracket dimensions.'
    return;
end


% Numbering of nodes and elements goes bottom to top and then left to right

% Nodes
xcoords = 0:h:W;
ycoords = flip(0:h:H);
[X,Y] = meshgrid(xcoords,ycoords);
inode = 0;
nnode_x = round(W/h) + 1;
nnode_y = round(H/h) + 1;
nnode_S = round(S/h);
nnodes = (nnode_x*nnode_y) - 2*nnode_S^2;
node = zeros( nnodes, 2);
% Array with node numbers
N = zeros(size(X));
for col=1:nnode_x
    for row=1:nnode_y
        thisx = X(row,col); thisy = Y(row,col);
        if ~( (thisx > (W-S) || thisx < S) && thisy > (H-S))
            inode = inode + 1;
            node(inode,:) = [thisx, thisy];
            N(row,col) = inode;
        end
    end
end

% Elements
ielem = 0;
nelem_x = nnode_x-1;
nelem_y = nnode_y-1;
nelem_S = nnode_S;
nelem = nelem_x*nelem_y - 2*nelem_S^2;
elem = zeros( nelem, 4);

for col=1:nelem_x
    for row=1:nelem_y
        this_elem_nodes = [N(row,col) N(row+1,col) N(row+1,col+1) N(row,col+1)];
        % Create element only if all nodes are non-zero (i.e., not in the 
        % cutout regions)
        if prod(this_elem_nodes)~=0
            ielem = ielem+1;
            elem(ielem,:) = this_elem_nodes;
        end
    end
end


% Store mesh in global array
FE.dim = 2;
FE.n_elem = nelem;
FE.n_node = nnodes;
FE.coords = node';
FE.elem_node = elem';



            