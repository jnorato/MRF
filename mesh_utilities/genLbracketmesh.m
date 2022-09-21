function genLbracketmesh()
% This function generates an L-bracket mesh of specified element size and
% dimensions
%
% W is the side of the L-bracket, S is the width of the cutout.
% h is the element side.


global FE

W = FE.mesh_input.L_side;
S = FE.mesh_input.L_cutout;
h = FE.mesh_input.L_element_size;

% Check that W and S are multiples of h
if mod(W,h) > 0 || mod(S,h) > 0
    disp 'Element size is not an exact divisor of L-bracket dimensions.'
    return;
end


% Numbering of nodes and elements goes bottom to top and then left to right

% Nodes
xcoords = 0:h:W;
ycoords = 0:h:W;
inode = 0;
nnode_W = round(W/h) + 1;
nnode_S = round(S/h);
nnodes = nnode_W^2 - nnode_S^2;
node = zeros( nnodes, 2);
for i=1:length(xcoords)
    for j=1:length(ycoords)
        thisx = xcoords(i); thisy = ycoords(j);
        if ~(thisx > (W-S) && thisy > (W-S))
            inode = inode + 1;
            node(inode,:) = [thisx, thisy];
        end
    end
end

% Elements
ielem = 0;
nelem_W = nnode_W-1;
nelem_S = nnode_S;
nelem = nelem_W^2 - nelem_S^2;
elem = zeros( nelem, 4);

% Vertical part of L-bracket + first column of elements of horizontal part
for col=1:nelem_W - nelem_S + 1
    for row=1:nelem_W
        if ~(col > nelem_W - nelem_S && row > nelem_W - nelem_S)
            ielem = ielem + 1;
            bl = nnode_W*(col-1) + row; 
            br = nnode_W*col + row;
            ur = br + 1;
            ul = bl + 1;
            elem(ielem,:) = [ bl br ur ul];
        end
    end
end

shift = nnode_W*(nnode_W-nnode_S);
% Rest of elements of horizontal part of L-bracket
for col=1:nelem_S - 1
    for row=1:nelem_W-nelem_S
        ielem = ielem + 1;
        bl = (nnode_W-nnode_S)*(col-1) + row + shift;
        br = (nnode_W-nnode_S)*col + row + shift;
        ur = br + 1;
        ul = bl + 1;
        elem(ielem,:) = [ bl br ur ul];     
    end
end

% Store mesh in global array
FE.dim = 2;
FE.n_elem = nelem;
FE.n_node = nnodes;
FE.coords = node';
FE.elem_node = elem';



            