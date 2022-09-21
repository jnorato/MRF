function [] = FE_assemble_stiffness_matrix()
%
% FE_ASSEMBLE assembles the global stiffness matrix, partitions it by 
% prescribed and free dofs, and saves the known partitions in the FE structure.

%% Reads: 
%        FE.
%           iK
%           jK
%           sK
%           fixeddofs_ind
%           freedofs_ind
%       OPT.
%           penalized_elem_dens

%% Writes:
%        FE.
%           Kpp 
%           Kpf
%           Kff

%% Declare global variables
global FE OPT

%% assemble and partition the global stiffness matrix

% Retrieve the penalized stiffness
    penalized_rho_e = permute(repmat(...
        OPT.pen_rho_e(:),...
            [1,FE.n_edof,FE.n_edof]),[2,3,1]);
        
% Ersatz material: 
penalized_Ke = penalized_rho_e .* FE.Ke;
FE.sK_penal = penalized_Ke(:);

% assemble the penalized global stiffness matrix (the sparse functions
% accumulates values with repeated indices, which allows to assemble the
% global stiffness matrix simply with this line).
K = sparse(FE.iK,FE.jK,FE.sK_penal);
K = (K+K')/2;

% partition the stiffness matrix and return these partitions to FE
FE.Kpp = K(FE.fixeddofs_ind,FE.fixeddofs_ind);
FE.Kfp = K(FE.freedofs_ind,FE.fixeddofs_ind);
% note: by symmetry Kpf = Kfp', so we don't store Kpf. Tall and thin
% matrices are stored more efficiently as sparse matrices, and since we
% generally have more free dofs than fixed, we choose to store the rows as
% free dofs to save on memory.
FE.Kff = K(FE.freedofs_ind,FE.freedofs_ind);

