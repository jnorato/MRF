function [c,grad_c] = compute_compliance()
%
% This function computes the mean compliance and its sensitivities
% based on the last finite element analysis
%
global FE OPT

% compute the compliance 
    c = full(dot(FE.U,FE.P));

% compute the design sensitivity
    Ke = FE.Ke;
    Ue = permute(repmat(...
        FE.U(FE.edofMat).',...
            [1,1,FE.n_edof]), [1,3,2]);
    Ue_trans = permute(Ue, [2,1,3]);

    Dc_Dpenalized_elem_dens = reshape(sum(sum( ...
        -Ue_trans.*Ke.*Ue, ...
            1),2),[1,FE.n_elem]);   

    grad_c = OPT.H' * (Dc_Dpenalized_elem_dens' .* OPT.dpen_rho_e); 
% save these values in the OPT structure
    OPT.compliance = c;
    OPT.grad_compliance = grad_c;