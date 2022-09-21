function [maxgrf,grad_maxgrf] = compute_maxGRF()
%
% This function computes the volume fraction and its sensitivities
% based on the last geometry projection
%
global OPT

% compute the volume fraction
    elgrf = 4*OPT.filt_rho_e.*(1 - OPT.filt_rho_e);
    p=24;
    [maxgrf, dSdx] = smooth_max(elgrf,p,'average');


% compute the design sensitivity
    grad_maxgrf = OPT.H' * (4*dSdx.*(1-2*OPT.filt_rho_e));
    
% output
    OPT.maxGRF = maxgrf;
    OPT.grad_maxGRF = grad_maxgrf;