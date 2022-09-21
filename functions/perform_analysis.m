function [] = perform_analysis()
%
% Filter and penalize densities, solve the finite
% element problem for the displacements and reaction forces, and then
% evaluate the relevant functions.
%
global OPT FE

% Filter densities
OPT.filt_rho_e = OPT.H * OPT.dv;

% Penalize densities
switch OPT.parameters.penalization_scheme
    case {'SIMP', 'RAMP'}
        [OPT.pen_rho_e, OPT.dpen_rho_e] = penalize(OPT.filt_rho_e, ...
            OPT.parameters.penalization_param, ...
            OPT.parameters.penalization_scheme);
    case {'modified_SIMP', 'modified_RAMP'}
        [OPT.pen_rho_e, OPT.dpen_rho_e] = penalize(OPT.filt_rho_e, ...
            OPT.parameters.penalization_param, ...
            OPT.parameters.penalization_scheme, ...
            FE.material.rho_min);
    otherwise
        warning('Unrecognized penalization scheme.');
        return;
end

% Perform FE analysis
FE_analysis();

% Evaluate objective and constraint functions
evaluate_relevant_functions();

end

