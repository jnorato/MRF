function [g, geq, gradg, gradgeq] = nonlcon(dv)
% [g, geq, gradg, gradgeq] = nonlcon(dv) returns the costraints

    global  OPT
    
    OPT.dv_old = OPT.dv;
    OPT.dv = dv(:);
    
    % Perform analysis only if design has changed
    tol = 1e-12;
    if max(OPT.dv - OPT.dv_old) > tol
        perform_analysis();
    end

    n_con = OPT.functions.n_func-1; % number of constraints
    g = zeros(n_con,1);
    gradg = zeros(OPT.n_dv,n_con);
    for i = 1:n_con
        g(i) = OPT.functions.f{i+1}.value;
        g = g - OPT.functions.constraint_limit(i);
        gradg(:,i) = OPT.functions.f{i+1}.grad;
    end
    geq = [];
    gradgeq = [];