function [f, gradf] = obj(dv)
    global  OPT
    
    OPT.dv_old = OPT.dv; % save the previous design
    OPT.dv = dv(:); % update the design
    
   
    % Update or perform the analysis only if design has changed
    tol = 1e-12;
    
    if max(OPT.dv - OPT.dv_old) > tol
        perform_analysis();
    end

    f = OPT.functions.f{1}.value;
    gradf = OPT.functions.f{1}.grad;
    