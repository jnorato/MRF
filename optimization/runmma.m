function history = runmma(x0,obj,nonlcon)
%
% Perform the optimization using MMA
%
global OPT FE


% Initialize history object
history.x = [];
history.fval = [];
history.fconsval = [];
history.grf = [];
history.norm_dfdx = [];

% Initialize lower and upper bounds vectors
switch OPT.parameters.penalization_scheme
    case {'SIMP', 'RAMP'}
        lb = FE.material.rho_min*ones(size(OPT.dv)); 
    case {'modified_SIMP', 'modified_RAMP'} 
        % lb = (FE.material.rho_min^2)*ones(size(OPT.dv));
        lb = (FE.material.rho_min/100)*ones(size(OPT.dv));
    otherwise
        warning('Unrecognized penalization scheme.');
        return;
end
ub = ones(size(OPT.dv)); 

%
ncons = OPT.functions.n_func - 1;  % Number of optimization constraints
ndv = OPT.n_dv; % Number of design variables

% Initialize vectors that store current and previous two design iterates
x = x0;
xold1 = x0; 
xold2 = x0;

% Initialize move limits 
ml_step = OPT.options.move_limit * abs(ub - lb);  % Compute move limits once

% Initialize lower and upper asymptotes
low = lb;
upp = ub;

% These are the MMA constants (Svanberg, 1998 DACAMM Course)
c = OPT.parameters.mma.c;
d = OPT.parameters.mma.d;
a0 = OPT.parameters.mma.a0;
a = OPT.parameters.mma.a;

% Evaluate the initial design and print values to screen 
iter = 0;
[f0val , df0dx] = obj(x);
[fval, ~, dfdx, ~] = nonlcon(x);
dfdx = dfdx';
fprintf('It. %i, Obj= %-12.5e, ConsViol = %-12.5e\n', ...
    iter, f0val, max(max(fval, zeros(ncons,1))));

%%%
% Save initial design to history
history.fval = [history.fval; f0val];
history.fconsval = [history.fconsval; fval];
history.x = [history.x x(:)];
if OPT.stress_needed
    history.gstress = [OPT.approx_h_max];
    history.true_stress_max = [OPT.true_stress_max];
end

%%%
% Plot initial design 
plotfun(iter);
          
%%%% Initialize stopping values
obj_change = 10*OPT.options.obj_tol;
OPT.grf = 4*dot(OPT.dv.*(1 - OPT.dv), FE.elem_vol)/sum(FE.elem_vol);

% Parameters for continuation on aggregation and rectifier functions
if OPT.parameters.continuation && OPT.stress_needed
    maxB = OPT.parameters.maxB;
    deltaB = OPT.parameters.deltaB;
    maxK = OPT.parameters.maxK;
    deltaK = OPT.parameters.deltaK;
    cont_started = false;
    cont_ended = false;
end



%
% ******* MAIN MMA LOOP STARTS *******
% 

while iter < OPT.options.max_iter
    
    iter = iter+1;

    % Impose move limits by modifying lower and upper bounds passed to MMA
    % Double move limits if design is feasible
    if OPT.stress_needed 
        if all(OPT.true_stress_max <= 1.2*OPT.parameters.slimit)
            step_fac = 2.0;
        else
            step_fac = 1.0;
        end
    else
            step_fac = 1.0;
    end
    mlb = max(lb, x - step_fac*ml_step);
    mub = min(ub, x + step_fac*ml_step); 

    %** continuation for aggregation and rectifier parameters
    if OPT.parameters.continuation && OPT.stress_needed
        if all(OPT.true_stress_max <= 2*OPT.parameters.slimit)
            OPT.parameters.aggregation_parameter = min(maxB, ...
                OPT.parameters.aggregation_parameter + deltaB);
            OPT.parameters.rectifier_parameter = min(maxK, ...
                OPT.parameters.rectifier_parameter + deltaK);
            if ~cont_started
                disp('**********');
                fprintf('Started aggregation continuation; GRF = %-12.5e.\n', OPT.grf);
                disp('**********');
                cont_started = true;
            end
        end

        if OPT.parameters.aggregation_parameter == maxB && ~cont_ended
            disp('**********');
            fprintf('Ended aggregation continuation at beta = %12.5e.\n', ...
                    OPT.parameters.aggregation_parameter);
            disp('**********');    
            cont_ended = true;
        end
    end
    
    
    % If requested, adaptive constraint scaling is done here. The reason it
    % is done here as opposed to in the function that evaluates the
    % response, is not to mess up the finite difference check, since the
    % ACS is inconsistent.
    if strcmpi(OPT.parameters.aggregation_type, 'p-norm') && OPT.stress_needed
        if OPT.parameters.ACS.use
            if iter <= 3
                cc = 1;
            else
                c_old3 = OPT.parameters.ACS.c(iter-3);
                c_old2 = OPT.parameters.ACS.c(iter-2);
                c_old1 = OPT.parameters.ACS.c(iter-1);
                if (c_old3 - c_old2)*(c_old2 - c_old1) < 0 
                    alpha = OPT.parameters.ACS.alpha_osc;
                else
                    alpha = OPT.parameters.ACS.alpha_no_osc;
                end
                cc = OPT.true_h_max/OPT.approx_h_max; % This works!
                cc = alpha*cc + (1-alpha)*c_old1;
            end
            OPT.parameters.ACS.c = [OPT.parameters.ACS.c cc];
        else
            cc = 1;
        end
        fval(1,1) = fval(1,1) - 1/cc + 1;
    end
    
    % Function scaling
    f0val = OPT.functions.objective_scale * f0val;
    df0dx = df0dx * OPT.functions.objective_scale;
    fval = OPT.functions.constraint_scale' .* fval;
    dfdx = OPT.functions.constraint_scale' .* dfdx;


    %%%% Solve MMA subproblem for current design x
    switch OPT.parameters.mma.version
        case  2007
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                 mmasub(ncons,ndv,iter,x,mlb,mub,xold1, ...
                     xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
        case 1999
            [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
                mma1999(ncons,ndv,iter,x,mlb,mub,xold1,xold2, ...
                f0val,df0dx,0*df0dx,fval,dfdx,0*dfdx,low,upp,a0,a,c,d);
        otherwise
            disp('Error: unknown version of MMA');
            exit;
    end    
    

    %%%% Updated design vectors of previous and current iterations
    xold2 = xold1;
    xold1 = x;
    x  = xmma;
    
    % Update function values and gradients
    % Note that OPT.dv gets updated inside these functions
    [f0val , df0dx] = obj(x);
    [fval, ~, dfdx, ~] = nonlcon(x);
    dfdx = dfdx';  
    
    % Compute norm of KKT residual vector
    [residu,kktnorm,residumax] = ...
    kktcheck(ncons,ndv,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
           lb,ub,df0dx,fval,dfdx,a0,a,c,d);
    
    % Produce output to screen
    if ~OPT.stress_needed
        fprintf('It. %i, Obj= %-12.5e, ConsViol = %-12.5e, KKT-norm = %-12.5e, Obj change = %-12.5e\n', ...
            iter, f0val, max(max(fval, zeros(ncons,1))), kktnorm, obj_change);
    else
        fprintf('It. %i, Obj= %-12.5e, ConsViol = %-12.5e, KKT-norm = %-12.5e, Obj change = %-12.5e, True max. stress = [', ...
            iter, f0val, max(max(fval, zeros(ncons,1))), kktnorm, obj_change);      
        fprintf('%-12.5e', OPT.true_stress_max(1:end));
        fprintf(']\n');
    end
    
%     % Save design to .mat file
    [folder, baseFileName, ~] = fileparts(OPT.options.mat_output_path);
    mat_filename = fullfile(folder, strcat(baseFileName, '.mat'));
    save(mat_filename, 'OPT');
    
    % Write to vtk file if requested.  
    if strcmpi(OPT.options.write_to_vtk, 'all')
        writevtk(OPT.options.vtk_output_path, 'dens', iter);
    end    
    
    % Update history
    history.fval = [history.fval; f0val];
    history.fconsval = [history.fconsval; fval];
    history.x = [history.x x(:)];
    if OPT.stress_needed
        history.gstress = [history.gstress OPT.approx_h_max];
        history.true_stress_max = [history.true_stress_max; OPT.true_stress_max];
    end
    
    % Compute change in objective function
    % Check only after first iteration
        
    if iter > 1
        obj_change = ((history.fval(iter)-history.fval(iter-1))/history.fval(iter-1));
        if abs(obj_change) < OPT.options.obj_tol
            fprintf('Objective function convergence tolerance satisfied.\n');
        end
    end
    
    % Compute gray region fraction, which serves as an indication of
    % convergence to 0-1 design.
    OPT.grf = 4*dot(OPT.dv.*(1 - OPT.dv), FE.elem_vol)/sum(FE.elem_vol);
    history.grf = [history.grf OPT.grf];
    
    % Store largest magnitude of stress constraint sensitivities
    if OPT.stress_needed
        history.norm_dfdx = [history.norm_dfdx max(abs(dfdx))];
    end
        
    % Plot current design
    plotfun(iter);
    
    % >>>>
    % Plot sensitivities
    % cap_up = 0.001; cap_low = -0.001;
    % plotfield2d(5, dfdx, cap_up, cap_low, true)
    % <<<
    
    % =================
    % Stopping criteria
    if abs(obj_change) < OPT.options.obj_tol && ...
            OPT.grf < OPT.options.max_GRF
        if (OPT.parameters.continuation && cont_ended) || ~OPT.parameters.continuation
            fprintf('Satisfied stopping criteria.\n');
            break;
        end
    end
    % =================    
%
end

% In the case of some 2d runs, plot the  outline
if FE.dim == 2 && OPT.options.plot
    plot_outline(1);
end

% In the case of 2d runs, plot the stress for each load case
offset_fig_n = 2; % Two other figures have been plotted before (dens and hist)
if FE.dim == 2 && OPT.options.plot
    for il=1:FE.nloads
        if ~ismember('maximum stress violation',OPT.functions.constraints) % If no stress constraints
            compute_max_stress_violation();
            lim = max(FE.svm(:,il));
        else
            lim = OPT.parameters.slimit(il);
        end
        plotstress2d(il+offset_fig_n, lim, il, false);
    end
end

% Write vtk for final iteration if requested
if strcmp(OPT.options.write_to_vtk, 'all') || ...
        strcmp(OPT.options.write_to_vtk, 'last')
    if ~OPT.stress_needed % Compute stresses regardless of whether or not there is a stress constraint
        compute_max_stress_violation();
    end
    OPT.write_stress_to_vtk = true;
    writevtk(OPT.options.vtk_output_path, 'dens', iter);
end 

% Store history array in OPT structure so that it can be recovered later
OPT.history = history;

% ============================================


    function plotfun(iter)
        % Note that this function has a slightly different format than its
        % equivalent for fmincon.
        
        if OPT.options.plot == true
            figure(1)
            plot_density(1)
            axis equal
            xlim([FE.coord_min(1), FE.coord_max(1)])
            ylim([FE.coord_min(2), FE.coord_max(2)])
            if FE.dim == 2
                view(2)
            else
                zlim([FE.coord_min(3), FE.coord_max(3)])
                view([50,22])
            end
            drawnow;
            
        end
    end
end