function FE_solve(type)
%
% This function solves the system of linear equations arising from the
% finite element discretization of Eq. (17).  It stores the displacement 
% and the reaction forces in FE.U and FE.P.
% 
% type is either 'primal' or 'adjoint' 

global FE
p = FE.fixeddofs_ind;
f = FE.freedofs_ind;

primal = strcmpi(type, 'primal');
adjoint = strcmpi(type, 'adjoint');
if ~primal && ~adjoint
    warning('Unrecognized FE analysis type.');
end

% save the system RHS
if primal
    FE.rhs = FE.P(f,:) - FE.Kfp * FE.U(p,:);
elseif adjoint
    FE.prhs = FE.dJdu(f,:);
    FE.lambda = zeros(FE.n_global_dof,FE.nloads);
end


if strcmpi(FE.analysis.solver.type, 'direct')
    if FE.analysis.solver.use_gpu == true
       warning('GPU solver selected, but only available for iterative solver, solving on CPU.')
       FE.analysis.solver.use_gpu = false;
    end
    % Note: the direct solution works for multiple RHSs
    if primal
        FE.R = chol(FE.Kff); % Choleski factorization; store the factor for sensitivity analysis
        FE.U(f,:) = FE.R\(FE.R'\FE.rhs);
    elseif adjoint
        % Use Choleski factorization from prior primal analysis
        FE.lambda(f,:) = FE.R\(FE.R'\FE.prhs); %%%%%
    end
elseif strcmpi(FE.analysis.solver.type, 'iterative')
    tol = FE.analysis.solver.tol;
    maxit = FE.analysis.solver.maxit;
    % check if the user has specified use of the gpu
    if ~isfield(FE.analysis.solver,'use_gpu')
        FE.analysis.solver.use_gpu = false;
    end
    
    if FE.analysis.solver.use_gpu == true
        %% gpu solver
        ME.identifier = [];
        try
            gpu = gpuDevice(1);
            gpu.wait
            gpu.reset;
            A = gpuArray(FE.Kff); % copy FE.Kff to gpu
            for iload=1:FE.nloads                
                if primal
                    b = FE.rhs(:,iload);
                    x = FE.U(f,iload); % use last solution as initial guess
                elseif adjoint
                    b = FE.prhs(:,iload);
                    x = FE.lambda(f,iload);
                end

                M1 = diag(diag(FE.Kff)); % Jacobi preconditioner
                M2 = [];

                x = gather(pcg(...
                    A, ...
                    b, ... % need to loop over columns
                    tol, maxit, ...
                    M1,M2, ... % preconditioner(s)
                    x ... % use last solution as initial guess. 
                    )); 
                gpu.wait;
                if primal
                    FE.U(f,iload) = x;
                elseif adjoint
                    FE.lambda(f,iload) = x;
                end
            end
            gpu.reset;
        catch ME
            % something went wrong, display it and revert to cpu solver
            disp(ME.identifier);
            FE.analysis.solver.use_gpu = false;
        end
    elseif FE.analysis.solver.use_gpu == false
        %% cpu solver
        
        % If not running on GPUs, compute and store preconditioner only for
        % primal analysis and reuse for adjoint analysis.
        if primal
            ME.identifier = [];
            try
                % For some reason modified ichol does not work for this
                % problem (keeps encountering nonpositive pivots), so just
                % use plain vanilla ichol
                FE.L = ichol(FE.Kff);
                FE.ichol_exists = true;
            catch ME    
            end 

            if (strcmp(ME.identifier,'MATLAB:ichol:Breakdown'))
              warning ('ichol encountered nonpositive pivot, using Jacobi preconditioner.');
              FE.ichol_exists = false;

                % Use Jacobi preconditioner
                FE.M = diag(diag(FE.Kff));
             
            end
        end

        if primal
            if FE.ichol_exists
                for iload=1:FE.nloads
                    % Use solution of previous iteration as initial guess
                    FE.U(f,iload) = pcg(FE.Kff, FE.rhs(:,iload), tol, ...
                        maxit, FE.L,FE.L', FE.U(f,iload));
                end
            else
                for iload=1:FE.nloads
                    FE.U(f,iload) = pcg(FE.Kff, FE.prhs(:,iload), tol, ...
                        maxit, FE.M, [], FE.U(f,iload));
                end
            end
        elseif adjoint
             if FE.ichol_exists
                 for iload=1:FE.nloads
                    FE.lambda(f,iload) = pcg(FE.Kff, FE.prhs(:,iload), tol, ...
                        maxit, FE.L,FE.L.', FE.lambda(f,iload));
                 end
             else
                for iload=1:FE.nloads
                    FE.lambda(f,iload) = pcg(FE.Kff, FE.prhs(:,iload), tol, ...
                        maxit, FE.M, [], FE.lambda(f,iload));
                end
            end           
        end

    end
end

% Compute reaction forces
if primal
    FE.P(p,:) = FE.Kpp*FE.U(p,:) + FE.Kfp' * FE.U(f,:);
end
