function inputs_Vframe2d()
%% Input file 
%
% *** THIS SCRIPT HAS TO BE CUSTOMIZED BY THE USER ***
%
% This script sets up all the parameters for the optimization, and also
% indicates what files contain the mesh, boundary conditions, and initial
% design.
%
%
% It is recommended that for different problems you make copies of this
% file in the input_files subfolder, so that to switch from one problem to
% another you need only change the run statement in the inputs.m file in
% the root folder.
%


% ** Do not modify this line **
global FE OPT 

% Set this flag to true if you want plotting during the optimization
plot_cond = true; 

% Set this flag to true if you want copying of output files
OPT.options.save_outputs = true;

%% =======================================================================
%% Mesh information
%
% First of all, an important clarification: in this code, we refer to mesh
% as exclusively the nodal coordinates and the element nodal connectivity
% array. You have to write/modify a separate matlab script (setup_bcs) to
% set up the loads and displacement boundary conditions).
%
% This code provides three options to populate a mesh, indicated in the 
% FE.mesh_input.type field: 'generate','read-home-made', 'read-gmsh' and
% '2DLbracket':
%
% 1- 'generate': Generate a rectangular/parallelepiped mesh on the fly by 
%                providing the dimensions and element size (using the 
%                generate_mesh function).
% 2- 'read-home-made': Load mesh from Matlab .mat file 
%                      (which you can create before running this code using 
%                      the makemesh function).
% 3- 'read-gmsh':  Read quadrilateral or hexahedral mesh generated by Gmsh 
%                  and exported to the Matlab format. For this code, we
%                  tested version 4.0.6 of Gmsh (but it is possible that
%                  earlier versions work; if so, try at your own risk).
% 4- '2DLbracket': Generate a 2D L-bracket mesh on the fly by providing the
%                  dimensions 
% 5- 'V-frame':    Generate rectangular mesh with inverted 'V' cutout by
%                  providing dimensions

FE.mesh_input.type = 'V-frame';

% If mesh input type is 'generate', you must specify the dimensions of
% the rectangle/cuboid and the number of elements along each direction:
% FE.mesh_input.box_dimensions = [20 5]; % ***
% FE.mesh_input.elements_per_side = [200 50]; % ***
% FE.mesh_input.elements_per_side = [16 4]; % ***

% If mesh input type is 'read-home-made', you must provide a
% mesh file name including extension (*.mat).
% 
% NOTE: all folders are relative to the root folder where the GPTO_b.m
% script is located.
%
% FE.mesh_input.mesh_filename = 'input_files/cantilever2d/2drectangle.mat';

% If mesh input type is 'read-gmsh', you must provide a
% mesh file name including extension (*.m). To produce this file, you must
% first generate a transfinite mesh (only quad elements in 2-d, only hexa
% elements in 3-d) in Gmsh, and then export it with Matlab format
% (including the .m extension).
% FE.mesh_input.gmsh_filename = 'input_files/cantilever2d/cantilever2d.m';

% If mesh input type is '2DLBracket', you must provide the dimension of
% the side of the square box of the L-shape, the % length of the side of 
% the square cutout on the top-left side, and the element size.
%
% FE.mesh_input.L_side = 100;
% FE.mesh_input.L_cutout = 60;
% FE.mesh_input.L_element_size = 2;

% If mesh input type is 'V-frame', you must specify the dimensions of
% the rectangle, the number of elements along each direction, and the
% width and height of the V-shaped cutout.
FE.mesh_input.box_dimensions = [120 60]; % ***
FE.mesh_input.cutout_dimensions = [110 30]; % ***
FE.mesh_input.elements_per_side = [240 80]; % ***


%% =======================================================================
%% Boundary conditions

    % Here, you must specify the path to a Matlab script file that sets up the
    % boundary conditions (which you must modify according to the problem)
    FE.mesh_input.bcs_file = 'input_files/Vframe2d/MRF/setup_bcs_Vframe2d.m';


%% =======================================================================
%% Material information
%
    % Specify the Young's modulus and Poisson ratio 
    FE.material.E = 1; % Young's modulus of the design material
    FE.material.nu = 0.3; % Poisson ratio of the design material

    FE.material.rho_min = 1e-3; % Minimum density in void region
    FE.material.nu_void = 0.3;  % Poisson ratio of the void material

%% =======================================================================    
%% Finite element solver
    FE.analysis.solver.type = 'direct'; % 'direct';  % Options: 'direct' or 'iterative'
    FE.analysis.solver.tol = 1e-8;          % only for iterative
    FE.analysis.solver.maxit = 1e4;         % only for iterative
    FE.analysis.solver.use_gpu = false; 
    % NOTE: gpus can only be used if:
    % (a) FE.analysis.solver.type = 'iterative'
    % (b) the system has a compatible nvidia gpu (and drivers)
    %     You can query your system's gpu with the matlab command
    %     >> gpuDevice()
    % The gpu solver may be slower than an iterative solver on the cpu for
    % smaller problems because of the cost to transfer data to the gpu. 
    
%% =======================================================================        
%% Optimization problem definition
% functions:
    % Name of objective function
    OPT.functions.objective = 'volume fraction';
    OPT.functions.objective_scale = 1.0;
    % Names of inequality (<=) constraints
    OPT.functions.constraints = {'maximum stress violation'};
    % Inequality constraint (upper) limits vector: should have the
    % constraint limit for each one of the constraints.
    % For stress constrains, use 0, since the stress_limit is incorporated
    % through the OPT.parameters.slimit variable.
    OPT.functions.constraint_limit = [0.0]; % ***
    OPT.functions.constraint_scale = [1];

%% =======================================================================        
%% Penalization and filtering parameters 
    % Penalization scheme 
    OPT.parameters.penalization_scheme = 'modified_SIMP'; % Options: 'SIMP', 'modified_SIMP', 'RAMP', 'modified_RAMP' 
    % Parameter to be used for the penalization
    OPT.parameters.penalization_param = 3;
    % Filtering radius given as a factor of largest element side
    OPT.parameters.filter_radius_factor = 2.5;
    
%% =======================================================================        
%% Stress-based optimization parameters    
%
% Required and suggested settings:
% For p-norm + ACS approach: 
%   Required: aggregation_type = 'p-norm', ACS.use = true, 
%             aggregation_function = 'p-norm', continuation = false;
% For MRF:
%   Required: aggregation_type = 'mrf', ACS.use = false, continuation =
%             true, aggregation_function = 'KS', rectifier_function =
%             'shifted_KS'
%              
%
    % Relaxation power 
    OPT.parameters.relaxation_param = 2.5;
    % Type of aggregation function. Options: 'p-norm', 'mrf'
    OPT.parameters.aggregation_type = 'mrf';
    % Whether or not to do continuation on the aggregation and rectifier
    % parameters
    OPT.parameters.continuation = true;
    % Aggregation parameter (p for p-norm, k for KS)
    OPT.parameters.aggregation_parameter_init = 6;
    OPT.parameters.aggregation_parameter = 16;
    OPT.parameters.aggregation_parameter_delta = 0.2;
    % Global aggregation function; options: 'p-norm', 'p-mean', 'KS',
    % 'LKS', 'avgKS', 'softmax', 'average'.
    OPT.parameters.aggregation_function = 'LKS';
    % Rectifier function parameter
    OPT.parameters.rectifier_parameter_init = 2*OPT.parameters.aggregation_parameter_init;
    OPT.parameters.rectifier_parameter = 2*OPT.parameters.aggregation_parameter;
    OPT.parameters.rectifier_parameter_delta = 2*OPT.parameters.aggregation_parameter_delta;
    % Rectifier function; options: 'shiftedKS' (recommended), 'softmax',
    % 'softplus', 'LKS', 'intH', 'ELU'
    OPT.parameters.rectifier_function = 'shiftedKS';
    % In the case of shifted KS, specify shift
    OPT.parameters.rectifier_eps = 1e-3;
    % Stress limit for each element
    OPT.parameters.slimit = 3.3;
    % Adaptive constraint scaling
    OPT.parameters.ACS.use = false;
    OPT.parameters.ACS.alpha_osc = 0.8;
    OPT.parameters.ACS.alpha_no_osc = 1;
    OPT.parameters.ACS.c = [];  

%% =======================================================================        
%% Density for initial design
    OPT.parameters.init_dens = 0.5;

%% =======================================================================        
%% Optimization parameters
    % Version of MMA; options are 1999 and 2007
    OPT.parameters.mma.version = 1999;  
    ncons = length(OPT.functions.constraints);
    OPT.parameters.mma.c = 1000*ones(ncons,1); % default = 1000
    OPT.parameters.mma.d = ones(ncons,1);
    OPT.parameters.mma.a0 = 1;
    OPT.parameters.mma.a = zeros(ncons, 1);        
    % Whether plots should be produced or not 
    OPT.options.plot = plot_cond; 
    % Write to a vkt file; options are 'none', 'last' (only write last 
    % iteration)and 'all' (write all iterations).  
    OPT.options.write_to_vtk = 'last';
    % Recall that paths are relative to the root folder 
    OPT.options.vtk_output_path = 'output_files/Vframe2d/MRF/';
    % Path to write .mat file with densities for current design
    OPT.options.mat_output_path = 'output_files/Vframe2d/MRF/';
    % Path to save files if requested
    OPT.options.outputs_path = strcat('output_files/Vframe2d/MRF/');    
    % Move limits as a fraction of the range between bounds 
    OPT.options.move_limit = 0.02;  
    % Maximum number of iterations 
    OPT.options.max_iter = 500; 
    % Maximum gray region fraction
    OPT.options.max_GRF = 0.1;
    % Minimum change in objective function to keep iteration
    OPT.options.obj_tol = 1e-4;  % *** 

%% =======================================================================        
%% Sensitivities finite difference check
%
% These options allow you to run a finite difference check on the cost
% function, and/or the constraint.
%
% Please note that if you do a finite difference check, the code will stop
% right after the check and not continue with the optimization (typically
% you want to do one or the other but not both).
%
    % Whether or not to perform sensitivities finite difference check
    OPT.make_fd_check = false;
    % Step size for finite difference
    OPT.fd_step_size = 1e-5;
    % Whether or not to check cost function sensitivities
    OPT.check_cost_sens = true;
    % Whether or not to check constraint sensitivities
    OPT.check_cons_sens = true;
    
    
