function setup_bcs_Vframe2d()
%% Input file 
%
% *** THIS SCRIPT HAS TO BE CUSTOMIZED BY THE USER ***
%
% This script sets up the displacement boundary conditions and the forces
% for the analysis. 
%
% Important note: you must make sure you do not simultaneously impose 
% displacement boundary conditions and forces on the same degree of
% freedom.

% ** Do not modify this line **
global FE

coord_x = FE.coords(1,:);
coord_y = FE.coords(2,:);
if FE.dim == 3
    coord_z = FE.coords(3,:);
end

%% ============================
%% Compute predefined node sets
% compute_predefined_node_sets({'T_edge','TR_pt'})
% for an overview of this function, use: help compute_predefined_node_sets

% Apply force over 1/15 of the width around midpoint of top edge
W = FE.mesh_input.box_dimensions(1);
H = FE.mesh_input.box_dimensions(2);
VW = FE.mesh_input.cutout_dimensions(1);
VH = FE.mesh_input.cutout_dimensions(2);
d = (W-VW)/2;
lw = W/15;

tol = FE.max_elem_side/1000;
load_nodes = find( (FE.coords(2,:) > (H-tol)) & (abs(FE.coords(1,:) - W/2) < (lw + tol)) );
fixed_nodes = find( (FE.coords(2,:) < tol) & (FE.coords(1,:) < (d+tol)) );
roller_nodes = find( (FE.coords(2,:) < tol) & (FE.coords(1,:) > (W-d-tol)) );
compute_predefined_node_sets({'BL_pt'}); %<<
BL_pt = FE.node_set.BL_pt; %<<
%% ============================        

%% ============================        
%% Number of load cases
% Note: in this implementation, load cases can have different forces, but
% displacement boundary conditions must be the same for all load cases.
FE.nloads = 1; 

%% Applied forces

% Load case 1 
net_mag1 = -8;  % Force magnitude (net over all nodes where applied)
load_dir1 = 2;   % Force direction 
    
load_region1 = load_nodes;
load_mag1 = net_mag1/length(load_region1);

% Here, we build the array with all the loads.  If you have multiple
% applied loads, the load_mat array must contain all the loads as follows:
%  - There is one row per each load on a degree of freedom
%  - Column 1 has the node id where the load is applied
%  - Column 2 has the direction (1 -> x, 2 -> y, 3 -> z)
%  - Column 3 has the load magnitude.
%  - Column 4 has the load case number.
%
load_mat1 = zeros(length(load_region1),4);
load_mat1(:,1) = load_region1;
load_mat1(:,2) = load_dir1;
load_mat1(:,3) = load_mag1;
load_mat1(:,4) = 1;
n_force_dofs1 = size(load_mat1,1);


%% ============================        
%% Combine load cases if needed
load_mat = load_mat1;
n_force_dofs = n_force_dofs1;

%% ============================ 
%% Displacement boundary conditions
%
% NOTE: only one set of displacement boundary conditions is supported.
% I.e., there can be different forces in different load cases, but only one
% set of supports.
%
disp_region1 = fixed_nodes;
%>> disp_dirs_f1 = ones(1, length(disp_region1)); 
disp_dirs_f2 = 2*ones(1, length(disp_region1)); 
disp_mag1 = zeros(1, length(disp_region1));
disp_region2 = roller_nodes;
disp_dirs_r2 = 2*ones(1, length(disp_region2));
disp_mag2 = zeros(1, length(disp_region2));
disp_region3 = BL_pt; %<<
disp_dirs_bl1 = 1; %<<
disp_mag3 = 0; %<<



% Combine displacement BC regions
%>> disp_region = [disp_region1 disp_region1 disp_region2];
%>> disp_dirs = [disp_dirs_f1 disp_dirs_f2 disp_dirs_r2];
%>> disp_mag = [disp_mag1 disp_mag1 disp_mag2];
disp_region = [disp_region1 disp_region2 disp_region3];
disp_dirs = [disp_dirs_f2 disp_dirs_r2 disp_dirs_bl1];
disp_mag = [disp_mag1 disp_mag2 disp_mag3];

% In this example we are constraining both the x- and y-directions along
% the fixed nodes on the bottom left side, and the y-directions on the
% bottom right side.

% Here, we build the array with all the displacement BCs. 
% The disp_mat array must contain all the loads as follows:
%  - There is one row per each load on a degree of freedom
%  - Column 1 has the node id where the displacement BC is applied
%  - Column 2 has the direction (1 -> x, 2 -> y, 3 -> z)
%  - Column 3 has the displacement magnitude.
% 
disp_mat = zeros(length(disp_region),3);
for idisp=1:length(disp_region)
    disp_mat(idisp, 1) = disp_region(idisp);
    disp_mat(idisp, 2) = disp_dirs(idisp);
    disp_mat(idisp, 3) = disp_mag(idisp);
end

% *** Do not modify the code below ***
%
% Write displacement boundary conditions and forces to the global FE
% structure.
%
% Note: you must assign values for all of the variables below.
%
FE.BC.n_pre_force_dofs = n_force_dofs; % # of prescribed force dofs
FE.BC.force_node =  load_mat(:,1)';
FE.BC.force_dof = load_mat(:,2)';
FE.BC.force_value = load_mat(:,3)';
FE.BC.force_id = load_mat(:,4);
FE.BC.n_pre_disp_dofs = size(disp_mat,1); % # of prescribed displacement dofs
FE.BC.disp_node = disp_mat(:,1)';
FE.BC.disp_dof = disp_mat(:,2)';
FE.BC.disp_value = disp_mat(:,3)';
