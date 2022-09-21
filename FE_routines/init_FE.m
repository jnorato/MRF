function init_FE()
%
% Initialize the Finite Element structure
%
global FE

switch FE.mesh_input.type
    case 'generate'
        generate_mesh();
    case 'read-home-made'
        load(FE.mesh_input.mesh_filename);
    case 'read-gmsh'
        read_gmsh();
    case '2DLbracket'
        genLbracketmesh();
    case 'double2DLb'
        gendoubleLbracketmesh();    
    case 'V-frame'
        genVframemesh();         
    otherwise
        error('Unidentified mesh input type.');
end

% Plot mesh (usually uncommented, this is to help debugging mesh
% generation)
% plot_2dmesh(10);

% Compute element volumes and centroidal coordinates
FE_compute_element_info();

% Setup boundary conditions
run(FE.mesh_input.bcs_file);

% initialize the fixed/free partitioning scheme:
FE_init_partitioning();

% assemble the boundary conditions
FE_assemble_BC();

% compute elastic coefficients
FE_compute_constitutive_matrices();

% compute the element stiffness matrices
FE_compute_element_stiffness();