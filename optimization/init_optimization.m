function init_optimization()
global OPT FE 

% Initialize functions to compute
% Concatenate list of functions to be computed
f_list = {OPT.functions.objective, OPT.functions.constraints{:}};

% here we list all the functions that are available to compute as f{i}

f{1}.name = 'compliance';
f{1}.function = 'compute_compliance';

f{2}.name = 'volume fraction';
f{2}.function = 'compute_volume_fraction';

f{3}.name = 'maximum stress violation';
f{3}.function = 'compute_max_stress_violation';

f{4}.name = 'maxGRF';
f{4}.function = 'compute_maxGRF';

% compare all functions available with the ones specified in inputs.m
n = length(f);
m = length(f_list);
for j = 1:m
    for i = 1:n
        if strcmpi( f{i}.name, f_list{j})
            OPT.functions.f{j} = f{i};
        end
    end
end

OPT.functions.n_func =  numel(OPT.functions.f);

% Determine if stress will be computed
OPT.stress_needed = false;
OPT.write_stress_to_vtk = false;
for i = 1:OPT.functions.n_func
    if strcmpi(OPT.functions.f{i}.name, 'maximum stress violation')
        OPT.stress_needed = true;
        OPT.write_stress_to_vtk = true;
        break;
    end
end

OPT.n_dv = FE.n_elem;   % Number of design variables
OPT.dv = OPT.parameters.init_dens*ones(OPT.n_dv,1); % Initial design

% Uncomment the following line to override the initial design with a
% random design, which is useful to test the correctness of the
% sensitivities calculations taking into account the filter.
% rng(1); OPT.dv = rand(OPT.n_dv,1);

% =================================
% Construct and store filter matrix

% For each element, find all neighboring elements whose centroids are at a
% distance <= the filter radius from the element's centroid.
% Use Matlab's rangesearch function, which uses an efficient kd-tree
% algorithm.
if OPT.parameters.filter_radius_factor < 1
    disp 'Error: filter radius factor must be at least 1.';
    return;
end
OPT.filter_radius = OPT.parameters.filter_radius_factor*FE.max_elem_side;
% Use a small tolerance to avoid including elements with zero weight in the
% filter matrix.
tol = FE.max_elem_side/1000;
search_dist = OPT.filter_radius - tol;
[neighbors_id,dist] = rangesearch(FE.centroids',FE.centroids',search_dist);

% Number of neighbors for each element
n_neigh = cell2mat(cellfun(@size,neighbors_id,'uni',false));
max_neigh = max(n_neigh(:,2));
iH = zeros(FE.n_elem * max_neigh,1);
jH = zeros(size(iH));
sH = zeros(size(iH));
idx=0;
for iel=1:FE.n_elem
    init_idx = idx + 1;
    final_idx = idx + n_neigh(iel,2);
    iH(init_idx:final_idx) = iel;
    jH(init_idx:final_idx) = neighbors_id{iel};
    num = 1 - dist{iel}/OPT.filter_radius;
    den = sum(num);
    sH(init_idx:final_idx) = num/den;
    idx = final_idx;
end
OPT.H = sparse(iH(1:idx), jH(1:idx), sH(1:idx));

% =================================
% Setup continuation on aggregation and rectifier parameters if needed
if OPT.parameters.continuation && OPT.stress_needed
    OPT.parameters.minB = OPT.parameters.aggregation_parameter_init;
    OPT.parameters.maxB = OPT.parameters.aggregation_parameter;
    OPT.parameters.deltaB = OPT.parameters.aggregation_parameter_delta;
    OPT.parameters.minK = OPT.parameters.rectifier_parameter_init;
    OPT.parameters.maxK = OPT.parameters.rectifier_parameter;
    OPT.parameters.deltaK = OPT.parameters.rectifier_parameter_delta;
    OPT.parameters.aggregation_parameter = OPT.parameters.minB;
    OPT.parameters.rectifier_parameter = OPT.parameters.minK;
end

    


