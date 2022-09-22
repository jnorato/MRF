function plot_density(fig)
%
% Create figure with plot of densities
%
if nargin < 1
    f = gcf; fig = f.Number;
end
%PLOT_DENSITY plots the density field into the specified figure
%
global FE OPT

% Use this commented code for a faster plot in the case of 'generate' and
% 'read-home-made' meshes.
%
% switch FE.mesh_input.type
%     case {'read-gmsh', '2DLbracket', 'double2DLb', 'V-frame'}
%         plot_density_cells(fig);
%     case {'generate', 'read-home-made'}
%         plot_density_levelsets(fig);
%     otherwise
%         warning('Unrecognized mesh type.');
% end

plot_density_cells(fig);


% title_string = sprintf('density, %s = %f',...
%     OPT.functions.objective,OPT.functions.f{1}.value);
% title(title_string);
% xlabel('x'); ylabel('y'); zlabel('z');
axis off;
title '';

