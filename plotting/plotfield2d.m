function plotfield2d(fig, x, cap_up, cap_low, meshflag)
%
% Plot sensitivities for 2d problems
%
% fig       is the figure number to be created
% x         is the field (with as many components as
%           elements in the mesh)
% cap_up, cap_low
%           are the largest and smallest values to cap the color scale 
%           (values above/below these values will be plotted gray/white).
% meshflag  is true if element edges are to be plotted


global FE 

F = FE.elem_node'; % matrix of faces to be sent to patch function
V = FE.coords'; % vertex list to be sent to patch function

figure(fig); cla; hold on;

% Plot mesh if requested
p = patch('Faces',F,'Vertices',V);
p.FaceColor = 'none';
if meshflag
    p.EdgeColor = 'k';
else
    p.EdgeColor = 'none';
end
axis equal;

% Plot field
newmap = jet;
newmap(256,:) = [0.5 0.5 0.5]; % Gray for overflow
newmap(1,:) = [1 1 1]; % White for underflow
colormap(newmap);
p.CData = x;
p.FaceColor = 'flat';
caxis([cap_low cap_up]);
c = colorbar;
% Add the ticks for the caps in the colorbar
t = c.Ticks;
if t(length(t)) < cap_up
    c.Ticks = [t cap_up];
end
t = c.Ticks;
if t(1) < cap_low
    c.Ticks = [cap_low t];
end
c.FontSize = 14;
axis off;
