function plotstress2d(fig, cap, load, meshflag)
%
% Plot stresses for 2d problems
%
% fig       is the figure number to be created
% cap       is the largest value to cap the color scale (values above this will
%           be plotted gray.
% load      is the load case for the stresses
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


% Plot stresses
newmap = jet;
newmap(256,:) = [0.5 0.5 0.5]; % Gray for overflow
colormap(newmap);
s = FE.svm(:,load);
p.CData = s;
p.FaceColor = 'flat';
caxis([0 cap]);
c = colorbar;
% Add the tick for the cap in the colorbar
t = c.Ticks;
if t(length(t)) < cap
    c.Ticks = [t cap];
end
c.FontSize = 14;
axis off;






