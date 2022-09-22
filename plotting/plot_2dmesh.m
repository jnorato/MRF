function plot_2dmesh(fig)
%
% Plot a 2d mesh
%
global FE 

F = FE.elem_node.'; % matrix of faces to be sent to patch function
V = FE.coords'; % vertex list to be sent to patch function

figure(fig); cla; hold on   
p = patch('Faces',F,'Vertices',V);
p.FaceColor = 'w';
p.EdgeColor = 'k';
axis equal;
    
end
