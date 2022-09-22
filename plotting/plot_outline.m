function plot_outline(fig)
% This function plots the outline of the Lbracket in 2D

global FE 

figure(fig);
mesh_type_known = true;

switch FE.mesh_input.type
    case '2DLbracket'
        a = FE.mesh_input.L_side;
        b = FE.mesh_input.L_cutout;
        x = [0 a a (a-b) (a-b) 0 0];
        y = [0 0 (a-b) (a-b) a a 0];
    case 'double2DLb'
        a = FE.mesh_input.L_width;
        b = FE.mesh_input.L_height;
        c = FE.mesh_input.L_cutout;
        x = [0 a a (a-c) (a-c) c c 0 0];
        y = [0 0 (b-c) (b-c) b b (b-c) (b-c) 0];   
    case 'generate'
        if FE.dim == 2
            a = FE.mesh_input.box_dimensions(1);
            b = FE.mesh_input.box_dimensions(2);
            x = [0 a a 0 0];
            y = [0 0 b b 0];
        end
    case 'V-frame'
        W = FE.mesh_input.box_dimensions(1);
        H = FE.mesh_input.box_dimensions(2);
        VW = FE.mesh_input.cutout_dimensions(1);
        VH = FE.mesh_input.cutout_dimensions(2);
        d = (W-VW)/2;
        x = [0 d W/2 (W-d) W W 0 0];
        y = [0 0 VH 0 0 H H 0];
    otherwise
        mesh_type_known = false;
        disp('Ignoring plot outline--unrecognized shape.');
end

if mesh_type_known
    line(x,y,'Color', 'k');
end