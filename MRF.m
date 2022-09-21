% =========================================================================
% 
% MRF
% 
% A Matlab code for stress-constrained topology optimization using the
% maximum rectifier approach.
% Version 1.0 -- June 2021 
%
% Julian Norato and Hollis Smith
% Department of Mechanical Engineering
% University of Connecticut
%
% A significant portion of this code is a derivative of the GPTO code
% written by the authors.
%
%
% Disclaimer
% ==========
% This software is provided by the contributors "as-is" with no explicit or
% implied warranty of any kind. In no event shall the University of
% Connecticut or the contributors be held liable for damages incurred by
% the use of this software.
%
% License
% =======
% This software is released under the Creative Commons CC BY 4.0
% license. As such, you are allowed to copy and redistribute the material 
% in any medium or format, and to remix, transform, and build upon the 
% material, as long as you give appropriate credit, provide a link to the 
% % license, and indicate if changes were made. You may do so in any 
% reasonable manner, but not in any way that suggests the licensor 
% endorses you or your use.
%
% Acknowledgments
% ===============
% This research was supported in part by the Air Force Research Laboratory, 
% Aerospace Systems Directorate, through the Air Force Oﬃce of Scientiﬁc 
% Research Summer Faculty Fellowship Program, Contract Numbers FA8750-15-3-6003 
% and FA9550-15-0001. The ﬁrst author also gratefully acknowledges partial 
% support to conduct this work from the US National Science Foundation, 
% award CMMI-1751211. 
%
% AFRL Distribution Statement A: Approved for Public Release; Distribution 
% is Unlimited. PA# AFRL-2022-3456
%
%
% GCMMA-MMA-code is redistributed under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 3 of 
% the License, or (at your option) any later version.
% =========================================================================


clear all; close all; clc;
%% source folders containing scripts not in this folder
addpath(genpath('FE_routines'))
addpath(genpath('functions'))
addpath(genpath('mesh_utilities'))
addpath(genpath('optimization'))
addpath(genpath('utilities'))
addpath(genpath('plotting'))


global OPT FE

%% Start timer
tic;

%% Start diary to save output on command window to file
diaryname = 'outfile.txt';
diary(diaryname);

%% Initialization
get_inputs();

init_FE();
tic;
init_optimization();
toc;

%% Analysis
perform_analysis(); 

%% Finite difference check of sensitivities
% (If requested)
if OPT.make_fd_check
    run_finite_difference_check();
    return;  % End code here
end

% Create output folder if needed
if OPT.options.save_outputs 
    out_folder = OPT.options.outputs_path;
    % Check if folder exists; if not, create:
    if ~exist(out_folder, 'dir')
           mkdir(out_folder);
    end
end
%% Optimization
OPT.history = runmma(OPT.dv,@(x)obj(x),@(x)nonlcon(x));

%%  Additional postprocessing
% Compute and report compliance regardless of whether or not it is an 
% optimization function.
% This allows us to compare the compliance of stress-constrained designs
[c,~] = compute_compliance();
fprintf('Compliance of final design = %-12.5e\n', c);
OPT.c = c;

%
% Report gray region fraction, which serves as an indication of
% convergence to 0-1 design.
fprintf('Gray region fraction of final design = %-12.5e\n', OPT.grf);

fprintf('Largest true max of relaxed stress of final design = %-12.5e\n', ...
    OPT.true_stress_max);

%  Save data structures to .mat file for future recovery
[folder, baseFileName, ~] = fileparts(OPT.options.mat_output_path);
mat_filename = fullfile(folder, strcat(baseFileName, '.mat'));
save(mat_filename, 'OPT');
    
%% Plot History
if OPT.options.plot == true
    plot_history(2);
end

%% Report time
toc

%% Turn off diary logging
diary off;

% ================================
%% Copy outputs to selected folder
% Update folder name as needed (as well as figure formats)
if OPT.options.save_outputs 
    if OPT.options.plot
        % Density plot
        saveas(1, strcat(out_folder, '/dens.fig'));
        saveas(1, strcat(out_folder, '/dens.pdf'));
        % History plot
        saveas(2, strcat(out_folder, '/hist.fig'));
        saveas(2, strcat(out_folder, '/hist.png'));
        if FE.dim == 2 
            offset_fig_n = 2; 
            for il=1:FE.nloads
                saveas(offset_fig_n + il, strcat(out_folder, '/stress', string(il), '.fig'));
                saveas(offset_fig_n + il, strcat(out_folder, '/stress', string(il), '.pdf'));
            end
        end
    end
    % OPT data structure
    save( strcat(out_folder, '/OPT.mat'), 'OPT');
    % Diary
    copyfile(diaryname, out_folder);
    % Input files
    copyfile(strcat('./', OPT.input_file_name), out_folder);
    copyfile(FE.mesh_input.bcs_file, out_folder);
end
% ================================


 