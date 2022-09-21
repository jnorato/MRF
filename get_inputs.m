function get_inputs ()
%
%% This script has a single line, in which the user must specify the 
%% location of the Matlab master input script.
%%

global OPT

% Comment out everything and uncomment the line corresponding to the
% problem you want to run.

%% Example 1 -- L-bracket
% ======================
%
% MRF
% ===
OPT.input_file_name = 'input_files/Lbracket2d/MRF/h2_B12_LKS_r25/inputs_Lbracket2d.m';  % h_e = 2, beta_max = 12
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h2_B16_LKS_r25/inputs_Lbracket2d.m';  % h_e = 2, beta_max = 16
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h2_B20_LKS_r25/inputs_Lbracket2d.m';  % h_e = 2, beta_max = 20
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h1_B12_LKS_r25/inputs_Lbracket2d.m';  % h_e = 1, beta_max = 12
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h1_B16_LKS_r25/inputs_Lbracket2d.m';  % h_e = 1, beta_max = 16
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h1_B20_LKS_r25/inputs_Lbracket2d.m';  % h_e = 1, beta_max = 20
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h05_B12_LKS_r25/inputs_Lbracket2d.m';  % h_e = 0.05, beta_max = 12
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h05_B16_LKS_r25/inputs_Lbracket2d.m';  % h_e = 0.05, beta_max = 16
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h05_B20_LKS_r25/inputs_Lbracket2d.m';  % h_e = 0.05, beta_max = 20
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h016_B16_LKS_r4/inputs_Lbracket2d.m';  % h_e = 0.016, beta_max = 16
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h020_B16_LKS_r3/inputs_Lbracket2d.m';  % h_e = 0.020, beta_max = 16
% OPT.input_file_name = 'input_files/Lbracket2d/MRF/h025_B16_LKS_r3/inputs_Lbracket2d.m';  % h_e = 0.025, beta_max = 16
%
% ACS
% ===
% OPT.input_file_name =
% 'input_files/Lbracket2d/ACS/h2_P6_r25/inputs_Lbracket2d.m'; % h_e = 2, P = 6
% 'input_files/Lbracket2d/ACS/h2_P8_r25/inputs_Lbracket2d.m'; % h_e = 2, P = 8
% 'input_files/Lbracket2d/ACS/h2_P10_r25/inputs_Lbracket2d.m'; % h_e = 2, P = 10
% 'input_files/Lbracket2d/ACS/h1_P6_r25/inputs_Lbracket2d.m'; % h_e = 1, P = 6
% 'input_files/Lbracket2d/ACS/h1_P8_r25/inputs_Lbracket2d.m'; % h_e = 1, P = 8
% 'input_files/Lbracket2d/ACS/h1_P10_r25/inputs_Lbracket2d.m'; % h_e = 1, P = 10
% 'input_files/Lbracket2d/ACS/h05_P6_r25/inputs_Lbracket2d.m'; % h_e = 0.5, P = 6
% 'input_files/Lbracket2d/ACS/h05_P8_r25/inputs_Lbracket2d.m'; % h_e = 0.5, P = 8
% 'input_files/Lbracket2d/ACS/h05_P10_r25/inputs_Lbracket2d.m'; % h_e = 0.5, P = 10


%% Example 2 -- Portal frame
% ==========================
%
% OPT.input_file_name = 'input_files/Vframe2d/MRF/inputs_Vframe2d.m'; % MRF
% OPT.input_file_name = 'input_files/Vframe2d/ACS/inputs_Vframe2d.m'; % ACS
% OPT.input_file_name = 'input_files/Vframe2d/COMP/inputs_Vframe2d.m'; % Compliance

%% Example 3 -- Cracked plate
% ===========================
%
% OPT.input_file_name = 'input_files/cracked_plate2d/MRF/inputs_cracked_plate2d.m'; MRF
% OPT.input_file_name = 'input_files/cracked_plate2d/ACS/inputs_cracked_plate2d.m'; ACS

%% Example 4 -- Double L-bracket
% ==============================
%
% OPT.input_file_name = 'input_files/doubleLbracket2d/MRF/inputs_doubleLbracket2d.m'; % MRF
% OPT.input_file_name = 'input_files/doubleLbracket2d/ACS/inputs_doubleLbracket2d.m'; % ACS

%% Example 5 -- 3D-cantilever
% ===========================
%
% OPT.input_file_name = 'input_files/cantilever3d/MRF/inputs_cantilever3d.m'; % MRF
% OPT.input_file_name = 'input_files/cantilever3d/ACS/inputs_cantilever3d.m'; % ACS
% OPT.input_file_name = 'input_files/cantilever3d/COMP/inputs_cantilever3d.m'; % COMP


% ==========================
% Do not uncomment this line
run(OPT.input_file_name);
 
