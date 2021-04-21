%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Breno Spinelli Coelho                                                %
% Adviser: Prof. Dr. Fernando Fernandes Paiva                                  %
% Institute of Physics of São Carlos - University of São Paulo                 %
% CIERMag-In vivo Magnetic Resonance Images and Spectroscopy Center            %
% Applied Physics Master Degree - Computational Physics                        %
% The aim of this script is to simulate MR signal influenced by diffusion and  %
% perfusion and to estimate the parameters f, D and D*                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==============================================================================
% ==============================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.
%   16/06/2020: ver 1.2 - comments, error catching session, 
% ===============================================================================
% ===============================================================================

clear;
clc;
close all;
% delete(gcp);

%% --------------Global Inputs-------------- %%
simula = 'reticulado';                                                  %'voxel' to simulate a voxel
                                                                        %'reticulado' to simulate a reticulation
                                                                        
method  = 'Tikh';                                                       %Tikhonov method for NNLS
weight = 'talwar';                                                      %Weight for robust LLS
SNR = 20:5:50;                                                          %Signal-to-noise ratio
x0 = [0.15 1e-2 1e-3];                                                  %Initial guess
b1 = [0 5 10 15 20 30 35 40 45 50 55 60 70 80 ...                       %b-values
    90 100 150 200 250 300 350 400 450 500 1000];
b2 = [0 10 20 30 40 50 75 100 200 400 600 1000];          
b3 = [0 5 10 15 30 60 120 250 500 1000];
b4 = [0 20 40 80 200 400 700 1000];

%% --------------Voxel Inputs-------------- %%
s1 = struct('f', 0.1:0.05:0.3, 'D_diff', 3e-3, 'D_perf', 4e-2);         %Variation of the perfusion fraction parameter
s2 = struct('f', 0.2, 'D_diff', 1e-3:0.5e-3:7e-3, 'D_perf', 4e-2);      %Variation of the diffusion coefficient [mm²/s]
s3 = struct('f', 0.2, 'D_diff', 3e-3, 'D_perf', 1e-2:0.5e-2:7e-2);      %Variation of the pseudodiffusion coefficient [mm²/s] 
s_sim = [s1 s2 s3];                                                     %Array with the latter structures
num = 1000;                                                             %Number of iterations for running estimations                                     
b_seq = {b1, b2, b3, b4};                                               %b-values array [s/mm²]

%% -----------Reticulation Inputs---------- % %
tecido =  [0.10, 1e-3, 1e-2;                                            %Tissue 1 |
           0.20, 2e-3, 2e-2;                                            %Tissue 2 |-> rows are tissues; columns 1, 2 and 3 are f, D and D* respectively
           0.30, 3e-3, 3e-2];                                           %Tissue 3 |
Nc = 160;                                                               %Reticulation number of columns
Nr = 60;                                                                %Reticulation number of rows
d_t = 10*[5, 5;                                                         %Tissue 1 |
       5, 5;                                                            %Tissue 2 |-> rows are tissues; columns 1 and 2 are c and r tissue dimensions respectively
       5, 5];                                                           %Tissue 3 |
   
ROI = [2, 3, 4, 5];                                                     %ROI dimensions
fit = 'nlls_lm';                                                        %Fitting method (nnls, nlls_lm, nlls_trr, lls, llsr, nlls2)
b = b2;                                                                 %b-values array [s/mm²]

%% Error Handling

met_nnls = struct('Tikh', 1, 'tsvd', 2, 'dsvd', 3, 'mtsvd', 4);
met_fit = struct('lls', 1, 'lls_robust', 2, 'nlls_lm', 3, 'nlls_trr', 4, 'bayes', 5, 'nnls', 6, 'nlls2', 7);

if (~isfield(met_nnls, method))
    error('"method" must be one of these: "Tikh", "tsvd", "dsvd" or "mtsvd" ');
elseif (length(x0) ~= 3)
    error('x0 must have 3 parameter guesses');
elseif (~isfield(met_fit, fit))
    error('"fit" must be one of these: "lls", "lls_robust", "nlls_lm", "nlls_trr", "bayes", "nnls", "nnls" or "nlls2"');
end

%% Simulation
tStart = tic;
switch simula
    
    case 'voxel'
        
        for i = 1 : 3 % Simulates each structure of parameters
            
           
    
            f = s(i).f;
            D_diff = s(i).D_diff;
            D_perf = s(i).D_perf;
            
            fprintf('Simulating structure %d...\n', i);
                
            if (i == 1) % s1

                for j = 1 : length(s1.f)
                    fprintf('   ... valor %d de %d \n', j, length(s1.f));

                    M1 = f(j) * exp(-b1 * D_perf) + (1 - f(j)) * exp(-b1 * D_diff);
                    M2 = f(j) * exp(-b2 * D_perf) + (1 - f(j)) * exp(-b2 * D_diff);
                    M3 = f(j) * exp(-b3 * D_perf) + (1 - f(j)) * exp(-b3 * D_diff);
                    M4 = f(j) * exp(-b4 * D_perf) + (1 - f(j)) * exp(-b4 * D_diff);

                    tbStart = tic;
                    fprintf('      ... distribuição b1\n');
                    [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1] = voxel(M1, SNR, b1, method, num, x0);
                    disp(['      ... b1 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b2\n');
                    [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2] = voxel(M2, SNR, b2, method, num, x0);
                    disp(['      ... b2 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b3\n');
                    [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3] = voxel(M3, SNR, b3, method, num, x0);
                    disp(['      ... b3 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b4\n');
                    [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4] = voxel(M4, SNR, b4, method, num, x0);
                    disp(['      ... b4 levou: ', num2str(toc(tbStart)), 's'])                   

                    R_b1_s1(j, :) = [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1];
                    R_b2_s1(j, :) = [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2];
                    R_b3_s1(j, :) = [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3];
                    R_b4_s1(j, :) = [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4];                       

                    save('workspace_simulation_fernando_1000_v3');                       
                end

            elseif (i == 2) % s2

                for j = 1 : length(s2.D_diff)
                    fprintf('   ... valor %d de %d...\n', j, length(s2.D_diff));

                    M1 = f * exp(-b1 * D_perf) + (1 - f) * exp(-b1 * D_diff(j));
                    M2 = f * exp(-b2 * D_perf) + (1 - f) * exp(-b2 * D_diff(j));
                    M3 = f * exp(-b3 * D_perf) + (1 - f) * exp(-b3 * D_diff(j));
                    M4 = f * exp(-b4 * D_perf) + (1 - f) * exp(-b4 * D_diff(j));
                    
                    tbStart = tic;
                    fprintf('      ... distribuição b1\n');
                    [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1] = voxel(M1, SNR, b1, method, num, x0);
                    disp(['      ... b1 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b2\n');
                    [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2] = voxel(M2, SNR, b2, method, num, x0);
                    disp(['      ... b2 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b3\n');
                    [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3] = voxel(M3, SNR, b3, method, num, x0);
                    disp(['      ... b3 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b4\n');
                    [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4] = voxel(M4, SNR, b4, method, num, x0);
                    disp(['      ... b4 levou: ', num2str(toc(tbStart)), 's'])                    

                    R_b1_s2(j, :) = [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1];
                    R_b2_s2(j, :) = [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2];
                    R_b3_s2(j, :) = [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3];
                    R_b4_s2(j, :) = [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4];                       

                    save('workspace_simulation_fernando_1000_v3');                     
                end

            else % s3

                for j = 1 : length(s3.D_perf)
                    fprintf('   ... valor %d de %d...\n', j, length(s3.D_perf));

                    M1 = f * exp(-b1 * D_perf(j)) + (1 - f) * exp(-b1 * D_diff);
                    M2 = f * exp(-b2 * D_perf(j)) + (1 - f) * exp(-b2 * D_diff);
                    M3 = f * exp(-b3 * D_perf(j)) + (1 - f) * exp(-b3 * D_diff);
                    M4 = f * exp(-b4 * D_perf(j)) + (1 - f) * exp(-b4 * D_diff);

                    tbStart = tic;
                    fprintf('      ... distribuição b1\n');
                    [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1] = voxel(M1, SNR, b1, method, num, x0);
                    disp(['      ... b1 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b2\n');
                    [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2] = voxel(M2, SNR, b2, method, num, x0);
                    disp(['      ... b2 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b3\n');
                    [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3] = voxel(M3, SNR, b3, method, num, x0);
                    disp(['      ... b3 levou: ', num2str(toc(tbStart)), 's'])
                    tbStart = tic;
                    fprintf('      ... distribuição b4\n');
                    [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4] = voxel(M4, SNR, b4, method, num, x0);
                    disp(['      ... b4 levou: ', num2str(toc(tbStart)), 's'])                    

                    R_b1_s3(j, :) = [lls_1, llsr_1, nlls2_1, nlls_1, lev_1, nnls_1];
                    R_b2_s3(j, :) = [lls_2, llsr_2, nlls2_2, nlls_2, lev_2, nnls_2];
                    R_b3_s3(j, :) = [lls_3, llsr_3, nlls2_3, nlls_3, lev_3, nnls_3];
                    R_b4_s3(j, :) = [lls_4, llsr_4, nlls2_4, nlls_4, lev_4, nnls_4];                       

                    save('workspace_simulation_fernando_1000_v3');                   
                end
            end                           
        end      
    
    case 'reticulado'
        
        path = 'C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\Imagens\Reticulado';
        odn = uigetdir(path, 'Select the output folder');
            
        [ntissue,~] = size(tecido);
        M = zeros(ntissue,length(b));
        
        for itissue = 1:ntissue
            f = tecido(itissue,1);
            D_diff = tecido(itissue,2);
            D_perf = tecido(itissue,3);
            M(itissue,:) = f*exp(-b*D_perf)+(1-f)*exp(-b*D_diff);
        end
       
        [Voxel_noise, Reticulation_data] = reticulado_v1_2(M, SNR, Nc, Nr, d_t, b, ROI, odn, method, x0, fit, weight, tecido);
        
        save('workspace_simulation_RETICULADO_v2');        
               
    otherwise 
        error('There is no such simulation option');
end 

warndlg({[simula ' simulation done!'];['Total time of execution: ', num2str(toc(tStart)), 's']}, 'Warning!');
disp(['Total time of execution: ', num2str(toc(tStart)), 's'])