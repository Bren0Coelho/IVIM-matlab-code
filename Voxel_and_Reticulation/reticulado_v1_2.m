  function [Voxel_noise, Reticulation_data] = reticulado_v1_2(M, SNR, Nc, Nr, d_t, b, ROI, odn, method, x0, fit, weight, tecido)    
%% Reticulado
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   09/11/2020: ver 1.0
%   10/08/2020: ver 1.1 - include parallel distribution of tasks with
%   creteJob and createTask
%   14/08/2020: ver 1.2 - include parallel distribution of tasks with
%   parfor in the number of ROIs
% ========================================================================
% ========================================================================
% This function simulates a reticulation of dimension Nr x Nc x nbdist of
% Voxels with signal M, adds noise to it and sends it to RegionOI function
% INPUT: M = simulated clean signal (array of double)
%        SNR = signal-to-noise ratio (array)
%        Nc = reticulation columns number
%        Nr = reticulation rows number
%        d_t = array of tissue dimensions [columns, rows] (double)
%        b = values of b (array of double)
%        ROI = ROI dimension (double)
%        odn = output directory 
%        method = method used in NNLS. It can be 'Tikh', 'tsvd', 'dsvd' or 'mtsvd' (char)             
%        x0 = [f, D, D*] initial guess (double)
%        fit = fitting algorithm (char)('nnls', 'nlls_lm', 'nlls_trr', 'lls', 'llsr', 'nlls2')
%        weight = weight for robust LLS (char)
%        tecido = 2D matrix with tissue [f, D, D*] (double)
%
% OUTPUT: Reticulation_data = struct with simulated signals and estimated parameters (cell).
%         Voxel_noise = 4D matrix with all reticulations for each b value
%         and each SNR
% ========================================================================
% ========================================================================

    nSNR = length(SNR);                             % number of SNR
    nbdist = length(b);                             % number of b-values
    nROI = length(ROI);                             % number of ROIs
    [ntissue,~] = size(tecido);                     % number of tissues
    Voxel = zeros(Nr, Nc, nbdist);                  % reticulation of voxels
    Voxel_noise = zeros(Nr, Nc, nbdist, nSNR);      % noisy reticulation of voxels
    Reticulation_data = cell(nROI, nSNR);           % simulation results are put here inside
    
    % Building up a reticulation
            
    t_ci = 1;
    t_ri = 1;        
    for i = 1:ntissue            
        
        S = M(i,:);
        t_cf = t_ci + d_t(i,1) - 1;            
        
        for c = t_ci:t_cf
            
            t_rf = d_t(i,2);
            
            for r = t_ri:t_rf
                
%                 Sc = complex(S, zeros(size(S)));
                Voxel(r, c, :) = S;                    
                
            end            
        end
        t_ci = t_cf + 1;            
    end
    
    shift_c = floor((Nc - sum(d_t(:,1)))/2);
    shift_r = floor((Nr - min(d_t(:,2)))/2);
    Voxel = circshift(Voxel, [shift_r shift_c]); % centers the tissues
    Obj = (Voxel(:,:,1) > 0);                    % tracking area

    clust = parcluster();                        % cluster object
    clust.NumWorkers = nROI;                     % number of workers to receive a task
    parpool(nROI);                               % number of parallel pools
    
    for j = 1:nSNR    
        
        % Adding noise
        for ibdist = 1:nbdist
%             Voxel_noise(:, :, ibdist, j) = abs(awgn(Voxel(:,:,ibdist), SNR(j), 'measured'));
            Voxel_noise(:, :, ibdist, j) = awgn(Voxel(:,:,ibdist), SNR(j), 'measured');
        end
        snr = SNR(j);
        Signal = Voxel_noise(:, :, :, j);
        tStart = tic;
        
        % Drawing ROI
        parfor iroi = 1:nROI

            disp(['Simulation of SNR = ', num2str(snr), ' and ROI ', num2str(ROI(iroi)), 'x', num2str(ROI(iroi))]);

            [D_diff_fit, D_perf_fit, f_fit] = RegionOI_v1_2(b, snr, ROI(iroi), Signal, odn, fit, method, weight, x0, Obj);
            
            Reticulation_data{iroi,j} = struct('Ddiff_fit', D_diff_fit, 'Dperf_fit', D_perf_fit, 'f_fit', f_fit);                        
        end
        disp(['Time of tracking: ', num2str(toc(tStart)), 's'])
        disp('----------------------------------')
    end
    delete(gcp);
 end