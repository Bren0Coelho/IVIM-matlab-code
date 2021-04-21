function [lls, llsr, nlls2, trr, lev, nnls] = voxel(M, SNR, b, method, num, x0)

%% Voxel
% [lls, llsr, nlls2, trr, lev, nnls] = voxel(M, SNR, b, method, num, x0)
% ==============================================================================
% ==============================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.
%   16/06/2020: ver 1.2 - comments
% ===============================================================================
% ===============================================================================
% Simulates corrupted DWI signals of a voxel
% INPUT: M = simulated clean signal (array of double)
%        SNR = signal-to-noise ratio (array)
%        b = values of b (array of double)
%        method = method used in NNLS. It can be 'Tikh', 'tsvd', 'dsvd' or 'mtsvd' (char)
%        num = number of iterations
%        X0 = initial guess
% OUTPUT: for each fitting method, it returns a struct with 'num' parameters 
%         estimation, average, standard deviation and corrupted signal
%         lls = struct with Linear Least Square data
%         llsr = struct with Robust Linear Least Square data
%         nlls2 = struct with segmented Non-lnear Least square data
%         nlls = struct with Trust-Region-Reflective data
%         lev = struct with Levenberg-Marquadt data
%         nnls = struct with Non-negative Least Square data
% - Every fitting method takes part in the simulation of the same corrupted signal
% - It simulates for a bunch of SNR values

    %% Final results matrixes
    
    % Average and standard deviation
    fit_mean_lls = zeros(length(SNR), 3);
    fit_std_lls = zeros(length(SNR), 3);
    
    fit_mean_llsr = zeros(length(SNR), 3);
    fit_std_llsr = zeros(length(SNR), 3);
     
    fit_mean_nlls2 = zeros(length(SNR), 3);
    fit_std_nlls2 = zeros(length(SNR), 3);
    
    fit_mean_trr = zeros(length(SNR), 3);
    fit_std_trr = zeros(length(SNR), 3);
     
    fit_mean_lev = zeros(length(SNR), 3);
    fit_std_lev = zeros(length(SNR), 3);
    
    fit_mean_nnls = zeros(length(SNR), 3);
    fit_std_nnls = zeros(length(SNR), 3);
    
    % Corrupted signal
    S = zeros(num, length(b), length(SNR));
        
    % Estimation parameters of all fitting methods for all iterations 
    fit_lls = zeros(num, 3, length(SNR));
    fit_llsr = zeros(num, 3, length(SNR));
    fit_nlls2 = zeros(num, 3, length(SNR));
    fit_trr = zeros(num, 3, length(SNR));
    fit_lev = zeros(num, 3, length(SNR));
    fit_nnls = zeros(num, 3, length(SNR));    
    
    %% For each SNR value, it runs 'num' times
    %parpool(2); %parallel pool
    
    for i = 1 : length(SNR)

        snr = SNR(i);
        fprintf('         ... SNR %d\n', snr);        
        
        tRepStart = tic;
        parfor j = 1 : num
        %for j = 1 : num
            
            %% Adds white Gaussian noise on the signal and stores it in S
            % - j: iteration (1 to 1000)
            % - i: SNR index(1 to 7)
            Mc = complex(M, zeros(size(M)));            
            M_noise = abs(awgn(Mc, snr, 'measured'));            
            S(j, :, i) = M_noise;

            %% Fitting
            
            [D_diff_fit_lls, D_perf_fit_lls, f_fit_lls] = LLS(M_noise, b);
            [D_diff_fit_llsr, D_perf_fit_llsr, f_fit_llsr] = LLS_Robust(M_noise, b, 'talwar');
            [D_diff_fit_nlls2, D_perf_fit_nlls2, f_fit_nlls2] = NLLS2(M_noise, b, x0);
            [D_diff_fit_nlls, D_perf_fit_nlls, f_fit_nlls] = NLLS_TRR(M_noise, b, x0);
            [D_diff_fit_lev, D_perf_fit_lev, f_fit_lev] = NLLS_LM(M_noise, b, x0);
            [D_diff_fit_nnls, D_perf_fit_nnls, f_fit_nnls] = NNLS(M_noise, b, method);
            

            %% 3D matrix of estimated parameters
            % - rows: iteration (1 to 1000)
            % - colums: estimated parametres
            % - j: SNR index (1 to 7)
            fit_lls(j,:,i) = [D_diff_fit_lls, D_perf_fit_lls, f_fit_lls];
            fit_llsr(j,:,i) = [D_diff_fit_llsr, D_perf_fit_llsr, f_fit_llsr];
            fit_nlls2(j,:,i) = [D_diff_fit_nlls2, D_perf_fit_nlls2, f_fit_nlls2];
            fit_trr(j,:,i) = [D_diff_fit_nlls, D_perf_fit_nlls, f_fit_nlls];
            fit_lev(j,:,i) = [D_diff_fit_lev, D_perf_fit_lev, f_fit_lev];
            fit_nnls(j,:,i) = [D_diff_fit_nnls, D_perf_fit_nnls, f_fit_nnls];

        end
        disp(['            ... repetitions take: ', num2str(toc(tRepStart)), 's'])
        
        %% 2D matrix of average and standard deviation
        % - row: SNR index (1 to 7)
        % - colum: average of 'num' estimated parameters        
        fit_mean_lls(i, 1 : 3) = mean(fit_lls(:, :, i));
        fit_std_lls(i, 1 : 3) = std(fit_lls(:, :, i));
        
        fit_mean_llsr(i, 1 : 3) = mean(fit_llsr(:, :, i));
        fit_std_llsr(i, 1 : 3) = std(fit_llsr(:, :, i));
        
        fit_mean_nlls2(i, 1 : 3) = mean(fit_nlls2(:, :, i));
        fit_std_nlls2(i, 1 : 3) = std(fit_nlls2(:, :, i));
        
        fit_mean_trr(i, 1 : 3) = mean(fit_trr(:, :, i));
        fit_std_trr(i, 1 : 3) = std(fit_trr(:, :, i));
         
        fit_mean_lev(i, 1 : 3) = mean(fit_lev(:, :, i));
        fit_std_lev(i, 1 : 3) = std(fit_lev(:, :, i));
         
        fit_mean_nnls(i, 1 : 3) = mean(fit_nnls(:, :, i));
        fit_std_nnls(i, 1 : 3) = std(fit_nnls(:, :, i));
       
    end
    
    %% Final structures    
    lls = struct('pontos', fit_lls, 'media', fit_mean_lls, 'dp', fit_std_lls, 'Sinal', S);
    llsr = struct('pontos', fit_llsr, 'media', fit_mean_llsr, 'dp', fit_std_llsr, 'Sinal', S);
    nlls2 = struct('pontos', fit_nlls2, 'media', fit_mean_nlls2, 'dp', fit_std_nlls2, 'Sinal', S);
    trr = struct('pontos', fit_trr, 'media', fit_mean_trr, 'dp', fit_std_trr, 'Sinal', S);
    lev = struct('pontos', fit_lev, 'media', fit_mean_lev, 'dp', fit_std_lev, 'Sinal', S);
    nnls = struct('pontos', fit_nnls, 'media', fit_mean_nnls, 'dp', fit_std_nnls, 'Sinal', S);
    
end