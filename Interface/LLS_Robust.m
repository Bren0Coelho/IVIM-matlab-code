function [D_diff_est, D_perf_est, f_est] = LLS_Robust(M_noise, b, weight)

%% Robust (Weighted) Linear Least Square.
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.   
% ========================================================================
% ========================================================================
% Robust Linear Regression (y = alpha*x + beta)
% INPUT: M_noise = simulated corrupted signal (array of double)
%        b = values of b (array of double)
%        weight = weight for LLS fitting
% OUTPUT: f_est = perfusion fraction
%         D_diff_est = diffusion coefficient [mm²/s]
%         D_perf_est = pseudodiffusion coefficient [mm²/s]
% - M = (1-f)exp(-b*D_diff) -> log(M/M_0) = log((1-f)exp(-b*D_diff))
% - log(M) = log(1-f) -b*D_diff
% Faltava o tratamento de outlier do D_perf

    M_diff_noise = M_noise(b > 200);

    y1 = log(M_diff_noise)';
    x1 = b(b > 200)';

    'stats:statrobustfit:IterationLimit';

    % Linear fitting
    result = robustfit(x1, y1, weight);

    f_est = 1 - exp(result(1));
    D_diff_est = -result(2);

    if D_diff_est < 0 || D_diff_est > 5e-2
        D_diff_est = 0;
        f_est = 0;
        D_perf_est = 0;
        return
    elseif f_est > 1 || f_est < 0
        D_diff_est = 0;
        f_est = 0;
        D_perf_est = 0;
        return
    end
  
    % Non-linear fitting 
    [~, D_perf_est, ~] = NLLS_TRR(M_noise, b, D_diff_est * 10, D_diff_est, f_est);
end
