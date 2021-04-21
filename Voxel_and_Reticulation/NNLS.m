function [D_diff_est, D_perf_est, f_est] = NNLS(M_noise, b, method)

%% Non-negative Least Square 
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.
% ========================================================================
% ========================================================================
% Finds a Dx1 vector x which minimizes ||(A; u*H)*x - (M_noise; u*H*x0)||
% INPUT: M_noise = simulated corrupted signal (array of double)
%        b = values of b (array of double)
%        method = method used in NNLS. It can be 'Tikh', 'tsvd', 'dsvd' or 'mtsvd' (char)
% OUTPUT: f = perfusion fraction
%         D_diff = diffusion coefficient [mm²/s]
%         D_perf = pseudodiffusion coefficient [mm²/s]
% - x stores the amplitudes of each exponential of A;
% - It assumes that we are dealing with a ill-posed problem and utilizes
% a regularization tool, the L-curve, from hich we take the regularization
% parameter u.

    %% General parameters
    
    M = 500;            % number of allowed diffusion components
    

    %% Calculate L-curve
    % Ref. Hansen 2005, The L-curve and its use in the numerical treatment
    % of inverse problems
    
    lj = 10 .^ (linspace(-5, 0, M));        % allowed diffusion components

    argexp = -b' * lj;
    A = exp(argexp);

    [U, sigma, ~] = csvd(A);
    [lambda_corner, ~, ~, ~] = l_curve(U, sigma, M_noise', method);

    x0 = zeros(M, 1);       % a priori estimate of x which is set to zero
                            % when no a priori information is available
    d = [M_noise'; x0];     % problem of the form min_x ||C x - d||^2_2

    Ident = eye(M, M);

    %% NNLS Fitting
    
    C = [A; lambda_corner * Ident];

    sk = lsqnonneg(C, d);

    [amps, locs, widths, proms] = findpeaks(sk, lj, ...
        'Annotate', 'extents', 'WidthReference', 'halfheight');

    j = 1;

    if length(amps) > 2
        tmp = [amps, locs', widths', proms];
        for i = 1 : 2        
            [~, idx] = max(tmp);
            pks(j, :) = tmp(idx(1), :);
            tmp(idx(1), :) = -Inf;
            j = j + 1;
        end
    else
        pks = [amps, locs', widths', proms];
    end

    if length(amps) == 1
        D_diff_est = locs;
        D_perf_est = 0;
        f_est = 0;
    else
        Dcoefs = pks(:, 2);
        D_diff_est = min(Dcoefs);
        D_perf_est = max(Dcoefs);

        tmp = find(sk > 0);
        idx = find(diff(tmp) > 1);
        if isempty(idx)
            f_est = widths(1) / sum(widths);
        else
            tmp = find(sk > 0);
            idx = find(diff(tmp) > 1);
            f_est = sum(sk(tmp(idx + 1) : tmp(end)));
        end
    end
    
    if D_diff_est < 0 || D_diff_est > 5e-2 || D_perf_est < 0 || D_perf_est > 5e-1
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

    % % Plot NNLS results
    % 
    % fit = A * sk;
    % 
    % figure
    % subplot(2, 1, 1);
    % plot(b, M_noise, 'o');
    % hold on
    % plot(b, fit, 'LineWidth', 2, 'MarkerSize', 15);
    % title('IVIM dataset');
    % xlabel('b-value (s/mm^2)');
    % ylabel('Signal Intensity (a.u.)');
    % hold off
    % 
    % subplot(2,1,2)
    % semilogx(1000 * lj, sk, 'LineWidth', 2)
    % hold on
    % title('D Spectrum');
    % xlabel('Diffusion Coefficient (10^{-3} mm^2/s)');
    % ylabel('Amplitude (a.u.)');
    % hold off
end