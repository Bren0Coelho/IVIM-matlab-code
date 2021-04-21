function [D_diff_est, D_perf_est, f_est] = NLLS2(M_noise, b, x0)

%% Nonlinear Segmented Least Square 
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.
%   06/17/2020: ver 1.2 - comments
% ========================================================================
% ========================================================================
% Non-linear regression of mono-exponential dacay (y = alpha*exp(x*w)) and 
%posterior biexponential decay (y = alpha*exp(x*w) + beta*exp(x*z))
% INPUT: M_noise = simulated corrupted signal (array of double)
%        b = values of b (array of double)
%        x0 = initial guess
% OUTPUT: f_est = perfusion fraction
%         D_diff_est = diffusion coefficient [mm²/s]
%         D_perf_est = pseudodiffusion coefficient [mm²/s]




    % Initial parameters and function definition
    fun = @(p, b) p(1) * exp(-b * p(2)) + p(3);

    M_diff_noise = M_noise(b > 200);

    % Mono-exponencial fitting
    y1 = log(M_diff_noise)';
    x1 = b(b > 200)';

    P_diff = regress(y1, [ones(size(x1)) x1], 1);

    f_est = 1 - exp(P_diff(1));
    D_diff_est = -P_diff(2);

    if (f_est < 1 && f_est > 0 && D_diff_est > 0 && D_diff_est < 1)            
        x0 =  [f_est, 10*D_diff_est, D_diff_est] ;            
    end

    y1 = M_diff_noise';

%     options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt',...
%         'Display', 'off', 'MaxFunctionEvaluations', 1e3,...
%         'MaxIterations', 1e3, 'FunctionTolerance', 1e-10,...
%         'StepTolerance', 1e-6, 'UseParallel', 1);
    
    options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt',...
        'Display', 'off', 'MaxFunEvals', 1e3,...
        'MaxIter', 1e3, 'TolFun', 1e-10,...
        'TolX', 1e-6);

    result = lsqcurvefit(fun, x0, x1, y1, [], [], options);

    % fitcurve = fun(result, b);
    % 
    % f1 = figure;
    % scatter(b, M_noise)
    % hold on
    % plot(b, fitcurve, 'r', 'LineWidth', 2)
    % hold off
    % close(f1)

    f_est = 1 - result(1);
    D_diff_est = result(2);

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

    % Bi-exponential fitting
    [~, D_perf_est, ~] = NLLS_TRR(M_noise, b, D_diff_est * 10, D_diff_est, f_est);
end
