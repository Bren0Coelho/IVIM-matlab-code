function  [D_diff_est, D_perf_est, f_est] = NLLS_TRR(M_noise, b, x0, D_diff_est_lls, f_est_lls)

%% Nonlinear Least Square - Trust Region Reflective
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   11/18/2019: ver 1.0
%   05/21/2020: ver 1.1 - adjusted for compatiility with Breno's code.
%   
% ========================================================================
% ========================================================================
% Non-linear least square with biexponential model (y = alpha*exp(x*w) + beta*exp(x*z))
% INPUT: M_noise = simulated corrupted signal (array of double)
%        b = values of b (array of double)
%        X0 = initial guess
%        D_diff_est_lls = diffusion coefficient from LLS fitting
%        f_est_lls = perfusion fraction from LLS fitting
% OUTPUT: f_est_nlls = perfusion fraction
%         D_diff_est_nlls = diffusion coefficient [mm²/s]
%         D_perf_est_nlls = pseudodiffusion coefficient [mm²/s]

    %% Definition of a function to evaluate square error summation
    if (nargin == 5)         
        fun = @(p, b) f_est_lls * exp(-b * p(1)) + ...
            (1 - f_est_lls) * exp(-b * D_diff_est_lls);      
    else
        fun = @(p, b) p(1) * exp(-b * p(2)) + ...
            (1 - p(1)) * exp(-b * p(3));
        
        M_diff_noise = M_noise(b > 200);
                
        y1 = log(M_diff_noise)';
        x1 = b(b > 200)';
        
        P_diff = regress(y1, [ones(size(x1)) x1], 1);

        f_est = 1 - exp(P_diff(1));
        D_diff_est = -P_diff(2);

        if (f_est < 1 && f_est > 0 && D_diff_est > 0 && D_diff_est < 1)            
            x0 =  [f_est, 10*D_diff_est, D_diff_est] ;            
        end           
    end
    
    % Fitting
    options = optimoptions('lsqcurvefit', 'Algorithm', 'trust-region-reflective',...
        'Display', 'off', 'MaxFunEvals', 1e3,...
        'MaxIter', 1e3, 'TolFun', 1e-10,...
        'TolX', 1e-6);
     
    result = lsqcurvefit(fun, x0, b, M_noise, [], [], options);
    
    if (length(result) > 1)
        D_perf_est = result(2);
        f_est = result(1);
        D_diff_est = result(3);
        
        if D_diff_est < 0 || D_diff_est > 5e-2
            D_diff_est = 0;
            f_est = 0;
            D_perf_est = 0;
            return
        elseif D_perf_est < 0 || D_perf_est > 5e-1
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
    else
        D_perf_est = result;
        f_est = [];
        D_diff_est = [];
        
        if D_perf_est < 0 || D_perf_est > 5e-1
            D_diff_est = 0;
            return
        end
    end
end