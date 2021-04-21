function [D_diff_fit, D_perf_fit, f_fit] = fitting(fit, method, M_noise, b, weight, x0)

%% Fitting
% Calls the methods of fitting according to the input parameter fit;
% INPUT: fit = fitting method to be called
%        method = method used in NNLS. It can be 'Tikh', 'tsvd', 'dsvd' or 'mtsvd' (char)
%        M_noise = simulated corrupted signal (array of double)
%        b = values of b (array of double)
%        i = thresholf value to segmented LLS regression
%        weight = weight for LLS fitting
%        x0 = initial guess
%
% OUTPUT: f_fit = perfusion fraction
%         D_diff_fit = diffusion coefficient [mm²/s]
%         D_perf_fit = pseudodiffusion coefficient [mm²/s]


%% Call

    switch fit
        
        case 'lls'
            
            % Linear Least Square
            
            [D_diff_fit, D_perf_fit,f_fit] = LLS(M_noise, b);
            
        case 'llsr'
            
            % Robust Linear Least Square
            
            [D_diff_fit, D_perf_fit,f_fit] = LLS_Robust(M_noise, b, weight);

        case 'nlls_lm'
            
            % Levenberg-Marquadt
            
            [D_diff_fit, D_perf_fit, f_fit] = NLLS_LM(M_noise, b, x0);

        case 'nlls_trr'
            
            % Trust-Region-Reflective
            
            [D_diff_fit, D_perf_fit, f_fit] = NLLS_TRR(M_noise, b, x0);
            
        case 'nlls2'
            
            % Segmented Nonlinear Least Square
            
            [D_diff_fit, D_perf_fit, f_fit] = NLLS2(M_noise, b, x0);
        
        case 'nnls'
            
            % Non-negative Least Square
            
            [D_diff_fit, D_perf_fit, f_fit] = NNLS(M_noise, b, method);
            
        case 'bayes'
            
            % Bayes Inference
            
%             lim = [0 0 0;1 5e-2 5e-1]; % order: f,D,D*
%             n = 1e4;
%             rician = false;
%             prior = {'flat','lognorm','lognorm'};
%             burns = 1;
%             meanonly = true;
%             
%             out =  bayes(M_noise,f,D,Dstar,b,lim,n,rician,prior,burns,meanonly);
%             D_diff_fit = out.D.mean;
%             D_perf_fit = out.Dstar.mean;
%             f_fit = out.f.mean;

              [D_diff_fit, D_perf_fit, f_fit] = bayes2(M_noise, b);
              
        otherwise
            
            error(['fitting method: ', fit, '. There is no such method to simulate. Options: lls, llsr, nlls2, nlls_lm, nlls_trr and nnls!']);          
            
    end
end