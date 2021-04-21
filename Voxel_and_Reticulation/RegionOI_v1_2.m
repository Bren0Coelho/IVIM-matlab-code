function [D_diff_fit, D_perf_fit, f_fit] = RegionOI_v1_2(b, SNR, ROI, Voxel_noise, ~, fit, method, weight, x0, Obj)
%% ROI
% =======================================================================
% =======================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   09/23/2020: ver 1.0
%   10/07/2020: ver 1.1 - new ROI's initial position; three parametric maps
%   in one figure; mean value and standard deviation estimation; new step
%   of ROI's motion; outlier treatment
%   10/14/2020: ver 1.2 - New ROI tracking with no differentiation by
%   tissue
% ========================================================================
% ========================================================================
% This function creates a set of ROIs to analyze reticulation information,
% creates a mask based on the ROIs with coordinates of reticulation to 
% assess IVIM parameters and calculate them
% INPUT: b = b-values array (double)
%        SNR = momentarily SNR value(double)
%        ROI = ROI dimension (double)
%        Voxel_noise = reticulation noisy data (3D double matrix)
%        odn = output directory 
%        fit = fitting algorithm (char)('nnls', 'nlls_lm', 'nlls_trr', 'lls', 'llsr', 'nlls2')
%        method = method used in NNLS. It can be 'Tikh', 'tsvd', 'dsvd' or 'mtsvd' (char)
%        weight = weight for robust LLS (char)
%        x0 = [f, D*, D] initial guess (double)
%        Obj = 2D logical reticulation mask -- 1 when in the tissue, 0 when out
%
% OUTPUT: D_diff_fit = parametric map of estimated D (double)
%         D_perf_fit = parametric map of estimated D* (double)
%         f_fit = parametric map of estimated f (double)
% =======================================================================
% =======================================================================

    MR_b0 = Voxel_noise(:,:,1);                                                         % image for b = 0 s/mm² 
    D_diff_fit = zeros(size(MR_b0));                                                    % parametric map of D
    D_perf_fit = zeros(size(MR_b0));                                                    % parametric map of D*
    f_fit = zeros(size(MR_b0));                                                         % parametric map of f
    text1 = ['Acquisition for b = ', num2str(b(1)),...
        ' s/mm^{2}, SNR = ', num2str(SNR), ' and ROI dimension of ', num2str(ROI)];     % title of image
    H = figure;           
    MR_gray = mat2gray(MR_b0);
    MRI = image(MR_gray, 'CDataMapping','scaled');                                      % image
    colormap(gray);
    step = ROI;                                                                         % tracking step
    title(text1);
    set(gcf, 'Position', get(0, 'Screensize'));
    [R,C] = find(Obj);                                                                  % find tracking area
    C = unique(C);
    ci = C(1); cf = C(end);                                                             % initial and final columns of tracking area
    ri = R(1);                                                                          % initial row of tracking area
    nbdist = length(b);                                                                 % number of b-values
    
    % Tracking
    for c = ci:step:cf
        
        % Initial and final position of ROI
        roi_ci = c; roi_cf = c + ROI - 1;
        roi_ri = ri; roi_rf = roi_ri + ROI - 1;
        
        while true          
                       
            M_noise = zeros(1,nbdist);
            roi = imrect(gca, [roi_ci roi_ri ROI-1 ROI-1]); %ROI           
            set(H, 'PaperUnits', 'centimeters');
            set(H, 'PaperPosition', [0 0 60 30]);
            set(gca, 'FontSize', 12);
%             ofn = fullfile(odn, ['Reticulado_SNR' num2str(SNR) '_ROI' num2str(ROI) '_c' num2str(c) '_r' num2str(roi_ri) '.jpg']);
%             saveas(H, ofn);
            BW = createMask(roi, MRI);                   
            BW(roi_ri:roi_rf,roi_ci:roi_cf) = 1;
            [s,q] = find(BW & Obj);
            delete(roi);
                   
            % Make IVIM signal
            for z = 1:nbdist
                M_noise(z) = sum(sum(Voxel_noise(roi_ri:roi_rf,roi_ci:roi_cf,z)));
            end

            % Fitting
            [D_diff_fit(s(1):s(end), q(1):q(end)), D_perf_fit(s(1):s(end), q(1):q(end)), f_fit(s(1):s(end), q(1):q(end))] =...
                fitting(fit, method, M_noise./(ROI*ROI), b, weight, x0); 

            % outliers treatment
            idx = (D_diff_fit > 5e-2 | D_diff_fit < 0); D_diff_fit(idx) = 0;
            idx = (D_perf_fit > 5e-1 | D_perf_fit < 0); D_perf_fit(idx) = 0;
            idx = (f_fit > 1 | f_fit < 0); f_fit(idx) = 0;
            
            
            % Boundaries treatment
            if (~all(sum(Obj(roi_ri:roi_rf+1, roi_ci:roi_cf),2)))
                break;
            end
            
            roi_ri = roi_rf + 1;
            roi_rf = roi_ri + step - 1;
            
        end
    end    
    close(H);
end