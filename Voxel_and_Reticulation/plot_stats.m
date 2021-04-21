%% plot_estat
% Plots statistics results of simulation done in voxel.m:
% - Scatter 3D D_diff por D_perf por f
% - Histogram of the estimated parameters
% - Plots error curve, standard deviation and outliers versus SNR for each
% model method

close all; clc;

%% Input
nest = 3;   % Number of simulated structures with parameters
nbdist = 4; % Number of simulated sequences of b

methods = {'LLS'; 'LLSR'; 'NLLS2'; 'TRR'; 'LEV'; 'NNLS'};
nmethods = length(methods);

%% Loads workspace (variables and simulation results)
[fn, dn] = uigetfile('*.mat', 'Select the workspace file');
odn = uigetdir(dn, 'Select the output folder');

load(fullfile(dn, fn));

%% Scatter and histogram
for ibdist = 1 : nbdist
    
    for imethods = 1 : nmethods        
        
        for iest = 2 : 2
            s = eval(matlab.lang.makeValidName(['s', num2str(iest)]));
            nfs = length(s.f);
            nDdiff = length(s.D_diff);
            nDperf = length(s.D_perf);

            data = eval(matlab.lang.makeValidName(['R_b', num2str(ibdist),...
                '_s', num2str(iest)]));

            b = eval(matlab.lang.makeValidName(['b', num2str(ibdist)]));
            nSNR = length(SNR);
            
            iparams = 1;
            H = figure(2);
            for iDdiff = 1 : nDdiff
                curDdiff = s.D_diff(iDdiff);
                
                for iDperf = 1 : nDperf
                    curDperf = s.D_perf(iDperf);
                    
                    for ifs = 1 : nfs
                        curf = s.f(ifs);
                        
                        erroDdiff = zeros(1, nSNR);
                        erroDperf = zeros(1, nSNR);
                        errof = zeros(1, nSNR);
                        
                        dpDdiff = zeros(1, nSNR);
                        dpDperf = zeros(1, nSNR);
                        dpf = zeros(1, nSNR);
                        
                        outDdiff = zeros(1, nSNR);
                        outDperf = zeros(1, nSNR);
                        outf = zeros(1, nSNR);
                        
                        for iSNR = 1 : nSNR
                            
                            curdata = data(iparams, imethods).pontos;
                            
                            curDdiffdata = data(iparams, imethods).pontos(:, 1, iSNR);
                            curDperfdata = data(iparams, imethods).pontos(:, 2, iSNR);
                            curfdata = data(iparams, imethods).pontos(:, 3, iSNR);
                            
                            vcurDdiffdata = curDdiffdata(curDdiffdata > 0 & curDdiffdata < 5e-2);
                            vcurDperfdata = curDperfdata(curDperfdata > 0 & curDperfdata < 5e-1);
                            vcurfdata = curfdata(curfdata > 0 & curfdata < 1);
                            
                            outDdiff(iSNR) = length(curDdiffdata) - length(vcurDdiffdata);
                            outDperf(iSNR) = length(curDperfdata) - length(vcurDperfdata);
                            outf(iSNR) = length(curfdata) - length(vcurfdata);
            
%                             % Scatter
%                             H = figure(1);
%                             subplot(2, 2, 4)
%                             scatter3(curDdiffdata, curDperfdata, curfdata, 'filled');
%                             title(['Scatter for SNR = ', num2str(SNR(iSNR)), ' - Method: ', methods{imethods}]);
%                             xlabel('D [mm²/s]');
%                             ylabel('D* [mm²/s]');
%                             zlabel('f');
% 
%                             % Histogram
%                                 %D_diff    
%                             subplot(2, 2, 1)
%                             histogram(vcurDdiffdata);
%                             hold on;
%                             line([curDdiff, curDdiff], ylim, 'LineWidth', 2, 'Color', 'r');
%                             line([mean(vcurDdiffdata), mean(vcurDdiffdata)], ylim, 'LineWidth', 2, 'Color', 'g');
%                             title(['Histogram of D for SNR = ', num2str(SNR(iSNR)), ' - ', methods{imethods}]);
%                             xlabel('D [mm²/s]');
%                             ylabel('Occurrence');
% 
%                                 %D_perf
%                             subplot(2, 2, 2)
%                             histogram(vcurDperfdata);
%                             hold on;
%                             line([curDperf, curDperf], ylim, 'LineWidth', 2, 'Color', 'r');
%                             line([mean(vcurDperfdata), mean(vcurDperfdata)], ylim, 'LineWidth', 2, 'Color', 'g');
%                             title(['Histogram of D^{*} for SNR = ', num2str(SNR(iSNR)), ' - ', methods{imethods}]);
%                             xlabel('D* [mm²/s]');
%                             ylabel('Occurrence');
%                             legend('data', 'simulated', 'mean', 'Location', 'best');
% 
%                                 %f
%                             subplot(2, 2, 3)
%                             histogram(vcurfdata);
%                             hold on;
%                             line([curf, curf], ylim, 'LineWidth', 2, 'Color', 'r');
%                             line([mean(vcurfdata), mean(vcurfdata)], ylim, 'LineWidth', 2, 'Color', 'g');
%                             title(['Histogram of f for SNR = ', num2str(SNR(iSNR)), ' - ', methods{imethods}]);
%                             xlabel('f');
%                             ylabel('Occurrence');
%                             
%                             set(gcf, 'Position', get(0, 'Screensize'));
%                             ofn = fullfile(odn, ['s' num2str(iest) '_b' num2str(ibdist)...
%                                 '_' methods{imethods} '_f' num2str(curf)...
%                                 '_Ddiff' num2str(curDdiff) '_Dperf' num2str(curDperf)...
%                                 '_SNR' num2str(SNR(iSNR)) '.jpg']);
%                             saveas(H, ofn, 'jpg')
%                             close(H)
                            
                            % Percentual error and standard deviation                            
                            erroDdiff(iSNR) = abs(mean (100 * (vcurDdiffdata - curDdiff) / curDdiff));
                            erroDperf(iSNR) = abs(mean (100 * (vcurDperfdata - curDperf) / curDperf));
                            errof(iSNR) = abs(mean (100 * (vcurfdata - curf) / curf));
                            
                            dpDdiff(iSNR) = std(vcurDdiffdata);
                            dpDperf(iSNR) = std(vcurDperfdata);
                            dpf(iSNR) = std(vcurfdata);                           
                        end
                        
                        % Graphics
                                                
                            %D_diff            
                        subplot(3, 3, 1);
                        semilogy(SNR, erroDdiff);
                        xlim([SNR(1) SNR(end)]);
                        title(['Relative error % of D - ', methods{imethods}]);     
                        xlabel('SNR'); ylabel('Error %');          
                        set(gca, 'FontSize', 12);
                        grid;
                        hold on

                        subplot(3, 3, 2);
                        plot(SNR, dpDdiff);
                        hold on;
                        line(xlim, [curDdiff, curDdiff], 'LineWidth', 2, 'Color', 'r');
                        xlim([SNR(1) SNR(end)]);
                        title(['SD of D - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Standard Deviation');
                        set(gca, 'FontSize', 12);
                        grid;                       
                        
                        subplot(3, 3, 3);
                        plot(SNR, outDdiff);
                        xlim([SNR(1) SNR(end)]); ylim([0 size(curdata, 1)]);
                        title(['Outliers of D - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Outliers');
                        set(gca, 'FontSize', 12);
                        grid;
                        hold on
                        
                            %D_perf            
                        subplot(3, 3, 4);
                        semilogy(SNR, erroDperf);
                        xlim([SNR(1) SNR(end)]);
                        title(['Relative error % of D^{*} - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Error %');
                        set(gca, 'FontSize', 12);
                        grid;
                        hold on
                        
                        subplot(3, 3, 5);
                        plot(SNR, dpDperf);
                        hold on;
                        line(xlim, [curDperf, curDperf], 'LineWidth', 2, 'Color', 'r');
                        xlim([SNR(1) SNR(end)]);
                        title(['SD of D^{*} - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Standard Deviation');
                        set(gca, 'FontSize', 12);
                        grid;
                        
                        
                        subplot(3, 3, 6);
                        plot(SNR, outDperf);
                        xlim([SNR(1) SNR(end)]); ylim([0 size(curdata, 1)]);
                        title(['Outliers of D^{*} - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Outliers');
                        set(gca, 'FontSize', 12);
                        grid;
                        hold on
            
                            %f
                        subplot(3, 3, 7);
                        semilogy(SNR, errof);
                        xlim([SNR(1) SNR(end)]);
                        title(['Relative error % of f - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Error %');
                        set(gca, 'FontSize', 12);
                        grid;
                        hold on
                        
                        subplot(3, 3, 8);
                        plot(SNR, dpf);
                        hold on;
                        line(xlim, [curf, curf], 'LineWidth', 2, 'Color', 'r');
                        xlim([SNR(1) SNR(end)]);
                        title(['SD of f - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Standard Deviation');
                        set(gca, 'FontSize', 12);
                        grid;                        
                        
                        subplot(3, 3, 9);
                        plot(SNR, outf);
                        xlim([SNR(1) SNR(end)]); ylim([0 size(curdata, 1)]);
                        title(['Outliers of f - ', methods{imethods}]);
                        xlabel('SNR'); ylabel('Outliers');
                        grid;
                        hold on
                        
                        set(H, 'PaperUnits', 'centimeters');
                        set(H, 'PaperPosition', [0 0 60 30]);
                        set(gca, 'FontSize', 12);
                        ofn = fullfile(odn, ['s' num2str(iest) '_b' num2str(ibdist)...
                            '_' methods{imethods} '_f' num2str(curf)...
                            '_Ddiff' num2str(curDdiff) '_Dperf' num2str(curDperf)...
                            '_Erro_DP.jpg']);                       
                        
                    end
                    iparams = iparams + 1;                    
                end
            end            
            saveas(H, ofn);
            close(H)
        end        
    end
end