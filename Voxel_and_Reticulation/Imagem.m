%% Imagem
% ==============================================================================
% ==============================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   07/08/2020: ver 1.0
%
% ===============================================================================
% ===============================================================================
% This function plots D, D*, f and fD* estimations as an N x N image
% ===============================================================================
% ===============================================================================



%% Input
clear;
load('workspace_simulation_RETICULADO_v3_NOISE');
close all;
clc;
nSNR = length(SNR);
nROI = length(ROI);
rDim = 4;

avrg_Ddiff = zeros(nROI, ntissue, nSNR);
dp_Ddiff = zeros(nROI, ntissue, nSNR);
avrg_Ddiff2 = zeros(nROI, ntissue, nSNR);
dp_Ddiff2 = zeros(nROI, ntissue, nSNR);

avrg_Dperf = zeros(nROI, ntissue, nSNR);
dp_Dperf = zeros(nROI, ntissue, nSNR);
avrg_Dperf2 = zeros(nROI, ntissue, nSNR);
dp_Dperf2 = zeros(nROI, ntissue, nSNR);

avrg_f = zeros(nROI, ntissue, nSNR);
dp_f = zeros(nROI, ntissue, nSNR);
avrg_f2 = zeros(nROI, ntissue, nSNR);
dp_f2 = zeros(nROI, ntissue, nSNR);

%% Loads workspace (variables and simulation results)

for isnr = 1:1
    
    path = ['C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\Imagens\Reticulado\SNR ', num2str(SNR(isnr))];        
    
    for iroi = 1:nROI            
           
            D_diff_show = Reticulation_data{iroi, isnr}.Ddiff_fit;

            D_perf_show = Reticulation_data{iroi, isnr}.Dperf_fit;

            f_show = Reticulation_data{iroi, isnr}.f_fit;
%            
%             %% Show image
            
            H = figure();
            ax1 = subplot(3,1,1);
            image(f_show,'CDataMapping','scaled');
            c = colorbar;
            c.Label.String = 'f';
            c.Label.FontSize = 14;
            set(ax1, 'XTick', [], 'YTick', []);
            set(gca, 'FontSize', 12);
            title(['Simulated parameters: [f_1,D_1,D^{*}_1] = [', num2str(tecido(1,1)), ',',  num2str(tecido(1,2)), ',', num2str(tecido(1,3)), ']', ...
                '; [f_2,D_2,D^{*}_2] = [', num2str(tecido(2,1)), ',',  num2str(tecido(2,2)), ',', num2str(tecido(2,3)), ']', ...
                '; [f_3,D_3,D^{*}_3] = [', num2str(tecido(3,1)), ',',  num2str(tecido(3,2)), ',', num2str(tecido(3,3)), ']', ...
                ' SNR = ', num2str(SNR(isnr)), ' ROI ', num2str(ROI(iroi)), ...
                'x', num2str(ROI(iroi))])
           
            ax2 = subplot(3,1,2);
            image(D_diff_show,'CDataMapping','scaled');
            c = colorbar;
            c.Label.String = 'D (mm^{2}/s)';
            c.Label.FontSize = 14;
            set(gca, 'FontSize', 12);
            set(ax2, 'XTick', [], 'YTick', []);            

            ax3 = subplot(3,1,3);
            image(D_perf_show,'CDataMapping','scaled');
            c = colorbar;
            c.Label.String = 'D^{*} (mm^{2}/s)';
            c.Label.FontSize = 14;
            set(gca, 'FontSize', 12);
            set(ax3, 'XTick', [], 'YTick', []);            
            
%             set(gcf, 'Position', get(0, 'Screensize'));
%             set(H, 'PaperUnits', 'centimeters');
%             set(H, 'PaperPosition', [0 0 60 30]);                                    
%             ofn = fullfile(path, ['SNR = ',num2str(SNR(isnr)), '_ROI', num2str(ROI(iroi)), 'x', num2str(ROI(iroi)) '_NOISE.jpg']);
%             saveas(H, ofn);            
%             close(H);            

        % Mean and standard deviation
        Obj = (D_diff_show ~= 0);
        [R,C] = find(Obj,1);
        for itissue = 1:ntissue

            avrg_Ddiff(iroi,itissue,isnr) = mean2(D_diff_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));
            dp_Ddiff(iroi,itissue,isnr) = std2(D_diff_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));

            avrg_Dperf(iroi,itissue,isnr) = mean2(D_perf_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));
            dp_Dperf(iroi,itissue,isnr) = std2(D_perf_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));

            avrg_f(iroi,itissue,isnr) = mean2(f_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));
            dp_f(iroi,itissue,isnr) = std2(f_show(R:R+d_t(itissue,2)-1,C:C+d_t(itissue,1)-1));
            
            % Region with no image background influence
            avrg_Ddiff2(iroi,itissue,isnr) = mean2(D_diff_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));
            dp_Ddiff2(iroi,itissue,isnr) = std2(D_diff_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));

            avrg_Dperf2(iroi,itissue,isnr) = mean2(D_perf_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));
            dp_Dperf2(iroi,itissue,isnr) = std2(D_perf_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));

            avrg_f2(iroi,itissue,isnr) = mean2(f_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));
            dp_f2(iroi,itissue,isnr) = std2(f_show(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim));
                     
                        
            I1 = figure(SNR(isnr));
            ax1 = subplot(ntissue,2,2*itissue-1);            
            plot(ROI, 100*abs(avrg_Ddiff(:,itissue,isnr)-tecido(itissue,2))./tecido(itissue,2), 'bo', ...
                ROI, 100*abs(avrg_Ddiff2(:,itissue,isnr)-tecido(itissue,2))./tecido(itissue,2), 'ro', 'Linewidth', 2);            
            grid on;
            ylabel(ax1, ['RE (%) Tissue ', num2str(itissue)], 'FontSize', 12);                                    
            ax2 = subplot(ntissue,2, 2*itissue);
            plot(ROI, dp_Ddiff(:,itissue,isnr), 'bo', ...
                ROI, dp_Ddiff2(:,itissue,isnr), 'ro', 'Linewidth', 2);            
%             line(xlim, [avrg_Ddiff(:,itissue,isnr), avrg_Ddiff(:,itissue,isnr)], 'LineWidth', 2, 'Color', 'r');
            ylabel(ax2, ['SD Tissue ', num2str(itissue)], 'FontSize', 12);                                                
            grid on;
            set(ax1, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            set(ax2, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            if itissue == ntissue                
                xlabel(ax1, 'ROI square dimension');                
                xlabel(ax2, 'ROI square dimension');
            elseif itissue == 1
                legend('Whole', 'Region', 'Location', 'eastoutside');
                title(ax1, ['RE of mean estimated D for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                           
                title(ax2, ['SD of estimated D for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                           
            end            
                        
            I2 = figure(SNR(isnr)+1);
            ax1 = subplot(ntissue,2,2*itissue-1);           
            plot(ROI, 100*abs(avrg_Dperf(:,itissue,isnr)-tecido(itissue,3))./tecido(itissue,3), 'bo', ...
                ROI, 100*abs(avrg_Dperf2(:,itissue,isnr)-tecido(itissue,3))./tecido(itissue,3), 'ro', 'Linewidth', 2);            
            grid on;
            ylabel(ax1, ['RE (%) Tissue ', num2str(itissue)], 'FontSize', 12);  
            ax2 = subplot(ntissue,2, 2*itissue);
            plot(ROI, dp_Dperf(:,itissue,isnr), 'bo', ...
                ROI, dp_Dperf2(:,itissue,isnr), 'ro', 'Linewidth', 2);            
%             line(xlim, [avrg_Dperf(:,itissue,isnr), avrg_Dperf(:,itissue,isnr)], 'LineWidth', 2, 'Color', 'r');
            ylabel(ax2, ['SD Tissue ', num2str(itissue)], 'FontSize', 12);                                    
            grid on;
            set(ax1, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            set(ax2, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            if itissue == ntissue
                xlabel(ax1, 'ROI square dimension');                
                xlabel(ax2, 'ROI square dimension');
            elseif itissue == 1
                legend('Whole', 'Region', 'Location', 'eastoutside');
                title(ax1, ['RE of mean estimated D* for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                           
                title(ax2, ['SD of estimated D* for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                           
            end            
            
            I3 = figure(SNR(isnr)+2);
            ax1 = subplot(ntissue,2,2*itissue-1);             
            plot(ROI, 100*abs(avrg_f(:,itissue,isnr)-tecido(itissue,1))./tecido(itissue,1), 'bo', ...
                ROI, 100*abs(avrg_f2(:,itissue,isnr)-tecido(itissue,1))./tecido(itissue,1), 'ro', 'Linewidth', 2);            
            grid on;
            ylabel(ax1, ['RE (%) Tissue ', num2str(itissue)], 'FontSize', 12);  
            ax2 = subplot(ntissue,2, 2*itissue);
            plot(ROI, dp_f(:,itissue,isnr), 'bo', ...
                ROI, dp_f2(:,itissue,isnr), 'ro', 'Linewidth', 2);            
%             line(xlim, [avrg_f(:,itissue,isnr), avrg_f(:,itissue,isnr)], 'LineWidth', 2, 'Color', 'r');
            ylabel(ax1,['RE (%) Tissue ', num2str(itissue)], 'FontSize', 12);
            ylabel(ax2,['SD Tissue ', num2str(itissue)], 'FontSize', 12);
            grid on;
            set(ax1, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            set(ax2, 'XTickLabel', ROI, 'XTick', 2:numel(ROI)+1, 'FontSize', 12);
            if itissue == ntissue
                xlabel(ax1, 'ROI square dimension');                
                xlabel(ax2, 'ROI square dimension');
            elseif itissue == 1
                legend('Whole', 'Region', 'Location', 'eastoutside');
                title(ax1, ['RE of mean estimated f for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                           
                title(ax2, ['SD of estimated f for SNR = ', num2str(SNR(isnr))], 'FontSize', 12);                                                      
            end
            
%             set(I1, 'PaperUnits', 'centimeters');
%             set(I1, 'PaperPosition', [0 0 60 30]);                                    
%             ofn = fullfile(path, ['SNR = ',num2str(SNR(isnr)), '_Ddiff_AD.jpg']);
%             saveas(I1, ofn);    
%             
%             set(I2, 'PaperUnits', 'centimeters');
%             set(I2, 'PaperPosition', [0 0 60 30]);                                    
%             ofn = fullfile(path, ['SNR = ',num2str(SNR(isnr)), '_Dperf_AD.jpg']);
%             saveas(I2, ofn);           
%             
%             set(I3, 'PaperUnits', 'centimeters');
%             set(I3, 'PaperPosition', [0 0 60 30]);                                    
%             ofn = fullfile(path, ['SNR = ',num2str(SNR(isnr)), '_f_AD.jpg']);
%             saveas(I3, ofn);      
            
            C = C + d_t(itissue,1);
                        
        end       
    end
    %close all;
end







