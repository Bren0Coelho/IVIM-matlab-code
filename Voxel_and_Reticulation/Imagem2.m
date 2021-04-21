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
load('workspace_simulation_RETICULADO_v2');
close all;
clc;
nSNR = length(SNR);
nROI = length(ROI);
start = 1; final = 0;
start2 = 1; final2 = 0;

tec_start_r = (Nr - d_t(itissue,2))/2 + 1;
tec_start_c = (Nc - sum(d_t(:,1)))/2 + 1;

tec_end_r = tec_start_r + d_t(itissue,2) - 1;
tec_end_c = tec_start_c + d_t(itissue,1) - 1;

%% Loads workspace (variables and simulation results)

for itissue = 1:ntissue   
    
    for isnr = 1:nSNR
        
       path = ['C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\Imagens\Reticulado\SNR ', num2str(SNR(isnr))];
       
       for iroi = 1:nROI
                    
          D_diff_show = Reticulation_data{iroi, isnr}.Ddiff_fit;
          D_perf_show = Reticulation_data{iroi, isnr}.Dperf_fit;
          f_show = Reticulation_data{iroi, isnr}.f_fit;
          final = final+d_t(itissue,1)*d_t(itissue,2);
                        
          D_diff_show_lin = D_diff_show(tec_start_r:tec_end_r, tec_start_c:tec_end_c);
          Ddiff(start:final) = D_diff_show_lin(:);
          Origin_Ddiff(start:final) = repmat(ROI(iroi), size(Ddiff(start:final)));                               
         
          D_perf_show_lin = D_perf_show(tec_start_r:tec_end_r, tec_start_c:tec_end_c);
          Dperf(start:final) = D_perf_show_lin(:);
          Origin_Dperf(start:final) = repmat(ROI(iroi), size(Dperf(start:final)));    

          f_show_lin = f_show(tec_start_r:tec_end_r, tec_start_c:tec_end_c);
          F(start:final) = f_show_lin(:);
          Origin_f(start:final) = repmat(ROI(iroi), size(F(start:final)));
          
          start = final + 1;
          
          % Region with no image background influence
          
          final2 = final2+(d_t(itissue,1)-8)*(d_t(itissue,2)-4);
          
          D_diff_show_lin2 = D_diff_show(tec_start_r:tec_end_r-4, tec_start_c+4:tec_end_c-4);
          Ddiff2(start2:final2) = D_diff_show_lin2(:);
          Origin_Ddiff2(start2:final2) = repmat(ROI(iroi), size(Ddiff2(start2:final2)));                               
         
          D_perf_show_lin2 = D_perf_show(tec_start_r:tec_end_r-4, tec_start_c+4:tec_end_c-4);
          Dperf2(start2:final2) = D_perf_show_lin2(:);
          Origin_Dperf2(start2:final2) = repmat(ROI(iroi), size(Dperf2(start2:final2)));    

          f_show_lin2 = f_show(tec_start_r:tec_end_r-4, tec_start_c+4:tec_end_c-4);
          F2(start2:final2) = f_show_lin2(:);
          Origin_f2(start2:final2) = repmat(ROI(iroi), size(F2(start2:final2)));
          
          start2 = final2 + 1;

       end   
        
        I = figure(SNR(isnr) + itissue);       
        subplot(3,1,1);       
        boxplot(F,Origin_f);
        title(['Estimated parameters per ROI dimension SNR = ', num2str(SNR(isnr)), ' and tissue ', num2str(itissue)]);
        ylabel('f', 'FontSize', 12);       

        subplot(3,1,2);
        boxplot(Ddiff,Origin_Ddiff);
        ylabel('D [mm²/s]', 'FontSize', 12);

        subplot(3,1,3);
        boxplot(Dperf,Origin_Dperf);
        ylabel('D* [mm²/s]', 'FontSize', 12);
        xlabel('ROI square dimension', 'FontSize', 12)
%         set(gcf, 'Position', get(0, 'Screensize'));

        set(I, 'PaperUnits', 'centimeters');
        set(I, 'PaperPosition', [0 0 60 30]);                                    
        ofn = fullfile(path, ['BOX_SNR = ',num2str(SNR(isnr)), '_tissue', num2str(itissue),'.jpg']);
        saveas(I, ofn);
        
        I = figure(SNR(isnr) + itissue);       
        subplot(3,1,1);       
        boxplot(F2,Origin_f2);
        title(['Estimated parameters per ROI dimension of region SNR = ', num2str(SNR(isnr)), ' and tissue ', num2str(itissue)]);
        ylabel('f', 'FontSize', 12);       

        subplot(3,1,2);
        boxplot(Ddiff2,Origin_Ddiff2);
        ylabel('D [mm²/s]', 'FontSize', 12);

        subplot(3,1,3);
        boxplot(Dperf2,Origin_Dperf2);
        ylabel('D* [mm²/s]', 'FontSize', 12);
        xlabel('ROI square dimension', 'FontSize', 12)
%         set(gcf, 'Position', get(0, 'Screensize'));

        set(I, 'PaperUnits', 'centimeters');
        set(I, 'PaperPosition', [0 0 60 30]);                                    
        ofn = fullfile(path, ['BOX_Region_SNR = ',num2str(SNR(isnr)), '_tissue', num2str(itissue),'.jpg']);
        saveas(I, ofn); 
    end 
    tec_start_r = (Nr - d_t(itissue,2))/2 + 1;
    tec_start_c =  tec_end_c + 1;
    
    tec_end_r = tec_start_r + d_t(itissue,2) - 1;
    tec_end_c = tec_start_c + d_t(itissue,1) - 1;
end

close all;