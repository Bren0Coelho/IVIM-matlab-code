clear;
load('workspace_simulation_RETICULADO_v3_NOISE');
close all;
clc;
nSNR = length(SNR);
nROI = length(ROI);
rDim = 0;
region = 1;

estim = input('What is the estimated parameter?\n', 's');      % f, D_diff, D_perf 
clc;

for isnr = 1:nSNR
    
    for iroi = 1:nROI        
        
        dataDiff = Reticulation_data{iroi, isnr}.Ddiff_fit;
        dataPerf = Reticulation_data{iroi, isnr}.Dperf_fit;
        dataF = Reticulation_data{iroi, isnr}.f_fit;        
        
        if region == 0
            
            D_diff_show(:, iroi) = dataDiff(:);
            D_perf_show(:, iroi) = dataPerf(:);
            f_show(:, iroi) = dataF(:);
            
        else
            
            Obj = (dataDiff ~= 0);
            [R,C] = find(Obj,1);
            start2 = 1; final2 = 0;
            for itissue = 1:ntissue
                
               final2 = final2+(d_t(itissue,1)-2*rDim)*(d_t(itissue,2)-rDim);
                                
               dataDiff2 = dataDiff(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim);
               dataPerf2 = dataPerf(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim);
               dataF2 = dataF(R:R+d_t(itissue,2)-1-rDim,C+rDim:C+d_t(itissue,1)-1-rDim);
                    
               D_diff_show(start2:final2, iroi) = dataDiff2(:);
               D_perf_show(start2:final2, iroi) = dataPerf2(:);
               f_show(start2:final2, iroi) = dataF2(:);
               C = C + d_t(itissue,1);
                
               start2 = final2 + 1;
                
            end           
        end       
    end    
    alfa = 0.05;
    correct = 'dunn-sidak';
    mostra = 'off';
    switch estim
        case 'D_diff'
            [~, ~, stats1] = kruskalwallis(D_diff_show, ROI, 'off');
            [c(:,:, isnr),~,~,~] = multcompare(stats1, 'Alpha', alfa, 'CType', correct, 'Display', mostra);
        case 'D_perf'
            [~, ~, stats2] = kruskalwallis(D_perf_show, ROI, 'off');
            [c(:,:, isnr),~,~,~] = multcompare(stats2, 'Alpha', alfa, 'CType', correct, 'Display', mostra);             
        case 'f'
            [~, ~, stats3] = kruskalwallis(f_show, ROI, 'off');
            [c(:,:, isnr),~,~,~] = multcompare(stats3, 'Alpha', alfa, 'CType', correct, 'Display', mostra);                              
        otherwise
            error('No such parameter was estimated, only f, D_diff and D_perf');
    end

    Group1 = ROI(c(:,1, isnr));
    Group2 = ROI(c(:,2, isnr));
    lowerCI = c(:,3, isnr);
    upperCI = c(:,5, isnr);
    meanDiff = c(:,4, isnr);
    pValue = c(:,6, isnr);

    T = table(Group1', Group2', lowerCI, upperCI, meanDiff, pValue);

    T.Properties.Description = ['Multiple comparison test of parameter ', estim, ', SNR = ', num2str(SNR(isnr))] ;           
    disp(T.Properties.Description);            
    format shortG; disp(T);
    disp('______________________________________________________________________________________________________');


end