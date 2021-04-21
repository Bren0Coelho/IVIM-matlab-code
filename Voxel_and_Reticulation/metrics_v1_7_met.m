%% metrics_v1.7_met
% ==============================================================================
% ==============================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   06/16/2020: ver 1.0
%   06/19/2020: ver 1.1 - normalization before difference between simulated
%   and estimated parameters (Fernando)
%   06/20/2020: ver 1.2 - normalization after difference between simulated
%   and estimated parameters; new treatment of outliers -> positive
%   outliers receive maximum values, whereas negative values receive 0
%   06/21/2020: ver 1.3 - new treatment of outliers; they receive the
%   maximum value of assumed physiologic limits; calculates how
%   often each method reaches the minimum distance and plots a histogram;
%   06/23/2020: ver 1.4 - shows the distances depending on the b sequences
%   06/26/2020: ver 1.5 - shows how often each b sequence reaches the
%   minimum values
%   07/10/2020: ver 1.6 - the normalized distances are calculated with all
%   the methods at the same time; calculates p-value
%   07/10/2020: ver 1.7 - the normalized distances are calculated with all
%   the methods at the same time with no comparison between b-value
%   distributions
%   07/28/2020: ver 1.8 - calculates a score of methods values performance
%   08/21/2020: ver 1.9 - shows tables of mutiple comparison tests between
%   methods
%  
% ===============================================================================
% ===============================================================================
% - Calculates statistics factors to compare simulated parameters to estimated parameters
% - Metrics: euclidean distance between simulated parameter and estimated
% parameter
% - Steps:  1) normalize all estimated data points of each parameter and the 
%           simulated parameter: (p_est - p_min)/(p_max - p_min) and (p_sim - p_min)/(p_max - p_min)
%           2) calculates the differnece between estimated and simulated
%           parameters
%           3) estimates the euclidean distance between parameters;      
%           4) calcultates the summation of the distances per method;
%           5) plots one graphic per SNR per sequence showing distances vs. simulated
%           6) plots scores of methods performance 
%           7) plots histograms to show how often each method reaches the
%           minimum distance


close all; 
clc; 
clear;
%% Loads workspace (variables and simulation results)
[fn, dn] = uigetfile('*.mat', 'Select the workspace file');
odn = uigetdir(dn, 'Select the output folder');

load(fullfile(dn, fn));

%% Input
bseq = {'b1' 'b2' 'b3' 'b4'};   % b sequences
methods = {'LLS'; 'LLSR'; 'NLLS2'; 'TRR'; 'LEV'; 'NNLS'};

nest = length(s_sim);           % Number of simulated structures with parameters
nbdist = length(bseq);          % Number of simulated sequences of b 
nSNR = length(SNR);             % Number of SNR values
id_b = 4;                       % b sequence whose distances are to be plotted
nmethods = length(methods);
melhor_metodo = zeros(nest, 13, nSNR, nbdist);
p1 = zeros(nest, 13, nSNR, nbdist);
p2 = zeros(nest, 13, nSNR, nbdist);
p3 = zeros(nest, 13, nSNR, nbdist);
ranking = zeros(nest,nmethods,nSNR, nbdist);


%% Normalizes the parameters

   
for iest = 1 : nest
    
    s = eval(matlab.lang.makeValidName(['s', num2str(iest)]));
    nfs = length(s.f);
    nDdiff = length(s.D_diff);
    nDperf = length(s.D_perf);
    if nfs > 1
        distance = zeros(nfs, nmethods, nbdist);
        var = s.f;
        n = nfs;
        
    elseif nDdiff > 1
        distance = zeros(nDdiff, nmethods, nbdist);
        var = s.D_diff;
        n = nDdiff;

    else
        distance = zeros(nDperf, nmethods, nbdist);
        var = s.D_perf;
        n = nDperf;

    end
    
    

    for iSNR = 1 : nSNR

        for ibdist = 1 : nbdist
            
            data = eval(matlab.lang.makeValidName(['R_b', num2str(ibdist),...
                '_s', num2str(iest)]));  
            
            b = eval(matlab.lang.makeValidName(['b', num2str(ibdist)]));

            iparams = 1;
            curdata = zeros(nmethods*1000,3);
            delta_total = zeros(nmethods*1000,3,n,nbdist);
            
            for ifs = 1 : nfs
                curf = s.f(ifs);

                for iDdiff = 1 : nDdiff
                    curDdiff = s.D_diff(iDdiff);

                    for iDperf = 1 : nDperf
                        curDperf = s.D_perf(iDperf);
                        curparams = [curDdiff curDperf curf];
                        minimo_data = 1;

                        for imethods = 1 : nmethods
                            
%                             curdata = data(iparams, imethods).pontos(:,:,iSNR); 
                            
%-----------------------------10/07/2020: ver 1.6%-------------------------                            
                            curdata(1000*imethods-999:1000*imethods,:) = data(iparams, imethods).pontos(:,:,iSNR);
%--------------------------------------------------------------------------

%-----------------------------19/06/2020: ver 1.0%-------------------------
%                             curDdiffdata = data(iparams, imethods).pontos(:, 1, iSNR);
%                             curDperfdata = data(iparams, imethods).pontos(:, 2, iSNR);
%                             curfdata = data(iparams, imethods).pontos(:, 3, iSNR);
%                             
%                             norm_curdata = zeros(num,3);    % Norm of the estimated data
%                             norm_cur = zeros(num,3);        % Norm of simulated parameters
% 
%                                 for idata = 1:num
%                                     
%                                     if ((curDdiffdata(idata) <= 0 || curDdiffdata(idata) > 5e-2) || (curDperfdata(idata) <= 0 || curDperfdata(idata)...
%                                             > 5e-1) || (curfdata(idata) > 1 || curfdata(idata) <= 0))
%                                         norm_curdata(idata,1:3) = 1;
%                                         norm_cur(idata,1:3) = 0;
% 
%                                     else                                        
%                                         norm_curdata(idata,1) = (curDdiffdata(idata)-min(curDdiffdata(curDdiffdata>0)))./...
%                                             (max(curDdiffdata)-min(curDdiffdata(curDdiffdata>0)));
%                                         
%                                         norm_curdata(idata,2) = (curDperfdata(idata)-min(curDperfdata(curDperfdata>0)))./...
%                                             (max(curDperfdata)-min(curDperfdata(curDperfdata>0)));
%                                         
%                                         norm_curdata(idata,3) = (curfdata(idata)-min(curfdata(curfdata>0)))./...
%                                             (max(curfdata)-min(curfdata(curfdata>0)));
%                                         
%                                         norm_cur(idata,1) = (curDdiff-min(curDdiffdata(curDdiffdata>0)))./...
%                                             (max(curDdiffdata)-min(curDdiffdata(curDdiffdata>0)));
%                                         
%                                         norm_cur(idata,2) = (curDperf-min(curDperfdata(curDperfdata>0)))./...
%                                             (max(curDperfdata)-min(curDperfdata(curDperfdata>0)));
%                                         
%                                         norm_cur(idata,3) = (curf-min(curfdata(curfdata>0)))./...
%                                             (max(curfdata)-min(curfdata(curfdata>0)));                                                                          
%                                     end
%                                     
%                                 end
%                                 
%                                 
%                                 avr_distance(iparams, imethods) = mean(sqrt((norm_curdata(:,1)-norm_cur(:,1)).^2 + ...
%                                     (norm_curdata(:,2)-norm_cur(:,2)).^2 + (norm_curdata(:,3)-norm_cur(:,3)).^2));
%--------------------------------------------------------------------------
%---------------------------19/06/2020: ver 1.1----------------------------                           
%                         idx = (curdata(:,1) <= 0 | curdata(:,1) > 5e-2);
%                         curdata (idx,:) = 0;                            
% 
%                         idx = (curdata(:,2) <= 0 | curdata(:,2) > 5e-1);
%                         curdata (idx,:) = 0;
% 
%                         idx = (curdata(:,3) <= 0 | curdata(:,3) > 1);
%                         curdata (idx,:) = 0;

%-------------------------------------------------------------------------- 
                                    


%---------------------------20/06/2020: ver 1.2----------------------------
%                             
%                             idx1 = (curdata(:, 1) <= 0);
%                             curdata(idx1, :) = 0;
%                             
%                             idx2 = (curdata(:, 1) > 5e-2);
%                             curdata(idx2, :) = 5e-2;
%                             
%                             idx1 = (curdata(:, 2) <= 0);
%                             curdata (idx1, :) = 0;
%                             
%                             idx2 = (curdata(:, 2) > 5e-1);
%                             curdata (idx2, :) = 5e-1;
%                             
%                             idx1 = (curdata(:, 3) <= 0);
%                             curdata(idx1, :) = 0;
%                             
%                             idx2 = (curdata(:, 3) > 1);
%                             curdata (idx2, :) = 1;
%--------------------------------------------------------------------------  

%---------------------------21/06/2020: ver 1.3----------------------------                           
%                             idx1 = (curdata(:, 1) <= 0 | curdata(:, 1) > 5e-2);
%                             curdata (idx1, 1) = 5e-2;
%                             curdata (idx1, 2) = 5e-1;
%                             curdata (idx1, 3) = 1;
%                             
%                             idx2 = (curdata(:, 2) <= 0 | curdata(:, 2) > 5e-1);
%                             curdata (idx2, 1) = 5e-2;
%                             curdata (idx2, 2) = 5e-1;
%                             curdata (idx2, 3) = 1;
%                             
%                             idx3 = (curdata(:, 3) <= 0 | curdata(:, 3) > 1);
%                             curdata (idx3, 1) = 5e-2;
%                             curdata (idx3, 2) = 5e-1;
%                             curdata (idx3, 3) = 1;
%-------------------------------------------------------------------------- 


%                             delta(:,1) = curdata(:,1) - curparams(1);
%                             delta(:,2) = curdata(:,2) - curparams(2);
%                             delta(:,3) = curdata(:,3) - curparams(3);
%                             delta = abs(delta);
%                             mindelta = min(delta);
%                             maxdelta = max(delta);                  
%                             diffdelta = maxdelta - mindelta;
% 
%                             for tmp = 1 : 3
%                                 idx = (curdata(:,tmp) == 0);
%                                 delta(idx, tmp) = maxdelta(tmp);
%                             end
% 
%                             normdelta(:,1) = (delta(:,1) - mindelta(1)) ./ diffdelta(1);
%                             normdelta(:,2) = (delta(:,2) - mindelta(2)) ./ diffdelta(2);
%                             normdelta(:,3) = (delta(:,3) - mindelta(3)) ./ diffdelta(3);
% 
%                             
% %---------------------------20/06/2020: ver 1.2---------------------------                
% 
%                             distance(iparams, imethods, ibdist) = sum(sqrt(sum((normdelta(1000*imethods-999:1000*imethods,:)').^2)))...
%                                 /sum(sqrt(sum(ones(3,1000).^2)));
% %--------------------------------------------------------------------------
% 
% %---------------------------20/06/2020: ver 1.1----------------------------
%     %                             distance(iparams, imethods) = sum(sqrt(sum((normdelta').^2)));
% %--------------------------------------------------------------------------
% 
% %---------------------------21/06/2020: ver 1.3----------------------------
%                             if (minimo_data > distance(iparams, imethods, ibdist))
%                                 minimo_data = distance(iparams, imethods, ibdist);
%                                 minimo_method = imethods;
%                             end

%--------------------------------------------------------------------------

                        end
                        
%---------------------------10/07/2020: ver 1.6----------------------------                           




                        % Kruskal-Wallis test
                        x1 = reshape(curdata(:,1), num, nmethods);
                        x2 = reshape(curdata(:,2), num, nmethods);
                        x3 = reshape(curdata(:,3), num, nmethods);
                        [~, ~, stats1] = kruskalwallis(x1, methods, 'off');
                        [~, ~, stats2] = kruskalwallis(x2, methods, 'off');
                        [~, ~, stats3] = kruskalwallis(x3, methods, 'off');
                        
                        
                        
                        
                        % Post-hoc test
                        % The two first columns of c are the binomial
                        % combinations of methods represented by their
                        % positions in "methods", e.g., c(1, 1:2) = [1 2]
                        % represents a comparison between LLS and LLSR
                        
                        alfa = 0.05;
                        correct = 'dunn-sidak';
                        mostra = 'off';
                        [c1(:,:,iest, iparams, iSNR, ibdist),m1(:,:,iest, iparams, iSNR, ibdist),~,~] =...
                            multcompare(stats1, 'Alpha', alfa, 'CType', correct, 'Display', mostra);  



                        [c2(:,:,iest, iparams, iSNR, ibdist),m2(:,:,iest, iparams, iSNR, ibdist),~,~] =...
                            multcompare(stats2, 'Alpha', alfa, 'CType', correct, 'Display', mostra);



                        [c3(:,:,iest, iparams, iSNR, ibdist),m3(:,:,iest, iparams, iSNR, ibdist),~,~] =...
                            multcompare(stats3, 'Alpha', alfa, 'CType', correct, 'Display', mostra);
                            
                        
                        % Distance calculation
                        idx = (curdata(:,1) <= 0 | curdata(:,1) > 5e-2);
                        curdata (idx,:) = 0;                            

                        idx = (curdata(:,2) <= 0 | curdata(:,2) > 5e-1);
                        curdata (idx,:) = 0;

                        idx = (curdata(:,3) <= 0 | curdata(:,3) > 1);
                        curdata (idx,:) = 0;
                        
                        delta(:,1) = curdata(:,1) - curparams(1);
                        delta(:,2) = curdata(:,2) - curparams(2);
                        delta(:,3) = curdata(:,3) - curparams(3);                        
                        delta_abs = abs(delta);
                        mindelta = min(delta_abs);
                        maxdelta = max(delta_abs);                  
                        diffdelta = maxdelta - mindelta;

                        for tmp = 1 : 3
                            idx = (curdata(:,tmp) == 0);
                            delta_abs(idx, tmp) = maxdelta(tmp);
                        end

                        normdelta(:,1) = (delta_abs(:,1) - mindelta(1)) ./ diffdelta(1);
                        normdelta(:,2) = (delta_abs(:,2) - mindelta(2)) ./ diffdelta(2);
                        normdelta(:,3) = (delta_abs(:,3) - mindelta(3)) ./ diffdelta(3);

                        for imethods = 1:nmethods
                

                            distance(iparams, imethods, ibdist) = sum(sqrt(sum((normdelta(1000*imethods-999:1000*imethods,:)').^2)))...
                                /sum(sqrt(sum(ones(3,1000).^2)));



                            if (minimo_data > distance(iparams, imethods, ibdist))
                                minimo_data = distance(iparams, imethods, ibdist);
                                minimo_method = imethods;
                            end                           
 
                        end
%--------------------------------------------------------------------------
                        
                        melhor_metodo(iest, iparams, iSNR, ibdist) = minimo_method;                        
                        iparams = iparams + 1;
                    end
                end
            end
            
%---------------------------07/28/2020: ver 1.8----------------------------
            [~,index] = sort(distance(:,:,ibdist),2);
            [M,N] = size(index);
            for i = 1 : M
                for j = 1 : N
                    ranking(iest,index(i,j),iSNR,ibdist) = ranking(iest,index(i,j),iSNR,ibdist) + N - j + 1;                    
                end
            end
%--------------------------------------------------------------------------
        end
        
        
        
%         H = figure(iSNR*10);
%         subplot(3, 1, iest);
%         plot(var', distance(:, :, id_b), 'o--', 'LineWidth', 2);        
%         ylabel('Distance');
%         text = ['Normalized b', num2str(id_b), ' distances from simulated parameters and SNR = ', num2str(SNR(iSNR))];
%         if iest == 1
%             title(text);
%             xlabel('f');
%         elseif iest == 2
%             xlabel('D [mm²/s]');
%         else
%             xlabel('D* [mm²/s]');
%         end
%         ylim([0 1]);
%         hold on;
%         grid;
%         legend(methods, 'Location', 'bestoutside');
%         hold off;
%         set(H, 'Position', get(0, 'Screensize'));
%         text_save = ['Distance b', num2str(id_b), ' SNR' num2str(SNR(iSNR)) '.jpg'];
%         ofn = fullfile(odn, text_save);
%         saveas(H, ofn, 'jpg');

        
 
      

    end    

end
% close(H);
%---------------------------07/28/2020: ver 1.8----------------------------
% figure;
% bar(reshape(sum(sum(ranking),3),nmethods,nbdist)');
% grid on;
% title('Rank of methods per sequence');
% ylabel('Score');
% xlabel('Sequences');
% legend(methods, 'FontSize', 9, 'Location', 'bestoutside');
% set(gca, 'XTickLabel', bseq, 'XTick', 1:numel(bseq), 'FontSize', 12);
% 
% figure;
% bar(reshape(sum(sum(ranking),4),nmethods,nSNR)');
% grid on;
% title('Rank of methods per SNR value');
% ylabel('Score');
% xlabel('SNR');
% legend(methods, 'FontSize', 9, 'Location', 'bestoutside');
% set(gca, 'XTickLabel', SNR, 'XTick', 1:numel(SNR), 'FontSize', 12);
% 
% figure;
% bar(sum(sum(sum(ranking),3),4)); 
% title('Rank of methods');
% ylabel('Score');
% xlabel('Methods');
% grid on;
% set(gca, 'XTickLabel', methods, 'XTick', 1:numel(methods), 'FontSize', 12);
%--------------------------------------------------------------------------


%---------------------------21/06/2020: ver 1.3----------------------------

% I = figure('Name','How often each method reaches the minimum distance - per sequence');
% set(I, 'Position', get(0, 'Screensize'));
% for i = 1:length(bseq)
%     subplot(2,2,i);
%     histogram(melhor_metodo(:,:,:,i), 'BinLimits', [.5 length(methods)+.5]);
%     title(bseq(i));
%     ylabel('Occurrence');
%     xlabel('Methods');
%     set(gca, 'XTickLabel', methods, 'XTick', 1:numel(methods));    
% end
% 
% ofn = fullfile(odn, 'How often each method reaches the minimum distance_per sequence');
% saveas(I, ofn, 'jpg');
% %--------------------------------------------------------------------------
%  
% 
% 
% %---------------------------26/06/2020: ver 1.5----------------------------
% 
% J = figure; histogram(melhor_metodo, 'BinLimits', [.5 length(methods)+.5]);
% set(J, 'Position', get(0, 'Screensize'));
% title('How often each method reaches the minimum distance - Total');
% ylabel('Occurrence');
% xlabel('Methods');
% set(gca, 'XTickLabel', methods, 'XTick', 1:numel(methods));
% ofn = fullfile(odn, 'How often each method reaches the minimum distance_Total');
% saveas(J, ofn, 'jpg');

% %--------------------------------------------------------------------------


%---------------------------21/08/2020: ver 1.9----------------------------
show_p = input('Would you like to see any p-value results for significance difference? [yes/no]\n', 's');

if isequal(show_p, 'yes')
    statistics = true;
    while statistics    
        
        snr = input('What is the SNR value?\n');                       % 20, 25, 30, 35, 40, 45, 50
        i_snr = find(SNR == snr);
        
        sequence = input('What is the b-value sequence?\n', 's');      % b1, b2, b3, b4       
        
        i_seq = 1;
        while (~isequal(bseq{i_seq},sequence))            
            i_seq = i_seq + 1;            
        end

        estim = input('What is the estimated parameter?\n', 's');      % f, D_diff, D_perf 
        
        switch estim
            case 'f'
                c = c3;               
            case 'D_diff'
                c = c1;
            case 'D_perf'
                c = c2;
            otherwise
                error('No such parameter was estimated, only f, D_diff and D_perf');

        end
        
        str = input('What is the structure?\n');                       % 1, 2, 3
        
        switch str
            case 1
                n = length(s1.f);
                s = s1.f;               
            case 2
                n = length(s2.D_diff);
                s = s2.D_diff;
            case 3
                n = length(s3.D_perf);
                s = s3.D_perf;
            otherwise
                error('No such structure was simulated, only 1, 2 and 3');
             
        end
                
        
        clc; 
        for param = 1:n
            Group1 = methods(c(:,1, str, param, i_snr, i_seq));
            Group2 = methods(c(:,2, str, param, i_snr, i_seq));
            lowerCI = c(:,3, str, param, i_snr, i_seq);
            upperCI = c(:,5, str, param, i_snr, i_seq);
            meanDiff = c(:,4, str, param, i_snr, i_seq);
            pValue = c(:,6, str, param, i_snr, i_seq);

            T = table(Group1, Group2, lowerCI, upperCI, meanDiff, pValue);
            
            T.Properties.Description = ['Multiple comparison test of parameter ', estim, ' for sequence ', bseq{i_seq}, ', SNR = ', num2str(snr), ', structure ', ...
                num2str(str), ' and parameter ', num2str(s(param))] ;           
            disp(T.Properties.Description);            
            format shortG; disp(T);
            disp('______________________________________________________________________________________________________');
            
        end
        
        show_p = input('Would you like to see any other p-value results for significance difference? [yes/no]\n', 's');
        if (isequal(show_p, 'no'))
            statistics = false;
        end
    end
end
%--------------------------------------------------------------------------