%% metrics_v1.7_b
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
%   the b-value distributions at the same time
%   07/28/2020: ver 1.8 - calculates a score of sequences values performance
%   08/21/2020: ver 1.9 - shows tables of mutiple comparison tests between
%   b-value sequences
%   
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
%           4) calcultates the summation of the distances per distribution;
%           5) plots one graphic per SNR per method showing distances vs. simulated
%           6) plots scores of sequences performance  
%           7) plots histograms to show how often each sequence reaches the
%           minimum distance.

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
id_met = 4;                     % Method whose distances are to be plotted
nmethods = length(methods);
melhor_seq = zeros(nest, 13, nSNR, nmethods);
p1 = zeros(nest, 13, nSNR, nmethods);
p2 = zeros(nest, 13, nSNR, nmethods);
p3 = zeros(nest, 13, nSNR, nmethods);
ranking = zeros(nest,nbdist,nSNR,nmethods);

%% Normalizes the parameters

   
for iest = 1 : nest
    
    s = eval(matlab.lang.makeValidName(['s', num2str(iest)]));
    nfs = length(s.f);
    nDdiff = length(s.D_diff);
    nDperf = length(s.D_perf);
    if nfs > 1
        distance = zeros(nfs, nbdist, nmethods);
        var = s.f;
        n = nfs;

    elseif nDdiff > 1
        distance = zeros(nfs, nbdist, nmethods);
        var = s.D_diff;
        n = nDdiff;

    else
        distance = zeros(nfs, nbdist, nmethods);
        var = s.D_perf;
        n = nDperf;

    end

    for iSNR = 1 : nSNR

        for imethods = 1 : nmethods   
            
            

            iparams = 1;
            curdata = zeros(nbdist*1000,3);
            
            for ifs = 1 : nfs
                curf = s.f(ifs);

                for iDdiff = 1 : nDdiff
                    curDdiff = s.D_diff(iDdiff);

                    for iDperf = 1 : nDperf
                        curDperf = s.D_perf(iDperf);
                        curparams = [curDdiff curDperf curf];
                        minimo_data = 1;

                        for ibdist = 1 : nbdist
                            
                            data = eval(matlab.lang.makeValidName(['R_b', num2str(ibdist),...
                            '_s', num2str(iest)]));  

                            b = eval(matlab.lang.makeValidName(['b', num2str(ibdist)]));
                            
%                             curdata = data(iparams, imethods).pontos(:,:,iSNR); 
                            
%-----------------------------10/07/2020: ver 1.6%-------------------------                            
                            curdata(1000*ibdist-999:1000*ibdist,:) = data(iparams, imethods).pontos(:,:,iSNR);

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
                        x1 = reshape(curdata(:,1), num, nbdist);
                        x2 = reshape(curdata(:,2), num, nbdist);
                        x3 = reshape(curdata(:,3), num, nbdist);
                        [~, ~, stats1] = kruskalwallis(x1, bseq, 'off');
                        [~, ~, stats2] = kruskalwallis(x2, bseq, 'off');
                        [~, ~, stats3] = kruskalwallis(x3, bseq, 'off');
                        
                        % Post-hoc test
                        % The two first columns of c are the binomial
                        % combinations of methods represented by their
                        % positions in "bseq", e.g., c(1, 1:2) = [1 2]
                        % represents a comparison between b1 and b2
                        
                        % Post-hoc test                          
                            
                        [c1(:,:,iest, iparams, iSNR, imethods),m1(:,:,iest, iparams, iSNR, imethods),~,~] =...
                            multcompare(stats1, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');  


                        [c2(:,:,iest, iparams, iSNR, imethods),m2(:,:,iest, iparams, iSNR, imethods),~,~] =...
                            multcompare(stats2, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');


                        [c3(:,:,iest, iparams, iSNR, imethods),m3(:,:,iest, iparams, iSNR, imethods),~,~] =...
                            multcompare(stats3, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');

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

                        for ibdist = 1:nbdist
                
                            distance(iparams, ibdist, imethods) = sum(sqrt(sum((normdelta(1000*ibdist-999:1000*ibdist,:)').^2)))...
                                /sum(sqrt(sum(ones(3,1000).^2)));
                            
                            if (minimo_data > distance(iparams, ibdist, imethods))
                                minimo_data = distance(iparams, ibdist, imethods);
                                minimo_seq = ibdist;
                            end
                        end
%--------------------------------------------------------------------------

                        melhor_seq(iest, iparams, iSNR, imethods) = minimo_seq;                        
                        iparams = iparams + 1;
                    end
                end
            end
%---------------------------07/28/2020: ver 1.8----------------------------
            [~,index] = sort(distance(:,:,imethods),2);
            [M,N] = size(index);
            for i = 1 : M
                for j = 1 : N
                    ranking(iest,index(i,j),iSNR,imethods) = ranking(iest,index(i,j),iSNR,imethods) + N - j + 1;                    
                end
            end
%--------------------------------------------------------------------------
        end
        
        
%         H = figure(iSNR*10); 
%         subplot(3, 1, iest);        
%         plot(var', distance(:, :, id_met), 'o--', 'LineWidth', 2);        
%         ylabel('Distance');
%         text = ['Normalized ', methods{id_met}, ' distances from simulated parameters and SNR = ', num2str(SNR(iSNR))];
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
%         legend(bseq, 'Location', 'bestoutside');
%         hold off;
%         set(H, 'Position', get(0, 'Screensize'));
%         ofn = fullfile(odn, ['Distance NNLS SNR' num2str(SNR(iSNR)) '.jpg']);
%         saveas(H, ofn, 'jpg');     
   

    end    

end
% close(H);

%---------------------------07/28/2020: ver 1.8----------------------------
% figure;
% %total summation of points per method: 10*(nf+nD+nD*)*nSNR
% bar(reshape(sum(sum(ranking),3),nbdist,nmethods)'); 
% grid on;
% title('Rank of sequences per method');
% ylabel('Score');
% xlabel('Methods');
% legend(bseq, 'FontSize', 10, 'Location', 'bestoutside');
% set(gca, 'XTickLabel', methods, 'XTick', 1:numel(methods), 'FontSize', 12);
% 
% figure;
% bar(reshape(sum(sum(ranking),4),nbdist,nSNR)');
% grid on;
% title('Rank of sequences per SNR value');
% ylabel('Score');
% xlabel('SNR');
% legend(bseq, 'FontSize', 10, 'Location', 'bestoutside');
% set(gca, 'XTickLabel', SNR, 'XTick', 1:numel(SNR), 'FontSize', 12);
% 
% figure;
% bar(sum(sum(sum(ranking),3),4)); 
% title('Rank of sequences');
% ylabel('Score');
% xlabel('Sequences');
% grid on;
% set(gca, 'XTickLabel', bseq, 'XTick', 1:numel(bseq), 'FontSize', 12);
% %--------------------------------------------------------------------------
% 
% %---------------------------26/06/2020: ver 1.5----------------------------
% 
% for iSNR = 1:nSNR
%     J = figure('Name', ['How often each sequence reaches the minimum distance; SNR = ', num2str(SNR(iSNR))]);
%     set(J, 'Position', get(0, 'Screensize'));
%     for i = 1:length(methods)
%         subplot(3,2,i);
%         histogram(melhor_seq(:,:,iSNR,i), 'BinLimits', [.5 length(bseq)+.5]);
%         title(methods(i));
%         ylabel('Occurrence');
%         xlabel('Sequence');
%         set(gca, 'XTickLabel', bseq, 'XTick', 1:numel(bseq));                
%     end    
%     ofn = fullfile(odn, ['How often each sequence reaches the minimum distance; SNR = ', num2str(SNR(iSNR))]);
%     saveas(J, ofn, 'jpg');
% %     
% end
% % % 
% % % 
% K = figure('Name', 'How often each sequence reaches the minimum distance per method');
% set(K, 'Position', get(0, 'Screensize'));
% for i = 1:length(methods)
%     subplot(3,2,i);
%     histogram(melhor_seq(:,:,:,i), 'BinLimits', [.5 length(bseq)+.5]);
%     title(methods(i));
%     ylabel('Occurrence');
%     xlabel('Sequence');
%     set(gca, 'XTickLabel', bseq, 'XTick', 1:numel(bseq));    
% end
% ofn = fullfile(odn, 'How often each sequence reaches the minimum distance per method');
% saveas(K, ofn, 'jpg');
% % 
% % 
% L = figure; histogram(melhor_seq, 'BinLimits', [.5 length(bseq)+.5]);
% set(L, 'Position', get(0, 'Screensize'));
% title('How often each sequence reaches the minimum distance - Total');
% ylabel('Occurrence');
% xlabel('Sequence');
% set(gca, 'XTickLabel', bseq, 'XTick', 1:numel(bseq));
% ofn = fullfile(odn, 'How often each sequence reaches the minimum distance_Total');
% saveas(L, ofn, 'jpg');
% %--------------------------------------------------------------------------

%---------------------------21/08/2020: ver 1.9----------------------------

%%
show_p = input('Would you like to see any p-value results for significance difference? [yes/no]\n', 's');

if isequal(show_p, 'yes')
    statistics = true;
    while statistics    
        
        snr = input('What is the SNR value?\n');                       % 20, 25, 30, 35, 40, 45, 50
        met = input('What is the method?\n', 's');                     % LLS, LLSR, NLLS2, TRR, LEV, NNLS
        str = input('What is the structure?\n');                       % 1, 2, 3
        estim = input('What is the estimated parameter?\n', 's');      % f, D_diff, D_perf
        
        i_met = 1;
        while (~isequal(methods{i_met},met))            
            i_met = i_met + 1;            
        end
        i_snr = find(SNR == snr);
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
            Group1 = bseq(c(:,1, str, param, i_snr, i_met));
            Group2 = bseq(c(:,2, str, param, i_snr, i_met));
            lowerCI = c(:,3, str, param, i_snr, i_met);
            upperCI = c(:,5, str, param, i_snr, i_met);
            meanDiff = c(:,4, str, param, i_snr, i_met);
            pValue = c(:,6, str, param, i_snr, i_met);

            T = table(Group1', Group2', lowerCI, upperCI, meanDiff, pValue);
            
            T.Properties.Description = ['Multiple comparison test of parameter ', estim, ' for method ', met, ', SNR = ', num2str(snr), ', structure ', ...
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
