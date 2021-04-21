%% metrics_v1.7_params
% ==============================================================================
% ==============================================================================
% Version Control
%
% Fernando Paiva @ CIERMag/IFSC/USP
%   06/16/2020: ver 1.0
%   
%
% ===============================================================================
% ===============================================================================
% - Calculates statistics factors to compare simulated parameters to 
% estimated parameters
% - Steps:  calculates p-value between distributions of estimated parameters
% to verify whether struct variations are likely to influence estimation of
% constant parameters (alpha = 5 %)

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
nmethods = length(methods);
p1 = zeros(nest, nbdist, nSNR, nmethods);
p2 = zeros(nest, nbdist, nSNR, nmethods);
p3 = zeros(nest, nbdist, nSNR, nmethods);

%% Calculates p-values
   
for iest = 1 : nest
    
        
    s = eval(matlab.lang.makeValidName(['s', num2str(iest)]));
    nfs = length(s.f);
    nDdiff = length(s.D_diff);
    nDperf = length(s.D_perf);
    if nfs > 1
        distance = zeros(nfs, nbdist, nmethods);
        group = s.f;
        n = nfs;

    elseif nDdiff > 1
        distance = zeros(nfs, nbdist, nmethods);
        group = s.D_diff;
        n = nDdiff;

    else
        distance = zeros(nfs, nbdist, nmethods);
        group = s.D_perf;
        n = nDperf;

    end

    for iSNR = 1 : nSNR

        for imethods = 1 : nmethods
            
            for ibdist = 1 : nbdist
                data = eval(matlab.lang.makeValidName(['R_b', num2str(ibdist),...
                '_s', num2str(iest)]));  

                b = eval(matlab.lang.makeValidName(['b', num2str(ibdist)]));
                curdata = zeros(n*num,3);

                for ifs = 1 : nfs
                    curf = s.f(ifs);

                    for iDdiff = 1 : nDdiff
                        curDdiff = s.D_diff(iDdiff);

                        for iDperf = 1 : nDperf
                            curDperf = s.D_perf(iDperf);
                            curparams = [curDdiff curDperf curf];
                            
                            for iparams = 1 : n

                                curdata(num*iparams-999:num*iparams,:) = data(iparams, imethods).pontos(:,:,iSNR); % [D D* f]

                            end
                            
                            % p1 (D)  - one must check iest 1 and 3,
                            % because struct 2 is the one that varies. So
                            % iest 2 will be < 0.05 at some point
                            % (probably for high SNR values)
                            
                            % p2 (D*) - one must check iest 1 and 2
                            % because struct 3 is the one that varies. So
                            % iest 3 will be < 0.05 at some point
                            % (probably for high SNR values)
                            
                            % p3 (f)  - one must check iest 2 and 3
                            % because struct 1 is the one that varies. So
                            % iest 1 will be < 0.05 at some point
                            % (probably for high SNR values)
                            
                            % Kruskal-Wallis test and Post-hoc test
                            % The two first columns of c are the binomial
                            % combinations of methods represented by their
                            % positions in "methods", e.g., c(1, 1:2) = [1 2]
                            % represents a comparison between LLS and LLSR                            

                            
                            if iest == 1
                                
                                x1 = reshape(curdata(:,1), num, n);
                                x2 = reshape(curdata(:,2), num, n);

                                [~, ~, stats1] = kruskalwallis(x1, group, 'off');
                                [~, ~, stats2] = kruskalwallis(x2, group, 'off');
                                
                                [c1_s1(:,:,1,ibdist,iSNR,imethods),m1_s1(:,:,1,ibdist,iSNR,imethods),~,~] =...
                                    multcompare(stats1, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');  


                                [c2_s1(:,:,1,ibdist,iSNR,imethods),m2_s1(:,:,1,ibdist,iSNR,imethods),~,~] =...
                                    multcompare(stats2, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');                             
                                
                            elseif iest == 2
                                
                                x2 = reshape(curdata(:,2), num, n);    
                                x3 = reshape(curdata(:,3), num, n);
                                
                                [~, ~, stats2] = kruskalwallis(x2, group, 'off');
                                [~, ~, stats3] = kruskalwallis(x3, group, 'off');
                               
                                [c2(:,:,iest, ibdist, iSNR, imethods),m2(:,:,iest, ibdist, iSNR, imethods),~,~] =...
                                    multcompare(stats2, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');


                                [c3(:,:,iest, ibdist, iSNR, imethods),m3(:,:,iest, ibdist, iSNR, imethods),~,~] =...
                                    multcompare(stats3, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');                               
                                
                            else
                                
                                x1 = reshape(curdata(:,1), num, n);    
                                x3 = reshape(curdata(:,3), num, n);
                                
                                [~, ~, stats1] = kruskalwallis(x1, group, 'off');
                                [~, ~, stats3] = kruskalwallis(x3, group, 'off');

                                [c1(:,:,iest, ibdist, iSNR, imethods),m1(:,:,iest, ibdist, iSNR, imethods),~,~] =...
                                    multcompare(stats1, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');                                 


                                [c3(:,:,iest, ibdist, iSNR, imethods),m3(:,:,iest, ibdist, iSNR, imethods),~,~] =...
                                    multcompare(stats3, 'Alpha', 0.05, 'CType', 'dunn-sidak', 'Display', 'off');
                                
                            end                            
                        end
                    end
                end
            end
        end
    end    
end

show_p = input('Would you like to see any p-value results for significance difference? [yes/no]\n', 's');

if isequal(show_p, 'yes')
    statistics = true;
    while statistics    
        
        snr = input('What is the SNR value?\n');                       % 20, 25, 30, 35, 40, 45, 50
        
        str = input('What is the structure?\n');                       % 1, 2, 3
        while(str ~= 1 && str ~= 2 && str ~= 3)
            str = input('The structure must be 1, 2 or 3\n', 's');
        end                
        
        sequence = input('What is the sequence?\n', 's');              % b1, b2, b3, b4
        i_seq = 1;
        while (~isequal(bseq{i_seq},sequence))            
            i_seq = i_seq + 1;            
        end
        i_snr = find(SNR == snr);        

        estim = input('What is the estimated parameter?\n', 's');      % f, D_diff, D_perf        
        switch str
            case 1
                if isequal(estim, 'D_diff')
                    c = c1_s1;
                elseif isequal(estim, 'D_perf')
                    c = c2_s1;
                else
                    error('The estimated parameter must be D_diff or D_perf\n');                     
                end
                n = length(s1.f);
                s = s1.f;               
            case 2
                if isequal(estim, 'f')
                    c = c3;
                elseif isequal(estim, 'D_perf')
                    c = c2;
                else
                    error('The estimated parameter must be f or D_perf\n');
                end
                n = length(s2.D_diff);
                s = s2.D_diff;
            case 3
                if isequal(estim, 'f')
                    c = c3;
                elseif isequal(estim, 'D_diff')
                    c = c1;
                else
                    error('The estimated parameter must be f or D_diff\n');
                end
                n = length(s3.D_perf);
                s = s3.D_perf;
        end
        
        clc;
        
        for param = 1:nmethods
            
            Group1 = s(c(:,1, str, i_seq, i_snr, param));
            Group2 = s(c(:,2, str, i_seq, i_snr, param));
            lowerCI = c(:,3, str, i_seq, i_snr, param);
            upperCI = c(:,5, str, i_seq, i_snr, param);
            meanDiff = c(:,4, str, i_seq, i_snr, param);
            pValue = c(:,6, str, i_seq, i_snr, param);

            T = table(Group1', Group2', lowerCI, upperCI, meanDiff, pValue);
            
            T.Properties.Description = ['Multiple comparison test of parameter ', estim, ' for sequence ', sequence, ', SNR = ', num2str(snr), ', structure ', ...
                num2str(str), ' and method ', methods{param}] ;           
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
