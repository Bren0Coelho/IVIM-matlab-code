close all; clc;

%% Input
nest = 3;   % Number of simulated structures with parameters
nbdist = 4; % Number of simulated sequences of b

methods = {'LLS'; 'LLSR'; 'NLLS2'; 'TRR'; 'LEV'; 'NNLS'};
nmethods = length(methods);
outDdiff_total = zeros(3,1);
outDperf_total = zeros(3,1);
outf_total = zeros(3,1);
erroDdiff_total = zeros(3,1);
erroDperf_total = zeros(3,1);
errof_total = zeros(3,1);
dpDdiff_total = zeros(3,1);
dpDperf_total = zeros(3,1);
dpf_total = zeros(3,1);
d = 0;

%% Loads workspace (variables and simulation results)
[fn, dn] = uigetfile('*.mat', 'Select the workspace file');
% odn = uigetdir(dn, 'Select the output folder');

load(fullfile(dn, fn));

%% Scatter and histogram
for iest = 1 : nest
    
    for ibdist = 1 : nbdist
        
        s = eval(matlab.lang.makeValidName(['s', num2str(iest)]));
        nfs = length(s.f);
        nDdiff = length(s.D_diff);
        nDperf = length(s.D_perf);
                     
        data = eval(matlab.lang.makeValidName(['R_b', num2str(ibdist),...
            '_s', num2str(iest)]));
        
        b = eval(matlab.lang.makeValidName(['b', num2str(ibdist)]));
        nSNR = length(SNR);
        
        iparams = 1;
        
        for ifs = 1 : nfs
            curf = s.f(ifs);
            
            for iDdiff = 1 : nDdiff
                curDdiff = s.D_diff(iDdiff);
                
                for iDperf = 1 : nDperf
                    curDperf = s.D_perf(iDperf);
                    
                    for imethods = 1 : nmethods                 

                        
                        for iSNR = 1 : nSNR
                            
                            curdata = data(iparams, imethods).pontos;
                            
                            curDdiffdata = data(iparams, imethods).pontos(:, 1, iSNR);
                            curDperfdata = data(iparams, imethods).pontos(:, 2, iSNR);
                            curfdata = data(iparams, imethods).pontos(:, 3, iSNR);
                            
                            vcurDdiffdata = curDdiffdata(curDdiffdata > 0 & curDdiffdata < 5e-2);
                            vcurDperfdata = curDperfdata(curDperfdata > 0 & curDperfdata < 5e-1);
                            vcurfdata = curfdata(curfdata > 0 & curfdata < 1);
                            
                            outDdiff_total(iest) = outDdiff_total(iest) + length(curDdiffdata) - length(vcurDdiffdata);
                            outDperf_total(iest) = outDperf_total(iest) + length(curDperfdata) - length(vcurDperfdata);
                            outf_total(iest) = outf_total(iest) + length(curfdata) - length(vcurfdata);
            
                           
                            % Percentual error and standard deviation                            
                            erroDdiff_total(iest) = erroDdiff_total(iest) + abs(mean ((vcurDdiffdata - curDdiff) / curDdiff));
                            erroDperf_total(iest) = erroDperf_total(iest) + abs(mean ((vcurDperfdata - curDperf) / curDperf));
                            errof_total(iest) = errof_total(iest) + abs(mean ((vcurfdata - curf) / curf));
                            
                            dpDdiff_total(iest) = dpDdiff_total(iest) + std(vcurDdiffdata);
                            dpDperf_total(iest) = dpDperf_total(iest) + std(vcurDperfdata);
                            dpf_total(iest) = dpf_total(iest) + std(vcurfdata);
                            
                            d = d + 1;
                        end            
                    end
                    iparams = iparams + 1;
                end
            end
        end
    end
end

outlier = [outDdiff_total outDperf_total outf_total]./d;
erro = 100*[erroDdiff_total erroDperf_total errof_total]./d;
dp = [dpDdiff_total dpDperf_total dpf_total]./d;