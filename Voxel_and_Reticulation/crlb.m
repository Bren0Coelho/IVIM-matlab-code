%% crlb
%% Formulação retirada de Fundamentals of Statistical Signal Processing -- Estimation Theory p.35
%% Funcao Cramer-Rao Lower Bounds calcula as variancias maximas dos parametros de ruido difusao e perfusao ;
%% Entrada: b, Signal;
%% Saida: var_D_perf, var_D_diff, var_f_fit

close all;
clear;
clc;
    
b1 = [0 5 10 15 20 30 35 40 45 50 55 60 70 80 90 100 150 200 250 300 350 400 450 500 1000];
b2 = [0 10 20 30 40 50 75 100 200 400 600 1000];          
b3 = [0 5 10 15 30 60 120 250 500 1000];
b4 = [0 20 40 80 200 400 700 1000];
bseq = {b1 b2 b3 b4};
bs = {'b1' 'b2' 'b3' 'b4'};
var = zeros(length(bseq),3,3);
%% Escrever a funcao de probabilidade 

sigma_2 = 3e-3;

%     p = exp(-sum((Signal-(1-f)*exp(-b*D_diff)-f*exp(-b*D_perf)).^2)/(2*sigma_2))/((2*pi*sigma_2)^(N/2));
%     ln = log(p);
%     s = (1-f)*exp(-b*D_diff) + f*exp(-b*D_perf);
%     Signal = awgn(s,40, 'measured');
%% Escrever os logaritmos da função de probabilidade

for i = 1:4
    
    b = bseq{i};
    
    for j = 1:3
        if j == 1
            D_diff = 3e-3;
            D_perf = 4e-2;
            f = 0.1:0.05:0.3;                      
        elseif j ==2
            D_diff = 1e-3:0.5e-3:7e-3;
            D_perf = 4e-2;
            f = 0.2;     
        else
            D_diff = 3e-3;
            D_perf = 1e-2:0.5e-2:7e-2;
            f = 0.2;
        end
        
        for k = 1:length(D_diff)
            for l = 1:length(D_perf)
                for m = 1:length(f)        
                    
                    dSdf = exp(-b*D_perf(l))-exp(-b*D_diff(k));
                    dSdDperf = (-b).*f(m).*exp(-b*D_perf(l));
                    dSdDdiff = (-b).*(1-f(m)).*exp(-b*D_diff(k));


                    %% Escrever as derivadas dos logaritmos na matriz de Fisher

                    % Termos da diagonal principal

                    Fisher(1,1) = (-1/sigma_2)*sum(dSdDdiff.^2);
                    Fisher(2,2) = (-1/sigma_2)*sum(dSdDperf.^2);
                    Fisher(3,3) = (-1/sigma_2)*sum(dSdf.^2);

                    % Termos restantes

                    Fisher(1,2) = -(1/sigma_2)*sum(dSdDdiff.*dSdDperf);
                    Fisher(1,3) = -(1/sigma_2)*sum(dSdDdiff.*dSdf);   

                    Fisher(2,1) = -(1/sigma_2)*sum(dSdDdiff.*dSdDperf);
                    Fisher(2,3) = -(1/sigma_2)*sum(dSdDperf.*dSdf);

                    Fisher(3,1) = -(1/sigma_2)*sum(dSdDdiff.*dSdf);
                    Fisher(3,2) = -(1/sigma_2)*sum(dSdDperf.*dSdf);   


                    Fisher_inv = inv(-Fisher);             



                    %% Escrever as variancias de acordo com a matriz de Fisher

                    %% plot CRLB

                    var(i,j,1) = var(i,j,1) + Fisher_inv(1,1);
                    var(i,j,2) = var(i,j,2) + Fisher_inv(2,2);
                    var(i,j,3) = var(i,j,3) + Fisher_inv(3,3);
                end
            end
        end
    end
end

figure;
subplot(3,1,1)
plot([1 2 3 4]',var(:,:,1), 'x--', 'LineWidth', 2);
title('CRLB variance');
ylabel('\Sigma \sigma_{D}^2', 'FontSize', 15); 
set(gca, 'XTickLabel', bs, 'XTick', 1:numel(bs), 'Fontsize', 12);
legend({'struct 1', 'struct 2', 'struct 3'}, 'Location', 'bestoutside','FontSize', 12);
grid;

subplot(3,1,2)
plot([1 2 3 4]',var(:,:,2), 'x--', 'LineWidth', 2);
ylabel('\Sigma \sigma_{D*}^2', 'FontSize', 15);
set(gca, 'XTickLabel', bs, 'XTick', 1:numel(bs), 'Fontsize', 12);
grid;

subplot(3,1,3) 
plot([1 2 3 4]',var(:,:,3), 'x--', 'LineWidth', 2);
ylabel('\Sigma \sigma_{f}^2', 'FontSize', 15);
xlabel('Sequences', 'Fontsize', 15);
set(gca, 'XTickLabel', bs, 'XTick', 1:numel(bs), 'Fontsize', 12);
grid;
