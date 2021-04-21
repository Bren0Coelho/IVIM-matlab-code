%% plot_sinais
%% Plota gráficos de sinais com parâmetros calculados em voxel.m:
%% - Com ruído
%% - Sem ruído
%% - Fitting


%% Carrega variáveis e resultados da simulação

close all;
clear;
clc;

load('workspace_simulation_fernando_1000_v3');

%% Entradas


est = 1;        % estrutura que guarda o modelo analisado (s1, s2 ou s3)
                % em s1, f varia; D_diff e D_perf são fixos
                % em s2, D_diff varia; f e D_perf são fixos
                % em s3, D_perf varia; f e D_diff são fixos

i = 3;          % variável que guarda o índice do parâmetro a ser analisado nos gráficos
                % Se est = 1 --> 1 <= i <= 6
                % Se est = 2 ou 3 --> 1 <= i <= 13
                % Esta variável deve ser coerente com a variável est!!
                
data = R_b1_s1; % estrutura que guarda os dados analsados segundo modelo e a sequência dos valores de b
                % Se est = 1 --> data = R_b1_s1 ou R_b2_s1 ou R_b3_s1 ou R_b4_s1
                % Se est = 2 --> data = R_b1_s2 ou R_b2_s2 ou R_b3_s2 ou R_b4_s2
                % Se est = 3 --> data = R_b1_s3 ou R_b2_s3 ou R_b3_s3 ou R_b4_s3
                % Esta variável deve ser coerente com a variável est!!
                
b = b1;         % variável que guarda a sequência de b para uma eventual demonstração de sinais
                % Se data = R_b1... --> b = b1
                % Se data = R_b2... --> b = b2
                % Se data = R_b3... --> b = b3
                % Se data = R_b4... --> b = b4
                % Esta variável deve ser coerente com a variável data!!
                
snr = 1;

c = randi([1 1000], 1, 10); % Gera 10 valores entre 1 e 1000 que serão os índices da matriz de sinais a serem plotados

%% Mostrar sinais


for n = 1:length(c)
    
  sin = c(n);

  if (est == 1)         

      M = s1.f(i)*exp(-b*s1.D_perf)+(1-s1.f(i))*exp(-b*s1.D_diff);
  else
      
      if (est == 2)
          
          M = s2.f*exp(-b*s2.D_perf)+(1-s2.f)*exp(-b*s2.D_diff(i));
          
      else
          
          M = s3.f*exp(-b*s3.D_perf(i))+(1-s3.f)*exp(-b*s3.D_diff);
          
      end
  end
      

%   f_fit_lls = data(i, 1).pontos(sin,3,snr) ; D_diff_fit_lls = data(i, 1).pontos(sin,1,snr); D_perf_fit_lls = data(i, 1).pontos(sin,2,snr);
%   f_fit_nlls = data(i, 2).pontos(sin,3,snr) ; D_diff_fit_nlls = data(i, 2).pontos(sin,1,snr); D_perf_fit_nlls = data(i, 2).pontos(sin,2,snr);
  f_fit_lev = data(i, 3).pontos(sin,3,snr) ; D_diff_fit_lev = data(i, 3).pontos(sin,1,snr); D_perf_fit_lev = data(i, 3).pontos(sin,2,snr);
%   f_fit_nnls = data(i, 4).pontos(sin,3,snr) ; D_diff_fit_nnls = data(i, 4).pontos(sin,1,snr); D_perf_fit_nnls = data(i, 4).pontos(sin,2,snr);

%   M_fit_lls = f_fit_lls*exp(-b*D_perf_fit_lls)+(1-f_fit_lls)*exp(-b*D_diff_fit_lls);
  
%   M_fit_nlls = f_fit_nlls*exp(-b*D_perf_fit_nlls)+(1-f_fit_nlls)*exp(-b*D_diff_fit_nlls);
  M_fit_lev = f_fit_lev*exp(-b*D_perf_fit_lev)+(1-f_fit_lev)*exp(-b*D_diff_fit_lev);
%   M_fit_nnls = f_fit_nnls*exp(-b*D_perf_fit_nnls)+(1-f_fit_nnls)*exp(-b*D_diff_fit_nnls);

  figure(sin);
%   set(gcf, 'Position', get(0, 'Screensize'));
%   subplot(2,2,1);
%   plot(b, M_fit_lls, 'rx-', b, data(i, 1).Sinal(sin,:,snr), 'ko-', b, M, 'b^-');
%   title(['Sinais ruidoso, limpo e fitting LLS; SNR = ', num2str(SNR(snr))]);
%   xlabel('b [s²/mm]');
%   ylabel('ln(Sinal)'); 
%   legend('Fitting', 'Com ruído', 'Sem ruído');
%   grid;

%   subplot(2,2,2);
%   plot(b, log(M_fit_nlls), 'rx-', b, log(data(i, 2).Sinal(sin,:,snr)), 'ko-', b, log(M), 'b^-');
%   title(['Sinais ruidoso, limpo e fitting NLLS; SNR = ', num2str(SNR(snr))]);
%   xlabel('b [s²/mm]');
%   ylabel('ln(Sinal)'); 
%   legend('Fitting', 'Com ruído', 'Sem ruído');
%   grid;
% 
%   subplot(2,2,3);
  plot(b, data(i, 3).Sinal(sin,:,snr), 'ko-', b, M, 'b^-');
  title(['Sinais ruidoso, limpo e fitting LEV; SNR = ', num2str(SNR(snr))]);
  xlabel('b [s²/mm]');
  ylabel('Sinal'); 
  legend('Com ruído', 'Sem ruído');
  grid;
% 
%   subplot(2,2,4);
%   plot(b, log(M_fit_nnls), 'rx-', b, log(data(i, 4).Sinal(sin,:,snr)), 'ko-', b, log(M), 'b^-');
%   title(['Sinais ruidoso, limpo e fitting NNLS; SNR = ', num2str(SNR(snr))]);
%   xlabel('b [s²/mm]');
%   ylabel('ln(Sinal)'); 
%   legend('Fitting', 'Com ruído', 'Sem ruído');
%   grid;
  
%   saveas(H, ['Sinais_', num2str(est),'_', num2str(a), '_', num2str(sin)], 'png')
  
end

