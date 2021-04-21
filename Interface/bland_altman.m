function [data_mean,data_diff,md,sd] = bland_altman(data1,data2)
% Function to generate Bland Altman plots. Barry Greene, September 2008
% Bland, J.M., Altman, D.G. 'Statistical methods for assessing agreement ...
% between two methods of clinical measurement'(1986) Lancet, 1 (8476), pp. 307-310.
% Inputs: data1: Data from first instrument
%         data2: Data from second instument  
% Produces Bland Altman plot with mean difference and mean difference +/-
% 2*SD difference lines.

[m,n] = size(data1);
if(n>m)
    data1 = data1';
end
if(size(data1)~=size(data2))
    error('Data matrices must be the same size');
end
figure;
text_title = {'f cerebellum mean signal', 'D cerebellum mean signal', 'D* cerebellum mean signal', 'f cerebellum parametric mean',...
    'D cerebellum parametric mean', 'D* cerebellum parametric mean'};
for i = 1:size(data1,2)
    data_mean = mean([data1(:,i),data2(:,i)],2);            % Mean of values from each instrument 
    data_diff = data1(:,i) - data2(:,i);                    % Difference between data from each instrument
    [h,p,~,~] = ttest2(data1(:,i),data2(:,i));              % t-test
    disp(h);disp(p);
    md = mean(data_diff);                                   % Mean of difference between instruments 
    sd = std(data_diff);                                    % Std dev of difference between instruments    
    ax = subplot(2,3,i);    
    plot(data_mean,data_diff,'ob','MarkerSize',8,'LineWidth',2);   % Bland Altman plot
    xlim([min(data_mean) max(data_mean)]);    
    YL = get(ax, 'YLim');
    maxlim = max(abs(YL));
    set(ax, 'YLim', [-maxlim maxlim]);
    
    if i == 1
        hold on ;
    end; 
    
    line(xlim,[md, md], 'Color', 'black', 'LineStyle', '--');                                % Mean difference line
    line(xlim,[1.96*sd, 1.96*sd], 'Color', 'red', 'LineStyle', '--');                        % Mean -1.96*SD line  
    line(xlim,[-1.96*sd, -1.96*sd], 'Color', 'red', 'LineStyle', '--');                      % Mean +1.96*SD line    
    gtext(num2str(md));
    gtext(num2str(1.96*sd), 'Color', 'red');
    gtext(num2str(-1.96*sd), 'Color', 'red');
    gtext(['P = ', num2str(p)], 'Color', 'blue');
    
    if i == 1
        legend('Differnces', 'Mean', '\pm1.96SD', 'Location', 'best');        
    end

    grid on;
    title(text_title{i},'FontSize',10);
    xlabel('Mean of two measures','FontSize',10);
    ylabel('Difference between two measures','FontSize',10);
end
hold off