function [parametric_map_f, parametric_map_Ddiff, parametric_map_Dperf] = calculateMAP(img, bvalues, slice, fit, outlierList)
    
    parametric_map_f = zeros(size(img,1),size(img,2));
    parametric_map_Ddiff = parametric_map_f;
    parametric_map_Dperf = parametric_map_f;

%         clust = parcluster();                     % cluster object
%         clust.NumWorkers = 6;                     % number of workers to receive a task
%         parpool(6);                               % number of parallel pools

    % Creates waitbar to monitor fitting progress through the image
    h = waitbar(0,'1','Name','Voxel by Voxel fitting',...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');            
    setappdata(h,'canceling',0)
    cancel = false;

    % Voxelwise fitting
    for col = 1:size(img,2)
        for row = 1:size(img,1)             % parallel execution of column fitting                               
            Signal = img(row,col,:);           % IVIM signal          
            Signal = Signal(:)./Signal(1,1,1); % S(b)/S0
            excludeNaN = isnan(Signal);        % exclude NaN from array
            excludeInf = isinf(Signal);        % exclude Inf from array
            Signal(excludeNaN) = 0;
            Signal(excludeInf) = 0;
            [parametric_map_Ddiff(row,col), parametric_map_Dperf(row,col), parametric_map_f(row,col)] =...
                fitting(fit, 'Tikh', Signal, bvalues, 'talwar', [0.15 1e-2 1e-3]);                 
        end            

        % Check for Cancel button press
        if getappdata(h,'canceling')
            delete(h);
            cancel = true;
            break;
        end

        % Report current estimate in the waitbar's message field
        waitbar(col/size(img,2), h, sprintf('%d columns out of %d were tracked in the image',col, size(img,2)))                       

    end
%         delete(gcp)     % Shut down workers
    delete(h)       % DELETE the waitbar; don't try to CLOSE it       

    if ~cancel

        % outliers treatment
        f_min = outlierList(1);f_max = outlierList(2);
        D_min = outlierList(3);D_max = outlierList(4);
        Dperf_min = outlierList(5);Dperf_max = outlierList(6);
        idx = (parametric_map_Ddiff > D_max | parametric_map_Ddiff < D_min); parametric_map_Ddiff(idx) = 0;
        idx = (parametric_map_Dperf > Dperf_max | parametric_map_Dperf < Dperf_min); parametric_map_Dperf(idx) = 0;
        idx = (parametric_map_f > f_max | parametric_map_f < f_min); parametric_map_f(idx) = 0;

        % show maps
        figure('Name', 'Parametric map');
        subplot(2,2,1);
        image(parametric_map_f,'CDataMapping','scaled');
        c = colorbar;
        c.Label.String = 'f';
        c.Label.FontSize = 14;
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'FontSize', 12);

        subplot(2,2,2);
        image(parametric_map_Ddiff,'CDataMapping','scaled');
        c = colorbar;
        c.Label.String = 'D (mm^{2}/s)';
        c.Label.FontSize = 14;
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'FontSize', 12);

        subplot(2,2,3);
        image(parametric_map_Dperf,'CDataMapping','scaled');
        c = colorbar;
        c.Label.String = 'D^{*} (mm^{2}/s)';
        c.Label.FontSize = 14;
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'FontSize', 12);

        subplot(2,2,4);
        image(parametric_map_Dperf.*parametric_map_f,'CDataMapping','scaled');
        c = colorbar;
        c.Label.String = 'fD^{*} (mm^{2}/s)';
        c.Label.FontSize = 14;
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'FontSize', 12);        

        % Save pixel values?
        choice = questdlg('Would you like to save parametric data?', ...
                          'Save pixel values', ...
                          'Yes','No','No');
        % Handle response
        if strcmp(choice,'Yes')
            filedir = uigetdir('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE', 'Select directory');
            filename = uiputfile('*.mat','Save Variable As');
            save([filedir '\' filename], 'parametric_map_f', 'parametric_map_Ddiff', 'parametric_map_Dperf',...
                'f_min', 'f_max', 'D_min', 'D_max', 'Dperf_min', 'Dperf_max', 'fit', 'slice');                          
        end
    end
end
