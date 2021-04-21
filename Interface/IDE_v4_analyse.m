function varargout = IDE_v4_analyse(varargin)
% IDE_V4_ANALYSE MATLAB code for IDE_v4_analyse.fig
%      IDE_V4_ANALYSE, by itself, creates a new IDE_V4_ANALYSE or raises the existing
%      singleton*.
%
%      H = IDE_V4_ANALYSE returns the handle to a new IDE_V4_ANALYSE or the handle to
%      the existing singleton*.
%
%      IDE_V4_ANALYSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IDE_V4_ANALYSE.M with the given input arguments.
%
%      IDE_V4_ANALYSE('Property','Value',...) creates a new IDE_V4_ANALYSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IDE_v4_analyse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IDE_v4_analyse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IDE_v4_analyse

% Last Modified by GUIDE v2.5 29-Mar-2021 21:23:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IDE_v4_analyse_OpeningFcn, ...
                   'gui_OutputFcn',  @IDE_v4_analyse_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IDE_v4_analyse is made visible.
function IDE_v4_analyse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IDE_v4_analyse (see VARARGIN)

% Choose default command line output for IDE_v4_analyse
handles.output = hObject;

parametric_map_Dperf = getappdata(0,'parametric_map_Dperf');
parametric_map_f = getappdata(0,'parametric_map_f');
[im, c] = showMAP(parametric_map_Dperf.*parametric_map_f, 'fD* [mm²/s]');
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', '< fD* <');


global ROI;  
ROI = {'ROIs'; 'Labels'; 'Shape'}; % Those are the line titles -> ROI = {'ROIs', roi1, roi2, ...;
%                                                                        'Labels', Thalamus, Cortex, ...; 
%                                                                        'Shape', Ellipse, Freehand, ...}

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IDE_v4_analyse wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IDE_v4_analyse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
findobj('tag', 'IDE_v4');
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in ROI_shape.
function ROI_shape_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROI_shape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_shape


% --- Executes during object creation, after setting all properties.
function ROI_shape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in create_roi.
function create_roi_Callback(hObject, eventdata, handles)
% hObject    handle to create_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ROI;

% Selecting ROI shape
ent_ROI_shape = get(handles.ROI_shape, 'String');
num_ROI_shape = get(handles.ROI_shape, 'Value');
text_ROI_shape = ent_ROI_shape(num_ROI_shape);
shape = cell2mat(text_ROI_shape); % shape of the ROI

% Drawing ROI
switch shape
    case 'Freehand'
        roi = imfreehand;
    case 'Polygon'
        roi = impoly;
    case 'Square'
        roi = imrect;
    case 'Ellipse'
        roi = imellipse;
end

% Label and color of ROIs
prompt = {'Enter ROI label', 'Enter label color (red, blue, cyan, green, yellow, magenta)'};
dlg_title = 'ROI label';
answer = inputdlg(prompt,dlg_title,[1 40; 1 40]); 
ROIlab = gtext(answer{1}, 'Color', answer{2});
setColor(roi,answer{2});
ROI{1,end+1} = roi;   
ROI{2,end} = ROIlab;
ROI{3,end} = shape;


% --- Executes on button press in import_roi.
function import_roi_Callback(hObject, eventdata, handles)
% hObject    handle to import_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
[filename, ~, ~] = uigetfile('*.txt', 'Select ROI Protocol file');
fileID = fopen(filename, 'r');    
nROI = 1;
fgetl(fileID); % skip title line
while ~feof(fileID)
    line = fgetl(fileID);
    if ismember('Label', line)            
        label = line(10:end);
    elseif ismember('Shape', line)
        shape = line(8:end);  % get ROI sahpe
        nROI = nROI + 1;      % increment number of ROIs            
        line = fgetl(fileID); % get position and dimension of ROIs            
        roiPos = [];
        lin_roiPos = 1;
        while ~isempty(line) 
            roiPos(lin_roiPos,:) = str2num(line); % store position and
            line = fgetl(fileID); 
            lin_roiPos = lin_roiPos + 1;
        end
        switch shape
            case 'Freehand'
                roi = imfreehand(gca, roiPos);
            case 'Polygon'
                roi = impoly(gca, roiPos);
            case 'Square'
                roi = imrect(gca, roiPos);
            case 'Ellipse'
                roi = imellipse(gca, roiPos);
        end
        prompt = 'Enter label color (red, blue, cyan, green, yellow, black)';
        dlg_title = 'ROI color';
        answer = inputdlg(prompt,dlg_title,[1 40]);            
        ROIlab = gtext(label, 'Color', answer{1});
        setColor(roi,answer{1});
        ROI{1,nROI} = roi;            
        ROI{2,nROI} = ROIlab;
        ROI{3,nROI} = shape;
    else
        continue; % ignore those lines                         
    end         
end
fclose(fileID); 


% --- Executes on button press in export_roi.
function export_roi_Callback(hObject, eventdata, handles)
% hObject    handle to export_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No ROI can be exported'}, 'Error', 'error');
elseif size(ROI,2) == 1
    msgbox({'There is no ROI on the image', 'Please, draw a ROI'}, 'Error', 'error');   
else   
    filename = uiputfile('*.txt','Save ROI protocol');
    fileID = fopen(filename, 'wt+');
    title = 'ROI PROTOCOL';
    fprintf(fileID, '%s\n', title);
    fprintf(fileID, '\n%s\n', '--------- POSITION & DIMENSION ---------');    
    nROI = size(ROI, 2)-1;
    for iROI = 2:nROI+1
        roiPos = round(ROI{1,iROI}.getPosition);
        fprintf(fileID, '\n%s', ['Label ', num2str(iROI-1), ': ', ROI{2, iROI}.String]);       
        fprintf(fileID, '\n%s\n', ['Shape: ' ROI{3,iROI}]);
        for i = 1:size(roiPos,1);            
            for j = 1:size(roiPos,2)
                fprintf(fileID, '%.0f ', roiPos(i,j));
            end
            fprintf(fileID, '%s\n', ' ');
        end
    end
    fprintf(fileID, '\n\n%s ', 'END OF FILE');
    fclose(fileID);
end

% --- Executes on button press in radiobutton_fDperf.
function radiobutton_f_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_fDperf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
parametric_map_f = getappdata(0,'parametric_map_f');
nROI = size(ROI,2)-1;
if size(ROI,2) > 1
    for iROI = 1:nROI
        roiPos{iROI} = ROI{1,iROI+1}.getPosition;        
        shape{iROI} = ROI{3,iROI+1};
        color{iROI,:} = ROI{1,iROI+1}.getColor;
        labelPos{iROI,:} = ROI{2,iROI+1}.Position;
        label{iROI} = ROI{2,iROI+1}.String;
    end
            
    [im, c] = showMAP(parametric_map_f, 'f');
           
    for iROI = 1:nROI        
        switch shape{iROI}
            case 'Freehand'
                roi = imfreehand(gca, roiPos{iROI});
            case 'Polygon'
                roi = impoly(gca, roiPos{iROI});                
            case 'Square'
                roi = imrect(gca, roiPos{iROI});
            case 'Ellipse'
                roi = imellipse(gca, roiPos{iROI});
        end
        setColor(roi, color{iROI,:});
        ROIlab = text('Position', labelPos{iROI,:}, 'Color', color{iROI,:}, 'String', label{iROI});
        ROI{1,iROI+1} = roi;
        ROI{2,iROI+1} = ROIlab;
        ROI{3,iROI+1} = shape{iROI};
    end
else    
    [im, c] = showMAP(parametric_map_f, 'f');
end
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', '< f <');


% --- Executes on button press in radiobutton_D.
function radiobutton_D_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
parametric_map_Ddiff = getappdata(0,'parametric_map_Ddiff');
nROI = size(ROI,2)-1;
if size(ROI,2) > 1
    for iROI = 1:nROI
        roiPos{iROI} = ROI{1,iROI+1}.getPosition;        
        shape{iROI} = ROI{3,iROI+1};
        color{iROI,:} = ROI{1,iROI+1}.getColor;
        labelPos{iROI,:} = ROI{2,iROI+1}.Position;
        label{iROI} = ROI{2,iROI+1}.String;
    end
               
    [im, c] = showMAP(parametric_map_Ddiff, 'D [mm²/s]');   
        
    for iROI = 1:nROI        
        switch shape{iROI}
            case 'Freehand'
                roi = imfreehand(gca, roiPos{iROI});
            case 'Polygon'
                roi = impoly(gca, roiPos{iROI});                
            case 'Square'
                roi = imrect(gca, roiPos{iROI});
            case 'Ellipse'
                roi = imellipse(gca, roiPos{iROI});
        end
        setColor(roi, color{iROI,:});
        ROIlab = text('Position', labelPos{iROI,:}, 'Color', color{iROI,:}, 'String', label{iROI});
        ROI{1,iROI+1} = roi;
        ROI{2,iROI+1} = ROIlab;
        ROI{3,iROI+1} = shape{iROI};
    end
else   
    [im, c] = showMAP(parametric_map_Ddiff, 'D [mm²/s]');
end
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', '< D <');



% --- Executes on button press in radiobutton_Dperf.
function radiobutton_Dperf_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Dperf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
parametric_map_Dperf = getappdata(0,'parametric_map_Dperf');
nROI = size(ROI,2)-1;
if size(ROI,2) > 1
    for iROI = 1:nROI
        roiPos{iROI} = ROI{1,iROI+1}.getPosition;        
        shape{iROI} = ROI{3,iROI+1};
        color{iROI,:} = ROI{1,iROI+1}.getColor;
        labelPos{iROI,:} = ROI{2,iROI+1}.Position;
        label{iROI} = ROI{2,iROI+1}.String;
    end
            
    [im, c] = showMAP(parametric_map_Dperf, 'D* [mm²/s]'); 
        
    for iROI = 1:nROI        
        switch shape{iROI}
            case 'Freehand'
                roi = imfreehand(gca, roiPos{iROI});
            case 'Polygon'
                roi = impoly(gca, roiPos{iROI});                
            case 'Square'
                roi = imrect(gca, roiPos{iROI});
            case 'Ellipse'
                roi = imellipse(gca, roiPos{iROI});
        end
        setColor(roi, color{iROI,:});
        ROIlab = text('Position', labelPos{iROI,:}, 'Color', color{iROI,:}, 'String', label{iROI});
        ROI{1,iROI+1} = roi;
        ROI{2,iROI+1} = ROIlab;
        ROI{3,iROI+1} = shape{iROI};
    end
else    
    [im, c] = showMAP(parametric_map_Dperf, 'D* [mm²/s]');
end
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', '< D* <');



% --- Executes on button press in radiobutton_fDperf.
function radiobutton_fDperf_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_fDperf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
parametric_map_Dperf = getappdata(0,'parametric_map_Dperf');
parametric_map_f = getappdata(0,'parametric_map_f');
nROI = size(ROI,2)-1;
if size(ROI,2) > 1
    for iROI = 1:nROI
        roiPos{iROI} = ROI{1,iROI+1}.getPosition;        
        shape{iROI} = ROI{3,iROI+1};
        color{iROI,:} = ROI{1,iROI+1}.getColor;
        labelPos{iROI,:} = ROI{2,iROI+1}.Position;
        label{iROI} = ROI{2,iROI+1}.String;
    end
            
    [im, c] = showMAP(parametric_map_Dperf.*parametric_map_f, 'fD* [mm²/s]');   
        
    for iROI = 1:nROI        
        switch shape{iROI}
            case 'Freehand'
                roi = imfreehand(gca, roiPos{iROI});
            case 'Polygon'
                roi = impoly(gca, roiPos{iROI});                
            case 'Square'
                roi = imrect(gca, roiPos{iROI});
            case 'Ellipse'
                roi = imellipse(gca, roiPos{iROI});
        end
        setColor(roi, color{iROI,:});
        ROIlab = text('Position', labelPos{iROI,:}, 'Color', color{iROI,:}, 'String', label{iROI});
        ROI{1,iROI+1} = roi;
        ROI{2,iROI+1} = ROIlab;
        ROI{3,iROI+1} = shape{iROI};
    end
else    
    [im, c] = showMAP(parametric_map_Dperf.*parametric_map_f, 'fD* [mm²/s]');    
end
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', '< fD* <');




% --- Executes on button press in histogram.
function histogram_Callback(hObject, eventdata, handles)
% hObject    handle to histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;

if size(ROI,2) > 1
    selected_map = get(handles.choose_map, 'SelectedObject');
    map = get(selected_map, 'String');

    switch map
        case 'f'
            parametric_map = getappdata(0,'parametric_map_f');

        case 'D'
            parametric_map = getappdata(0,'parametric_map_Ddiff');

        case 'D*'
            parametric_map = getappdata(0,'parametric_map_Dperf');

        case 'fD*'
            parametric_map_Dperf = getappdata(0,'parametric_map_Ddiff');
            parametric_map_f = getappdata(0,'parametric_map_f');
            parametric_map = parametric_map_f.*parametric_map_Dperf;
    end

    nROI = size(ROI,2)-1;               % total number of ROIs

    % Number of beams
    nbeams = str2double(get(handles.n_beams, 'String'));
    figure;
    textLegend = cell(1,2*nROI);
    for iROI = 1:nROI
        % Create logical mask for ROI pixels
        roi = ROI{1,iROI+1};
        roiPos = roi.getPosition';
        ROI{2, iROI+1}.Position(1:2) = roiPos(1:2)-3; % ROI label is char, there is no position
        BW = createMask(roi);                         % Makes a mask from ROI

        Obj = (parametric_map(:,:,1) > 0);
        [s,q] = find(BW & Obj);                                   % Select Signal within the ROI

        % Calculate mean parameter
        ROIdim = length(s);
        values = zeros(1,ROIdim);
        for i = 1:ROIdim
            values(i) = parametric_map(s(i),q(i));            
        end
        
%         figure('Name', ['Histogram of ' ROI{2,iROI+1}.String]); 
        
        histogram(values, 'BinLimits', [min(values) max(values)], 'BinWidth', abs(min(values)-max(values))/nbeams);
        if strcmp(map,'D')
            xlabel('D [mm²/s]', 'FontSize', 14);
        elseif strcmp(map,'D*')
            xlabel('D* [mm²/s]', 'FontSize', 14);
        elseif strcmp(map,'fD*')
            xlabel('fD* [mm²/s]', 'FontSize', 14);
        else
            xlabel('f', 'FontSize', 14);
        end                
        ylabel(['number of ocurrencies of ' map], 'FontSize', 14);    
%         title(['Histogram of ' map ' in ' ROI{2,iROI+1}.String])
        title(['Histogram of ' map]);
        set(gca, 'FontSize', 14);
        y_limit = ylim;
        if iROI > 1
            if y_limit(2) > h.YData(2)
                h.YData(2) = y_limit(2);
            end
        end
        textLegend{2*iROI - 1} = ROI{2,iROI+1}.String;        
        textLegend{2*iROI} = [ROI{2,iROI+1}.String ' mean'];
        h = line([mean(values), mean(values)], ylim, 'LineWidth', 2, 'Color', roi.getColor);
        annotation('textbox', [0.66,0.7,0.1,0.1],'String', {[ROI{2,iROI+1}.String], ['mean = ', num2str(mean(values))], ['std = \pm', num2str(std(values))]}, 'FontSize', 14);
        hold on;
    end
    legend(textLegend, 'FontSize', 12, 'Location', 'best');
    hold off;
else
    msgbox({'There is no ROI on the map', 'Please, draw a ROI or import it'}, 'Error', 'error');
end

% --- Executes on button press in mean_parameters.
function mean_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to mean_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;

if size(ROI,2) == 1
    msgbox({'There is no ROI on the map', 'Please, draw a ROI or import it'}, 'Error', 'error');
else
    selected_map = get(handles.choose_map, 'SelectedObject');
    map = get(selected_map, 'String');

    switch map
        case 'f'
            parametric_map = getappdata(0,'parametric_map_f');
            line = 1;
        case 'D'
            parametric_map = getappdata(0,'parametric_map_Ddiff');
            line = 2;
        case 'D*'
            parametric_map = getappdata(0,'parametric_map_Dperf');
            line = 3;
        case 'fD*'
            parametric_map_Dperf = getappdata(0,'parametric_map_Ddiff');
            parametric_map_f = getappdata(0,'parametric_map_f');
            parametric_map = parametric_map_f.*parametric_map_Dperf;
            line = 4;
    end

    nROI = size(ROI,2)-1;               % total number of ROIs
    avrg_parameter = cell(line,nROI);
    name_column = cell(1,nROI);

    for iROI = 1:nROI

        % Create logical mask for ROI pixels
        roi = ROI{1,iROI+1};
        roiPos = roi.getPosition';
        ROI{2, iROI+1}.Position(1:2) = roiPos(1:2)-3; % ROI label is char, there is no position
        BW = createMask(roi);                         % Makes a mask from ROI

        Obj = (parametric_map(:,:,1) > 0);
        [s,q] = find(BW & Obj);                       % Select Signal within the ROI

        % Calculate mean parameter
        ROIdim = length(s);
        values = zeros(1,ROIdim);
        for i = 1:ROIdim
            values(i) = parametric_map(s(i),q(i));          
        end
        str = [num2str(mean(values)), ' ', char(177), ' ', num2str(std(values))]; % mean +/- std
%         str = num2str(mean(values)); 
        avrg_parameter(line, iROI) = {str};
        set(handles.parameters, 'data', avrg_parameter);        
        name_column{iROI} = ROI{2,iROI+1}.String; % set the name of the column: 1, 2, 3, ...
    end
    set(handles.parameters, 'ColumnName', name_column'); % set column names
end
        
        


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
ROI = {'ROIs'; 'Labels'; 'Shape'};
filedir = uigetdir('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE', 'Select directory');
filename = uigetfile([filedir '\*.mat'], 'Select data');
load([filedir '\' filename], 'parametric_map_Dperf', 'parametric_map_f', 'parametric_map_Ddiff');
setappdata(0,'parametric_map_Dperf', parametric_map_Dperf);
setappdata(0,'parametric_map_f', parametric_map_f);
setappdata(0,'parametric_map_Ddiff', parametric_map_Ddiff);

selected_map = get(handles.choose_map, 'SelectedObject');
map = get(selected_map, 'String');

switch map
    case 'f'
        parametric_map = parametric_map_f;
        label = 'f';        
    case 'D'
        parametric_map = parametric_map_Ddiff;
        label = 'D [mm²/s]';       
    case 'D*'
        parametric_map = parametric_map_Dperf;
        label = 'D* [mm²/s]';        
    case 'fD*'
        parametric_map = parametric_map_f.*parametric_map_Dperf;
        label = 'fD* [mm²/s]';        
end
[im, c] = showMAP(parametric_map, label);
set(handles.lower_limit, 'UserData', {c, im});
set(handles.upper_limit, 'UserData', {c, im});
set(handles.map_shown, 'String', ['<' map '<']);



% --- Executes on button press in change_label.
function change_label_Callback(hObject, eventdata, handles)
% hObject    handle to change_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in delete_rois.
function delete_rois_Callback(hObject, eventdata, handles)
% hObject    handle to delete_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if size(ROI, 2) == 1
    msgbox({'There is no ROI on the image', 'Please, draw a ROI'}, 'Error', 'error');
elseif size(ROI, 2) == 2
    roi = ROI{1,2};
    label = ROI{2,2};
    roi.delete;
    delete(label);
    ROI = {'ROIs'; 'Label'; 'Shape'};
else
    toDelete = get(handles.uibuttongroup2, 'SelectedObject');
    toDelete = get(toDelete, 'String');
    if strcmp(toDelete, 'All')
        nROI = size(ROI,2)-1;
        for iROI = 2:nROI+1
            roi = ROI{1,iROI};
            label = ROI{2,iROI};
            roi.delete;
            delete(label);
        end
        ROI = {'ROIs'; 'Label'; 'Shape'};
    else
        prompt = {'Enter label of the ROI to be deleted'};
        dlg_title = 'Deleted ROIs';
        answer = inputdlg(prompt,dlg_title,[1 40]);
        nROI = size(ROI,2)-1;
        for iROI = 1:nROI
            if strcmp(answer{1}, ROI{2,iROI+1}.String)
                roi = ROI{1,iROI+1};
                label = ROI{2,iROI+1};
                roi.delete;
                delete(label);
                ROI(:, iROI+1) = [];                
                break;
            end
        end
    end
end


function n_beams_Callback(hObject, eventdata, handles)
% hObject    handle to n_beams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_beams as text
%        str2double(get(hObject,'String')) returns contents of n_beams as a double


% --- Executes during object creation, after setting all properties.
function n_beams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_beams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clone_roi.
function clone_roi_Callback(hObject, eventdata, handles)
% hObject    handle to clone_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;

if strcmp(get(handles.axes1, 'Visible'), 'off') % If there is no image
    
    msgbox({'No .PAR image selected', 'No ROI can be cloned'}, 'Error', 'error');
    
elseif size(ROI,2) == 1 % If there is no ROI
    
    msgbox({'No ROI was detected', 'No ROI can be cloned'}, 'Error', 'error');
    
elseif size(ROI,2) == 2 % If there is only one ROI
    
    % New label
    prompt = {'Write a new label'};
    dlg_title = 'Clone ROI';
    answer = inputdlg(prompt,dlg_title,[1 40]);
    
    % Clone
    roi = ROI{1,2};
    shape = ROI{3,2};
    color = roi.getColor;
    switch shape
        case 'Freehand'                    
            roi = imfreehand(gca, roi.getPosition + 50);
        case 'Ellipse'
            roiPos = roi.getPosition;
            roi = imellipse(gca, [roiPos(1)+50 roiPos(2)+50 roiPos(3) roiPos(4)]);
        case 'Polygon'                    
            roi = impoly(gca, roi.getPosition + 50);                    
        case 'Square'
            roiPos = roi.getPosition;
            roi = imrect(gca, [roiPos(1)+50 roiPos(2)+50 roiPos(3) roiPos(4)]);
    end
    setColor(roi, color);
    ROIlab = gtext(answer{1}, 'Color', color);
    
    % ROI record update
    ROI{1, 3} = roi;
    ROI{2, 3} = ROIlab;
    ROI{3, 3} = shape;    
    
else % If there are more ROIs
    prompt = {'Enter label of ROI to be cloned', 'Write a new label'};
    dlg_title = 'Clone ROI';
    answer = inputdlg(prompt,dlg_title,[1 40; 1 40]);
    nROI = size(ROI,2)-1;
    found = false;
    for iROI = 1:nROI
        if strcmp(answer{1}, ROI{2,iROI+1}.String)        
            roi = ROI{1,iROI+1};
            shape = ROI{3,iROI+1};
            color = roi.getColor;
            switch shape
                case 'Freehand'                    
                    roi = imfreehand(gca, roi.getPosition + 50);
                case 'Ellipse'
                    roiPos = roi.getPosition;
                    roi = imellipse(gca, [roiPos(1)+50 roiPos(2)+50 roiPos(3) roiPos(4)]);
                case 'Polygon'                    
                    roi = impoly(gca, roi.getPosition + 50);                    
                case 'Square'
                    roiPos = roi.getPosition;
                    roi = imrect(gca, [roiPos(1)+50 roiPos(2)+50 roiPos(3) roiPos(4)]);
            end
            setColor(roi, color);
            ROIlab = gtext(answer{2}, 'Color', color);
            found = true;
            break;
        end
    end
    if found
        ROI{1, end+1} = roi;
        ROI{2, end} = ROIlab;
        ROI{3, end} = shape;
    else
        msgbox({'No ROI with such label was detected', 'No ROI can be cloned'}, 'Error', 'error');
    end
end



function lower_limit_Callback(hObject, eventdata, handles)
% hObject    handle to lower_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lower_limit as text
%        str2double(get(hObject,'String')) returns contents of lower_limit as a double
data = get(hObject, 'UserData');
c = data{1};
im = data{2};
low = str2double(get(hObject, 'String'));
if low < min(im.CData(:))
    msgbox({'Lower limit is lower than the minimum value in image'}, 'Error', 'error');
else
    c.Limits(1) = low;
    idx = im.CData < low;
    im.CData(idx) = low;
end


% --- Executes during object creation, after setting all properties.
function lower_limit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lower_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upper_limit_Callback(hObject, eventdata, handles)
% hObject    handle to upper_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upper_limit as text
%        str2double(get(hObject,'String')) returns contents of upper_limit as a double
data = get(hObject, 'UserData');
c = data{1};
im = data{2};
upp = str2double(get(hObject, 'String'));
if upp > max(im.CData(:))
    msgbox({'Upper limit is higher than the minimum value in image'}, 'Error', 'error');
else
    c.Limits(2) = upp;
    idx = im.CData > upp;
    im.CData(idx) = upp;
end


% --- Executes during object creation, after setting all properties.
function upper_limit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upper_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function map_shown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to map_shown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in choose_map.
function choose_map_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in choose_map 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)   


% --------------------------------------------------------------------
function save_map_Callback(hObject, eventdata, handles)
% hObject    handle to save_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selected_map = get(handles.choose_map, 'SelectedObject');
map = get(selected_map, 'String');        

old_ax = handles.axes1; 
figure;
new_ax = gca;
copyobj(get(old_ax, 'children'), new_ax);

c = colorbar;
switch map
    case 'f'
        c.Label.String = 'f';
    case 'D'
        c.Label.String = 'D [mm²/s]';
    case 'D*'
        c.Label.String = 'D* [mm²/s]';
    case 'fD*'
        c.Label.String = 'fD* [mm²/s]';
end
c.Label.FontSize = 14;
axis ij image;
set(gca, 'XTick', [], 'YTick', []);
