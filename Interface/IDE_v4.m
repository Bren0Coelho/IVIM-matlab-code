function varargout = IDE_v4(varargin)
% IDE_V4 MATLAB code for IDE_v4.fig
%      IDE_V4, by itself, creates a new IDE_V4 or raises the existing
%      singleton*.
%
%      H = IDE_V4 returns the handle to a new IDE_V4 or the handle to
%      the existing singleton*.
%
%      IDE_V4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IDE_V4.M with the given input arguments.
%
%      IDE_V4('Property','Value',...) creates a new IDE_V4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IDE_v4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IDE_v4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% =======================================================================
% =======================================================================
% Version Control
%
% Breno Spinelli Coelho @ CIERMag/IFSC/USP
%   Last Modified by GUIDE v2.5 25-Feb-2021 10:23:55
% ========================================================================
% ========================================================================
% IDE_v4 is a simulation interface to estimate diffusion and perfusion
% parameters (f, D and D*) within regions of interest, SNR values from MRI 
% .PAR/.REC v4.2 archives and parametric maps of IVIM parameters
%
%
% The use of the interface or code is allowed under conditions of
% acknowledgement: COELHO, B S; PAIVA, F F; Desenvolvimento e implementação
% de métodos de processamento de imagens por RM para obtenção de parâmetros
% relacionados à perfusão e à difusão. 2021 - Instituto de Física de São
% Carlos, Universidade de São Paulo



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IDE_v4_OpeningFcn, ...
                   'gui_OutputFcn',  @IDE_v4_OutputFcn, ...
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


% --- Executes just before IDE_v4 is made visible.
function IDE_v4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IDE_v4 (see VARARGIN)

% Choose default command line output for IDE_v4
handles.output = hObject;
clc;

% That is the global variable ROI (cell), which stores the ROIs as objects
% (1st line), texts (2nd line) and chars (3rd line)
global ROI;  
ROI = {'ROIs'; 'Labels'; 'Shape'}; % Those are the line titles -> ROI = {'ROIs', roi1, roi2, ...;
%                                                                        'Labels', Thalamus, Cortex, ...; 
%                                                                        'Shape', Ellipse, Freehand, ...}

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IDE_v4 wait for user response (see UIRESUME)
% uiwait(handles.IDE_v4);


% --- Outputs from this function are returned to the command line.
function varargout = IDE_v4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider of slice selection movement.
function slice_selection_Callback(hObject, eventdata, handles)
% hObject    handle to slice_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global ROI;
ROI = {'ROIs'; 'Labels'; 'Shape'};
% Slice selection through the slider
data = get(hObject, 'UserData');
MR_b0 = data{1};                          % image when b = 0 s/mm²
header = data{2};                         % archive header
slice = floor(get(hObject, 'Value') + 1); % slice 
set(handles.slice_selected, 'String', num2str(slice)); 

% Catching index of b-value
nbvalues = header.tbl(:, header.tblcols.diff_b_fac); 
bvalues = unique(nbvalues);                          % b-values
border = floor(get(handles.b_selection, 'Value')+1); % index of b-value whose image is shown on display

% Change image on the screen
MRI = mat2gray(MR_b0(:,:,slice,:)); % grayscale
handFig = montage(MRI(:,:,border), 'size', [1 1]);
set(handles.axes1, 'Visible', 'On');
set(gca, 'XTick', [], 'YTick', []);
set(handles.Image, 'UserData', handFig);

% set handles data
set(handles.roi, 'UserData', {MR_b0, header, handFig});                  % ROI data
set(handles.histogram, 'UserData', MRI(:,:,border));                     % Histogram data
set(handles.snr_roi, 'UserData', MR_b0(:,:,slice,border));               % SNR data
set(handles.param_map, 'UserData', {MR_b0(:,:,slice,:), bvalues, slice});% Parametric map data
set(handles.complement, 'UserData', MR_b0(:,:,slice,border));            % Complementary image data
set(handles.adjust_histogram, 'UserData', handFig);                      % Histogram adjustment data




% --- Executes during object creation, after setting all properties.
function slice_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider of b-value selection movement.
function b_selection_Callback(hObject, eventdata, handles)
% hObject    handle to b_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Slice selection through the slider
global ROI;
data = get(hObject, 'UserData');
MR_b0 = data{1};                                          % image when b = 0 s/mm²
header = data{2};                                         % archive header
slice = floor(get(handles.slice_selection, 'Value') + 1); % slice 

nbvalues = header.tbl(:, header.tblcols.diff_b_fac);
bvalues = unique(nbvalues);                               % b-values
border = floor(get(hObject, 'Value')+1);                  % index of b-value whose image is shown on display
bvalue = bvalues(border);                                 % b-value onto display
nROI = size(ROI,2)-1;

% Keep ROIs on image while b-values change
if size(ROI,2) > 1
    for iROI = 1:nROI
        roiPos{iROI} = ROI{1,iROI+1}.getPosition;        
        shape{iROI} = ROI{3,iROI+1};
        color{iROI,:} = ROI{1,iROI+1}.getColor;
        labelPos{iROI,:} = ROI{2,iROI+1}.Position;
        label{iROI} = ROI{2,iROI+1}.String;
    end
            
    % Change image on the screen
    MRI = mat2gray(MR_b0(:,:,slice,:));                       % grayscale
    handFig = montage(MRI(:,:,border), 'size', [1 1]);        
    set(handles.axes1, 'Visible', 'On');
    set(gca, 'XTick', [], 'YTick', []);
    set(handles.Image, 'UserData', handFig);
        
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
    % Change image on the screen
    MRI = mat2gray(MR_b0(:,:,slice,:));                       % grayscale
    handFig = montage(MRI(:,:,border), 'size', [1 1]);        
    set(handles.axes1, 'Visible', 'On');
    set(gca, 'XTick', [], 'YTick', []);
    set(handles.Image, 'UserData', handFig);
end

% set handles data
set(handles.b_selected, 'String', num2str(bvalue));             % selected b-value
set(handles.roi, 'UserData', {MR_b0, header, handFig});         % ROI data
set(handles.histogram, 'UserData', MRI(:,:,border));            % Histogram data    
set(handles.snr_roi, 'UserData', MR_b0(:,:,slice,border));      % SNR data
set(handles.complement, 'UserData', MR_b0(:,:,slice,border));   % Complementary image data
set(handles.adjust_histogram, 'UserData', handFig);             % Histogram adjustment data 



% --- Executes during object creation, after setting all properties.
function b_selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function slice_selected_Callback(hObject, eventdata, handles)
% hObject    handle to slice_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slice_selected as text
%        str2double(get(hObject,'String')) returns contents of slice_selected as a double


% --- Executes during object creation, after setting all properties.
function slice_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b_selected_Callback(hObject, eventdata, handles)
% hObject    handle to b_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_selected as text
%        str2double(get(hObject,'String')) returns contents of b_selected as a double

% --- Executes during object creation, after setting all properties.
function b_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ROI_ROI_shape.
function ROI_shape_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_ROI_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROI_ROI_shape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROI_ROI_shape


% --- Executes during object creation, after setting all properties.
function ROI_shape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROI_ROI_shape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button Create ROI press.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get image data

global ROI; % global variable ROI
if strcmp(get(handles.axes1, 'Visible'), 'off') % verifies if there is any image on display
    msgbox({'No .PAR image selected', 'No ROI can be drawn'}, 'Error', 'error');
else
    
    % catching UserData
    data = get(hObject, 'UserData'); 
    img = data{1};                                              % image pixel values
    header = data{2};                                           % archive header
    slice = str2double(get(handles.slice_selected, 'String'));  % slice selected

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

    % Set IVIM signal
    set(handles.IVIM_signal, 'UserData', {header, img(:,:,slice,:), cell2mat(text_ROI_shape)});
end




% --- Executes on button press in IVIM_signal.
function IVIM_signal_Callback(hObject, eventdata, handles)
% hObject    handle to IVIM_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off') % verifies if there is any image on display
    
    msgbox({'No .PAR image selected', 'No signal can be caught'}, 'Error', 'error');
    
elseif size(ROI,2) == 1 % when no ROIs are drawn
    
    msgbox({'There is no ROI on the image', 'Please, draw a ROI'}, 'Error', 'error');
    
else
    
    % Catching data from handle
    data = get(hObject, 'UserData');
    header = data{1};
    img = data{2};

    nbvalues = header.diffvalues;
    nROI = size(ROI,2)-1;               % total number of ROIs
    MR_b0 = img(:,:,1);                 % image when b = 0 s/mm²
    M_noise = zeros(nROI, nbvalues);    % Signal array

    % State b-values from header
    bvalues = header.tbl(:, header.tblcols.diff_b_fac);
    b = unique(bvalues)';


    for iROI = 1:nROI

        % Create logical mask for ROI pixels
        roi = ROI{1,iROI+1};
        roiPos = roi.getPosition';
        ROI{2, iROI+1}.Position(1:2) = roiPos(1:2)-3; % ROI label is char, there is no position
        BW = createMask(roi);                         % Makes a mask from ROI
        Obj = (mat2gray(img(:,:,1)) > 0);
        [s,q] = find(BW & Obj);                       % Select Signal within the ROI

        % Create IVIM signal
        ROIdim = length(s);
        Signal = zeros(ROIdim, nbvalues);
        for z = 1:nbvalues
            MR_b = img(:,:,z);
            for i = 1:ROIdim
                Sb = MR_b(s(i),q(i));
                S0 = MR_b0(s(i),q(i));
                Signal(i,z) = Sb/S0;
            end
        end
        M_noise(iROI,:) = mean(Signal);
        
    end

    % legend text
    for iROI = 1:nROI
        text{iROI} = ROI{2,iROI+1}.String;
    end
    
    % plot    
    figure('Name', 'IVIM signal'); plot(b, M_noise, 'o', 'Linewidth', 2);grid;legend(text, 'Location', 'best');
    title('IVIM mean signal in the ROIs');
    ylabel('S/S_{0}');
    xlabel('b [s/mm²]');
    
    % set data to estimations
    set(handles.estimate, 'UserData', {M_noise, b, header});
end

% --- Executes on button press in b_values.
function b_values_Callback(hObject, eventdata, handles)
% hObject    handle to b_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No b-value can be shown'}, 'Error', 'error');
else
    b = get(hObject, 'UserData');
    msgbox(num2str(b'), 'b-values [s/mm²]', 'help');
end

% --- Executes on button press in Histogram.
function histogram_Callback(hObject, eventdata, handles)
% hObject    handle to histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No histogram can be estimated'}, 'Error', 'error');
else
    figure('Name', 'Histogram of image'); imhist(get(hObject, 'UserData'));
    xlabel('Intensity[0 - 1]');
    ylabel('Number of pixels');
end

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in File > Open.
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;

% Catch file name
[filename, ~, ~] = uigetfile('Dados de IVIM\*.PAR', 'Select PAR file');
set(hObject, 'UserData', filename);

% Read file
[img, header] = readrec_V4_2(['Dados de IVIM\' filename]);
img = squeeze(img); % removes singleton dimensions

% Pixel Values
FP = double(img);                                % floating point value
RS = header.tbl(1,header.tblcols.rescale_slope); % rescale slope
RI = header.tbl(1,header.tblcols.rescale_int);   % rescale intercept
SS = header.tbl(1,header.tblcols.scale_slope);   % scale slope
DV = FP.*(RS*SS);                                % displayed value on console 
PV = (DV-RI)./RS;                                % pixel value in REC file

% Set general informations
te = header.tbl(1, end-18); % Echo time
nbvalues = header.tbl(:, header.tblcols.diff_b_fac); 
bvalues = unique(nbvalues);% b-values
info = {header.patname ; header.fov(1); header.fov(2); header.fov(3); header.nslices; header.tr; te; header.datetime; header.diffvalues};
set(handles.general_information, 'data', info); % set information to table

% Slice panel
set(handles.slice_selected, 'UserData', {PV, header});
set(handles.slice_selected, 'String', '1');

% Show image
MR_gray = mat2gray(PV(:,:,1,1));
handFig = montage(MR_gray, 'size', [1 1], 'Indices', 1);
set(handles.axes1, 'Visible', 'On');
set(gca, 'XTick', [], 'YTick', []);

% Select Slice
theRangeSlices = header.nslices - 1;
stepsSlice = [1/theRangeSlices, 1/theRangeSlices];
set(handles.slice_selection, 'UserData', {PV, header});
set(handles.slice_selection, 'Value', 0);
set(handles.slice_selection, 'SliderStep', stepsSlice);
set(handles.slice_selection, 'Max', header.nslices-1);
set(handles.slice_selected, 'String', '1');  % show nº of slice

% Select b-value
theRangeB = header.diffvalues-1;
stepsB = [1/theRangeB, 1/theRangeB];
set(handles.b_selection, 'UserData', {PV, header});
set(handles.b_selection, 'Value', 0);
set(handles.b_selection, 'SliderStep', stepsB);
set(handles.b_selection, 'Max', header.diffvalues-1);
set(handles.b_selected, 'String', '0');

% SNR calculation
set(handles.snr_roi, 'UserData', PV(:,:,1,1));

% set data to handles
set(handles.b_values, 'UserData', bvalues);
set(handles.histogram, 'UserData', MR_gray);
set(handles.roi, 'UserData', {PV, header});
set(handles.import_roi, 'UserData', {header, PV});
set(handles.snr_roi, 'UserData', PV(:,:,1,1));
set(handles.all_header, 'UserData',header); % shows all of the informations in header
set(handles.param_map, 'UserData', {PV(:,:,1,:), bvalues, 1});
set(handles.complement, 'UserData', PV(:,:,1,1));
set(handles.adjust_histogram, 'UserData', handFig);

ROI = {'ROIs'; 'Labels'; 'Shape'};


% --------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ROI fitting.
function estimate_Callback(hObject, eventdata, handles)
% hObject    handle to estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No parameter can be estimated'}, 'Error', 'error');
elseif size(ROI,2) == 1
    msgbox({'There is no ROI on the image', 'Please, draw a ROI'}, 'Error', 'error');
else
    
    % Catching data from UserData
    data = get(hObject, 'UserData');
    Signal = data{1}; % IVIM signal
    b = data{2};
    header = data{3};

    nROI = size(Signal,1);
    nbvalues = length(b);

    Signal_fit = zeros(nROI,nbvalues);
    D_diff_fit = zeros(nROI,1);
    D_perf_fit = zeros(nROI,1);
    f_fit = zeros(nROI,1);

    fit = get(handles.Method, 'UserData');

    % Fitting per ROI
    for iROI = 1:nROI
        [D_diff_fit(iROI), D_perf_fit(iROI), f_fit(iROI)] = fitting(fit, 'Tikh', Signal(iROI,:), b, 'talwar', [0.15 1e-2 1e-3]);
        Signal_fit(iROI,:) = f_fit(iROI)*exp(-b'*D_perf_fit(iROI)) + (1-f_fit(iROI))*exp(-b'*D_diff_fit(iROI));       
        output(:,iROI) = {f_fit(iROI); D_diff_fit(iROI); D_perf_fit(iROI); f_fit(iROI)*D_perf_fit(iROI)};
        set(handles.parameters, 'data', output);        
        name_column{iROI} = ROI{2,iROI+1}.String; % set the name of the column: 1, 2, 3, ...
        r_2 = max(0, 1-sum((Signal(iROI,:) - Signal_fit(iROI,:)).^2)/sum((Signal(iROI,:) - mean(Signal(iROI,:))).^2)); digits(2); r_2 = vpa(r_2);       
        R2{iROI} = ['Fit ', num2str(iROI),' R²: ', char(r_2)];
    end
    
    set(handles.parameters, 'ColumnName', name_column'); % set column names 
    textLegend = cell(1,2*nROI); 

    % legend text
    for iROI = 1:nROI
        textLegend{iROI} = ROI{2,iROI+1}.String;
        textLegend{iROI+nROI} = [ROI{2,iROI+1}.String ' fit'];
    end

    % plot
    F = figure('Name', 'Fitting curves'); plot(b, Signal, 'o', b, Signal_fit, 'Linewidth', 2); grid;
    annotation('textbox', [0.15,0.3,0.1,0.1],'String', R2);
    ylabel('S(b)/S_{0}');
    xlabel('b[s/mm²]');
    title('Fitting curves and IVIM mean signal per ROI');
    legend(textLegend, 'Location', 'best');
    set(handles.Data, 'UserData', {Signal, b, header, output'});
    filename = regexprep(get(handles.Open, 'UserData'), '.PAR', '_ESTIMATION.jpg');
    saveas(F, fullfile('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE', filename));
end

% --- Executes on selection change in Fitting method pop-up.
function Method_Callback(hObject, eventdata, handles)
% hObject    handle to Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Method
ent_metodo = get(hObject, 'String');
num_metodo = get(hObject, 'Value');
texto_metodo = cell2mat(ent_metodo(num_metodo));

switch texto_metodo
    
    case 'Levenberg-Marquadt' % Default method
        
        fit = 'nlls_lm';
        
    case 'Trust-Region-Reflective'
        
        fit = 'nlls_trr';
                
    case 'Segmented Nonlinear Least Square'
        
        fit = 'nlls2';
               
    case 'Nonnegative Least Square'
        
        fit = 'nnls';
                
    case 'Linear Least Square'
        
        fit = 'lls';

    case 'Robust Linear Least Square'
        
        fit = 'llsr';     
        
end
set(hObject, 'UserData', fit);


% --- Executes during object creation, after setting all properties.
function Method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in Export > Data.
function Data_Callback(hObject, eventdata, handles)
% hObject    handle to Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
data = get(hObject, 'UserData');
if size(data,2) == 4 
    
    Signal = data{1};
    b = data{2};
    header = data{3};    
    output = cell2mat(data{4});
    
    slice = get(handles.slice_selected, 'String');
    te = header.tbl(1, header.tblcols.echo_time); % Echo time
    nROI = size(Signal,1);
    SigLen = size(Signal,2);
    filename = regexprep(get(handles.Open, 'UserData'), '.PAR', '.txt');
    fileID = fopen(filename, 'wt+');

    switch get(handles.Method, 'UserData')

         case 'nlls_lm' 

            fit = 'Levenberg-Marquadt';

        case 'nlls_trr'

            fit = 'Trust-Region-Reflective';

        case 'nlls2'

            fit = 'Segmented Nonlinear Least Square';

        case 'nnls'

            fit = 'Nonnegative Least Square';

        case 'lls'

            fit = 'Linear Least Square';

        case 'llsr'

            fit = 'Robust Linear Least Square';   

    end

    pMode = get(handles.uibuttongroup1, 'SelectedObject');
    pMode = get(pMode, 'String');

    title = 'IVIM SIGNAL ACQUIRED FROM IMAGE, B-VALUES UTILIZED AND ESTIMATED PARAMETERS';
    top = '--------------GENERAL INFORMATION--------------';
    fprintf(fileID, '%s\n', title);
    fprintf(fileID, '\n%s\n\n', top);
    info1 = {'Patient Name: ', 'FOV(A/P) [mm]: ', 'FOV(F/H) [mm]: ', 'FOV(R/L) [mm]: ', 'Nº Slices: ', 'Slice: ', 'TR [ms]: ', 'TE [ms]: ',...
        'Examination date/time: ', 'Nº b-values: ', 'Fitting method: ', 'SNR: ', 'Noise scale: '};
    info2 = {header.patname; num2str(header.fov(1)); num2str(header.fov(2)); num2str(header.fov(3)); num2str(header.nslices); slice;...
        header.tr; te; header.datetime; num2str(header.diffvalues); fit; get(handles.snr_estimated, 'String'); pMode};


    for i = 1:numel(info2)
        fprintf(fileID, '%s', info1{i});
        fprintf(fileID, '%s \n', info2{i});
    end

    fprintf(fileID, '\n%s\n', '--------------SIGNALS & B-VALUES--------------');

    for iROI = 1:nROI
        fprintf(fileID, '\n%s ', [ROI{2,iROI+1}.String ': ']);
        for iSigLen = 1:SigLen
            fprintf(fileID, '%3.2f ', Signal(iROI,iSigLen));
        end    
    end

    fprintf(fileID, '\n%s ', 'b: '); 

    for k = 1:length(b)
        fprintf(fileID, '%d ', b(k));    
    end

    fprintf(fileID, '\n\n%s\n\n', '--------------ESTIMATED PARAMETERS--------------');

    param = 'Format -> Label: f D [mm²/s] D* [mm²/s]';

    fprintf(fileID, '%s\n', param);
    for iROI = 1:nROI
        fprintf(fileID, '\n%s ', [ROI{2,iROI+1}.String ': ']);
        for j = 1:3
            if (j == 1)
                format = '%5.3f  ';
            else
                format = '%7.6f  ';
            end
            fprintf(fileID, format, output(iROI,j));
        end    
    end
    fclose(fileID);
else
    msgbox({'There are still some data to be generated', 'No data can be exported'}, 'Error', 'error');
end


% --- Executes on button press in Export > Image.
function Image_Callback(hObject, eventdata, handles)
% hObject    handle to Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save image with ROIs
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No image can be exported'}, 'Error', 'error');
else
    old_ax = handles.axes1; 
    handFig = figure;
    new_ax = gca;
    copyobj(get(old_ax, 'children'), new_ax);
    axis ij image;
    colormap(gray);
    set(gca, 'XTick', [], 'YTick', []);
    filenameImage = regexprep(get(handles.Open, 'UserData'), '.PAR', '.jpg');
    saveas(new_ax, fullfile('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE', filenameImage));
    close(handFig);
end

% --- Executes on button press in Estimate.
function snr_estimate_Callback(hObject, eventdata, handles)
% hObject    handle to snr_estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No SNR can be estimated'}, 'Error', 'error');
else
    data = get(hObject, 'UserData');
    sig = data{1};
    noise = data{2};
    pMode = get(handles.uibuttongroup1, 'SelectedObject');
    pMode = get(pMode, 'String');

    % Signal power and noise power
    sigPower = sum(abs(sig(:)).^2)/length(sig(:));
    noisePower = sum(abs(noise(:)).^2)/length(noise(:));

    % Decibel
    if(strcmp(pMode,'dB'))
      sigPower = 10*log10(sigPower);
      noisePower = 10*log10(noisePower);
      reqSNR = sigPower-noisePower;
    else
      reqSNR = sigPower/noisePower;    
    end

    set(handles.snr_estimated, 'String', num2str(round(reqSNR)));
end

% --- Executes on button press in SNR ROI.
function snr_roi_Callback(hObject, eventdata, handles)
% hObject    handle to snr_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'on')
    uiwait(msgbox({'First, draw a ROI on the background.', 'Then, draw a ROI on the tissue'}, 'Warning!', 'help'));
    MRI = get(hObject, 'UserData');
    h1 = imfreehand; BW1 = createMask(h1); % First background ROI
    h2 = imfreehand; BW2 = createMask(h2); % Second tissue ROI
    [s1,q1] = find(BW1);
    [s2,q2] = find(BW2);
    noise = MRI(s1,q1); % Select Signal within the ROI
    sig = MRI(s2,q2);   % Select Signal within the ROI
    set(handles.snr_estimate, 'UserData', {sig, noise});
    h1.delete;
    h2.delete;
else
    msgbox({'No .PAR image selected', 'No SNR can be estimated'}, 'Error', 'error');
end



% --- Executes on button press in snr_roi.
function snr_estimated_Callback(hObject, eventdata, handles)
% hObject    handle to snr_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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


% --- Executes on button press in import_roi.
function import_roi_Callback(hObject, eventdata, handles)
% hObject    handle to import_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No ROI can be imported'}, 'Error', 'error');
else
    data = get(hObject, 'UserData');
    header = data{1};
    img = data{2};
    slice = str2double(get(handles.slice_selected, 'String'));
    [filename, ~, ~] = uigetfile('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE\ROIs\*.txt', 'Select ROI Protocol file');
    fileID = fopen(['C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE\ROIs\' filename], 'r');    
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
            prompt = 'Enter label color (red, blue, cyan, green, yellow)';
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
    set(handles.IVIM_signal, 'UserData', {header, img(:,:,slice,:)});
    ROI{1,1} = 'ROIs';            
    ROI{2,1} = 'Label';
    ROI{3,1} = 'Shape';  
end


% --- Executes on button press in delete_rois.
function delete_rois_Callback(hObject, eventdata, handles)
% hObject    handle to delete_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No ROI can be deleted'}, 'Error', 'error');
elseif size(ROI, 2) == 1
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


% --- Executes on button Change label press.
function change_label_Callback(hObject, eventdata, handles)
% hObject    handle to change_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ROI;
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No label can be changed'}, 'Error', 'error');
elseif size(ROI, 2) == 1
     msgbox({'There is no ROI on the image', 'Please, draw a ROI'}, 'Error', 'error');
else
    prompt = {'Enter label to be changed', 'Write a new label'};
    dlg_title = 'Change ROI label';
    answer = inputdlg(prompt,dlg_title,[1 40; 1 40]);
    nROI = size(ROI,2)-1;
    for iROI = 2:nROI+1
        if strcmp(answer{1}, ROI{2,iROI}.String)        
            roi = ROI{1,iROI}; 
            label = ROI{2,iROI};        
            delete(label);
            ROIlab = gtext(answer{2}, 'Color', roi.getColor);
            ROI{2, iROI} = ROIlab;
            break;
        end
    end
end


% --- Executes on button press in all_header.
function all_header_Callback(hObject, eventdata, handles)
% hObject    handle to all_header (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No information can be shown'}, 'Error', 'error');
else
    msgbox('See command window for more information', 'Warning', 'help');
    fprintf('ALL INFORMATIONS FROM HEADER\n');
    disp(get(hObject, 'UserData'));
end


% --- Executes on button press in param_map.
function param_map_Callback(hObject, eventdata, handles)
% hObject    handle to param_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No parameter can be estimated'}, 'Error', 'error');
else
    
    data = get(hObject, 'UserData');
    img = data{1};
    bvalues = data{2};
    slice = data{3};
    fit = get(handles.Method, 'UserData');
    
    f_min = str2double(get(handles.outlier_fmin, 'String'));f_max = str2double(get(handles.outlier_fmax, 'String'));
    D_min = str2double(get(handles.outlier_Dmin, 'String'));D_max = str2double(get(handles.outlier_Dmax, 'String'));
    Dperf_min = str2double(get(handles.outlier_Dperf_min, 'String'));Dperf_max = str2double(get(handles.outlier_Dperf_max, 'String'));
    outlierList = [f_min, f_max, D_min, D_max, Dperf_min, Dperf_max];
    
    if ismember(1, isnan(outlierList))
        choice = questdlg({'It seems you did not define some outliers'; 'Would you like to proceed?'},'Outliers Undifined','Yes','No','No');        
        if strcmp(choice,'Yes')
            if ~isempty(fit)                
                [parametric_map_f, parametric_map_Ddiff, parametric_map_Dperf] = calculateMAP(img, bvalues, slice, fit, outlierList); 
                
                % Send data to analyse
                set(handles.analyse_map, 'UserData', {parametric_map_Ddiff, parametric_map_f, parametric_map_Dperf});
            else
                msgbox({'No fitting method was chosen', 'Please, choose a fitting method'}, 'Error', 'error');
            end       
        end
    else
        [parametric_map_f, parametric_map_Ddiff, parametric_map_Dperf] = calculateMAP(img, bvalues, slice, fit, outlierList);
        
        % Send data to analyse
        set(handles.analyse_map, 'UserData', {parametric_map_Ddiff, parametric_map_f, parametric_map_Dperf});
    end    
end


% --- Executes on button press in complement.
function complement_Callback(hObject, eventdata, handles)
% hObject    handle to complement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No complement can be shown'}, 'Error', 'error');
else
    MR_gray_0 = mat2gray(get(hObject, 'UserData'));
    MR_gray = imcomplement(MR_gray_0);
    handFig = montage(MR_gray, 'size', [1 1], 'Indices', 1);
    set(handles.axes1, 'Visible', 'On');
    set(gca, 'XTick', [], 'YTick', []);
    set(handles.histogram, 'UserData', MR_gray);
    set(handles.adjust_histogram, 'UserData', handFig);
end

% --- Executes on button press in adjust_histogram.
function adjust_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to adjust_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(handles.axes1, 'Visible'), 'off')
    msgbox({'No .PAR image selected', 'No adjustment can be made'}, 'Error', 'error');
else
    handFig = get(hObject, 'UserData');
    h = imcontrast(handFig);
    uiwait(h);
    set(handles.Image, 'UserData', handFig);
    set(handles.histogram, 'UserData', handFig.CData);
    set(handles.complement, 'UserData', handFig.CData);
end



function outlier_fmin_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_fmin as text
%        str2double(get(hObject,'String')) returns contents of outlier_fmin as a double


% --- Executes during object creation, after setting all properties.
function outlier_fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlier_fmax_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_fmax as text
%        str2double(get(hObject,'String')) returns contents of outlier_fmax as a double


% --- Executes during object creation, after setting all properties.
function outlier_fmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_fmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlier_Dmin_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_Dmin as text
%        str2double(get(hObject,'String')) returns contents of outlier_Dmin as a double


% --- Executes during object creation, after setting all properties.
function outlier_Dmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlier_Dmax_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_Dmax as text
%        str2double(get(hObject,'String')) returns contents of outlier_Dmax as a double


% --- Executes during object creation, after setting all properties.
function outlier_Dmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlier_Dperf_min_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_Dperf_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_Dperf_min as text
%        str2double(get(hObject,'String')) returns contents of outlier_Dperf_min as a double


% --- Executes during object creation, after setting all properties.
function outlier_Dperf_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_Dperf_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function outlier_Dperf_max_Callback(hObject, eventdata, handles)
% hObject    handle to outlier_Dperf_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outlier_Dperf_max as text
%        str2double(get(hObject,'String')) returns contents of outlier_Dperf_max as a double


% --- Executes during object creation, after setting all properties.
function outlier_Dperf_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outlier_Dperf_max (see GCBO)
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


% --------------------------------------------------------------------
function analyse_Callback(hObject, eventdata, handles)
% hObject    handle to analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function analyse_map_Callback(hObject, eventdata, handles)
% hObject    handle to analyse_map (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(hObject, 'UserData');

if isempty(data) % If there is no parametric data from current .PAR archive, choose data
    filedir = uigetdir('C:\Users\breno\Desktop\Mestrado\MATLAB\Fernando\IVIM IDE', 'Select directory');
    filename = uigetfile([filedir '\*.mat'], 'Select data');
    load([filedir '\' filename], 'parametric_map_Dperf', 'parametric_map_f', 'parametric_map_Ddiff');
else
    parametric_map_Ddiff = data{1};
    parametric_map_f = data{2};
    parametric_map_Dperf = data{3};    
end

% Send data to 
setappdata(0, 'parametric_map_f', parametric_map_f);
setappdata(0, 'parametric_map_Ddiff', parametric_map_Ddiff);
setappdata(0, 'parametric_map_Dperf', parametric_map_Dperf);

IDE_v4_analyse;
