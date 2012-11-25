function varargout = TrafficDensity(varargin)
% TRAFFICDENSITY M-file for TrafficDensity.fig
%      TRAFFICDENSITY, by itself, creates a new TRAFFICDENSITY or raises the existing
%      singleton*.
%
%      H = TRAFFICDENSITY returns the handle to a new TRAFFICDENSITY or the handle to
%      the existing singleton*.
%
%      TRAFFICDENSITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRAFFICDENSITY.M with the given input arguments.
%
%      TRAFFICDENSITY('Property','Value',...) creates a new TRAFFICDENSITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrafficDensity_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrafficDensity_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrafficDensity

% Last Modified by GUIDE v2.5 25-Nov-2012 20:53:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrafficDensity_OpeningFcn, ...
                   'gui_OutputFcn',  @TrafficDensity_OutputFcn, ...
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


% --- Executes just before TrafficDensity is made visible.
function TrafficDensity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrafficDensity (see VARARGIN)

% Choose default command line output for TrafficDensity
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TrafficDensity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrafficDensity_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnOpen.
function btnOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,FilePath] = uigetfile('*.avi','Select the video');
set(handles.btnVideoOpen, 'Enable', 'off');
set(handles.edtVideoGen, 'String', 'Není k dispozici');

if ~isequal(FileName,0)                                     % If correct select file
    set(handles.edtOpen, 'String', FileName);
    set(handles.btnParamGen, 'Enable', 'on');
    [non name] = fileparts(FileName);                       % separe filen name and extension
    handles.FilePath = FilePath;
    handles.FileName = name; 
    if exist(fullfile(FilePath, [name '.mat']),'file') ~= 0         % Exists parameters of the road?
        if exist(fullfile(FilePath, [name '.png']),'file') ~= 0 % Exists background ?
            set(handles.btnParamShow, 'Enable', 'on');
            set(handles.btnVideoGen, 'Enable', 'on');
            set(handles.edtParam, 'String', 'Úspìšnì naèteny');
            handles.ImageBcg = imread(fullfile(FilePath, [name '.png']));
            load(fullfile(FilePath, [name '.mat']));
            handles.ImageParam = roadLane;
        else
            OpenConfFail(handles);
        end
    else
        OpenConfFail(handles);
    end
else
    set(handles.edtOpen, 'String', 'Žádné');
    set(handles.btnParamGen, 'Enable', 'off');
    OpenConfFail(handles);
end

guidata(hObject, handles);

function OpenConfFail (handles)
    set(handles.btnParamShow, 'Enable', 'off');
    set(handles.btnVideoGen, 'Enable', 'off');
    set(handles.edtParam, 'String', 'Nejsou k dispozici');

function edtOpen_Callback(hObject, eventdata, handles)
% hObject    handle to edtOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtOpen as text
%        str2double(get(hObject,'String')) returns contents of edtOpen as a double


% --- Executes during object creation, after setting all properties.
function edtOpen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtParam_Callback(hObject, eventdata, handles)
% hObject    handle to edtParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtParam as text
%        str2double(get(hObject,'String')) returns contents of edtParam as a double


% --- Executes during object creation, after setting all properties.
function edtParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnParamGen.
function btnParamGen_Callback(hObject, eventdata, handles)
% hObject    handle to btnParamGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
trafficObj = mmreader(fullfile(handles.FilePath, [handles.FileName '.avi'])); %nactu video
bcg = uint8(get_background(trafficObj,30));
imwrite(bcg,fullfile(handles.FilePath, [handles.FileName '.png']),'png');
roadLane = GetTrafficLane(bcg,false);
save(fullfile(handles.FilePath, [handles.FileName '.mat']),'roadLane');
set(handles.btnParamShow, 'Enable', 'on');
handles.ImageBcg = bcg;
handles.ImageParam = roadLane;
guidata(hObject, handles);

% --- Executes on button press in btnParamShow.
function btnParamShow_Callback(hObject, eventdata, handles)
% hObject    handle to btnParamShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(1);
imshow(handles.ImageBcg);
hold on;
title(['Jízdní pruhy pro: ' handles.FileName '.avi']);
for i=1:3
   plot(handles.ImageParam.left(:,1,i),handles.ImageParam.left(:,2,i),'LineWidth',3); 
   plot(handles.ImageParam.right(:,1,i),handles.ImageParam.right(:,2,i),'LineWidth',3); 
end



function edtVideoGen_Callback(hObject, eventdata, handles)
% hObject    handle to edtVideoGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtVideoGen as text
%        str2double(get(hObject,'String')) returns contents of edtVideoGen as a double


% --- Executes during object creation, after setting all properties.
function edtVideoGen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtVideoGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnVideoGen.
function btnVideoGen_Callback(hObject, eventdata, handles)
% hObject    handle to btnVideoGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnVideoOpen.
function btnVideoOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btnVideoOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
