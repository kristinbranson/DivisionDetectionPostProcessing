function varargout = PlayDivision(varargin)
% PLAYDIVISION MATLAB code for PlayDivision.fig
%      PLAYDIVISION, by itself, creates a new PLAYDIVISION or raises the existing
%      singleton*.
%
%      H = PLAYDIVISION returns the handle to a new PLAYDIVISION or the handle to
%      the existing singleton*.
%
%      PLAYDIVISION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYDIVISION.M with the given input arguments.
%
%      PLAYDIVISION('Property','Value',...) creates a new PLAYDIVISION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlayDivision_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlayDivision_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlayDivision

% Last Modified by GUIDE v2.5 23-May-2017 10:55:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlayDivision_OpeningFcn, ...
                   'gui_OutputFcn',  @PlayDivision_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}),
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PlayDivision is made visible.
function PlayDivision_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlayDivision (see VARARGIN)

persistent detectionmatfile;
if isempty(detectmatfile) || ~exist(detectionmatfile,'file'),
  detectionmatfile = '';
end

[handles.detections,handles.detectionmatfile,handles.boxrad] = ...
  myparse(varargin,'detections',[],...
  'detectionmatfile','','boxrad',[100,100,100,20]);
if isempty(handles.detections),
  if isempty(handles.detectionmatfile),
    [fn,pn] = uigetfile('*.mat','Select detection mat file',detectionmatfile);
    if ~ischar(fn),
      delete(handles.figure1);
      return;
    end
    handles.detectionmatfile = fullfile(pn,fn);
  end
  handles.detections = load(detectionmatfile);
end

if isempty(detectionmatfile),
  detectionmatfile = handles.detectionmatfile;
end

if ~isfield(handles.detections,'t0'),
  handles.detections.t0 = 1;
end

[~,detectioni] = max(handles.detections.scores);
handles = SetDetectionIndex(handles,detectioni);

% Choose default command line output for PlayDivision
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlayDivision wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles = SetDetectionIndex(handles,detectioni)

handles.detectioni = detectioni;
handles.loc = handles.detections.locs(handles.detectioni,:);

% read in video info
v0 = max(1,floor(handles.loc)-handles.boxrad);
v1 = min([rawsz,predsz(3),ntimepoints],ceil(handles.loc)+handles.boxrad);
ncurr = v1 - v0 + 1;
switch handles.detections.rawfiletype,
  case 'klb',
    handles.rawvid = zeros(ncurr,'uint16');
    for t = v0(end):v1(end),
      handles.rawvid(:,:,:,t) = readKLBroi(handles.detections.inputdatafiles{t+handles.detections.t0-1},[v0(1:3);v1(1:3)]);
    end
  otherwise
    error('Not implemented');
end

switch handles.detections.predfiletype,
  case 'klb'
    handles.predvid = zeros(ncurr,'uint8');
    for t = v0(end):v1(end),
      handles.predvid(:,:,:,t) = readKLBroi(handles.detections.preddatafiles{t+t0-1},[v0(1:3);v1(1:3)]);
    end
  case 'h5'
    handles.predvid = h5read(handles.detections.preddatafile,handles.detections.preddatasetname,v0+[0,0,0,handles.detections.t0-1],v1-v0+1);
end

% --- Outputs from this function are returned to the command line.
function varargout = PlayDivision_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_z_Callback(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_time_Callback(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_time as text
%        str2double(get(hObject,'String')) returns contents of edit_time as a double


% --- Executes during object creation, after setting all properties.
function edit_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_z as text
%        str2double(get(hObject,'String')) returns contents of edit_z as a double


% --- Executes during object creation, after setting all properties.
function edit_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
