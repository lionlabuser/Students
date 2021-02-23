function varargout = spm_fnirs_viewer_timeseries(varargin)
% GUI for plotting (i) intensity, (ii) standard deviation, and (iii)
% frequency amplitude of optical density and hemoglobin changes 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_fnirs_viewer_timeseries_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_fnirs_viewer_timeseries_OutputFcn, ...
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


% --- Executes just before spm_fnirs_viewer_timeseries is made visible.
function spm_fnirs_viewer_timeseries_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for spm_fnirs_viewer_timeseries
handles.output = hObject;

%--------------------------------------------------------------------------
% read time series of raw and (filtered) data 
y = varargin{1,1}; 
P = varargin{1,2}; 
fs = P.fs; 

fy = []; 
try, fy = varargin{1,3}; end 

if size(varargin, 2) == 4, 
    chs = varargin{1,4}; 
else
    chs = find(P.mask ~= 0); 
end 

%--------------------------------------------------------------------------
lcolor{1} = [0 0 1]; 
if isstruct(y),
    Y{1} = reshape(spm_vec(y), P.ns, P.nch, []);
else
    Y{1} = y;
end

if isempty(fy) 
    str_list{1,1} = 'OD (wave1)';
    str_list{2,1} = 'OD (wave2)';
    str_list{3,1} = 'HbO';
    str_list{4,1} = 'HbR';
    str_list{5,1} = 'HbT';
    set(handles.radio_filtered, 'enable', 'off', 'value', 0);

else 
    str_list{1,1} = 'HbO';
    str_list{2,1} = 'HbR';
    str_list{3,1} = 'HbT';
    set(handles.radio_filtered, 'enable', 'on', 'value', 1);
    
    try % down-sampled data
        fs(2) = K.D.nfs;
        ns = K.D.ns;
    catch % full data
        fs(2) = P.fs;
        ns = P.ns;
    end
    
    if isstruct(fy);
        Y{2} = reshape(spm_vec(fy), [ns P.nch n]);
    else
        Y{2} = fy;
    end
    lcolor{2} = [1 0 0];
end

%--------------------------------------------------------------------------
% plot time series 
axes(handles.axes_signal); 
etime = P.ns/P.fs; 
for i = 1:size(Y, 2)
    time{i} = linspace(0, etime, size(Y{i},1)); 
    plot(time{i}, Y{i}(:, chs(1), 1), 'color', lcolor{i}); 
    xlabel('Time [s]'); 
    hold on; 
end 
axis tight;
hold off; 

handles.Y = Y; 
handles.P = P; 
handles.fs = fs; 
handles.time = time; 
handles.chs = chs; 
handles.lcolor = lcolor; 

%--------------------------------------------------------------------------
% set-up GUI
nsteps = size(chs, 2); 
set(handles.slider_ch, 'sliderstep', [1/(nsteps-1), 1/(nsteps-1)], 'max', nsteps, 'min', 1, 'value', 1);
set(handles.edit_info, 'string', ['Ch ' num2str(chs(1))]);
set(handles.popup_wave, 'value', 1, 'string', str_list); 
set(handles.popup_signal, 'value', 1);
set(handles.radio_raw, 'value', 1); 
set(handles.edit_fname, 'string', P.fname.nirs); 
if P.mask(chs(1)) == 1 
    set(handles.radio_ch, 'value', 1); 
elseif P.mask(chs(1)) == -1 % bad channel
    set(handles.radio_ch, 'value', 0); 
end

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = spm_fnirs_viewer_timeseries_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_ch_Callback(hObject, eventdata, handles)
num_s = round(get(handles.slider_ch, 'value')); set(handles.slider_ch, 'value', num_s); 

ch = handles.chs(num_s); 
wav = get(handles.popup_wave, 'value');
stype = get(handles.popup_signal, 'value');
h = [get(handles.radio_raw, 'value') get(handles.radio_filtered, 'value')]; 

for i = 1:2
    if h(i) == 1
        y = handles.Y{i}(:, ch, wav);
        time = handles.time{i}; 
        plot_signal(time, y, stype, handles, i); 
        hold on; 
    end
end
hold off;

set(handles.edit_info, 'string', ['Ch ' num2str(ch)]);
if handles.P.mask(ch) == 1 
    set(handles.radio_ch, 'value', 1); 
elseif handles.P.mask(ch) == -1 % bad channel
    set(handles.radio_ch, 'value', 0); 
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_ch_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_info_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_info_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_wave.
function popup_wave_Callback(hObject, eventdata, handles)
num_s = get(handles.slider_ch, 'value');
ch = handles.chs(num_s); 
wav = get(handles.popup_wave, 'value'); 

h = [get(handles.radio_raw, 'value') get(handles.radio_filtered, 'value')]; 

axes(handles.axes_signal); 
for i = 1:2 
    if h(i) == 1
        plot(handles.time{i}, handles.Y{i}(:,ch,wav), 'color', handles.lcolor{i}); 
        xlabel('Time [s]');
        hold on; 
    end 
end
hold off;
axis tight; 

% reset ui
set(handles.popup_signal, 'value', 1); % light intensity 

% --- Executes during object creation, after setting all properties.
function popup_wave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_signal.
function popup_signal_Callback(hObject, eventdata, handles)
num_s = get(handles.slider_ch, 'value'); 
ch = handles.chs(num_s); 
wav = get(handles.popup_wave, 'value');
stype = get(handles.popup_signal, 'value');

if stype == 2
    prompt = {'moving window length [sec]', 'threshold factor-motion detection'};
    def = {'1', '3'};
    answer = inputdlg(prompt, '', 1, def);
    if isempty(answer), return; end 
    L = str2num(answer{1,1}); handles.L = L;
    th = str2num(answer{2,1}); handles.th = th;
end 

h = [get(handles.radio_raw, 'value') get(handles.radio_filtered, 'value')]; 

for i = 1:2
    if h(i) == 1
        y = handles.Y{i}(:, ch, wav);
        time = handles.time{i}; 
        plot_signal(time, y, stype, handles, i); 
        hold on; 
    end
end
hold off;

guidata(hObject, handles);
        

% --- Executes during object creation, after setting all properties.
function popup_signal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in radio_raw.
function radio_raw_Callback(hObject, eventdata, handles)
num_s = get(handles.slider_ch, 'value'); 
ch = handles.chs(num_s); 
wav = get(handles.popup_wave, 'value');
stype = get(handles.popup_signal, 'value');

h = [get(handles.radio_raw, 'value') get(handles.radio_filtered, 'value')]; 
axes(handles.axes_signal); 
cla; 

for i = 1:2
    if h(i) == 1
        y = handles.Y{i}(:, ch, wav);
        time = handles.time{i}; 
        plot_signal(time, y, stype, handles, i); 
        hold on; 
    end
end
hold off;

% --- Executes on button press in radio_filtered.
function radio_filtered_Callback(hObject, eventdata, handles)
num_s = get(handles.slider_ch, 'value'); 
ch = handles.chs(num_s); 
wav = get(handles.popup_wave, 'value');
stype = get(handles.popup_signal, 'value');

h = [get(handles.radio_raw, 'value') get(handles.radio_filtered, 'value')]; 
axes(handles.axes_signal); 
cla; 

for i = 1:2
    if h(i) == 1
        y = handles.Y{i}(:, ch, wav);
        time = handles.time{i}; 
        plot_signal(time, y, stype, handles, i); 
        hold on; 
    end
end
hold off;

% --- Executes on button press in radio_ch.
function radio_ch_Callback(hObject, eventdata, handles)
num_s = get(handles.slider_ch, 'value'); 
ch = handles.chs(num_s); 
P = handles.P; 

if get(handles.radio_ch, 'value') 
    P.mask(ch) = 1; 
else
    P.mask(ch) = -1;
end
save(P.fname.nirs, 'P', '-append', spm_get_defaults('mat.format'));

handles.P = P; 
guidata(hObject, handles);


function edit_fname_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_fname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_signal(x, y, stype, handles, i)

fs = handles.fs(i);

axes(handles.axes_signal); 
switch stype
    case 1 % intensity 
        plot(x, y, 'color', handles.lcolor{i});
        xlabel('Time [s]'); 
        axis tight;
        
    case 2 % standard deviation 
        y = spm_fnirs_MovStd(y, round(handles.L*fs/2));
        indx_n = find(isnan(y) == 0);
        th = ones(size(y, 1), 1) .* (handles.th*mean(y(indx_n)));
        plot(x, y, 'color', handles.lcolor{i}); hold on;
        plot(x, th, 'color', [0 0 0]);
        xlabel('Time [s]'); 
        axis tight;
        
    case 3 % frequency amplitude 
        ns = size(y, 1);
        nfft = 2^nextpow2(ns);
        y = fft(y, nfft)./ns;
        
        nf = nfft/2+1;
        f = fs ./ 2 * linspace(0,1,nf)';
        fy = 2*abs(y(1:nf));
        
        plot(f, fy, 'color', handles.lcolor{i});
        xlabel('Frequency [Hz]'); 
        axis([min(f) max(f) min(fy) mean(fy) *10]);
end % switch

