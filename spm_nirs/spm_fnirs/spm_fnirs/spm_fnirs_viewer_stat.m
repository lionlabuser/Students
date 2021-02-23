function varargout = spm_fnirs_viewer_stat(varargin)
% Specify contrast vector, compute inference SPM, and then visualize
% thresholded SPM on a surface of rendered brain template 
% FORMAT spm_fnirs_viewer_stat(SPM) 
%
% SPM       structure array of estimated GLM parameters 
%
%--------------------------------------------------------------------------
% note: 
% In this code, 
% (i) contrast vector is specified, using spm_conman.m 
% (ii) inference SPM is computed, using spm_fnirs_contrasts.m 
% (iii) height threshold is computed, using spm_uc.m 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @spm_fnirs_viewer_stat_OpeningFcn, ...
    'gui_OutputFcn',  @spm_fnirs_viewer_stat_OutputFcn, ...
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


function spm_fnirs_viewer_stat_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

SPM = varargin{1,1};

load(SPM.xY.VY, 'P');
load(P.fname.pos); clear P;

% set slider ui
nsteps = 6;
view = 2; % intial view of rendered brain (dorsal view) 
set(handles.slider_view, 'sliderstep', [1/(nsteps-1), 1/(nsteps-1)], 'max', nsteps, 'min', 1, 'value', view);

% display rendered brain
cmap = load('Split.mat');
brain = R.rend{view}.ren; 

axes(handles.axes_image);
imagesc(brain, [0 2]);
colormap(cmap.split);
axis image
axis off

handles.SPM = SPM;
handles.R = R;
handles.Z = cell(1, nsteps);
handles.u = []; % height threshold 

set(handles.listbox_pixel, 'string', {}); 

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = spm_fnirs_viewer_stat_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


function slider_view_Callback(hObject, eventdata, handles)
R = handles.R;
Z = handles.Z;
u = handles.u; 

view = round(get(handles.slider_view, 'value')); set(handles.slider_view, 'value', view);
brain = R.rend{view}.ren .* 64;

%----------------------------------------------------------------------
% overlay statisitc image with rendered brain
if ~isempty(Z{view})
    if ~isempty(u),
        indx_m = find(Z{view} > u);
    else
        indx_m = find(isnan(Z{view}) == 0);
    end
    min_Z = min(Z{view}(indx_m)); max_Z = max(Z{view}(indx_m));
    brain(indx_m) = spm_fnirs_adjust_colorscale(Z{view}(indx_m), min_Z, max_Z);
end

%----------------------------------------------------------------------
% display results on GUI
axes(handles.axes_image);
image(brain);
axis off
axis image

cmap = load('Split.mat');
colormap(cmap.split),

if ~isempty(Z{view}) && ~isempty(indx_m)
    hc = colorbar;
    set(hc, 'YLim', [65 128], 'YTick', linspace(65, 128, 5)');
    y_tick = linspace(min_Z, max_Z, 5)';
    y_tick = round(10.*y_tick)./10;
    
    set(hc, 'YTickLabel', num2str(y_tick));
    set(hc, 'FontSize', 12);
end

guidata(hObject, handles);


function slider_view_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function push_contrast_Callback(hObject, eventdata, handles)
SPM = handles.SPM;
R = handles.R;

%----------------------------------------------------------------------
% specify contrast vector, and compute inference SPM 
cd(SPM.swd); 

try, xCon = SPM.xCon; catch, xCon = {}; end
[ic,xCon] = spm_conman(SPM,'T&F',Inf,...
    '    Select contrasts...',' for conjunction',1);

SPM.xCon = xCon;
SPM = spm_fnirs_contrasts(SPM, ic);

load(SPM.xCon(ic).Vspm);

view = get(handles.slider_view, 'value');
brain = R.rend{view}.ren .* 64;

%----------------------------------------------------------------------
% represent inference SPM on the rendered brain surface 
if ~isempty(Z{view})
    indx_m = find(isnan(Z{view}) == 0);
    min_Z = min(Z{view}(indx_m)); max_Z = max(Z{view}(indx_m));
    brain(indx_m) = spm_fnirs_adjust_colorscale(Z{view}(indx_m), min_Z, max_Z);
end

%----------------------------------------------------------------------
% display inference SPM 
axes(handles.axes_image);
image(brain);
axis off
axis image

cmap = load('Split.mat');
colormap(cmap.split)

if ~isempty(Z{view})
    hc = colorbar;
    set(hc, 'YLim', [65 128], 'YTick', linspace(65, 128, 5)');
    y_tick = linspace(min_Z, max_Z, 5)';
    y_tick = round(10.*y_tick)./10;
    
    set(hc, 'YTickLabel', num2str(y_tick));
    set(hc, 'FontSize', 12);
end

% Update handles structure
handles.Z = Z;
handles.SPM = SPM;
handles.ic = ic;
handles.u = []; 

guidata(hObject, handles);

function push_activation_Callback(hObject, eventdata, handles)
SPM = handles.SPM;
R = handles.R;
Z = handles.Z;
ic = handles.ic;

view = get(handles.slider_view, 'value');
%----------------------------------------------------------------------
% compute thresholded SPM (given p-value), using random field theory 
brain = R.rend{view}.ren .* 64;

if ~isempty(Z{view})
    %----------------------------------------------------------------------
    % calculate height thresholds
    thresDesc = spm_input('p value adjustment to control', 1, 'b', 'FWE|none', [], 1);
    
    df = [SPM.xCon(ic).eidf, SPM.xX.erdf];
    
    switch thresDesc
        case 'FWE' % Family-wise false positive rate
            %----------------------------------------------------------------------
            u = spm_input('p value (FWE)','+0','r',0.05,1,[0,1]);
            u = spm_uc(u,df, SPM.xCon(ic).STAT , SPM.xVol.R{view},1,SPM.xVol.S{view});
            
        case 'none' % No adjustment
            %----------------------------------------------------------------------
            u = spm_input(['threshold {',SPM.xCon(ic).STAT,' or p value}'],'+0','r',0.001,1);
            
            if u <= 1
                thresDesc = ['p<' num2str(u) ' (unc.)'];
                u = spm_u(u,df, SPM.xCon(ic).STAT);
            end
    end
    
    indx_over = find(Z{view} > u);
    min_Z = min(Z{view}(indx_over)); max_Z = max(Z{view}(indx_over));
    brain(indx_over) = spm_fnirs_adjust_colorscale(Z{view}(indx_over), min_Z, max_Z);
end

%----------------------------------------------------------------------
% display thresholded SPM (brain activation) 
axes(handles.axes_image);
image(brain);
axis off
axis image

cmap = load('Split.mat'); 
colormap(cmap.split),

if ~isempty(Z{view}) && ~isempty(indx_over)
    hc = colorbar;
    set(hc, 'YLim', [65 128], 'YTick', linspace(65, 128, 5)');
    y_tick = linspace(min_Z, max_Z, 5)';
    y_tick = round(10.*y_tick)./10;
    
    set(hc, 'YTickLabel', num2str(y_tick));
    set(hc, 'FontSize', 12);
end

handles.u = u; 
guidata(hObject, handles);


function listbox_pixel_Callback(hObject, eventdata, handles)

function listbox_pixel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function push_pixel_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------------
% Identify MNI coordinates and statistic value of a specified pixel 
R = handles.R;
Z = handles.Z; 
ic = handles.ic; 
SPM = handles.SPM;

%----------------------------------------------------------------------
% identify (x,y) positions of point 
[x, y] = ginput(1); 
x = round(x); y = round(y); 

view = get(handles.slider_view, 'value');

%----------------------------------------------------------------------
% identify MNI poitions of point 
xyz = R.rend{view}.xyz; 
xyz = reshape(xyz, [3 size(R.rend{view}.ren)]);
xyz_p = xyz(:, y, x); 

str{1,1} = sprintf('xyz [mm]: %3.2f %3.2f %3.2f \n', xyz_p); 

%----------------------------------------------------------------------
% identify statistic value 
str{2,1} = sprintf('%c stat: %3.2f \n', handles.SPM.xCon(handles.ic).STAT, Z{view}(y,x)); 

%----------------------------------------------------------------------
[label acc] = nfri_anatomlabel_spm(xyz_p', 'test', 8, 6);
label = label{1,1}; acc = acc{1}; 
nlabel = size(label, 1);
for i = 1:nlabel 
    str{2+i, 1} = sprintf('%s : %3.2f \n', label{i,1}, acc(i,1)); 
end
set(handles.listbox_pixel, 'string', str); 


function push_save_Callback(hObject, eventdata, handles)
%----------------------------------------------------------------------
% save the current axes as png file format
[fname fdir] = uiputfile('*.png', 'Save the current image as');
if fname == 0, return; end 

sname = fullfile(fdir, fname); 

F = getframe(gca);
imwrite(F.cdata, sname);

sname = sprintf('%s_with_GUI.png', fname(1:end-4)); 
F = getframe(gcf);
imwrite(F.cdata, sname); 

