function varargout = spm_fnirs(varargin)
% SPM for fNIRS toolbox: (startup function)
%
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% 
% $Id$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spm_fnirs_OpeningFcn, ...
                   'gui_OutputFcn',  @spm_fnirs_OutputFcn, ...
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


% --- Executes just before spm_fnirs is made visible.
function spm_fnirs_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for spm_fnirs
handles.output = hObject;

mdir = fileparts(which('spm_fnirs'));
sdir = fullfile(mdir, 'nfri_functions');
addpath(genpath(sdir));
addpath(fullfile(mdir, 'canonical')); 
addpath(genpath(fullfile(mdir, 'external'))); 
addpath(fullfile(mdir, 'gen_conimg')); 

clc; 
fprintf('\n');
disp( 'Statistical Parametric Mapping for Functional Near-Infrared Spectroscopy.'); 
fprintf('\n');

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = spm_fnirs_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
% Convert fNIRS data to optical density and hemoglobin changes 
%--------------------------------------------------------------------------
function push_convert_Callback(hObject, eventdata, handles)
spm_fnirs_convert_ui;

%--------------------------------------------------------------------------
% Transform optode/channel positions from subject space to MNI space 
%--------------------------------------------------------------------------
function push_spatial_Callback(hObject, eventdata, handles)
 spm_fnirs_spatialpreproc_ui;

%--------------------------------------------------------------------------
% Apply temporal filters to fNIRS data 
%--------------------------------------------------------------------------
function push_temporal_Callback(hObject, eventdata, handles)
spm_fnirs_temporalpreproc_ui; 

%--------------------------------------------------------------------------
% Specify general linear model (GLM) of fNIRS for the 1st level analysis 
%--------------------------------------------------------------------------
function push_specify_Callback(hObject, eventdata, handles)
spm_fnirs_specify1st_ui; 

%--------------------------------------------------------------------------
% Estimate GLM parameters for each channel 
%--------------------------------------------------------------------------
function push_estimate_Callback(hObject, eventdata, handles)
 spm_fnirs_estimate_ui
 
%--------------------------------------------------------------------------
% Interpolate GLM parameters and compute thresholded SPM 
%--------------------------------------------------------------------------
function push_results_Callback(hObject, eventdata, handles)
spm_fnirs_results_ui; 

