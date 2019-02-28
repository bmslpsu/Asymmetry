function [] = Experiment_Asymmetry_Control(Fn,spatFreq)
% Experiment_Asymmetry_Control_Contrast: runs an experiment using LED arena and fly panel controller
% This code is written for Panel Controller v3 and NiDAQ seesion mode
% NOTES:    1. The pattern can be of HIGH or LOW contrast or INTERPOLATED MOTION (requires three different sets of patterns)
%           2. User must maunally set the root folder in the fucntion
%   INPUTS:
%       Fn          :  	fly number
%       spatFreq    :   spatial frequency (0=Random)
%       steps       :   gain factor for interpolated motion (steps=1 for no interpolation)
%---------------------------------------------------------------------------------------------------------------------------------
daqreset
imaqreset
%% Set Directories & Controller Parameters %%
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = ['D:\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(spatFreq) '\'];
viddir = [rootdir 'Vid\'];
% validSpatFreq = 7.5*[3,4,8];
switch spatFreq
    case 0
        patID = 6;
        posY  = 1;
    case 22.5
        patID = 8;
        posY  = 3;
    case 30
        patID = 8;
        posY  = 4;
    case 60
        patID = 8;
        posY  = 6;
    otherwise
        error('Invalid Spatial Frequency!!')
end
%% EXPERIMENTAL PARAMETERS %%
%---------------------------------------------------------------------------------------------------------------------------------
n_tracktime = 10;  	% seconds for each EXPERIMENT
n_resttime = 2;    	% seconds for each REST
n_pause = 0.2;      % pause between commands
n_AI = 6;        	% # of analog input channels
n_rep = 5;          % number of cycles through velocities for each fly

%% SETUP DATA AQUISITION: NiDAQ (session mode)
%---------------------------------------------------------------------------------------------------------------------------------
% Find DAQ
devices = daq.getDevices;
s = daq.createSession('ni');
addAnalogInputChannel(s,devices.ID, 1:n_AI, 'Voltage'); % Add Analog Input Channels
    % AI1 = AO_0: trigger
    % AI2 = DAC1: x-position of stimulus
    % AI3 = DAC2: y-position of stimulus
    % AI4 = L: left wing beat amplitude
    % AI5 = R: right wing beat amplitude
    % AI6 = FREQ: frequency
    
% Add analog output channels
ch.AO = addAnalogOutputChannel(s,devices.ID,'ao0','Voltage');
ch.AO.Name = 'Trigger';

% Set Sampling Rate
s.Rate = 5000;                   % samples per second
s.IsContinuous = false;          % continuous data collection until stopped
% s.DurationInSeconds = 10;

% Setup Sampling
s.Rate = 5000; % samples per second
s.IsContinuous = false;	% continuous data collection until stopped


FrameRate = 200; % camera frame rate
nFrame = FrameRate * n_tracktime; % # of frames to log

t = 0:1/s.Rate:n_tracktime;
TriggerSignal = (square(2*pi*FrameRate*t,90) + 1)';

disp('DAQ Setup Done...')
%% SETUP CAMERA INPUT %%
%---------------------------------------------------------------------------------------------------------------------------------
adaptorName = 'gige';
deviceID = 1;
vidFormat = 'Mono8';
vid = videoinput(adaptorName, deviceID, vidFormat);

% Configure vidobj properties.
set(vid, 'ErrorFcn', @imaqcallback);
set(vid, 'LoggingMode', 'memory');
set(vid,'FrameGrabInterval',1);
ROI.x = 500;
ROI.y = 250;
ROI.xoff = (round(659 - ROI.x)/2);
ROI.yoff = (round(494 - ROI.y)/2);
vid.ROIPosition = [ROI.xoff ROI.yoff ROI.x ROI.y];
set(vid, 'FramesPerTrigger', 1); % frames to log for each trigger
vid.TriggerRepeat = nFrame-1; % # triggers

% Configure vidobj source properties.
srcObj1 = get(vid, 'Source');
srcObj1.Gamma = 0.386367797851563;
srcObj1.GainRaw = 964;
srcObj1.ExposureTimeAbs = 4000; % 200 Hz frame rate
% srcObj1(1).ExposureMode = 'Timed'; % exposure time controlled by pulse width

% Trigger config
triggerconfig(vid, 'hardware','DeviceSpecific','DeviceSpecific')
set(srcObj1(1),'LineSelector','Line1');
set(srcObj1(1),'TriggerActivation','RisingEdge');
set(srcObj1(1),'TriggerMode','on');
set(srcObj1(1),'TriggerSelector','FrameStart');

disp('VID Setup Done...')
%% Set Variable to Control Pattern Velocity %%
%---------------------------------------------------------------------------------------------------------------------------------
% Cycle through [pattern_ID , xpos , ypos , velocity mode nexptime marker] 
Speed = (30:30:150)'; % [deg/s] speed
Vel = [Speed;-Speed]; % [deg/s] CW & CCW speeds
nVel = length(Vel);   % # of velocities
Gain = Vel/3.75;      % Gains corresponding to each velocity
Gain_all = repmat(Gain,n_rep,1); % repeat gains for n_rep
Gain_rand = Gain_all(randperm(size(Gain_all,1)),:); % reshuffle randomly

%% START EXPERIMENT AND DATA COLLECTION %%
%---------------------------------------------------------------------------------------------------------------------------------
tic
for kk = 37:n_rep*nVel
    disp('-------------------------------------------------------')
    gain = Gain_rand(kk);               % gain corresponding to velocity for trial
    disp(['Trial:  ' num2str(kk)])      % prints counter to command line
   	preview(vid);       % open video preview window
	start(vid)          % start video buffer
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)               % set pattern to "Pattern_Fourier_bar_barwidth=8"
    Panel_com('set_position',[1, 1]); pause(n_pause)            % set start position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)                 % 0=open,1=closed,2=fgen,3=vmode,4=pmode)
    Panel_com('set_funcX_freq', 50); pause(n_pause)             % default X update rate
	Panel_com('set_funcY_freq', 50); pause(n_pause)           	% default Y update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)     % [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                         % start closed-loop tracking
    pause(n_resttime)
    Panel_com('stop');
    %-----------------------------------------------------------------------------------------------------------------------------
    % SETUP EXPERIMENT
    pause(1)
    disp(['move ground ' num2str(gain)])
    Panel_com('stop');pause(n_pause);
    Panel_com('set_pattern_id', patID);pause(n_pause)               % pattern = 
    Panel_com('set_position',[1, posY]);pause(n_pause)              % set starting position (xpos,ypos) [ypos = spatFreq]
    Panel_com('set_funcX_freq', 50);pause(n_pause);                 % default X update rate
    Panel_com('set_funcY_freq', 50);pause(n_pause);                 % default Y update rate
    Panel_com('set_mode', [0,0]);pause(n_pause)                     % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    Panel_com('send_gain_bias',[gain,0,0,0]);pause(n_pause)         % open-loop
    %-----------------------------------------------------------------------------------------------------------------------------
    % RUN EXPERIMENT AND COLLECT DATA
    queueOutputData(s,TriggerSignal) % set trigger AO signal
    Panel_com('start')  % run trial
    [data,t_p] = s.startForeground;
	stop(vid) % stop video buffer
    Panel_com('stop')
    %-----------------------------------------------------------------------------------------------------------------------------
    % GET DATA AND SAVE TO .mat FILE
    % Create named structure array to store data
	[vidData,t_v] = getdata(vid, vid.FramesAcquired); % get video data
      
    % CLOSED LOOP BAR TRACKING
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)               % set pattern to "Pattern_Fourier_bar_barwidth=8"
    Panel_com('set_position',[1, 1]); pause(n_pause)            % set start position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)                 % 0=open,1=closed,2=fgen,3=vmode,4=pmode)
    Panel_com('set_funcX_freq', 50); pause(n_pause)             % default X update rate
	Panel_com('set_funcY_freq', 50); pause(n_pause)           	% default Y update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)     % [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                         % start closed-loop tracking
    
    % Save data
    disp('Saving...')
    save([rootdir 'Fly_' num2str(Fn) '_Trial_' num2str(kk) '_Vel_' num2str(3.75*Gain_rand(kk)) '_SpatFreq_' ...
        num2str(spatFreq) '.mat'],'-v7.3','data','t_p');
    save([viddir 'Fly_' num2str(Fn) '_Trial_' num2str(kk) '_Vel_' num2str(3.75*Gain_rand(kk)) '_SpatFreq_' ...
        num2str(spatFreq) '.mat'],'-v7.3','vidData','t_v');
    %-----------------------------------------------------------------------------------------------------------------------------
end
toc
%---------------------------------------------------------------------------------------------------------------------------------
disp('Done');
daqreset
imaqreset
clear
end
