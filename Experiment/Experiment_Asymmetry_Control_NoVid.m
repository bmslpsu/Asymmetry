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
spatFreqFolder = ['D:\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(spatFreq) '\'];
% spatFreqFolder = ['E:\Experiment_Asymmetry_Control_V2\LowContrast' num2str(spatFreq)];
% spatFreqFolder = ['E:\Experiment_Asymmetry_Control_V2\InterpolatedMotion' num2str(spatFreq)];
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
n_resttime = 5;    	% seconds for each REST
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
    
% Set Sampling Rate
s.Rate = 1000;                  % samples per second
s.IsContinuous = false;          % continuous data collection until stopped
s.DurationInSeconds = 10;
disp('DAQ Setup Done...')
%% SETUP CAMERA INPUT %%
%---------------------------------------------------------------------------------------------------------------------------------
[vid] = RigidFlyCamSettingsBC_acA640_120gm_V1();
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
for kk = 1:n_rep*nVel
    disp('-------------------------------------------------------')
    gain = Gain_rand(kk);               % gain corresponding to velocity for trial
    disp(['Trial:  ' num2str(kk)])      % prints counter to command line
    preview(vid);
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
    Panel_com('set_pattern_id', patID);pause(n_pause)               % pattern = (Pattern_spatFreq_22_30_60 or Pattern_Random_Ground_48)
    Panel_com('set_position',[1, posY]);pause(n_pause)              % set starting position (xpos,ypos) [ypos = spatFreq]
    Panel_com('set_funcX_freq', 50);pause(n_pause);                 % default X update rate
    Panel_com('set_funcY_freq', 50);pause(n_pause);                 % default Y update rate
    Panel_com('set_mode', [0,0]);pause(n_pause)                     % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    Panel_com('send_gain_bias',[gain,0,0,0]);pause(n_pause)  	% open-loop
    %-----------------------------------------------------------------------------------------------------------------------------
    % RUN EXPERIMENT AND COLLECT DATA
    % Open log file and start data aquisition
%     fid1 = fopen('log.bin','w');
%     lh = addlistener(s,'DataAvailable',@(src, event)logData(src, event, fid1));
%     startBackground(s); % start data collection for WBA 
    Panel_com('start')  % run trial
    [raw_data,time] = s.startForeground;
%     pause(n_tracktime)
    Panel_com('stop')
    s.stop; % stop data colection
%     delete(lh); fclose(fid1); % reset data aquisition
    %-----------------------------------------------------------------------------------------------------------------------------
    % GET DATA AND SAVE TO .mat FILE
%     [raw_data,~]    = fread(fopen('log.bin','r'),[n_AI+1,inf],'double');   % get DAQ data
    % Create named structure array to store data
    dataDAQ.Time    = time;
    dataDAQ.Trig    = raw_data(1,:)';
    dataDAQ.Xpos    = raw_data(2,:)';
    dataDAQ.Ypos    = raw_data(3,:)';
    dataDAQ.LWing   = raw_data(4,:)';
    dataDAQ.RWing   = raw_data(5,:)';
    dataDAQ.WBF     = raw_data(6,:)';
    % Save data
    disp('Saving...')
    
    save([spatFreqFolder 'fly_' num2str(Fn) '_trial_' num2str(kk) '_Vel_' num2str(3.75*Gain_rand(kk)) '_SpatFreq_' ...
        num2str(spatFreq) '.mat'],'-v7.3','dataDAQ');
    %-----------------------------------------------------------------------------------------------------------------------------
end
toc
%---------------------------------------------------------------------------------------------------------------------------------
disp('Done');
daqreset
imaqreset
clear
end
