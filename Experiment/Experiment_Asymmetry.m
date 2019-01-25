function [] = Experiment_Asymmetry(Fn)
%% STAF_EXPERIMENT runs a experiment using the LED arena and fly panel
% Controller to generate a figure/ground STAF.
% Fn is the fly number
% STAF Experiment Template written by Ben Cellini
% This code is written for Panel Controller v3 and NiDAQ seesion mode
%---------------------------------------------------------------------------------------------------------------------------------
%rootdir = uigetdir({}, 'Select folder to save data'); % define directory to save file
rootdir = 'C:\BC\Rigid_data\Experiment_Asymmetry\';
viddir = 'C:\BC\Rigid_data\Experiment_Asymmetry\vid\';
file_name = ['fly_' num2str(Fn)]; % create file name

% EXPERIMENTAL PARAMETERS
n_tracktime = 10;           % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 5;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_trial = 10;               % # trials per fly
n_AI = 6;                   % # of analog input channels

tic;
%---------------------------------------------------------------------------------------------------------------------------------
% Set up data acquisition on NiDAQ (session mode)
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
s.IsContinuous = true;          % continuous data collection until stopped

% Set file name
fname = [rootdir '\' file_name '.daq'];

%---------------------------------------------------------------------------------------------------------------------------------
% Setup camera input
[vid] = RigidFlyCamSettingsBC_acA640_120gm_V1();
%---------------------------------------------------------------------------------------------------------------------------------
% Random sequence of -1 and 1  
Rdir = [1 1 -1 1 -1 1 -1 -1 1 -1 ];
%---------------------------------------------------------------------------------------------------------------------------------
% EXPERIMENT LOOP
for ii = 1:10      
   
    disp('trial')
    display(num2str(ii));              %prints counter to command line
    preview(vid);
    %-----------------------------------------------------------------------------------------------------------------------------
    % CLOSED LOOP BAR TRACKING for n_rest_pat;
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)       % set output to p_rest_pat (Pattern_Fourier_bar_barwidth=8)
    Panel_com('set_position',[1, 1]); pause(n_pause)    % set starting position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)      	% closed loop tracking [xpos,ypos] (NOTE: 0=open, 1=closed)
    Panel_com('set_posfunc_id',[1,0]); pause(n_pause)
    Panel_com('set_funcX_freq', 50); pause(n_pause)
    Panel_com('set_posfunc_id',[2,0]); pause(n_pause)
    Panel_com('set_funcY_freq', 50); pause(n_pause)
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)      % [xgain,xoffset,ygain,yoffset]
    Panel_com('start');
    pause(n_resttime)
    Panel_com('stop');
    
    %-----------------------------------------------------------------------------------------------------------------------------
    pause(1)    
    
    disp('move ground in random direction')
    Panel_com('set_pattern_id', 2); pause(n_pause)                	% set output to p_rest_pat (Pattern_Random_Ground_48)
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_position',[randi([1, 96]), 1]); pause(n_pause)  	% set starting position (xpos,ypos)
    Panel_com('set_posFunc_id',[1 0]); pause(n_pause)               % arg1 = channel (x=1); arg2 = funcxid = 1
    Panel_com('set_funcX_freq', 50); pause(n_pause)
    Panel_com('set_posFunc_id',[2 0]); pause(n_pause)               % arg1 = channel (y=2); arg2 = funcyid = 2
    Panel_com('set_funcY_freq', 50); pause(n_pause)
    Panel_com('set_mode', [0,0]); pause(n_pause)                    % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    Panel_com('send_gain_bias',[Rdir(ii)*12,0,0,0]); pause(n_pause)	% open loop
    
    % Open log file and start data collection
    fid1 = fopen('log.bin','w');
    lh = addlistener(s,'DataAvailable',@(src, event)logData(src, event, fid1));
    startBackground(s);
    start(vid)
    Panel_com('start')
    pause(n_tracktime);
    Panel_com('stop')
    s.stop;
    stoppreview(vid);
    stop(vid)
    delete(lh);
    fclose(fid1);
    
    [raw_data,count] = fread(fopen('log.bin','r'),[n_AI+1,inf],'double');
    
    t_p = raw_data(1,:);
    data = raw_data(2:n_AI+1,:);
    
    [vidData,t_v] = getdata(vid);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Save data according to wavelength and direction
    disp('Saving...')
    
    if Rdir(ii) == 1
        save([rootdir 'fly_' num2str(Fn) '_trial_' num2str(ii) '_CW'],'-v7.3','data','t_p');
        save([viddir  'fly_' num2str(Fn) '_trial_' num2str(ii) '_CW'],'-v7.3','vidData','t_v');
    elseif Rdir(ii) == -1
        save([rootdir 'fly_' num2str(Fn) '_trial_' num2str(ii) '_CCW'],'-v7.3','data','t_p');
        save([viddir  'fly_' num2str(Fn) '_trial_' num2str(ii) '_CCW'],'-v7.3','vidData','t_v');

    end
    
    pause(1)

end

delete(vid)
disp('Done');

clear all;
daqreset;
toc;


