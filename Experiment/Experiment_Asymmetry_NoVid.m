function [] = Experiment_Asymmetry_NoVid(Fn)
%% Experiment_Asymmetry_NoVid: runs an experiment using the LED arena and fly panel
% Fn is the fly number
% This code is written for Panel Controller v3 and NiDAQ seesion mode
%---------------------------------------------------------------------------------------------------------------------------------
daqreset
imaqreset

%% Set directories & experimental paramters %%
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'E:\Experiment_Asymmetry_Verification\';

% EXPERIMENTAL PARAMETERS
n_tracktime = 10;           % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 5;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_trial = 20;               % # trials per fly
n_AI = 6;                   % # of analog input channels

%% Set up data acquisition on NiDAQ (session mode) %%
%---------------------------------------------------------------------------------------------------------------------------------
devices = daq.getDevices; % get DAQ ID
s = daq.createSession('ni'); % create session

% Add analog input channels
ch.AI = addAnalogInputChannel(s,devices.ID, 1:n_AI, 'Voltage');
% Assign AI channel names
chNames = {'Trigger','Pat x-pos','Pat y-pos','L-WBA','R-WBA','WB-Freq'};
for kk = 1:length(ch.AI)
    ch.AI(kk).Name = chNames{kk};
        % AI1 = AO_0: trigger
        % AI2 = DAC1: x-position of stimulus
        % AI3 = DAC2: y-position of stimulus
        % AI4 = L: left wing beat amplitude
        % AI5 = R: right wing beat amplitude
        % AI6 = FREQ: frequency
end

% Add analog output channels
ch.AO = addAnalogOutputChannel(s,devices.ID,'ao0','Voltage');
ch.AO.Name = 'Trigger';

% Setup Sampling
s.Rate = 5000; % samples per second
s.IsContinuous = false;	% continuous data collection until stopped

FrameRate = 100; % camera frame rate [Hz]

t = 0:1/s.Rate:n_tracktime;
TriggerSignal = (square(2*pi*FrameRate*t,90) + 1)';

%% Set Experimental Sequence %%
%---------------------------------------------------------------------------------------------------------------------------------
% Random sequence of -1 and 1  
Rdir = [1 1 -1 1 -1 1 -1 -1 1 -1 ];
Rdir = repmat(Rdir,1,2);

%% EXPERIMENT LOOP %%
%---------------------------------------------------------------------------------------------------------------------------------
for ii = 1:n_trial      
    disp('Trial')
    disp(num2str(ii))
    
    % CLOSED LOOP BAR TRACKING %
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)           % set output to p_rest_pat (Pattern_Fourier_bar_barwidth=8)
    Panel_com('set_position',[1, 1]); pause(n_pause)        % set starting position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)             % closed loop tracking [xpos,ypos] (NOTE: 0=open, 1=closed)
    Panel_com('set_posfunc_id',[1,0]); pause(n_pause)       % no position function for x-channel
    Panel_com('set_funcX_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('set_posfunc_id',[2,0]); pause(n_pause)       % no position function for y-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)	% [xgain,xoffset,ygain,yoffset]
    Panel_com('start');                                     % start closed-loop rest
    pause(n_resttime)                                       % rest 5 seconds
    Panel_com('stop');                                      % stop closed-loop rest
       
	% EXPERIMENT SETUP %
    pause(1)
    disp(['Play Stimulus: dir = ' num2str(Rdir(ii))])
    Panel_com('set_pattern_id', 2); pause(n_pause)                	% set pattern to "Pattern_Random_Ground_48"
    Panel_com('set_position',[randi([1, 96]), 1]); pause(n_pause)  	% set starting position (xpos,ypos)
    Panel_com('set_posfunc_id',[1, 0]); pause(n_pause)            	% arg1 = channel (x=1,y=2); arg2 = funcid
	Panel_com('set_funcX_freq', 50); pause(n_pause)                 % 50Hz update rate for x-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)              	% 50Hz update rate for y-channel
    Panel_com('set_mode', [0,0]); pause(n_pause)                    % 0=open,1=closed,2=fgen,3=vmode,4=pmode
    Panel_com('send_gain_bias',[Rdir(ii)*12,0,0,0]); pause(n_pause)	% open loop at 45 deg/s
	
    % START EXPERIMENT & DATA COLLECTION %
    queueOutputData(s,TriggerSignal) % set trigger AO signal
    tic
        Panel_com('start') 	% start stimulus
        [data, t_p ] = s.startForeground; % data collection
        Panel_com('stop')  	% stop stimulus
    toc
 	
    % CLOSED LOOP BAR TRACKING %
    disp('rest');
    Panel_com('stop'); pause(n_pause)
    Panel_com('set_pattern_id', 1);pause(n_pause)           % set output to p_rest_pat (Pattern_Fourier_bar_barwidth=8)
    Panel_com('set_position',[1, 1]); pause(n_pause)        % set starting position (xpos,ypos)
    Panel_com('set_mode',[1,0]); pause(n_pause)             % closed loop tracking [xpos,ypos] (NOTE: 0=open, 1=closed)
    Panel_com('set_posfunc_id',[1,0]); pause(n_pause)       % no position function for x-channel
    Panel_com('set_funcX_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('set_posfunc_id',[2,0]); pause(n_pause)       % no position function for y-channel
    Panel_com('set_funcY_freq', 50); pause(n_pause)         % 50Hz update rate
    Panel_com('send_gain_bias',[-15,0,0,0]); pause(n_pause)	% [xgain,xoffset,ygain,yoffset]
	Panel_com('start');                                     % start closed-loop rest
	
    % Save data %
    disp('Saving...') ; disp('----------------------------------------------------------------------')
    if Rdir(ii) == 1
        fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_CW.mat'];
        save([rootdir   fname],'-v7.3','data','t_p');
    elseif Rdir(ii) == -1
        fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_CCW.mat'];
        save([rootdir   fname],'-v7.3','data','t_p');
    end
end
disp('Done');

daqreset
imaqreset
end