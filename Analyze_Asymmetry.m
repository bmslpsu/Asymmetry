
%% Analyze Asymmetry: Summary of this function goes here
%   Detailed explanation goes here
% INPUTS:
    % root: root directory
    % Amp: amplitude  of stimulus [deg]
    % HeadWings: head=1, wings=0
    % plotTime: show time domain plot (logical)
    % plotFreq: show frequency domain plots (logical)
% clear;close all;clc
clear
showplot.Time = 0;
showplot.Freq = 0;

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = 'H:\Experiment_Asymmetry\';
root.vid = [root.daq 'Vid\Angles\'];

% Select files
[FILES, PATH.vid] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.vid, 'MultiSelect','on');
FILES = FILES';

PATH.daq = uigetdir(root.daq);
PATH.daq = [PATH.daq '\'];

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% clearvars -except FILES PATH Amp
% close all
clc

nTrial = length(FILES); % total # of trials
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Vel = zeros(nTrial,1); WINGS.FileCells = cell(nTrial,6);% preallocate arrays
for jj = 1:nTrial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    WINGS.FileCells(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}}; % separate strings into six rows with one item in each cell
    Fly(jj,1) = str2double(temp{2}); % store fly #
    Trial(jj,1) = str2double(temp{4}); % store trial #
    
    if strcmp(temp{5},'CW')
        Vel(jj,1) =  45;
    elseif strcmp(temp{5},'CCW')
        Vel(jj,1) = -45;
    end
end

%% Set up indexing convention for files %%
% Normalize to start at fly#1 and increment by 1 for each fy (same for trials)
%---------------------------------------------------------------------------------------------------------------------------------
uFly = sort(unique(Fly)); % original # of each unique fly
nFly = length(uFly); % # flies
trialFly = cell(nFly,2);
for kk = 1:nFly
    trialFly{kk,1} = uFly(kk);  % fly # label
	trialFly{kk,2} = length(find(Fly==uFly(kk))); % # trials per fly
end
% Make indexing array for flies
newFly = 1:nFly;
idxFly = zeros(length(Fly),1);
pp = 1;
for kk = 1:nFly
    idxFly(pp:pp+trialFly{kk,2}-1,1) = newFly(kk); % fly index
    pp = pp+trialFly{kk,2};
end

% Make indexing array for velocity
uVel = sort(unique(Vel)); % find all unique velocity
nVel = length(uVel); % # of unique velocity
idxVel = zeros(nTrial,1);
for kk = 1:nVel
   idx = find(Vel == uVel(kk));
   idxVel(idx) = kk; % velocity index
end

fprintf('Total Flies''s: %i \n',nFly) ; fprintf('Total Trials''s: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear WINGS
wings.daq.EI = 10001;
wings.vid.EI = 1000;
% Preallocate data cells %
for kk = 1:nFly
    WINGS.daq.wba{kk,1} = cell(nVel,1);
    WINGS.vid.wba{kk,1} = cell(nVel,1);
end
WINGS.daq.time = 0:(1/1000):10;
WINGS.vid.time = 0:(1/100):10;
WINGS.vid.time = WINGS.vid.time(1:end-1);
figure (1) ; clf
figure (2) ; clf
% Save data in cells %
for kk = 1:nTrial
    clear data t_p lAngles rAngles% clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH.daq  FILES{kk}]); % load WBA
	load([PATH.vid   FILES{kk}]); % load VDID
	%-----------------------------------------------------------------------------------------------------------------------------
  	% Get wing data from daq %
   	Fc = 20; % cutoff frequency [Hz]
    Fs = 1000; % sampling frequency [Hz]
  	[b,a] = butter(2,Fc/(Fs/2)); % butterworth filter 
    wings.daq.L = filtfilt(b,a,data(4,1:wings.daq.EI))';
    wings.daq.R = filtfilt(b,a,data(5,1:wings.daq.EI))';
    wings.daq.wba = filtfilt(b,a,wings.daq.L - wings.daq.R);
    wings.daq.wba = wings.daq.wba - wings.daq.wba(1);
 	%-----------------------------------------------------------------------------------------------------------------------------
  	% Get wing data from vid %
   	Fc = 20; % cutoff frequency [Hz]
    Fs = 200; % sampling frequency [Hz]
  	[b,a] = butter(2,Fc/(Fs/2)); % butterworth filter 
    wings.vid.L = rAngles(1:wings.vid.EI)';
    wings.vid.R = lAngles(1:wings.vid.EI)';
    
    X = linspace(1, length(wings.vid.L), length(wings.vid.L));
    wings.vid.L = filtfilt(b,a,hampel(X, wings.vid.L, 50, 4));
	wings.vid.R = filtfilt(b,a,hampel(X, wings.vid.R, 50, 4));

    wings.vid.wba = filtfilt(b,a,wings.vid.L - wings.vid.R);
   	wings.vid.wba = wings.vid.wba - wings.vid.wba(1);
 	%-----------------------------------------------------------------------------------------------------------------------------
  	% Check wing beat frequency %
    wings.F = data(6,1:wings.daq.EI)';
    if min(wings.F)<1.70
        disp(['Low WBF:  Fly ' num2str(Fly(kk)) ' Trial ' num2str(Trial(kk))])
        continue
    end
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
    % Head
	WINGS.daq.wba 	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = wings.daq.wba;
	WINGS.vid.wba  	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = wings.vid.wba;
    %-----------------------------------------------------------------------------------------------------------------------------
    if showplot.Time
    figure (1) ; hold on ; box on
        xlabel('Time (s)')
        ylabel('V')
        if Vel(kk)>0
            plot(WINGS.daq.time,wings.daq.wba,'b')
        elseif Vel(kk)<0
            plot(WINGS.daq.time,wings.daq.wba,'r')
        end
    figure (2) ; hold on ; box on
        xlabel('Time (s)')
        ylabel('deg')
        if Vel(kk)>0
            plot(WINGS.vid.time,wings.vid.wba,'b')
        elseif Vel(kk)<0
            plot(WINGS.vid.time,wings.vid.wba,'r')
        end
    end
    
	if showplot.Freq
	figure (3) ; hold on ; box on
        xlabel('Time (s)')
        ylabel('WB Frequency (Hz)')
        if Vel(kk)>0
            plot(WINGS.daq.time,wings.F,'b')
        elseif Vel(kk)<0
            plot(WINGS.daq.time,wings.F,'r')
        end
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('DONE')

%% Averages & Figures %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
figure (4) ; clf ; hold on ; box on ; ylim([-5 5])
xlabel('Time (s)')
ylabel('V')
WINGS.daq.FlyMean = [];
WINGS.daq.GrandMean = [];
for kk = 1:nFly
    for jj = 1:nVel
        WINGS.daq.FlyMean.wba{kk,1}(:,jj) = mean(WINGS.daq.wba{kk}{jj},2);
        figure (4) ; hold on
        if jj==1
            h = plot(WINGS.daq.time,WINGS.daq.wba{kk}{jj},'r','LineWidth',1);
        elseif jj==2
            h = plot(WINGS.daq.time,WINGS.daq.wba{kk}{jj},'b','LineWidth',1);
        end
        
        for ii = 1:length(h)
            h(ii).Color(4) = 0.2;
        end
    end
    WINGS.daq.FlyMean.off{kk,1} = WINGS.daq.FlyMean.wba{kk,1}(:,1) + WINGS.daq.FlyMean.wba{kk,1}(:,2);
end

for kk = 1:nFly
    figure (4) ; hold on
        h1 = plot(WINGS.daq.time,WINGS.daq.FlyMean.wba{kk,1}(:,1),'r','LineWidth',3);
        h2 = plot(WINGS.daq.time,WINGS.daq.FlyMean.wba{kk,1}(:,2),'b','LineWidth',3);
        
        h1.Color(4) = 0.3;
        h2.Color(4) = 0.3;
end

% for kk = 1:nFly
%     figure (4) ; hold on
%      	plot(WINGS.daq.time,WINGS.daq.FlyMean.off{kk,1},'k','LineWidth',2)
% end

WINGS.daq.GrandMean.off = mean(cat(3,WINGS.daq.FlyMean.off{:}),3);
WINGS.daq.GrandMean.wba = mean((cat(3,WINGS.daq.FlyMean.wba{:})),3);

figure (4) ; hold on
	plot(WINGS.daq.time,WINGS.daq.GrandMean.wba(:,1),'r','LineWidth',7)
  	plot(WINGS.daq.time,WINGS.daq.GrandMean.wba(:,2),'b','LineWidth',7)
	plot(WINGS.daq.time,WINGS.daq.GrandMean.off,'k','LineWidth',10)
	plot(WINGS.daq.time,0*WINGS.daq.time,'--g','LineWidth',2)
    
%% ---------------------------------------------------------------------------------------------------------------------------------    
figure (5) ; clf ; hold on ; box on ; ylim([-50 50])
xlabel('Time (s)')
ylabel('deg')
WINGS.vid.FlyMean = [];
WINGS.vid.GrandMean = [];
for kk = 1:nFly
    for jj = 1:nVel
        WINGS.vid.FlyMean.wba{kk,1}(:,jj) = mean(WINGS.vid.wba{kk}{jj},2);
        figure (5) ; hold on
        if jj==1
            h = plot(WINGS.vid.time,WINGS.vid.wba{kk}{jj},'r','LineWidth',1);
        elseif jj==2
            h = plot(WINGS.vid.time,WINGS.vid.wba{kk}{jj},'b','LineWidth',1);
        end
        
        for ii = 1:length(h)
            h(ii).Color(4) = 0.2;
        end
    end
    WINGS.vid.FlyMean.off{kk,1} = WINGS.vid.FlyMean.wba{kk,1}(:,1) + WINGS.vid.FlyMean.wba{kk,1}(:,2);
end

for kk = 1:nFly
    figure (5) ; hold on
        h1 = plot(WINGS.vid.time,WINGS.vid.FlyMean.wba{kk,1}(:,1),'r','LineWidth',3);
        h2 = plot(WINGS.vid.time,WINGS.vid.FlyMean.wba{kk,1}(:,2),'b','LineWidth',3);
        
        h1.Color(4) = 0.3;
        h2.Color(4) = 0.3;
end

% for kk = 1:nFly
%     figure (5) ; hold on
%      	plot(WINGS.vid.time,WINGS.vid.FlyMean.off{kk,1},'k','LineWidth',2)
% end

WINGS.vid.GrandMean.off = mean(cat(3,WINGS.vid.FlyMean.off{:}),3);
WINGS.vid.GrandMean.wba = mean((cat(3,WINGS.vid.FlyMean.wba{:})),3);

figure (5) ; hold on
	plot(WINGS.vid.time,WINGS.vid.GrandMean.wba(:,1),'r','LineWidth',7)
  	plot(WINGS.vid.time,WINGS.vid.GrandMean.wba(:,2),'b','LineWidth',7)
	plot(WINGS.vid.time,WINGS.vid.GrandMean.off,'k','LineWidth',10)
	plot(WINGS.vid.time,0*WINGS.vid.time,'--g','LineWidth',2)
    
    
