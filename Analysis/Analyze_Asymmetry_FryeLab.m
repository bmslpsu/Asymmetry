function [] = Analyze_Asymmetry_FryeLab()
%% Analyze_Asymmetry_FryeLab: calculates WBA in response to CW & CCW
% ramps, compares DAQ & VIDEO measurments for new data
%   INPUTS:
%
%   OUTPUTS:
%
%---------------------------------------------------------------------------------------------------------------------------------
% User sets these variables %
showplot.Time = 0; % shows all WBA trials when loading data
showplot.Freq = 0; % shows all WBF trials when loading data
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
% root.daq = 'Q:\Box Sync\Research\Asymmetry\FryeLab\';
root.daq = 'H:\EXPERIMENTS\Experiment_Asymmetry_Frye\';
% root.vid = [root.daq 'Vid\Angles\'];

% Select VIDEO angle files & set DAQ file directory
[FILES, PATH.daq] = uigetfile({'*.mat', 'DAQ-files'}, 'Select VIDEO angle files', root.daq, 'MultiSelect','on');
FILES = FILES';
%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
nTrial = length(FILES); % total # of trials
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Vel = zeros(nTrial,1); WINGS.FileCells = cell(nTrial,6); % preallocate arrays
for jj = 1:nTrial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    WINGS.FileCells(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}};  % separate strings columns with one item in each cell
    Fly(jj,1) = str2double(temp{2});   % store fly #
    Trial(jj,1) = str2double(temp{4}); % store trial #
    
    % Store velocity
    if strcmp(temp{5},'CW')
        Vel(jj,1) =  45;
    elseif strcmp(temp{5},'CCW')
        Vel(jj,1) = -45;
    end
end
%% Set up indexing convention for files %%
% Normalize to start at fly#1 and increment by 1 for each fy (same for velcoity & trials)
%---------------------------------------------------------------------------------------------------------------------------------
uFly = sort(unique(Fly)); % # of each unique flies
nFly = length(uFly); % # of flies
trialFly = cell(nFly,2); % preallocate cell to store # of trials per fly
for kk = 1:nFly
    trialFly{kk,1} = uFly(kk);  % fly # label
	trialFly{kk,2} = length(find(Fly==uFly(kk))); % # trials per fly
end
% Make indexing array for flies
newFly = 1:nFly; % new fly #'s
idxFly = zeros(length(Fly),1);
pp = 1;
for kk = 1:nFly
    idxFly(pp:pp+trialFly{kk,2}-1,1) = newFly(kk); % fly index
    pp = pp+trialFly{kk,2};
end
% Make indexing array for velocities
uVel = sort(unique(Vel)); % find all unique velocity
nVel = length(uVel); % # of unique velocity
idxVel = zeros(nTrial,1);
for kk = 1:nVel
   idxVel(Vel == uVel(kk)) = kk; % velocity index
end

% Show user fly-trial stats
fprintf('Total Flies''s: %i \n',nFly) ; fprintf('Total Trials''s: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear WINGS
wings.daq.EI = 50000; % end index for DAQ files
% wings.vid.EI = 1000;  % end index for VID files
% Preallocate data cells
for kk = 1:nFly
    WINGS.daq.wba{kk,1} = cell(nVel,1);
    WINGS.vid.wba{kk,1} = cell(nVel,1);
end
WINGS.daq.time = (0:(1/5000):10)'; % time vector for DAQ
WINGS.daq.time = WINGS.daq.time(1:end-1); % time vector for DAQ
% WINGS.vid.time = (0:(1/100):10)';  % time vector for VID
% WINGS.vid.time = WINGS.vid.time(1:end-1);
% Save data in cells
for kk = 1:nTrial
    clear data t_p lAngles rAngles % clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH.daq  FILES{kk}],'data','t_p'); % load WBA
% 	load([PATH.vid   FILES{kk}],'lAngles','rAngles'); % load VID
	%-----------------------------------------------------------------------------------------------------------------------------
  	% Check wing beat frequency %
    wings.F = data(1:wings.daq.EI,7);
    if min(wings.F)<1.40
        disp(['Low WBF:  Fly ' num2str(Fly(kk)) ' Trial ' num2str(Trial(kk))])
        continue
    end
	%-----------------------------------------------------------------------------------------------------------------------------
  	% Get wing data from DAQ %
   	Fc = 20; % cutoff frequency [Hz]
    Fs = 5000; % sampling frequency [Hz]
  	[b,a] = butter(2,Fc/(Fs/2)); % butterworth filter 
    wings.daq.L = filtfilt(b,a,data(1:wings.daq.EI,5))'; % filter left wing
    wings.daq.R = filtfilt(b,a,data(1:wings.daq.EI,6))'; % filter right wing
    wings.daq.wba = filtfilt(b,a,wings.daq.L - wings.daq.R); % filter L-R wing
    wings.daq.wba = wings.daq.wba - wings.daq.wba(1); % normalize
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
	WINGS.daq.wba 	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = wings.daq.wba;
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
%% DAQ Figure %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (4) ; clf ; hold on ; box on ; ylim(2*[-1 1])
set(gcf,'Color','w')
xlabel('Time (s)')
ylabel('V')
WINGS.daq.FlyMean = [];
WINGS.daq.GrandMean = [];
for kk = 1:nFly
    for jj = 1:nVel
        WINGS.daq.FlyMean.wba{kk,1}(:,jj) = mean(WINGS.daq.wba{kk}{jj},2); % fly mean
        
        % plot all trials
        if jj==1
            h = plot(WINGS.daq.time,WINGS.daq.wba{kk}{jj},'r','LineWidth',1);
        elseif jj==2
            h = plot(WINGS.daq.time,WINGS.daq.wba{kk}{jj},'b','LineWidth',1);
        end
        
        for ii = 1:length(h)
            h(ii).Color(4) = 0.1;
        end
    end
    WINGS.daq.FlyMean.off{kk,1} = WINGS.daq.FlyMean.wba{kk,1}(:,1) + WINGS.daq.FlyMean.wba{kk,1}(:,2);
end

% plot fly  means
for kk = 1:nFly
    h1 = plot(WINGS.daq.time,WINGS.daq.FlyMean.wba{kk,1}(:,1),'r','LineWidth',3);
    h2 = plot(WINGS.daq.time,WINGS.daq.FlyMean.wba{kk,1}(:,2),'b','LineWidth',3);

    h1.Color(4) = 0.3;
    h2.Color(4) = 0.3;
end

for kk = 1:nFly
	plot(WINGS.daq.time,WINGS.daq.FlyMean.off{kk,1},'k','LineWidth',2)
end

% grand means
WINGS.daq.GrandMean.off = mean(cat(3,WINGS.daq.FlyMean.off{:}),3);
WINGS.daq.GrandMean.wba = mean((cat(3,WINGS.daq.FlyMean.wba{:})),3);
% grand STD
WINGS.daq.GrandSTD.off = std(cat(3,WINGS.daq.FlyMean.off{:}),0,3);
WINGS.daq.GrandSTD.wba = std((cat(3,WINGS.daq.FlyMean.wba{:})),0,3);

% plot grand means
[~,~] = PlotPatch(WINGS.daq.GrandMean.wba(:,1),WINGS.daq.GrandSTD.wba(:,1),WINGS.daq.time,1,nFly,'r',[0.5 0.5 0.5],0.8,5);
[~,~] = PlotPatch(WINGS.daq.GrandMean.wba(:,2),WINGS.daq.GrandSTD.wba(:,2),WINGS.daq.time,1,nFly,'b',[0.5 0.5 0.5],0.8,5);
[~,~] = PlotPatch(WINGS.daq.GrandMean.off,(WINGS.daq.GrandSTD.wba(:,1) + WINGS.daq.GrandSTD.wba(:,2)),WINGS.daq.time,1,nFly,...
    'k',[0.5 0.5 0.5],0.8,6);
plot(WINGS.daq.time,0*WINGS.daq.time,'-g','LineWidth',1)
end