function [] = Analyze_Asymmetry_Control(subRows,subCurrent)
%---------------------------------------------------------------------------------------------------------------------------------
% Analyze_Asymmetry_Control: calculates WBA in response to CW & CCW ramps, compares DAQ & VIDEO measurments
    % INPUTS:
        % -
    % OUTPUTS:
        % -
%---------------------------------------------------------------------------------------------------------------------------------
% User sets these variables %
showplot.Time = 0; % shows all WBA trials when loading data
showplot.Freq = 0; % shows all WBF trials when loading data
% subRows = 4;
% subCurrent = 2;
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
% root = 'H:\Experiment_Asymmetry_Control\HighContrast\';
root = 'H:\Experiment_Asymmetry_Control\LowContrast\';
% Select VIDEO angle files & set DAQ file directory
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'}, 'Select DAQ files', root, 'MultiSelect','on');
FILES = FILES';
%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
temp = textscan(PATH, '%s', 'delimiter','\'); % get spatial frequency (0 = Random)
spatFreq = temp{1}{4};
if strcmp(spatFreq,'0')
    spatFreq = 'Random';
else
    spatFreq = [spatFreq ' ' char(176)];
end

nTrial = length(FILES); % total # of trials
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Vel = zeros(nTrial,1); Dir = zeros(nTrial,1); % preallocate arrays
WINGS.FileCells = cell(nTrial,6);
for jj = 1:nTrial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    WINGS.FileCells(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}};  % separate strings columns with one item in each cell
    Fly(jj,1) = str2double(temp{2});            % store fly #
    Trial(jj,1) = str2double(temp{4});          % store trial #
    Vel(jj,1) = 3.75*abs(str2double(temp{6}));  % store velocity
    
    % Store direction
    if strcmp(temp{5},'CW')
        Dir(jj,1) = 2;
    elseif strcmp(temp{5},'CCW')
        Dir(jj,1) = 1;
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

nDir = length(unique(Dir)); % # directions

% Show user fly-trial stats
fprintf('Total Flies''s: %i \n',nFly) ; fprintf('Total Trials''s: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear WINGS
wings.daq.EI = 9900; % end index for DAQ files
% Preallocate data cells
for kk = 1:nFly
    WINGS.daq.wba{kk,1} = cell(nVel,2);
end
% WINGS.daq.time = (0:(1/1000):10)'; % time vector for DAQ
WINGS.daq.time = (linspace(0,10,wings.daq.EI))'; % time vector for DAQ
% Save data in cells
for kk = 1:nTrial
    clear WBAdata t_p  % clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH  FILES{kk}],'WBAdata','t_p'); % load WBA
	%-----------------------------------------------------------------------------------------------------------------------------
  	% Check wing beat frequency %
    wings.F = WBAdata(6,1:wings.daq.EI)';
    if min(wings.F)<1.70
        disp(['Low WBF:  Fly ' num2str(Fly(kk)) ' Trial ' num2str(Trial(kk))])
%         continue
    end
	%-----------------------------------------------------------------------------------------------------------------------------
  	% Get wing data from DAQ %
   	Fc = 20; % cutoff frequency [Hz]
    Fs = 1000; % sampling frequency [Hz]
  	[b,a] = butter(2,Fc/(Fs/2)); % butterworth filter 
    wings.daq.L = filtfilt(b,a,WBAdata(4,1:wings.daq.EI))'; % filter left wing
    wings.daq.R = filtfilt(b,a,WBAdata(5,1:wings.daq.EI))'; % filter right wing
    wings.daq.wba = filtfilt(b,a,wings.daq.L - wings.daq.R); % filter L-R wing
    wings.daq.wba = wings.daq.wba - wings.daq.wba(1); % normalize
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
	WINGS.daq.wba 	{idxFly(kk),1}{idxVel(kk),Dir(kk)}(:,end+1) = wings.daq.wba;
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

%% Response figure for varying speeds %%
%---------------------------------------------------------------------------------------------------------------------------------
WINGS.daq.FlyMean = [];
WINGS.daq.FlyMeanVel = [];
WINGS.daq.GrandMean = [];
for kk = 1:nFly
    WINGS.daq.FlyMean.wba{kk,1} = cellfun(@(x) mean(x,2),WINGS.daq.wba{kk},'UniformOutput',false);
end

for kk = 1:nFly
    for jj = 1:nVel
        for ii = 1:nDir
            WINGS.daq.FlyMeanVel.wba{jj,ii}(:,kk) = WINGS.daq.FlyMean.wba{kk,1}{jj,ii}; % temp to store fly mean data until mean is taken
        end
    end
end

WINGS.daq.GrandSTD.wba  = cellfun(@(x) std(x,0,2),WINGS.daq.FlyMeanVel.wba,'UniformOutput',false);
WINGS.daq.GrandMean.wba = cellfun(@(x) mean(x,2),WINGS.daq.FlyMeanVel.wba,'UniformOutput',false);

for jj = 1:nVel
    WINGS.daq.GrandMean.wba{jj,3} = WINGS.daq.GrandMean.wba{jj,2} + WINGS.daq.GrandMean.wba{jj,1}; 
end

F = figure (3);
set(F, 'Position', [100, 100, 1200, 200*subRows])
for jj = 1:nVel
    subidx = jj + nVel*(subCurrent-1);
    h = subplot(subRows,nVel,subidx) ; hold on ; box on
        h.FontSize = 12;
        
        if (subCurrent==1)
            title([num2str(uVel(jj)) ' ' char(176) '/s']);
        end
        
        if (jj==3) && (subRows==subCurrent)
            xlabel('Time (s)')
        end
        
        if jj==1
            ylabel({spatFreq;'V'})
        else
            yticklabels('')
        end
        
        if ~(subRows==subCurrent)
            xticklabels('')
        end
                
        xlim([0 10])
        xticks([0 10])
        ylim([-6 6])
        yticks([-5 5])
    
    [~,~] = PlotPatch(WINGS.daq.GrandMean.wba{jj,1},WINGS.daq.GrandSTD.wba{jj,1},WINGS.daq.time,1,nFly,'r',[0.5 0.5 0.5],0.8,5);
    [~,~] = PlotPatch(WINGS.daq.GrandMean.wba{jj,2},WINGS.daq.GrandSTD.wba{jj,2},WINGS.daq.time,1,nFly,'b',[0.5 0.5 0.5],0.8,5);
    [~,~] = PlotPatch(WINGS.daq.GrandMean.wba{jj,3},(WINGS.daq.GrandSTD.wba{jj,1}+ WINGS.daq.GrandSTD.wba{jj,2}),...
        WINGS.daq.time,1,nFly,'k',[0.5 0.5 0.5],0.8,5);
    plot(WINGS.daq.time,0*WINGS.daq.time,'-g','LineWidth',1)
end
end