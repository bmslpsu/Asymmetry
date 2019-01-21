
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

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = 'H:\Experiment_Asymmetry\';
root.vid = [root.daq 'Vid\Angles\'];

% Select files
[FILES, PATH.vid] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.vid, 'MultiSelect','on');
FILES = FILES';

PATH.daq = uigetdir(root.daq);
PATH.daq = [PATH.pat '\'];

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
clear WING
% Preallocate data cells %
WINGS.WBA.Pos  = cell(nFly,1);
WINGS.VID.Pos  = cell(nFly,1);
% Save data in cells %
for kk = 1:nTrial
    clear data t_p lAngles rAngles% clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH.daq  FILES{kk}]); % load WBA
	load([PATH.vid   FILES{kk}]); % load VDID
	%-----------------------------------------------------------------------------------------------------------------------------
	
    WING.WBA.Time = t_p;
    
    % Setup filter for wings %
   	Fc = 20; % cutoff frequency [Hz]
    Fs = 1/mean(diff(WING.WBA.Time)); % sampling frequency DAQ [Hz]
  	[b,a] = butter(2,Fc/(Fs/2)); % butterworth filter
    wings.L = data(:,4)
  	Fs = 200; % sampling frequency VID [Hz]
    
    
    head.Pos = filtfilt(b,a,hAngles); % filter head position [deg]
    head.Pos = head.Pos - mean(head.Pos); % subtract DC component

    %-----------------------------------------------------------------------------------------------------------------------------
    % Get pattern data from DAQ %
    pat.Time = t_p; % pattern time
    pat.Pos = data(:,2); % pattern position [deg]
    pat.Fs = 1/mean(diff(pat.Time)); % pattern sampling frequency [Hz]
    [pat.Pos] = FitPanel(pat.Pos,pat.Time,head.Time); % fit pattern posiion
    pat.Vel = diff(pat.Pos)./(1/head.Fs); % pattern velocity
    pat.Vel = [pat.Vel ; pat.Vel(end)]; % pattern velocity
 	%-----------------------------------------------------------------------------------------------------------------------------
   
	%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
    % Head
	WINGS.Time       {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Time;     % head time
	WINGS.Pos        {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Pos;      % head position
    WINGS.Vel        {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Vel;      % head velocity
    WINGS.VelMed     {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.VelMed;   % mean head velocity
    WINGS.VelSTD     {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.VelSTD;   % STD head velocity
    WINGS.Err.Pos    {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Pos;  % head position error
    WINGS.Err.Vel    {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Vel;  % head velocity error
    WINGS.ErrSum.Pos	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Pos;  % head position error sum
    WINGS.ErrSum.Vel	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Vel;  % head velocity error sum
    
    % Pattern
	PAT.Time  	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Time;   % pattern time
	PAT.Pos    	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Pos;    % pattern position
	PAT.Vel    	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Vel;    % pattern velocity
    %-----------------------------------------------------------------------------------------------------------------------------
    if showplot.Time
    rows = 4;
    figure (100)
        subplot(ceil(nTrial/rows),rows,kk) ; hold on
            title(['Fly ' num2str(Fly(kk)) ' Vel ' num2str(Vel(kk))])
            plot(head.Time,pat.Pos,'k')
            plot(head.Time,head.Pos,'b')
            plot(head.Time,head.Err.Pos,'r')
            box on
            hold off
            xlim([0 10])
%             ylim([-1000 1000])
   
    figure (101)
        subplot(ceil(nTrial/rows),rows,kk) ; hold on
            title(['Fly ' num2str(Fly(kk)) ' Vel ' num2str(Vel(kk))])
            plot(head.Time,pat.Vel,'k','LineWidth',3)
            plot(head.Time,head.Vel,'b')
            plot(head.Time,head.Err.Vel,'--r')
            box on
            hold off
            xlim([0 10])
%             ylim([-1000 1000])
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('DONE')

%% Saccade Detection %%
%---------------------------------------------------------------------------------------------------------------------------------
clear Sacs pks 
if showplot.Sacd
figure (34) ; clf
end
varNames = {'Duration','StartPos','PeakPos','EndPos','Amplitude','StartVel','PeakVel','EndVel','StartIndex',...
    'PeakIndex','EndIndex','StartTime','PeakTime','EndTime','MaxPreVelocity','Interval','Rate'};
pp = 1;
rows = 4;
count = 0; % rejected saccades
for kk = 1:nFly
    for jj = 1:nVel
        for ii = 1:size(WINGS.Pos{kk}{jj},2)
            SacdThresh = WINGS.VelMed{kk}{jj}(:,ii) + 3*WINGS.VelSTD{kk}{jj}(:,ii); % threshold for saccade detetcion (by trial)
            
            TimeData = WINGS.Time{kk}{jj}(:,ii); % get time data
            VelData = WINGS.Vel{kk}{jj}(:,ii); % get velocity data
            AbsVelData = abs(VelData); % absolute value fo velocity
            SacdVel = AbsVelData; SacdVel(SacdVel<SacdThresh) = 0; % set all data below threshold to 0
            PosData = WINGS.Pos{kk}{jj}(:,ii); % get position data
         	Fs = 1/mean(diff(TimeData));

             % Find local maxima
            [pks, pIdx] = findpeaks(SacdVel,'MINPEAKDISTANCE',25);
            step = 30; % window length (in samples)
            ta = 1*Fs; % adaptation time (in samples)
            [I] = find(pIdx > (ta) & pIdx < ((length(TimeData) - (step+1)))); % ignore saccades at beginning and end
            pIdx = pIdx(I);
            pks = pks(I);
            clear Sacs
            mm = 1;
            preWin = 15;
            if ~isempty(pks)
                for ww = 1:length(pIdx)
                   % DURATION: define interval as 1/4 peak amplitude
                    sIdx(mm) = find(AbsVelData(1:pIdx(ww)) <= pks(ww)/4,1,'last'); % saccade start index               
                    Eind = find(AbsVelData <= pks(ww)/4); % all values below 1/4 peak
                    Es = find(Eind > pIdx(ww),1,'first'); % first value after the start index is the end idex
                    if ~isempty(Es) % make sure data did not start above threshold
                        eIdx(mm) = Eind(Es); % saccade end index
                    else
                       break;
                       disp('here')
                    end

                    Sacs(mm,1) = 1000*(TimeData(eIdx(mm)) - TimeData(sIdx(mm))); % duration

                    % POSITION & AMPLITUDE: define as ~60 ms (4 samples) window at start and end of interval
                    Sacs(mm,2) =  PosData(sIdx(mm)); % start position
                    Sacs(mm,3) =  PosData(pIdx(ww)); % peak position
                    Sacs(mm,4) =  PosData(eIdx(mm)); % end position
                    Sacs(mm,5) = abs(Sacs(mm,4) - Sacs(mm,2)); % amplitude

                    % VELOCITY
                    Sacs(mm,6) = AbsVelData(sIdx(mm)); % start velocity
                    Sacs(mm,7) = pks(ww); % peak velocity
                    Sacs(mm,8) = AbsVelData(eIdx(mm)); % end velocity

                    % INDEX
                    Sacs(mm,9) = sIdx(mm);  % start index
                    Sacs(mm,10) = pIdx(ww); % peak index
                    Sacs(mm,11) = eIdx(mm); % end index

                    % TIME
                    Sacs(mm,12) = TimeData(sIdx(mm));  % start time
                    Sacs(mm,13) = TimeData(pIdx(mm));  % peak time
                    Sacs(mm,14) = TimeData(eIdx(mm));  % end time

                    % MAX PRE-SACCADE VELOCITY
                    if ww == 1
                        if pIdx(ww)<=preWin
                            Sacs(mm,15) = NaN;
                        else
                            Sacs(mm,15) = abs(median(VelData(pIdx(ww)-preWin:pIdx(ww))));
                        end
                    elseif (abs((pIdx(ww-1)-pIdx(ww))) > preWin) % make sure fly stops before next saccade
                        Sacs(mm,15) = abs(median(VelData(sIdx(mm)-preWin:sIdx(mm))));
                    else
                        Sacs(mm,15) = NaN;
                    end  

                    % INTER-SACCADE INTERVAL (ISI)
                    if ww == 1
                        Sacs(mm,16) = NaN;
                    else
                        Sacs(mm,16) = (pIdx(ww) - pIdx(ww-1));
                    end

                    % RATE
                    zeroSc = Sacs(:,1); % test column to determine # of saccades
                    Sacs(:,17) = size(zeroSc(~isnan(zeroSc)),1);
                    
                    if Sacs(mm,7)>1000
                        Sacs(mm,:) = nan;
                        count = count + 1;
                    end
                    
                    WINGS.SACD.ALL{kk,1}{jj,1}{ii,1} = splitvars(table(Sacs));
                    WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}.Properties.VariableNames = varNames;
                    mm = mm + 1;
                end
            else
                Sacs(1,1:16) = nan;
                Sacs(1,17) = 0;
                WINGS.SACD.ALL{kk,1}{jj,1}{ii,1} = splitvars(table(Sacs));
                WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}.Properties.VariableNames = varNames;
            end
            if showplot.Sacd
            figure (34)
                subplot(rows,ceil(nTrial/rows),pp) ; hold on ; title(['Fly ' num2str(uFly(kk)) ' Vel ' num2str(uVel(jj))])
                    plot( TimeData , AbsVelData ) % velocity data
                    plot( TimeData , SacdThresh*ones(length(TimeData),1) ,'k--') % threshold line
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,13)), ... % peaks
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,7)) ,'b*')
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,12)), ... % start 
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,6)) ,'g*')
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,14)), ... % end 
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,8)) ,'r*')
                    box on
                    hold off
                    
            figure (102)
                subplot(ceil(nTrial/rows),rows,pp) ; hold on ; title(['Fly ' num2str(kk) ' Vel ' num2str(uVel(jj))])
                    plot(WINGS.Time{kk}{jj}(:,ii) , WINGS.Pos{kk}{jj}(:,ii), 'k' )
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,13)), ... % peaks
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,3)) ,'r*')
                    xlim([0 10])
                    ylim([-15 15])
                    box on ; grid on
                    hold off
            end  
            pp = pp + 1; % subplot index
        end
    end
end

%% Saccade Intervals %%
%---------------------------------------------------------------------------------------------------------------------------------

for kk = 1:nFly
    for jj = 1:nVel
        for ii = 1:size(WINGS.Pos{kk}{jj},2)
            TimeData = WINGS.Time{kk}{jj}(:,ii); % get time data
            VelData = WINGS.Vel{kk}{jj}(:,ii); % get velocity data
            PosData = WINGS.Pos{kk}{jj}(:,ii); % get position data
         	Fs = 1/mean(diff(TimeData));

             % Find local maxima
            [pks, pIdx] = findpeaks(SacdVel,'MINPEAKDISTANCE',25);
            step = 30; % window length (in samples)
            ta = 1*Fs; % adaptation time (in samples)
%             [I] = find(pIdx > (ta) & pIdx < ((length(TimeData) - (step+1)))); % ignore saccades at beginning and end
%             pIdx = pIdx(I);
            clear Sacs
            mm = 1;
            preWin = 15;
            if ~isempty(pks)
                for ww = 1:length(pIdx)
                   % DURATION: define interval as 1/4 peak amplitude
                    sIdx(mm) = find(AbsVelData(1:pIdx(ww)) <= pks(ww)/4,1,'last'); % saccade start index               
                    Eind = find(AbsVelData <= pks(ww)/4); % all values below 1/4 peak
                    Es = find(Eind > pIdx(ww),1,'first'); % first value after the start index is the end idex
                    if ~isempty(Es) % make sure data did not start above threshold
                        eIdx(mm) = Eind(Es); % saccade end index
                    else
                       break;
                       disp('here')
                    end

                    Sacs(mm,1) = 1000*(TimeData(eIdx(mm)) - TimeData(sIdx(mm))); % duration

                    % POSITION & AMPLITUDE: define as ~60 ms (4 samples) window at start and end of interval
                    Sacs(mm,2) =  PosData(sIdx(mm)); % start position
                    Sacs(mm,3) =  PosData(pIdx(ww)); % peak position
                    Sacs(mm,4) =  PosData(eIdx(mm)); % end position
                    Sacs(mm,5) =  abs(Sacs(mm,4) - Sacs(mm,2)); % amplitude

                    % VELOCITY
                    Sacs(mm,6) = AbsVelData(sIdx(mm)); % start velocity
                    Sacs(mm,7) = pks(ww); % peak velocity
                    Sacs(mm,8) = AbsVelData(eIdx(mm)); % end velocity

                    % INDEX
                    Sacs(mm,9) = sIdx(mm);  % start index
                    Sacs(mm,10) = pIdx(ww); % peak index
                    Sacs(mm,11) = eIdx(mm); % end index

                    % TIME
                    Sacs(mm,12) = TimeData(sIdx(mm));  % start time
                    Sacs(mm,13) = TimeData(pIdx(mm));  % peak time
                    Sacs(mm,14) = TimeData(eIdx(mm));  % end time

                    % MAX PRE-SACCADE VELOCITY
                    if ww == 1
                        if pIdx(ww)<=preWin
                            Sacs(mm,15) = NaN;
                        else
                            Sacs(mm,15) = abs(median(VelData(pIdx(ww)-preWin:pIdx(ww))));
                        end
                    elseif (abs((pIdx(ww-1)-pIdx(ww))) > preWin) % make sure fly stops before next saccade
                        Sacs(mm,15) = abs(median(VelData(sIdx(mm)-preWin:sIdx(mm))));
                    else
                        Sacs(mm,15) = NaN;
                    end  

                    % INTER-SACCADE INTERVAL (ISI)
                    if ww == 1
                        Sacs(mm,16) = NaN;
                    else
                        Sacs(mm,16) = (pIdx(ww) - pIdx(ww-1));
                    end

                    % RATE
                    zeroSc = Sacs(:,1); % test column to determine # of saccades
                    Sacs(:,17) = size(zeroSc(~isnan(zeroSc)),1);
                    
                    if Sacs(mm,7)>1000
                        Sacs(mm,:) = nan;
                        count = count + 1;
                    end
                    
                    WINGS.SACD.ALL{kk,1}{jj,1}{ii,1} = splitvars(table(Sacs));
                    WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}.Properties.VariableNames = varNames;
                    mm = mm + 1;
                end
            else
                Sacs(1,1:16) = nan;
                Sacs(1,17) = 0;
                WINGS.SACD.ALL{kk,1}{jj,1}{ii,1} = splitvars(table(Sacs));
                WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}.Properties.VariableNames = varNames;
            end
            if showplot.Sacd
            figure (34)
                subplot(rows,ceil(nTrial/rows),pp) ; hold on ; title(['Fly ' num2str(uFly(kk)) ' Vel ' num2str(uVel(jj))])
                    plot( TimeData , AbsVelData ) % velocity data
                    plot( TimeData , SacdThresh*ones(length(TimeData),1) ,'k--') % threshold line
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,13)), ... % peaks
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,7)) ,'b*')
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,12)), ... % start 
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,6)) ,'g*')
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,14)), ... % end 
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,8)) ,'r*')
                    box on
                    hold off
                    
            figure (102)
                subplot(ceil(nTrial/rows),rows,pp) ; hold on ; title(['Fly ' num2str(kk) ' Vel ' num2str(uVel(jj))])
                    plot(WINGS.Time{kk}{jj}(:,ii) , WINGS.Pos{kk}{jj}(:,ii), 'k' )
                    plot(table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,13)), ... % peaks
                        table2array(WINGS.SACD.ALL{kk,1}{jj,1}{ii,1}(:,3)) ,'r*')
                    xlim([0 10])
                    ylim([-15 15])
                    box on ; grid on
                    hold off
            end  
            pp = pp + 1; % subplot index
        end
    end
end

























%% Saccade Means %%
%---------------------------------------------------------------------------------------------------------------------------------
WINGS.SACD.TrialMean = cell(nFly,1);
for kk = 1:nFly
    for jj = 1:nVel
        for ii = 1:size(WINGS.SACD.ALL{kk}{jj},1)
            WINGS.SACD.TrialMean{kk,1}{jj,1}(ii,:) = nanmean(table2array(WINGS.SACD.ALL{kk}{jj}{ii}),1); % fly mean
        end
    end
end

WINGS.SACD.HEAD.FlyMean = cell(nFly,1);
for kk = 1:nFly
    for jj = 1:nVel
%         HEAD.SACD.FlyMean{kk,1}(jj,:) = nanmean(HEAD.SACD.TrialMean{kk}{jj},1); % fly mean
        temp2 = [];
        for ii = 1:size(WINGS.SACD.ALL{kk}{jj},1)
            temp1 = table2array(WINGS.SACD.ALL{kk}{jj}{ii});
            temp2 = [temp2 ; temp1];
            WINGS.SACD.AllTrial{kk,1}{jj,1} = temp2;
            WINGS.SACD.FlyMean{kk,1}(jj,:) = nanmean(temp2,1);
        end
    end
end
WINGS.SACD.GrandMean = nanmean(cat(3,WINGS.SACD.FlyMean{:}),3);  % calculate grand mean at each frequency

% Saccade Mean Plots %%
%---------------------------------------------------------------------------------------------------------------------------------
% Duration, Amplitude, & Peak Velocity
figure (44) ; clf ; hold on
dataIdx = [1 5 7 17];
titleIdx = {'Duration','Amplitude','Peak Velocity','Rate'};
yaIdx = {'ms','deg','deg/s','# / 10s'};
for qq = 1:length(dataIdx)
	subplot(2,2,qq); hold on ; box on ; grid on ; grid minor
    title(titleIdx{qq})
    ylabel(yaIdx{qq}) ; xlabel('Velocity (deg/s)')
    xlim([-ceil(max(uVel)) ceil(max(uVel))])
    for kk = 1:nFly
        for jj = 1:nVel
            for ii = 1:size(WINGS.SACD.ALL{kk}{jj},1)
                sData = table2array(WINGS.SACD.ALL{kk}{jj}{ii}); % induvidual saccades
                plot(uVel(jj),sData(:,dataIdx(qq)),'g*')
                
%                 sData = HEAD.SACD.TrialMean{kk}{jj}; % trial mean
%                 plot(uFreq,sData(ii,dataIdx(qq)),'b*')
            end
        end
        plot(uVel , WINGS.SACD.FlyMean{kk}(:,dataIdx(qq)) ,'-*') % fly means
    end
    plot(uVel , WINGS.SACD.GrandMean(:,dataIdx(qq)) ,'k-o','LineWidth',4,'MarkerSize',5) % grand means
end
       
%% FUNCTION: FitPanel
function [x_out] = FitPanel(x_in,time,time_new)
%% ---------------------------------------------------------------------------------------------------------------------------------
% FitPanel: unwraps panel position and fits to straight line
    % INPUTS:
        % x_in: pattern position
        % time: pattern time 
        % time_new: new ouput time
	% OUTPUTS
        % x_out: fit data
debug = 0;
% x_in = data(:,2);
% time = pat.Time;
% time_new = head.Time;
% Fs = pat.Fs;
% t = 1:(1/Fs):max(t_out);
% t = pat.Time;
%---------------------------------------------------------------------------------------------------------------------------------      
% Fs = mean(diff(t_in))
% close all ; clc

x.deg   = 3.75*round((96/10)*(x_in)); % convert to deg
x.norm  = x.deg - x.deg(1); % normalize to start at 0
x.uw = x.norm;

% Find when panel position changes (transition)
pp = 1;
for jj = 2:length(x.norm)
    dP = x.norm(jj) - x.norm(jj-1);
    if abs(dP) > 1
%         mrkIdx(pp) = jj-1;
        mrkPos(pp) = x.norm(jj-1);
    	mkrTime(pp) = time(jj-1);
        pp = pp + 1;
    end
end
% indexMean = nan(1,length(mrkIdx));
% xPosMean = nan(1,length(mrkIdx));
% for jj = 2:length(mrkIdx)
%     indexMean(jj-1) = round(mean(mrkIdx(jj-1):1:mrkIdx(jj)));
%     xPosMean(jj-1)  = x.norm  (indexMean(jj-1));
%   	timeMean(jj-1)  = time   (indexMean(jj-1));
% end

dx.all  = diff(x.norm);
jump.idx = find(abs(dx.all)>180);
dir = median(diff(mrkPos));
for kk = 1:length(jump.idx)
    jump.mag = x.norm(jump.idx+1) - x.norm(jump.idx);
    x.uw((jump.idx(kk)+1):end) = x.uw((jump.idx(kk)+1):end) - jump.mag(kk) + dir;    
end

m = (x.uw(end) - x.uw(1))/(10 - 0);
% t = linspace(time_new(1),time_new(length(time_new)),length(time_new))';
t = linspace(0,10,length(time_new))';
% t = t - time_new(1);
x_out = m*t;

if debug
figure; clf ; hold on
%     idx = 1:length(x_in);
    plot(time,x.norm,'k')
    plot(mkrTime,mrkPos,'*r')
    plot(time,x.uw,'b')
    plot(time_new,x_out,'c')
    pause
end
%---------------------------------------------------------------------------------------------------------------------------------
end






