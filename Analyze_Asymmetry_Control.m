%%%% Analysis Script for Spatial Frequency Control Experiments (Rigid Tether-WBA) %%%%
%---------------------------------------------------------------------------------------------------------------------------------
%% Set directories %%
clear ; close all; clc
root = 'C:\Users\boc5244\Box Sync\Research\Asymmetry\Experiment_Asymmetry_Control';

[Files, dirpath] = uigetfile({'*.mat', 'MAT-files'},'Select fly trials', root, 'MultiSelect','on');

% Sort files naturally based on fly# & trial# (external function)
% Files = natsortfiles (Files);
nTrial = length(Files);
fprintf('Total Trials: %i \n',nTrial)
%---------------------------------------------------------------------------------------------------------------------------------
%% Parse filename data %%
clearvars -except root dirpath Files nTrial
clc
% Read in data from file names: [Fly=FLy#, Trial=trial# for each fly,Vel=velocity of stimulus]
% Preallocate file data arrays
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Vel = zeros(nTrial,1); Dir = zeros(nTrial,1); 
% Store file data in cells
for jj = 1:nTrial
    temp = textscan(char(Files(jj)), '%s', 'delimiter', '_.'); temp = temp{1}; % read individual strings into temp variable
    C(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6} temp{7}}; % separate strings into 7 rows with 1 word in each cell
    Fly(jj,1) = str2double(temp{2}); % store fly number
    Trial(jj,1) = str2double(temp{4}); % store trial number
    Vel(jj,1) = str2double(temp{6}); % store stimulus velocity (WBA gain) 
    % Store direction {1=CW,2=CCW}'
    if sign(Vel(jj,1))==1
        Dir(jj,1) = 1;
    elseif sign(Vel(jj,1))==-1
        Dir(jj,1) = 2;
    end
end

nFly   = length(unique(Fly)); fprintf('# Fly''s: %i \n',nFly) % # of flies
nSpeed = length(unique(abs(Vel))); % # of speeds (unsigned)
nDir   = length(unique(Dir)); % # of directions

% Set up new naming convention for files %%
% Renaming convention: normalize so we start at fly#1 and increment by 1 for each fy
normFly = Fly - min(Fly) + 1;   % Starts analysis at fly#1 even if file name specifies different #
nFly = length(unique(Fly));     % Number of flies being analyzed
uFly = unique(Fly);             % Fly number of each unique fly
trialFly = cell(nFly,2);
for kk=1:nFly
    trialFly{kk,1} = uFly(kk);  % Fly number
	trialFly{kk,2} = length(find(Fly==uFly(kk))); % Number of trials per fly
end

newFly = 1:nFly;
indexFly = zeros(length(Fly),1);
pp = 1;
for kk=1:nFly
    indexFly(pp:pp+trialFly{kk,2}-1,1) = newFly(kk);
    pp = pp+trialFly{kk,2};
end

pp = 1;
indexTrial = zeros(nTrial,1);
for kk=1:nFly
   indexTrial(pp:pp+trialFly{kk,2}-1) = 1:trialFly{kk,2};
   pp = pp+trialFly{kk,2};
end
%---------------------------------------------------------------------------------------------------------------------------------
%% Load Trial Data %%
clear DATA_ALL
clc
Fs = 1000;  % sampling rate
[b, a] = butter(2, 10/(Fs/2),'low'); % filter design
SI = 1; EI = 9900;  % constrain data to 0-10 seconds
S = 100;

% Initialize data cells
for ii = 1:nFly % # of flys
    for jj = 1:nDir % # of directions
        for kk = 1:nSpeed % # of speeds
            DATA_ALL{ii,jj}{kk,1} = [];
        end
    end
end
DATA_FLY = cell(nFly,1);
% Transfer trial data into cell array [DATA_ALL]
excluded_trials = 0; pp = 1; % keep track of rejected trials due to low frequency
nIdx_2 = 1;
for jj = 1:nTrial
    
    fname = [dirpath char(Files(jj))]; % create path for data files    
    load(fname); % load data files into workspace
    t_p = t_p'; WBAdata = WBAdata'; % put data in column vectors

    wba = WBAdata(:,4) - WBAdata(:, 5); % calculate wingbeat amplitude (L-R)
    Lwing = filtfilt(b, a, WBAdata(:,4)); % left amplitude
    Rwing = filtfilt(b, a, WBAdata(:,5)); % right amplitude
    wbaFilt = filtfilt(b, a, wba); % filter wingbeat amplitude (L-R)
    wbf = WBAdata(:,6); % frequency data
    
    if (mean(wbf) < 1.00) % reject trial if fly stopped
        wbaFilt = wbaFilt.*nan;
        excluded_trials = excluded_trials + 1;
        badTrial{pp,1} = ['Fly_' num2str(Fly(jj)) '_Trial_' num2str(Trial(jj))];
        pp = pp + 1;
    end

    wbaTrial = wbaFilt - median(wbaFilt(1:S)); % normalize by pre-stim WBA
    
    Vind = abs(Vel(jj)/8); % index from velocity
        
    DATA_ALL{indexFly(jj), Dir(jj)}{Vind,1}(end+1, 1:length(wbaTrial(SI:EI))) = wbaTrial(SI:EI)-wbaTrial(SI);
    
    DATA_FLY{indexFly(jj),1}(end+1, 1:length(wbaTrial(SI:EI))) = wbaTrial(SI:EI)-wbaTrial(SI); % all trials grouped by fly #
    
end

tt = linspace(0,10,EI);

fprintf('Bad Trials: %i \n',excluded_trials)
%---------------------------------------------------------------------------------------------------------------------------------
%% Induviudal fly mean & STD calculations %%
clc
% Calculate mean for each fly
clear flyMean
flyMean = cell(nFly,nSpeed);
for kk = 1:nFly
    for jj = 1:nDir % each fly
        flyMean{kk,jj}   = cellfun(@(x) nanmean(x,1), DATA_ALL{kk,jj}, 'UniformOutput', false);
    end
end

DATA_FLY_mean = cellfun(@(x) nanmean(x,2), DATA_FLY, 'UniformOutput', false);
DATA_MEAN_ALL = cell2mat(cellfun(@(x) nanmean(x,1), DATA_FLY_mean, 'UniformOutput', false));

%%
lim = SI:EI;
for mm = 1:nSpeed % each speed
    for jj = 1:nFly % each fly
        speed_ALL{mm,1}(jj,:) = flyMean{jj,1}{mm}(lim) + flyMean{jj,2}{mm}(lim);
        speed_ALL_CW{mm,1}(jj,:) = flyMean{jj,1}{mm}(lim);
        speed_ALL_CCW{mm,1}(jj,:) = flyMean{jj,2}{mm}(lim);
    end
end

% Mean for each speed across all individual flies
if (1==nFly) % For one fly only
    speed_mean = speed_ALL;
    speed_mean_CW = speed_ALL_CW;
    speed_mean_CCW = speed_ALL_CCW;
else
    speed_mean = cellfun(@nanmean, speed_ALL, 'UniformOutput', false);
    speed_std = cellfun(@nanstd, speed_ALL, 'UniformOutput', false);

    speed_mean_CW = cellfun(@nanmean, speed_ALL_CW, 'UniformOutput', false);
    speed_std_CW = cellfun(@nanstd, speed_ALL_CW, 'UniformOutput', false);

    speed_mean_CCW = cellfun(@nanmean, speed_ALL_CCW, 'UniformOutput', false);
    speed_std_CCW = cellfun(@nanstd, speed_ALL_CCW, 'UniformOutput', false);
end
%---------------------------------------------------------------------------------------------------------------------------------
%% T-TEST
clear TTEST meanTrial_ALL meanFly_ALL

TTEST = zeros(nFly,nDir);
for kk = 1:nFly
    [Hp,P] = ttest(DATA_FLY_mean{kk,1}, 0, 'tail', 'left','Alpha',0.01);
    TTEST(kk,1) = Hp;

    [Hp,P] = ttest(DATA_FLY_mean{kk,1}, 0, 'tail', 'right','Alpha',0.01);
    TTEST(kk,2) = Hp;
end

valCCW = sum(TTEST(:,1)==1);
valCW = sum(TTEST(:,2)==1);
valNone = nFly - valCW - valCCW;

% figure (10) ; clf ; hold on
%     histogram(DATA_MEAN_ALL,15,'BinWidth',0.25)
%     ylim([0 10])
%     xlim([-3 3])

% figure (11) ; clf
%     bar([valCW , valCCW , valNone]);
%     ylim([0 20])
%     xticks([1 2 3])
%     xticklabels({'CW','CCW','None'})

%% FIGURES %%
figure(2);
speedBank = {'30','60','90','120','150'};
for jj = 1:nSpeed    
    subplot(1, nSpeed, jj ) ; hold on
        title(speedBank{jj})
        plot(tt, speed_mean_CW{jj}','b','LineWidth',2)
        plot(tt, speed_mean_CCW{jj}','r','LineWidth',2)
        plot(tt, speed_mean{jj}','k','LineWidth',4) 
        plot([0 10],[0 0],'--g','linewidth',2)
        axis([0 10 -5 5])
        if (3==jj); xlabel('Time (s)'); end
end

%%
figure (3) ; hold on
speed_array = (8:8:40).*3.75;
limX = 10;
limY = 5;
rows = 3;
rN = 3;
for jj = 1:nSpeed
    subplot(rows, nSpeed + 2, jj + (rN-1)*(nSpeed+2))
        cla
        if (jj==1); ylabel('60') ; end
        if (rN==1) ; title(speed_array(jj)); end

        % CW
        uE = speed_mean_CW{jj} + speed_std_CW{jj}./sqrt(nFly);
        lE = speed_mean_CW{jj} - speed_std_CW{jj}./sqrt(nFly);
        yP=[lE,fliplr(uE)];
        xP=[tt,fliplr(tt)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        H.patch=patch(xP,yP,1,'facecolor',[0.5 0.5 0.5],...
            'edgecolor','none');
        hold on;
        plot(tt, speed_mean_CW{jj}','b','LineWidth',2)
        box on
        % CCW
        uE = speed_mean_CCW{jj} + speed_std_CCW{jj}./sqrt(nFly);
        lE = speed_mean_CCW{jj} - speed_std_CCW{jj}./sqrt(nFly);
        yP=[lE,fliplr(uE)];
        xP=[tt,fliplr(tt)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        H.patch=patch(xP,yP,1,'facecolor',[0.5 0.5 0.5],...
            'edgecolor','none');
        hold on;
        plot(tt, speed_mean_CCW{jj}','r','LineWidth',2)   

        % SUM
        uE = speed_mean{jj} + speed_std{jj}./sqrt(nFly);
        lE = speed_mean{jj} - speed_std{jj}./sqrt(nFly);
        yP=[lE,fliplr(uE)];
        xP=[tt,fliplr(tt)];
        xP(isnan(yP))=[];
        yP(isnan(yP))=[];
        H.patch=patch(xP,yP,1,'facecolor',[0.5 0.5 0.5],...
            'edgecolor','none');
        hold on;
        plot(tt, speed_mean{jj}','k','LineWidth',2)

        xlim([0 limX]);
        ylim([-limY limY]);
        line([0 limX],[0 0],'Color','g');
end

%%
hold off
figure (3)
    subplot(rows, nSpeed + 2, (nSpeed+1) + (rN-1)*(nSpeed+2)) ; cla ; hold off
        histogram(DATA_MEAN_ALL,15,'BinWidth',0.2)
        ylim([0 5])
        xlim([-3 3])

    subplot(rows, nSpeed + 2, (nSpeed+2) + (rN-1)*(nSpeed+2)) ; cla ; hold off
        bar([valCW , valCCW , valNone]);
        ylim([0 20])
        xticks([1 2 3])
        xticklabels({'CW','CCW','None'})

%%