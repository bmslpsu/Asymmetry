function [] = Analyze_Asymmetry_Verification_ALL()
%---------------------------------------------------------------------------------------------------------------------------------
% Analyze_Asymmetry_Verification_ALL:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root = [];
root{end+1}	= 'S:\Restricted\BC\EXPERIMENTS\Experiment_Asymmetry';
root{end+1} = 'S:\Restricted\BC\EXPERIMENTS\Experiment_Asymmetry_Verification';
root{end+1} = 'S:\Restricted\BC\EXPERIMENTS\Experiment_Asymmetry_Arena2';
root{end+1} = 'S:\Restricted\BC\EXPERIMENTS\Experiment_Asymmetry_George';
root{end+1} = 'S:\Restricted\BC\EXPERIMENTS\Experiment_Asymmetry_JMM_v2';
n.root = length(root);

% Select Files
FILES   = cell(n.root,1);
PATH    = FILES;
D       = FILES;
I       = FILES;
N       = FILES;
T       = FILES;
U       = FILES;
ALL     = [];
for ww = 1:n.root 
    [D{ww},I{ww},N{ww},U{ww},T{ww},FILES{ww},PATH{ww}] = GetFileData(root{ww});
    normI = I{ww};
    if ww==1
        normI{:,1} = normI{:,1};
    else
      normI{:,1} = normI{:,1} + ALL{end,1};
    end
    
    ALL = [ALL ; normI];
end
n.fly = max(ALL{:,1});
n.file = size(ALL,1);
n.cond = N{1}{1,3};

clear normI ww kk
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
WING.Set  	= cell(n.root,1);
WING.WBA    = cell(n.fly,N{1}{1,3});
WING.Time   = linspace(0,10,10*1000)'; % time vector
pp = 0;
for ww = 1:n.root % data set
    WING.Set{ww} = cell(N{ww}{1,1},N{ww}{1,3});
    for kk = 1:N{ww}{1,end} % files per data set
        pp = pp + 1;
        % Load head & DAQ data %
        load([PATH{ww} FILES{ww}{kk}],'data','t_p'); % load WBA
        if size(data,2)>size(data,1)
            data = data';
        end
        %-------------------------------------------------------------------------------------------------------------------------
        % Check wing beat frequency %
        wing.F = data(:,6);
        if min(wing.F)<1.40
            disp(['Low WBF:  Set ' num2str(ww) '/ Fly ' num2str(D{ww}{kk,1}) ' Trial ' num2str(D{ww}{kk,2})])
%             continue
        end
        %-------------------------------------------------------------------------------------------------------------------------
        % Get wing data from DAQ %
        Fc = 20; % cutoff frequency [Hz]
        wing.time = t_p; % wing time [s]
        Fs = 1/mean(diff(t_p)); % sampling frequency [Hz]
        [b,a] = butter(2,Fc/(Fs/2)); % butterworth filter 
        wing.L = filtfilt(b,a,data(:,4))'; % filter left wing
        wing.R = filtfilt(b,a,data(:,5))'; % filter right wing
        wing.wba = filtfilt(b,a,wing.L - wing.R); % filter L-R
        wing.wba = wing.wba - wing.wba(1); % normalize
        wing.wba = interp1(wing.time, wing.wba , WING.Time, 'linear');
        %-------------------------------------------------------------------------------------------------------------------------
        % Store data in cells %
        WING.WBA 	{ALL{pp,1},ALL{kk,3}}(:,end+1) = wing.wba;
        WING.Set 	{ww}{I{ww}{kk,1},I{ww}{kk,3}}(:,end+1) = wing.wba;
    end
end
disp('DONE')
% clear kk ww data t_p wing Fc Fs wing a b

n.trial = size(WING.WBA,1);

% Means
WING.FlyWBA = cellfun(@(x) median(x,2), WING.WBA, 'UniformOutput',false);
for kk = 1:n.fly
    WING.FlyWBA{kk,n.cond+1} = WING.FlyWBA{kk,1} + WING.FlyWBA{kk,2};
end

WING.FlyAllWBA = cell(1,n.cond+1);
for ii = 1:size(WING.FlyWBA,2)
    WING.FlyAllWBA{ii} = cat(2,WING.FlyWBA{:,ii});
end
WING.GrandWBA = cellfun(@(x) mean(x,2), WING.FlyAllWBA, 'UniformOutput',false);
WING.GrandSTD = cellfun(@(x) std(x,[],2), WING.FlyAllWBA, 'UniformOutput',false);
WING.GrandWBA = cat(2,WING.GrandWBA{:});
WING.GrandSTD = cat(2,WING.GrandSTD{:});

%% Grand Mean Figure %%
%---------------------------------------------------------------------------------------------------------------------------------
Fig1 = figure (1) ; clf
Fig1.Name = 'Mean Delta WBA';
Fig1.Color = 'w';
Fig1.Position = [100 100 800 500];
movegui(Fig1,'center')

ax1 = gca; hold on
ax1.FontSize = 12;
ax1.Title.String = 'Mean Delta WBA';
ax1.Title.Color = 'k';
ax1.Title.FontSize = 16;
ax1.YLabel.String = 'V';
ax1.YLabel.FontSize = 14;
ax1.YLim = 4*[-1 1];
ax1.XLabel.String = 'Time (s)';
ax1.XLabel.FontSize = ax1.YLabel.FontSize;
ax1.XLabel.Color = 'k';
ax1.XLim = [0 10];

cList = {'r','b','k'}; % {CCW,CW,Sum}

% % Trial
% for kk = 1:n.trial
%     for jj = 1:n.cond
%         for ii = 1:size(WING.WBA{kk,jj},2)
%             h.trial = plot(WING.Time,WING.WBA{kk,jj}(:,ii),'Color',cList{jj},'LineWidth',0.5);
%             h.trial.Color(4) = 0.1;
%         end
%     end
% end

TrialCount = [];
% Fly
for jj = 1:n.cond+1-1
    for kk = 1:size(WING.FlyAllWBA{jj},2)
%         h.fly = plot(WING.Time,WING.FlyAllWBA{jj}(:,kk),'Color',cList{jj},'LineWidth',2);
%         h.fly.Color(4) = 0.2;
        TrialCount = [TrialCount; kk jj 1 size(WING.WBA{kk,jj},2)];
        if TrialCount(end,4)<3
            disp('Here')
        end
    end
end

% Grand
for jj = 1:n.cond+1
    PlotPatch(WING.GrandWBA(:,jj),WING.GrandSTD(:,jj),WING.Time,2,n.fly,cList{jj},[0.5 0.5 0.5],0.5,8);
end
plot([0 ax1.XLim(2)],[0 0],'-g','LineWidth',1)
% plot(WING.Time,WING.GrandWBA(:,2),'.-c','LineWidth',1)

%% Population Asymmetry Distribution Figure %%
%---------------------------------------------------------------------------------------------------------------------------------
Fig2 = figure (2) ; clf
Fig2.Name = 'Mean Delta WBA';
Fig2.Color = 'w';
Fig2.Position = [100 100 800 500];

ax2 = gca; hold on
ax2.FontSize = 12;
ax2.Title.String = 'Mean Offset Fly Distribution';
ax2.Title.FontSize = 14;
ax2.XLabel.String = 'Offset (V)';
ax2.XLabel.FontSize = 14;
ax2.XLim = 12*[-1 1];
ax2.YLim = [0 0.5];
ax2.YTick = ax2.YLim;

mean_off = cellfun(@(x) mean(x,1), WING.FlyAllWBA, 'UniformOutput', false); % get all means for sum

set(Fig2, 'Renderer', 'painters', 'Units','inch', 'Position', [1,1, 7, 7]);
 
histogram(mean_off{1,3}, 10, 'BinWidth',1, 'Normalization', 'probability')
line([0,0],[0 1], 'Color','m')

[hh, p] = ttest(mean_off{1,3}, 0); % test if individuals are > 0 V
S = skewness(mean_off{1,3});

text(-4.5, 0.4,  ['p = ', num2str(round(p,3))]);
text(-4.5, 0.38, ['s = ', num2str(round(S,3))]);
text(-4.5, 0.36, ['h = ', num2str(hh)]);

movegui(Fig2,'center')

%% Induvidual Asymmetry Distribution Figure %%
%---------------------------------------------------------------------------------------------------------------------------------
% Make stats table for testing individuals
StatT = nan;
M = nan;
mm = 1;
for ii = 1:size(mean_off,1)
    for jj = 1:n.fly
        StatT(mm,1) = jj;                  	% animal
        StatT(mm,2) = 45*ii;            	% speed   
        StatT(mm,3) = 0;
        StatT(mm,4) = mean_off{ii,3}(jj);   % sum

        M(ii,jj) = mean_off{ii,3}(jj);

        mm = mm + 1;
    end
end

tbl = table(StatT(:,1),StatT(:,2),StatT(:,3),StatT(:,4),'VariableNames',{'Animal','Speed','SpatWav','Sum'});

B = varfun(@median,tbl, 'InputVariables','Sum','GroupingVariables','Animal');
B = table2array(B);

ydata = M; 
[r, c] = size(ydata);
xdata = repmat(1:c, r, 1);

% figure; boxplot(StatT(:,4),StatT(:,1))

Fig3 = figure (6); clf
set(Fig3, 'Renderer', 'painters', 'Units','inch', 'Position', [1,1, 12, 3]);

ax3 = gca; hold on
ax3.YLim = 6*[-1 1];
ax3.XLim = [0 n.fly+1];
ax3.XTick = 1:2:n.fly;
% ax.XTickLabels = [];

scatter(xdata(:), ydata(:), 'r', 'filled', 'jitter','on', 'jitterAmount', 0.3);
plot([xdata(1,:)-0.15; xdata(1,:) + 0.15], repmat(median(ydata, 1), 2, 1), 'k-', 'LineWidth',2)
plot([0 n.fly+1],[0 0],'-g')
errorbar(xdata(1,:), median(ydata,1), std(ydata,[],1)./sqrt(size(ydata,1)), 'k.')

% asymmetry index
for jj = 1:size(ydata,2)  
    [~, Pval] = ttest(ydata(:,jj));
    if Pval <= 0.05
        scatter(xdata(1,jj), 10, '*k')
    end
end

movegui(Fig3,'center')

%% Induvidual Trials Figure %%
%---------------------------------------------------------------------------------------------------------------------------------
% Individual trials per speed
ani = TrialCount(:,1);
speed = TrialCount(:,2);
dir = TrialCount(:,3);
tcount = TrialCount(:,4);

Fig4 = figure (4); clf
set(Fig4, 'Renderer', 'painters', 'Units','inch', 'Position', [1,1, 7, 7]);

boxplot(tcount, {ani speed dir},'labels', {ani speed dir}, 'plotstyle','compact');

ax4 = gca;
ax4.YLim = [0 15];
ax4.YTick = 1:ax4.YLim(end);
ax4.XTickLabels = '';

movegui(Fig4,'center')

%%
figure; scatter(StatT(:,1),StatT(:,3),'r.','jitter','on','jitterAmount',0.05);
hold on;
plot([B(:,1)-0.15; B(:,1) + 0.15], [(B(:,3));(B(:,3))], 'k-')

end