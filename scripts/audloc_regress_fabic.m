%% Behavioural analysis: auditory localiser

clear;
close all;
clc;

% Suppress figures
set(0,'DefaultFigureVisible','on')

% matfiles for group analysis
aud_resp = [];
aud_rt = [];
res=struct;

data_path = [pwd, '\data\experiment2\unis_audio\'];

% Subjects
subjID={'70';'71';'72';'73';'74';'75';'76';'77';'78';'79';'80';'81';'82';'83';'84';'85';'86';'87';...
    '88';'89';'90';'91';'93';'94';'95';'96';'97';'98';'99';'100';'101';'102';'103';'104';'106';'107';...
    '111';'112';'113';'114';'115'};

for iSubj = 1:length(subjID)
    
    % Select the data folder of the subject
    dataPath_audloc = [data_path, 'Log_audio_fab', subjID{iSubj}, '.mat'];

    load(dataPath_audloc)
    aud_tdata = LogFile;
    
    %% Mean responses 
    aud_subj_resp = grpstats(aud_tdata(aud_tdata.MissedResp==0,:), ...
        'Aud_Loc', {'mean', 'sem', 'std', 'numel'}, 'DataVars', 'Response');
    aud_subj_rt = grpstats(aud_tdata(aud_tdata.MissedResp==0,:), ...
        'Aud_Loc', {'mean', 'sem', 'std', 'numel'}, 'DataVars', 'RT');
    
    % Save current subjects data
    aud_resp = cat(1, aud_resp, [aud_subj_resp.mean_Response(1), aud_subj_resp.mean_Response(2), aud_subj_resp.mean_Response(3)]);
    aud_rt = cat(1, aud_rt, [aud_subj_rt.mean_RT(1), aud_subj_rt.mean_RT(2), aud_subj_rt.mean_RT(3)]);
    
    clear aud_subj_resp aud_subj_rt
    
    %% Regression
    X_temp = aud_tdata.Aud_Loc(aud_tdata.MissedResp==0,:);
    res(iSubj).X=X_temp;
    X(:,2)=X_temp;
    X(:,1)=ones(length(X_temp),1);   
    Y = aud_tdata.Response(aud_tdata.MissedResp==0,:);
    res(iSubj).Y=Y;
    % perform stats
    % b(1)=intercept, b(2)=slope
    % bint = confidence intervals for each coefficient estimate
    % regstats = R^2 statistic, F-statistic and its p-value, estimate of error variance
    [res(iSubj).b,res(iSubj).bint,~,~,res(iSubj).regstats] = regress(Y,X);
    res(iSubj).RMSE = sqrt(mean((Y-X(:,2)).^2));
    [res(iSubj).corr_rho,res(iSubj).corr_pval] = corr(X(:,2),Y,'Tail','right');
    coef_group(iSubj,:)=res(iSubj).b';
    RMSE_group(iSubj,:)=res(iSubj).RMSE;
    rho_group(iSubj,:)=res(iSubj).corr_rho;
    clear X Y X_temp
    
    disp(' ');
    disp([subjID{iSubj} ' completed']);
    disp(' ');
  
end
% Fisher-z transform rho correlation and compute across-participants' mean
rho_tot=atanh(rho_group);
rho_tot_mean_ftran=mean(rho_tot);
% transform back
rho_tot_mean=tanh(mean(rho_tot));

% effect size (z) and 95% CI
sigmaz=std(rho_tot)/sqrt(size(aud_resp,1));
l95=tanh(mean(rho_tot)-1.96*sigmaz);
u95=tanh(mean(rho_tot)+1.96*sigmaz);

%[h,p,ci,stats] = ttest(rho_tot);
%[h,p,ci,stats] = ttest(rho_tot, 0, 'Tail', 'Right')
%p1 = signrank(rho_tot);

regres_param=[coef_group,RMSE_group];

% save results
cd(data_path);
save('audloc_resp_group_attrel', 'aud_resp', 'aud_rt', 'res', 'regres_param');

aud_group_resp_mean = mean(aud_resp);
aud_group_resp_std = std(aud_resp);
aud_group_resp_sem = std(aud_resp)/sqrt(size(aud_resp,1));

aud_group_rt_mean = mean(aud_rt);
aud_group_rt_std = std(aud_rt);
aud_group_rt_sem = std(aud_rt)/sqrt(size(aud_rt,1));

%% Figure mean responses and RT

ni = 1; % number of figure rows (single participant)
nj = 2; % number of figure columns
loc = [-9 0 9];
% size
positionXY = [0, 0, 400, 150];
% colors
cols  = [0 0 0.8];
figure('color', [1 1 1], 'Position', positionXY);

% subplot loop (columns)
for j = 1:nj    
    % determine plot id and position
    subplot(ni,nj,j);    
    switch j
        case 1 % response       
            %plot([-9 0 9], [aud_resp(:,1), aud_resp(:,2), aud_resp(:,3)], 'MarkerEdgeColor', cols, 'Color', cols, 'LineWidth', 1); hold on
            %plot([-9 0 9], [aud_resp(:,1), aud_resp(:,2), aud_resp(:,3)], 'o', 'MarkerEdgeColor', cols, 'Color', cols, 'LineWidth', 1); hold on            
            plot(loc, aud_group_resp_mean, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', cols, 'MarkerFaceColor', cols); hold on
            %plot(loc, aud_group_resp_mean, 'MarkerSize', 7, 'MarkerEdgeColor', cols, 'MarkerFaceColor', cols, 'LineWidth', 3); hold on
            errorbar(loc, aud_group_resp_mean, aud_group_resp_sem, aud_group_resp_sem, 'Color', cols, 'LineWidth', 3);            
            % plot settings
            % give some space around extreme locations
            xl = [-loc(end)-6 loc(end)+6]; xlim(xl);
            yl = [min(loc)-6, max(loc)+6]; ylim(yl);            
            % adjust x and y ticks
            ticks = [-9 0 9];            
            set(gca, 'YTick', ticks);
            set(gca, 'XTick', ticks);
            set(gca,'FontName', 'Helvetica');
            set(gca,'FontSize', 15);            
        case 2 % RTs
            %plot([-9 0 9], [aud_rt(:,1), aud_rt(:,2), aud_rt(:,3)], 'o', 'MarkerEdgeColor', cols, 'Color', cols, 'LineWidth', 1); hold on            
            plot(loc, aud_group_rt_mean, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', cols, 'MarkerFaceColor', cols); hold on
            %plot(loc, aud_group_rt_mean, 'MarkerSize', 7, 'MarkerEdgeColor', cols, 'MarkerFaceColor', cols, 'LineWidth', 3); hold on
            errorbar(loc, aud_group_rt_mean, aud_group_rt_sem, aud_group_rt_sem, 'Color', cols, 'LineWidth', 3);            
            % plot settings
            % give some space around extreme locations
            xl = [-loc(end)-6 loc(end)+6]; xlim(xl);
            yl = [0 1800]; ylim(yl);            
            % adjust x and y ticks
            ticksX = [-9 0 9];
            ticksY = (200:400:1600);            
            set(gca, 'YTick', ticksY);
            set(gca, 'XTick', ticksX);
            set(gca, 'YTickLabel', (200:400:1600));
            set(gca,'FontName', 'Helvetica');
            set(gca,'FontSize', 15);            
    end
end

saveas(gcf, 'group_audloc_resp_av09.png');

%% Figure regression

coef_group_mean=mean(coef_group(:,:),1);
coef_group_sem=std(coef_group(:,:),0,1)/sqrt(length(subjID));

RMSE_group_mean=mean(RMSE_group(:,:),1);
RMSE_group_sem=std(RMSE_group(:,:),0,1)/sqrt(length(subjID));

cols = [0.6 0.6 0.6];

%% regression line group

positionXY = [0, 0, 170, 170];
figure('Position', positionXY);

% guideline
line([-16 16],[-16 16],'Color', [0.8 0.8 0.8], 'LineStyle', '-', 'LineWidth', 1);hold on;

xvect=[];
yvect=[];
for isubj=1:length(subjID)
    xvect=[xvect;res(isubj).X];
    yvect=[yvect;res(isubj).Y];
end
jit_val=6*rand(length(yvect),1);
yvect_plot=yvect+jit_val;
scatter(xvect,yvect_plot-mean(yvect_plot),10,'MarkerEdgeColor',cols,'LineWidth',0.01);

regressline=refline(coef_group_mean(2),coef_group_mean(1));hold on;
set(regressline,'Color',cols,'LineWidth',1);hold on;

% plot settings
% give some space around extreme locations
xl = [-16 16]; xlim(xl);
yl = [-16 16]; ylim(yl);

set(gca,'XTick',[-9 0 9],...
    'YTick',[-9 0 9],...
    'FontName', 'Helvetica', ...
    'TickLength', [0.02 0.02]);%, ...
%'TickDir','both');
saveas(gcf,'group_audloc_regres_av09.png');

%% regression line representative subj


for ind_subj=1:length(subjID)

positionXY = [0, 0, 170, 170];
figure('Position', positionXY);

% guideline
line([-16 16],[-16 16],'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);hold on;

xvect=[];
yvect=[];

xvect=[xvect;res(ind_subj).X];
yvect=[yvect;res(ind_subj).Y];

jit_val=2*rand(length(yvect),1);
jit_val_x=2*rand(length(yvect),1);
yvect_plot=yvect+jit_val;
scatter(xvect+jit_val_x,yvect_plot-mean(yvect_plot),10,'MarkerEdgeColor',cols,'LineWidth',0.01);

regressline=refline(coef_group(ind_subj,2),coef_group(ind_subj,1));hold on;
set(regressline,'Color',cols,'LineWidth',1.5);hold on;

% plot settings
% give some space around extreme locations
xl = [-16 16]; xlim(xl);
yl = [-16 16]; ylim(yl);

set(gca,'XTick',[-9 0 9],...
    'YTick',[-9 0 9],...
    'FontName', 'Helvetica', ...
    'TickLength', [0.02 0.02]);%, ...
%'TickDir','both');
figname=['group_audloc_regres_rep_subj',num2str(ind_subj),'.png'];
saveas(gcf,figname);
end
