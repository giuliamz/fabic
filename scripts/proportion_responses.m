%% Compute proportions of participants' responses and BCI predictions

close all;
clear;
clc;

%% general settings

% 9 AV locations in factorial design
ni = 3; % A loc
nj = 3; % V loc
nCond = 4; % 2 Action x 2 Resp
exp='exp1'; % exp1; exp2
mainPath = fullfile('<your path>',exp); % experiment folder
dataPath = fullfile(mainPath,'results'); % data folder
savePath = fullfile(dataPath,['results_' exp]); % folder to save data
% subject list
if strcmp(exp,'exp1')
    subjID={'05';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';...
    '38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';...
    '51';'52';'53';'54';'55';'56';'57';'58';'59'};
elseif strcmp(exp,'exp2')
    subjID = {'70';'71';'72';'73';'75';'76';'77';'79';'80';'81';...
    '83';'84';'85';'86';'87';'88';'89';'93';'94';'95';'98';'99';'101';...
    '102';'103';'104';'106';'107';'111';'112';'113';'114';'115'};
end
nSubj=length(subjID);
% initialize results matrix
hist_mat=nan(9,3,nCond,nSubj);

for isubj=1:nSubj
    
    % Load current subject dataset
    if strcmp(exp,'exp1')
        file = (fullfile(dataPath,['s' subjID{isubj} '_bci_simulations_avert_nosac_sep_best10']));
    elseif strcmp(exp,'exp2')
        file = (fullfile(dataPath,['s' subjID{isubj} '_bci_simulations_control_nosac_sep_best10']));
    end
    load(file, 'bciSimulations');
    
    data_behA=[];
    data_behV=[];
    data_predA=[];
    data_predV=[];
    
    for icond=1:18
        data_behA=[data_behA;bciSimulations(icond).freq_dataA];
        data_behV=[data_behV;bciSimulations(icond).freq_dataV];
        data_predA=[data_predA;bciSimulations(icond).freq_predA];
        data_predV=[data_predV;bciSimulations(icond).freq_predV];
    end
    hist_mat_bci(:,:,1,isubj)=data_predA(1:9,:);
    hist_mat_bci(:,:,2,isubj)=data_predA(10:18,:);
    hist_mat_bci(:,:,3,isubj)=data_predV(1:9,:);
    hist_mat_bci(:,:,4,isubj)=data_predV(10:18,:);
    
    hist_mat_beh(:,:,1,isubj)=data_behA(1:9,:);
    hist_mat_beh(:,:,2,isubj)=data_behA(10:18,:);
    hist_mat_beh(:,:,3,isubj)=data_behV(1:9,:);
    hist_mat_beh(:,:,4,isubj)=data_behV(10:18,:);
    
end

%% group means and std/sem

hist_mat_group_mean = mean(hist_mat_beh,4);
hist_mat_group_std = std(hist_mat_beh,0,4);
hist_mat_group_sem = hist_mat_group_std/sqrt(nSubj);

hist_mat_bci_group_mean = mean(hist_mat_bci,4);
hist_mat_bci_group_std = std(hist_mat_bci,0,4);
hist_mat_bci_group_sem = hist_mat_bci_group_std/sqrt(nSubj);

cols.Aud1  = [230 97 1]/255;
cols.Aud2  = [253 184 99]/255;
cols.Vis1 = [94 60 153]/255;
cols.Vis2 = [178 171 210]/255;

%% plot group distributions

positionXY = [0, 0, 500, 500];
loc = [1 2 3];

% auditory report
figure('color', [1 1 1], 'Position', positionXY);
% subplot loop (rows)
for i = 1:ni
    % subplot loop (columns)
    for j = 1:nj
        % determine plot id and position
        curplotid = (i-1)*nj+j;
        % determine plot id and position
        hsp = subplot(ni,nj,curplotid);
        plot(loc, hist_mat_group_mean(curplotid,:,1),'Color',cols.Aud1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,1),'--','Color',cols.Aud1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_group_mean(curplotid,:,2),'Color',cols.Aud2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,2),'--','Color',cols.Aud2,'LineWidth',1.2); hold on        
        % plot settings
        xl = [0.8 3.2]; xlim(xl);
        yl = [-0.05 1.05]; ylim(yl);
        set(gca,'FontName', 'Arial');
        set(gca,'FontSize', 12);
        ticksX = [];
        ticksY = [];
        set(gca, 'YTick', ticksY);
        set(gca, 'XTick', ticksX);
        h = gca; h.YAxis.Visible = 'off';
        set(gca,'LineWidth',1.2);
        box off
    end
end
saveas(gcf, fullfile(savePath,['prop_resp_BCI_Aud_' exp]), 'svg');

% visual report
figure('color', [1 1 1], 'Position', positionXY);
% subplot loop (rows)
for i = 1:ni
    % subplot loop (columns)
    for j = 1:nj
        % determine plot id and position
        curplotid = (i-1)*nj+j;
        % determine plot id and position
        hsp = subplot(ni,nj,curplotid);
        plot(loc, hist_mat_group_mean(curplotid,:,3),'Color',cols.Vis2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,3),'--','Color',cols.Vis2,'LineWidth',1.2); hold on
        plot(loc, hist_mat_group_mean(curplotid,:,4),'Color',cols.Vis1,'LineWidth',1.2); hold on
        plot(loc, hist_mat_bci_group_mean(curplotid,:,4),'--','Color',cols.Vis1,'LineWidth',1.2); hold on        
        % plot settings
        xl = [0.8 3.2]; xlim(xl);
        yl = [-0.05 1.05]; ylim(yl);
        set(gca,'FontName', 'Arial');
        set(gca,'FontSize', 12);
        ticksX = [];
        ticksY = [];
        set(gca, 'YTick', ticksY);
        set(gca, 'XTick', ticksX);
        h = gca; h.YAxis.Visible = 'off';
        set(gca,'LineWidth',1.2);
        box off
    end
end
saveas(gcf, fullfile(savePath,['prop_resp_BCI_Vis_' exp]), 'svg');
