%% Experiment 1: model comparison
% 2 (BCI vs FF model) x 2 (communicative vs non-communicative action) model space

% Note: the analysis presupposes that SPM12 is in the MATLAB path
% Download SPM12 from: https://www.fil.ion.ucl.ac.uk/spm/software/spm12

clear;
close all;
clc;

% Path to SPM12
addpath(genpath('<your path>'));
% Data path
dataPath='<your path>';
% Path to save data
savePath=fullfile(dataPath,'results','results_exp1');

% Subjects
subjID = {'05';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';...
    '38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';...
    '51';'52';'53';'54';'55';'56';'57';'58';'59'};
nsubj=length(subjID);

%% Load data

% Initialize (n subjects; 4 models)
bic_group=nan(nsubj,4);
param_group=cell(nsubj,1,4);

for iSubj = 1:nsubj
    % model 1: BCI no com action
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_bci_simulations_avert_nosac_pool_best10']),'bic','fm_bestParameters');
    bic_group(iSubj,1)=bic;
    param_group(iSubj,:,1)=fm_bestParameters;
    % model 2: BCI com action (1=com; 2=non-com)
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_bci_simulations_avert_nosac_sep_best10']),'bic','fm_bestParameters');
    bic_group(iSubj,2)=bic;
    param_group(iSubj,:,2)={[fm_bestParameters{1},fm_bestParameters{2}]};
    % model 3: FF no com action
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_fus_simulations_avert_nosac_pool_best10']),'bic','fm_bestParameters');
    bic_group(iSubj,3)=bic;
    param_group(iSubj,:,3)=fm_bestParameters;
    % model 4: FF com action (1=com; 2=non-com)
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_fus_simulations_avert_nosac_sep_best10']),'bic','fm_bestParameters');
    bic_group(iSubj,4)=bic;
    param_group(iSubj,:,4)={[fm_bestParameters{1},fm_bestParameters{2}]};
end

%% Random-effect analysis

% Run Bayesian model selection
[rfx.alpha,rfx.exp_r,rfx.xp,rfx.pxp,rfx.bor] = spm_BMS(bic_group);

% Plot
pxp_toplot = [rfx.pxp(1) rfx.pxp(2); rfx.pxp(3) rfx.pxp(4)];
map_vector=0:0.01:1;
map=repmat(map_vector',1,3);
positionXY = [0, 0, 350, 300];
figure('Position',positionXY);
imagesc(pxp_toplot);
colormap(map);
colorbar('LineWidth',1.2,'TickLength',0.02);
set(gca,'XTick',[1 2 3],'YTick',[1 2 3],...
    'FontName','Arial','FontSize', 12,...
    'TickLength',[0 0],'LineWidth', 1.2,...
    'YTickLabel',[],'XTickLabel',[]);
saveas(gcf, fullfile(savePath,'pxp_4models_exp1'), 'svg');

%% Model parameters

% Group statistics (mean and SEM)
param_m1=cell2mat(param_group(:,:,1));
param_m2=cell2mat(param_group(:,:,2));
param_m3=cell2mat(param_group(:,:,3));
param_m4=cell2mat(param_group(:,:,4));
% model 1: BCI no com action
param_group_m1(1,:)=mean(param_m1);
param_group_m1(2,:)=std(param_m1)/sqrt(nsubj);
% model 2: BCI com action (params [1-4]: com; params [5-8]: non-com)
param_group_m2(1,:)=mean(param_m2);
param_group_m2(2,:)=std(param_m2)/sqrt(nsubj);
% model 3: FF no com action
param_group_m3(1,:)=mean(param_m3);
param_group_m3(2,:)=std(param_m3)/sqrt(nsubj);
% model 4: FF com action (params [1-4]: com; params [5-8]: non-com)
param_group_m4(1,:)=mean(param_m4);
param_group_m4(2,:)=std(param_m4)/sqrt(nsubj);

% Prepare data for JASP analysis
conditionNames = {'pCcom','sigPcom','sigAcom','sigVcom','pCncom','sigPncom','sigAncom','sigVncom'};
T = array2table(param_m2, 'VariableNames', conditionNames);
writetable(T,fullfile(savePath,'params_winning_model_m2_avert_nosac.csv'));

% Compute difference Com - NCom for pcommon, sigmaA and sigmaV
pcommon_diff_Com_NCom_exp1=param_m2(:,1)-param_m2(:,5);
sigmaA_diff_Com_NCom_exp1=param_m2(:,3)-param_m2(:,7);
sigmaV_diff_Com_NCom_exp1=param_m2(:,4)-param_m2(:,8);
save(fullfile(savePath,"params_diff_Com_NCom_exp1"),"pcommon_diff_Com_NCom_exp1",...
    "sigmaA_diff_Com_NCom_exp1","sigmaV_diff_Com_NCom_exp1");

%% Plot group results

positionXY=[0, 0, 850, 300];
figure('Position',positionXY);
% colors
cols.com=[217 217 217]/255;
cols.nocom=[89 89 89]/255;
Cols=[cols.com;cols.nocom];
% nr columns
row = 1; % number of figure rows
col = 4; % number of figure columns

for c = 1:col
    hsp = subplot(1,col,c);
    switch c
        case 1
            plot_matrix_mean=[param_group_m2(1,1) param_group_m2(1,5)];
            plot_matrix_sem=[param_group_m2(2,1) param_group_m2(2,5)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.7,...
                    'FaceColor',Cols(i,:),...
                    'LineWidth',1.2)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.2)
            end
        case 2
            plot_matrix_mean=[param_group_m2(1,2) param_group_m2(1,6)];
            plot_matrix_sem=[param_group_m2(2,2) param_group_m2(2,6)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.7,...
                    'FaceColor',Cols(i,:),...
                    'LineWidth',1.2)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.2)
            end
        case 3
            plot_matrix_mean=[param_group_m2(1,3) param_group_m2(1,7)];
            plot_matrix_sem=[param_group_m2(2,3) param_group_m2(2,7)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.7,...
                    'FaceColor',Cols(i,:),...
                    'LineWidth',1.2)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.2)
            end
        case 4
            plot_matrix_mean=[param_group_m2(1,4) param_group_m2(1,8)];
            plot_matrix_sem=[param_group_m2(2,4) param_group_m2(2,8)];
            for i=1:2
                hold on
                bar(i,plot_matrix_mean(i),0.7,...
                    'FaceColor',Cols(i,:),...
                    'LineWidth',1.2)
                line([i i],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
                    plot_matrix_mean(i)-plot_matrix_sem(i)],...
                    'Color',[0 0 0],'LineWidth',1.2)
            end
    end    
    % set x and y axes and ticks
    xl = [0.3 2.7]; xlim(xl);
    if c == 1
        yl = [0 0.5]; ylim(yl);
        ticksY = (0:0.1:0.5);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', ([0 0.1 0.2 0.3 0.4 1]));
    elseif c == 2
        yl = [0 25]; ylim(yl);
        ticksY = (0:5:25);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:5:25));
    elseif c == 3
        yl = [0 10]; ylim(yl);
        ticksY = (0:2:10);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:2:10));
    elseif c == 4
        yl = [0 5]; ylim(yl);
        ticksY = (0:1:5);
        set(gca, 'YTick', ticksY);
        set(gca, 'YTickLabel', (0:1:5));
    end
    set(gca,'FontName', 'Arial');
    set(gca,'FontSize', 12);
    set(gca, 'XTickLabel', '');
    set(gca, 'XTick', []);
    set(gca,'TickLength', [0.02 0.02]);
    set(gca,'LineWidth',1.2)
end
saveas(gcf, fullfile(savePath,'params_winning_model_m2_exp1'), 'svg');