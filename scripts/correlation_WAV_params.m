%% Correlation analysis WAV and parameters winning model
% WAV Auditory responses [Com - NCom]
% Parameters winning model (BCI / Fixed Criterion Separated): pcommon / criterion [Com - NCom], sigmaA [Com - NCom]

clear;
close all;
clc;

%% General settings

% Experiment
exp='exp1'; % exp1; exp2
dataPath = fullfile('<your path>',exp); % data folder

% Subjects
nsubj=34;

% Model
model = 'fixCrit'; % bci fixCrit

% Load data
if strcmp(exp,'averted_nosaccades')
    load(fullfile(dataPath,['params_diff_' model '_Com_NCom_exp1.mat']));
    x1 = pcommon_diff_Com_NCom_exp1;
    x2 = sigmaA_diff_Com_NCom_exp1;
    load(fullfile(dataPath,'wav_aud_diff_Com_NCom_exp1.mat'));
    y = WAV_AUD_COM_NC';
elseif strcmp(exp,'control_nosaccades')
    load(fullfile(dataPath,['params_diff_' model '_Com_NCom_exp2.mat']));
    x1 = pcommon_diff_Com_NCom_exp2;
    x2 = sigmaA_diff_Com_NCom_exp2;
    load(fullfile(dataPath,'wav_aud_diff_Com_NCom_exp2.mat'));
    y = WAV_AUD_COM_NC';
end

%% Model 1: x = pcommon / criterion (Com - NCom)

% Fit linear model using fitlm (includes confidence intervals)
mdl1 = fitlm(x1, y);
rho1 = sqrt(mdl1.Rsquared.Ordinary);
pval1 = mdl1.Coefficients.pValue(2);

if strcmp(model,'fixCrit')
    min_x = -10;
    max_x = 10;
elseif strcmp(model,'bci')
    min_x = -0.3;
    max_x = 0.7;
end

% Sort x for plotting the fitted line and CI in order
xfit1 = linspace(min_x, max_x, nsubj)';
[yfit1, yci1] = predict(mdl1, xfit1);  % yci contains lower and upper CI

% Plot regression line
figure;plot(xfit1, yfit1, 'k-', 'LineWidth', 2.5);hold on;

% Plot confidence intervals
plot(xfit1, yci1(:,1), 'k--', 'LineWidth', 1.5);  % lower bound
plot(xfit1, yci1(:,2), 'k--', 'LineWidth', 1.5);  % upper bound

% Scatter plot
scatter(x1, y, 70, [217 217 217]/255, 'filled',...
    'LineWidth', 1.5, 'MarkerEdgeColor', [100 100 100]/255);

% Axis, labels and legend
box off
set(gca,'FontName', 'Arial');
set(gca,'FontSize', 12);
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2);
xlabel('X Com - NCom');
if strcmp(model,'fixCrit')
    xlim([-12 12]);
    set(gca,'XTick', min_x:2:max_x);
    ylim([-0.5 0.6]);
    set(gca,'YTick', -0.5:0.1:0.6);
elseif strcmp(model,'bci')
    set(gca,'XTick', min_x:0.1:max_x);
    ylim([-0.4 0.6]);
    xlim([-0.4 0.8]);
end
set(gca,'YLabel', text('String', 'W_{AV} (a.u.) Com - NCom'));

% Save figure
if strcmp(model,'fixCrit')
    saveas(gcf, fullfile(dataPath,['corr_WAV_criterion_' exp]), 'svg');
else
    saveas(gcf, fullfile(dataPath,['corr_WAV_pcommon_' exp]), 'svg');
end

%% Model 2: x = sigmaA (Com - NCom)

% Fit linear model using fitlm (includes confidence intervals)
mdl2 = fitlm(x2, y);
rho2 = sqrt(mdl2.Rsquared.Ordinary);
pval2 = mdl2.Coefficients.pValue(2);

if strcmp(model,'fixCrit')
    min_x = -10;
    max_x = 10;
elseif strcmp(model,'bci')
    min_x = -6;
    max_x = 6;
end

% Sort x for plotting the fitted line and CI in order
xfit2 = linspace(min_x, max_x, nsubj)';
[yfit2, yci2] = predict(mdl2, xfit2);  % yci contains lower and upper CI

% Plot regression line
figure;plot(xfit2, yfit2, 'k-', 'LineWidth', 2.5);hold on;

% Plot confidence intervals
plot(xfit2, yci2(:,1), 'k--', 'LineWidth', 1.5);  % lower bound
plot(xfit2, yci2(:,2), 'k--', 'LineWidth', 1.5);  % upper bound

% Scatter plot
scatter(x2, y, 70, [217 217 217]/255, 'filled',...
    'LineWidth', 1.5, 'MarkerEdgeColor', [100 100 100]/255);

% Axis, labels and legend
box off
set(gca,'FontName', 'Arial');
set(gca,'FontSize', 12);
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2);
xlabel('SigmaA Com - NCom');
if strcmp(model,'fixCrit')
    ylim([-0.4 0.6]);
    xlim([-12 12]);
    set(gca,'XTick', -10:2:10);
elseif strcmp(model,'bci')
    ylim([-0.4 0.6]);
    xlim([-7 7]);
    set(gca,'XTick', -6:1:6);
end

set(gca,'YLabel', text('String', 'W_{AV} (a.u.) Com - NCom'));

% Save figure
saveas(gcf, fullfile(dataPath,['corr_WAV_sigmaA_' model '_' exp]), 'svg');