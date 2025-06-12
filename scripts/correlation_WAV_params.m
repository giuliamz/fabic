%% Correlation analysis WAV and parameters winning model
% WAV Auditory responses [Com - NCom]
% Parameters winning model (BCI Separated): pcommon [Com - NCom], sigmaA [Com - NCom]

clear;
close all;
clc;

%% General settings

% Experiment
exp='exp1'; % exp1; exp2
dataPath = fullfile('<your path>',exp); % data folder

% Subjects
nsubj=34;

% Load data
if strcmp(exp,'exp1')
    load(fullfile(dataPath,'params_diff_Com_NCom_exp1.mat'));
    x1 = pcommon_diff_Com_NCom_exp1;
    x2 = sigmaA_diff_Com_NCom_exp1;
    load(fullfile(dataPath,'wav_aud_diff_Com_NCom_exp1.mat'));
    y = WAV_AUD_COM_NC';
elseif strcmp(exp,'exp2')
    load(fullfile(dataPath,'params_diff_Com_NCom_exp2.mat'));
    x1 = pcommon_diff_Com_NCom_exp2;
    x2 = sigmaA_diff_Com_NCom_exp2;
    load(fullfile(dataPath,'wav_aud_diff_Com_NCom_exp2.mat'));
    y = WAV_AUD_COM_NC';
end

%% Model 1: x = pcommon (Com - NCom)

% Fit linear model using fitlm (includes confidence intervals)
mdl1 = fitlm(x1, y);
rho1 = sqrt(mdl1.Rsquared.Ordinary);
pval1 = mdl1.Coefficients.pValue(2);

% Sort x for plotting the fitted line and CI in order
xfit1 = linspace(-0.3, 0.7, nsubj)';
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
xlabel('Pcommon Com - NCom');
ylim([-0.4 0.6]);
xlim([-0.4 0.8]);
set(gca,'XTick', -0.3:0.1:0.7);
set(gca,'YLabel', text('String', 'W_{AV} (a.u.) Com - NCom'));

% Save figure
saveas(gcf, fullfile(dataPath,['corr_WAV_pcommon_' exp]), 'svg');

%% Model 2: x = sigmaA (Com - NCom)

% Fit linear model using fitlm (includes confidence intervals)
mdl2 = fitlm(x2, y);
rho2 = sqrt(mdl2.Rsquared.Ordinary);
pval2 = mdl2.Coefficients.pValue(2);

% Sort x for plotting the fitted line and CI in order
xfit2 = linspace(-6, 6, nsubj)';
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
ylim([-0.4 0.6]);
xlim([-7 7]);
set(gca,'XTick', -6:1:6);
set(gca,'YLabel', text('String', 'W_{AV} (a.u.) Com - NCom'));

% Save figure
saveas(gcf, fullfile(dataPath,['corr_WAV_sigmaA_' exp]), 'svg');