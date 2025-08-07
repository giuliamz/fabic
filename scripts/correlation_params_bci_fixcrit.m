%% Correlation analysis pcommon and criterion parameters
% Fixed Criterion Separated model: criterion [Com - NCom]
% BCI Separated model: pcommon [Com - NCom]

clear;
close all;
clc;

%% General settings

% Experiment
exp='exp1'; % exp1; exp2
mainPath = fullfile('<your path>',exp);
dataPath = fullfile(mainPath,'results',['results_' exp]); % data folder

% Subjects
nsubj=34;

% Load data
if strcmp(exp,'exp1')
    load(fullfile(dataPath,'params_diff_fixCrit_Com_NCom_exp1.mat'));
    x = pcommon_diff_Com_NCom_exp1;
    load(fullfile(dataPath,'params_diff_bci_Com_NCom_exp1.mat'));
    y = pcommon_diff_Com_NCom_exp1;
elseif strcmp(exp,'exp2')
    load(fullfile(dataPath,'params_diff_fixCrit_Com_NCom_exp2.mat'));
    x = pcommon_diff_Com_NCom_exp2;
    load(fullfile(dataPath,'params_diff_bci_Com_NCom_exp2.mat'));
    y = pcommon_diff_Com_NCom_exp2;
end

%% Model: x = criterion (Com - NCom); y = pcommon (Com - NCom)

% Fit linear model using fitlm (includes confidence intervals)
mdl = fitlm(x, y);
rho = sqrt(mdl.Rsquared.Ordinary);
pval = mdl.Coefficients.pValue(2);

min_x = -10;
max_x = 10;

% Sort x for plotting the fitted line and CI in order
xfit = linspace(min_x, max_x, nsubj)';
[yfit, yci] = predict(mdl, xfit);  % yci contains lower and upper CI

% Plot regression line
figure;plot(xfit, yfit, 'k-', 'LineWidth', 2.5);hold on;

% Plot confidence intervals
plot(xfit, yci(:,1), 'k--', 'LineWidth', 1.5);  % lower bound
plot(xfit, yci(:,2), 'k--', 'LineWidth', 1.5);  % upper bound

% Scatter plot
scatter(x, y, 70, [217 217 217]/255, 'filled',...
    'LineWidth', 1.5, 'MarkerEdgeColor', [100 100 100]/255);

% Axis, labels and legend
box off
set(gca,'FontName', 'Arial');
set(gca,'FontSize', 12);
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1.2);
xlabel('Criterion Com - NCom');
set(gca,'XTick', min_x:2:max_x);
ylim([-0.6 0.7]);
xlim([-12 12]);
set(gca,'YLabel', text('String', 'pCommon Com - NCom'));

% Save figure
saveas(gcf, fullfile(dataPath,['corr_criterion_pcommon_' exp]), 'svg');