%% Compute WAV based on predictions of winning model (BCI / fixCrit Separated)

close all;
clear;
clc;

%% general settings

% which experiment?
exp='exp1'; % exp1; exp2
mainPath = fullfile('<your path>',exp);
dataPath = fullfile(mainPath,'results'); % data folder
savePath = fullfile(dataPath,['results_' exp]); % folder to save data

% which model?
model = 'bci'; % bci fixcrit

% use group AV congruent responses?
% 1 = yes: compute group AV congruent responses across participants and conditions
% 0 = no: use participant-specific and condition-specific AV congruent responses
AVcon_group = 1;

% index conditions in bciSimulations table (loaded below)
ind_Com = 1:9;
ind_NCom = 10:18;

% index locations in bciSimulations table (loaded below)
% AV congruent
ind_con = [5 9 1 9 1 5; 1 1 5 5 9 9];
% AV incongruent
ind_inc = [2 3 4 6 7 8];
% index by AV disparity
AVdisp = {[1 3 4 6];[2 5]};

% index group AV congruent locations across participants and conditions
ind_con_group = [2 3 1 3 1 2; 1 1 2 2 3 3];

% conditions: 2 Att x 2 Resp
nCond = 4;

% disparity: low high
nDisp = 2;

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

%% load and organize data

pred_resp=nan(length(ind_Com),nCond,nSubj);
for isub=1:nSubj
    % Load current subject dataset
    if strcmp(exp,'exp1')
        file = (fullfile(dataPath,['s' subjID{isub} '_' model '_simulations_avert_nosac_sep_best10']));
    elseif strcmp(exp,'exp2')
        file = (fullfile(dataPath,['s' subjID{isub} '_' model '_simulations_control_nosac_sep_best10']));
    end
    load(file, 'bciSimulations');
    % Extract and reorganize necessary info
    for icond=1:9
        % repA Com
        pred_resp(icond,1,isub) = bciSimulations(ind_Com(icond)).sA_resp;
        % repA NCom
        pred_resp(icond,2,isub) = bciSimulations(ind_NCom(icond)).sA_resp;
        % repV Com
        pred_resp(icond,3,isub) = bciSimulations(ind_Com(icond)).sV_resp;
        % repV NCom
        pred_resp(icond,4,isub) = bciSimulations(ind_NCom(icond)).sV_resp;
    end
end

% Compute group AV congruent responses across participants and conditions
pred_resp_group = mean(pred_resp,3);
AVloc_con = mean(pred_resp_group([1 5 9],:),2);

%% compute WAV

for isub=1:nSubj
    for icon = 1:nCond
        % WAV by AV incong loc, action, response
        if AVcon_group == 1
            % Option 1: use group AV congruent responses across participants and conditions
            for iloc = 1:numel(ind_inc)
                AVloc1 = AVloc_con(ind_con_group(1,iloc));
                AVloc2 = AVloc_con(ind_con_group(2,iloc));
                WAV(iloc,icon,isub) = (pred_resp(ind_inc(iloc),icon,isub) - AVloc1) / (AVloc2 - AVloc1);
            end
        else
            % Option 2: use participant-specific and condition-specific AV congruent responses
            for iloc = 1:numel(ind_inc)
                AVloc1 = pred_resp(ind_con(1,iloc),icon,isub);
                AVloc2 = pred_resp(ind_con(2,iloc),icon,isub);
                WAV(iloc,icon,isub) = (pred_resp(ind_inc(iloc),icon,isub) - AVloc1) / (AVloc2 - AVloc1);
            end
        end
        % WAV by AV disp, action, response
        for idisp = 1:nDisp
            WAV_cond8(idisp,icon,isub) = mean(WAV(AVdisp{idisp},icon,isub));
        end
        % WAV by action, response
        WAV_cond4(icon,isub) = mean(WAV(:,icon,isub));
    end
end

%% group means and std/sem

% WAV by AV incong loc, action, response
WAV_group_mean = mean(WAV,3);
WAV_group_std = std(WAV,0,3);
WAV_group_sem = WAV_group_std/sqrt(nSubj);

% WAV by AV disp, action, response
WAV_cond8_group_mean = mean(WAV_cond8,3);
WAV_cond8_group_std = std(WAV_cond8,0,3);
WAV_cond8_group_sem = WAV_cond8_group_std/sqrt(nSubj);

% WAV by action, response
WAV_cond4_group_mean = mean(WAV_cond4,2);
WAV_cond4_group_std = std(WAV_cond4,0,2);
WAV_cond4_group_sem = WAV_cond4_group_std/sqrt(nSubj);

%% plot WAV by AV disp, action, response

% plot matrices
plot_matrix_line=[WAV_cond8_group_mean(1) WAV_cond8_group_mean(3); WAV_cond8_group_mean(2) WAV_cond8_group_mean(4);
    WAV_cond8_group_mean(5) WAV_cond8_group_mean(7); WAV_cond8_group_mean(6) WAV_cond8_group_mean(8)];

plot_matrix_dot = [WAV_cond8_group_mean(1); WAV_cond8_group_mean(3); WAV_cond8_group_mean(2); WAV_cond8_group_mean(4);
    WAV_cond8_group_mean(5); WAV_cond8_group_mean(7); WAV_cond8_group_mean(6); WAV_cond8_group_mean(8)];

plot_matrix_sem = [WAV_cond8_group_sem(1); WAV_cond8_group_sem(3); WAV_cond8_group_sem(2); WAV_cond8_group_sem(4);
    WAV_cond8_group_sem(5); WAV_cond8_group_sem(7); WAV_cond8_group_sem(6); WAV_cond8_group_sem(8)];

% size
positionXY = [0, 0, 170, 350];
figure('color', [1 1 1], 'Position', positionXY);

% lines
conLineStyles={
    '-'
    '--'
    '-'
    '--'
    };

% colors
cols.Aud1 = [230 97 1]/255; % Aud Com
cols.Aud2 = [253 184 99]/255; % Aud NCom
cols.Vis1 = [94 60 153]/255; % Vis Com
cols.Vis2 = [178 171 210]/255; % Vis NCom

conCols=[cols.Aud1;cols.Aud2;cols.Aud1;cols.Aud2;...
    cols.Vis1;cols.Vis2;cols.Vis1;cols.Vis2];

% plot lines
for i=1:4
    hold on
    if i==1 || i==3
        k=0.05;
    else
        k=-0.05;
    end
    plot((1:2)+k,plot_matrix_line(i,:),...
        'Color','k',...
        'LineStyle',conLineStyles{i},...
        'LineWidth',1.5)
end

% plot mean and sem
a=repmat([1 2],1,4);
for i=1:8
    if ismember(i, [1,2,5,6])
        k=0.05;
    else
        k=-0.05;
    end
    hold on
    line([a(i)+k a(i)+k],[plot_matrix_dot(i)+plot_matrix_sem(i) ...
        plot_matrix_dot(i)-plot_matrix_sem(i)],...
        'Color',conCols(i,:),'LineWidth',1.5);hold on;
    plot(a(i)+k,plot_matrix_dot(i),'o',...
        'Color',conCols(i,:),...
        'LineWidth',1.5,'MarkerSize',6,...
        'MarkerEdgeColor',conCols(i,:),...
        'MarkerFaceColor',conCols(i,:));hold on;
end

% adjust x and y ticks
xl = [0.7 2.3]; xlim(xl);
yl = [0 1.05]; ylim(yl);
% adjust x and y ticks
set(gca,'FontName', 'Arial');
set(gca,'FontSize', 8.5);
set(gca,'YAxisLocation','left');
set(gca,'XTick', [1 2]);
set(gca,'XTickLabel', {'Com';'NCom'});
set(gca,'YLabel', text('String', '\it{W}_{\it{AV}} (a.u.)', 'FontSize', 11))
set(gca,'YTick', (0:0.1:1));
set(gca,'YTickLabel', (0:0.1:1));
set(gca,'TickLength', [0.01 0.01]);
set(gca,'LineWidth',1);
saveas(gcf, fullfile(savePath,['WAV_pred_' exp '_' model]), 'svg');