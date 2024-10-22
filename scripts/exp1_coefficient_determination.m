%% Experiment 1: compute coefficient of determination (R^2)

clear;
close all;
clc;

% Path to BCI scripts
addpath(genpath('<your path>'));
% Data path
dataPath='<your path>';

% Subjects
subjID = {'05';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';...
    '38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';...
    '51';'52';'53';'54';'55';'56';'57';'58';'59'};
nsubj=length(subjID);

%% logLike models of interest

% Initialize (n subjects; 4 models)
logLike_group=nan(nsubj,4);
for iSubj = 1:nsubj
    % model 1: BCI no com action
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_bci_simulations_avert_nosac_pool_best10']),'fm_bestlogLike');
    logLike_group(iSubj,1)=fm_bestlogLike{:};
    % model 2: BCI com action (1=com; 2=non-com)
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_bci_simulations_avert_nosac_sep_best10']),'fm_bestlogLike');
    logLike_group(iSubj,2)=fm_bestlogLike{:};
    % model 3: FF no com action
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_fus_simulations_avert_nosac_pool_best10']),'fm_bestlogLike');
    logLike_group(iSubj,3)=fm_bestlogLike{:};
    % model 4: FF com action (1=com; 2=non-com)
    load(fullfile(dataPath,'results',['s' subjID{iSubj} '_fus_simulations_avert_nosac_sep_best10']),'fm_bestlogLike');
    logLike_group(iSubj,4)=fm_bestlogLike{:};
end

%% logLike null model

logLike_null=nan(nsubj,4);
R2 = nan(nsubj,4);
maxR2 = nan(nsubj,1);

for iSubj = 1:nsubj
    
    % Load data file
    load(fullfile(dataPath,['Log_s' subjID{iSubj} '.mat']));
    bciData = LogFile;    
    % Get rid of missed responses
    bciData = renamevars(bciData,["Vis","Aud","Res"], ...
        ["V","A","Response"]);
    bciData = bciData(~isnan(bciData.Response),:);
    
    % Prepare table for model fitting
    % Put locA and locV
    bciData.Properties.VariableNames=["locV","locA","ResMod","Action","Response","RT"];
    
    % Put respA and respV
    bciData.respA = NaN(size(bciData,1),1);
    bciData.respV = NaN(size(bciData,1),1);
    for j = 1:size(bciData,1)
        % respA
        if bciData.ResMod(j)==1
            bciData.respA(j) = bciData.Response(j);
            bciData.respV(j) = NaN;
        % respV
        elseif bciData.ResMod(j)==2
            bciData.respV(j) = bciData.Response(j);
            bciData.respA(j) = NaN;
        end
    end
    
    % Keep only useful data (locaV, locA, respV, respA)
    bciData = bciData(:,{'locV','locA','respV','respA'});
    responseLoc = unique(bciData.locV);
    
    % Fit null model
    logLike_null(iSubj) = fitModelNull(bciData,responseLoc);
    
    %% Compute R2 for discrete responses based on Nagelkerke (1991)
    n = size(bciData,1);
    maxR2(iSubj,1) = 1-exp((2/n)*(-logLike_null(iSubj)));
    
    % we invert logLike_null and logLike_modelX because we have negative logLike
    R2(iSubj,1) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_group(iSubj,1)));
    R2(iSubj,2) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_group(iSubj,2)));
    R2(iSubj,3) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_group(iSubj,3)));
    R2(iSubj,4) = 1-exp((-2/n)*(logLike_null(iSubj)-logLike_group(iSubj,4)));   
    % normalization by maxR2
    R2(iSubj,1) = R2(iSubj,1)/maxR2(iSubj,1);
    R2(iSubj,2) = R2(iSubj,2)/maxR2(iSubj,1);
    R2(iSubj,3) = R2(iSubj,3)/maxR2(iSubj,1);
    R2(iSubj,4) = R2(iSubj,4)/maxR2(iSubj,1);
    
end % End of loop over subjects

meanR2 = mean(R2);
semR2 = std(R2)/sqrt(length(subjID));