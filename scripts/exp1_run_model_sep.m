%% Experiment 1 models estimation: 
% Forced Fusion Separated
% Bayesian Causal Inference Separated

%% Preparations
clc;
clear;

% Checking setup for parallel port
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    progrMonitor = false;
else
    isServer = false;
    progrMonitor = true;
    dbstop if error;
end

% Opening parallel pool
% if there is no parallel pool running, open one
currPool = gcp('nocreate');
if isempty(currPool)
    if isServer
        parpool('local',16);
    else
        parpool('local');
    end
end

% Add path to BCI scripts and fminsearchbnd function
addpath(genpath('<your path>'));
% Add path to consolidator function
addpath(genpath('<your path>'));
% Data path
dataPath='<your path>';

%% General settings

% Subjects
subjID = {'05';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';...
    '38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';...
    '51';'52';'53';'54';'55';'56';'57';'58';'59'};
% Model
model='bci'; % bci fus
% Number of minima to evaluate with fminsearch
nmin = 10;

%% Run modeling for each subject

for iSubj = 1:length(subjID)
    
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
    
    % Keep only useful data (Action, locaV, locA, respV, respA)
    bciData = bciData(:,{'Action','locV','locA','respV','respA'});
    
    %% Perform gridsearch on a grid of parameters
    
    % Response locations
    responseLoc = unique(bciData.locV);
    % Readout type: model averaging (1) vs. model selection (2)
    readout = 1;
    % For gridsearch, compute only loglike
    my_actions = 0;
    
    % Preparing the parameter space
    % Grid resolution
    n = 5;
    % Parameter settings
    if strcmp(model, 'fus')
        p_common=1;
    elseif strcmp(model, 'segA') || strcmp(model, 'segV') || strcmp(model, 'taskRel')
        p_common=0;
    else
        p_common=linspace(0.1,0.9,n);
    end
    sigV = linspace(0.1,30,n);
    sigA = linspace(0.1,30,n);
    sigP = linspace(0.1,30,n);
    % These parameters can be used as well, but not included in this analysis
    % kernelWidth = 1; 'kW'
    % muP = 0; % mean of spatial prior
    
    % Parameter names
    if strcmp(model, 'fus') || strcmp(model, 'taskRel')
        parameterNames = {'sigP','sigA','sigV'};
        gridVectors = {sigP,sigA,sigV};
    elseif strcmp(model, 'segA')
        parameterNames = {'sigP','sigA'};
        gridVectors = {sigP,sigA};
    elseif strcmp(model, 'segV')
        parameterNames = {'sigP','sigV'};
        gridVectors = {sigP,sigV};
    else
        parameterNames = {'p_common','sigP','sigA','sigV'};
        gridVectors = {p_common,sigP,sigA,sigV};
    end
    
    % Full factorial expansion of the specified parameter vectors
    nParameters = numel(gridVectors);
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectors{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinations = cat(2,coords{:});
    
    % This is actually relevant if you have two or more V reliabilities (it works with one nevertheless)
    sigNames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*')));
    otherParamIdx = cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*'));
    comLevels = unique(bciData.Action);
    nLevelsCom = numel(comLevels);
    
    % Initialise structures for subject results
    grid_logLike_com = cell(1,nLevelsCom);
    grid_10bestIdx = cell(1,nLevelsCom);
    grid_10bestParameters = cell(1,nLevelsCom);
    grid_10bestlogLike = cell(1,nLevelsCom);
    fm_10bestParameters = cell(1,nLevelsCom);
    fm_10bestErrorval = cell(1,nLevelsCom);
    fm_10bestLogLike = cell(1,nLevelsCom);
    fm_bestParameters = cell(1,nLevelsCom);
    fm_bestlogLike =  cell(1,nLevelsCom);
    bciSimulations = struct([]);
    
    % Initialise log likelihood output
    logLike_sum = 0;
    
    for iCom = 1:nLevelsCom
        
        % Select data for each action level (Communicative = 1; Non-Communicative = 2)
        actData = bciData((bciData.Action == comLevels(iCom)),:);
        
        % Select variance parameter
        actSig = paramCombinations(:,ismember(parameterNames,sigNames));
        
        % Generate the parameter setting for one condition
        actParameters = [paramCombinations(:,otherParamIdx),actSig];
        actParamNames = [parameterNames(otherParamIdx),'sigV'];
        
        % Pre-allocating array for data collection
        logLike_all = NaN(size(paramCombinations,1),1);
        
        % Starting the timer and printing details
        cStart = clock;
        fprintf('Performing gridsearch... \n');

        % Performing grid search on the specified parameter space
        switch model
            case 'bci'
                parfor i = 1:size(actParameters,1)
                    [logLike_all(i)] = bci_fitmodel(actParameters(i,:),actParamNames,actData,responseLoc,readout);
                end
            case 'fus'
                parfor i = 1:size(actParameters,1)
                    [logLike_all(i)] = fitModelFus(actParameters(i,:),actParamNames,actData,responseLoc);
                end
        end
        
        % Printing elapsed time
        cTemp = clock;
        fprintf('gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
            datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
               
        % Finding the n minima of the log-likelihood values and the corresponding parameter combination
        % Sort gridsearch log-likelihoods from the 1st best (global minimum)
        grid_logLike_com{iCom} = logLike_all;
        logLike_sort = sort(grid_logLike_com{iCom});
        
        % Get id number and parameter combination for the n best log-likelihoods
        for i = 1:nmin
            x=find(grid_logLike_com{iCom}==logLike_sort(i));
            grid_10bestIdx{iCom}(i,1) = x(1);
            grid_10bestParameters{iCom}{i,1} = actParameters(grid_10bestIdx{iCom}(i),:);
        end
        
        % Get n best log-likelihoods and save in subj-specific structure
        grid_10bestlogLike{iCom}(:,1) = logLike_sort(1:nmin);
        
        %% Perform fminsearch to refine the best parameters
        % Take the best parameters found in the gridsearch and use them as starting point in the fmincon

        % Options for fminsearchbnd
        opts = optimset('fminsearch');
        opts.Display = 'iter';
        opts.MaxFunEvals = 4000;
        opts.TolFun = 1.e-12;
        opts.MaxIter = 1000;
        
        % Set upper and lower bounds on parameters
        % 'p_common','sigP','sigA','sigV'
        if strcmp(model, 'fus') || strcmp(model, 'taskRel')
            LB = [0.001 0.001 0.001];
            UB = [30    30    30];
        elseif strcmp(model, 'segA') || strcmp(model, 'segV')
            LB = [0.001 0.001];
            UB = [30    30];
        else
            LB = [0 0.001 0.001 0.001 ];
            UB = [1 30    30    30    ];
        end
        
        % Creating anonymous function for input to fmincon       
        switch model
            case 'bci'
                fun = @(param) bci_fitmodel(param,actParamNames,actData,responseLoc,readout);
            case 'fus'
                fun = @(param) fitModelFus(param,actParamNames,actData,responseLoc);
        end
        
        % Performing fminsearch for each parameter combination
        for i = 1:nmin
            fprintf(['Performing fminsearch ' num2str(i) ' of subject ' num2str(iSubj) '... \n']);           
            [fm_10bestParameters{iCom}{i,1},fm_10bestErrorval{iCom}(i,1)] = fminsearchbnd(fun,grid_10bestParameters{iCom}{i,1},LB,UB,opts);
        end
        
        % Finishing timer and printing elapsed time
        fprintf('fminsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
            datestr(etime(clock,cTemp)/86400,'dd HH:MM:SS'));
        
        %% Evaluating model fit of best model
        
        % Save all simulation data and plot
        my_actions = 1;
        
        % Evaluate model fit of the n best parameters combinations
        switch model
            case 'bci'
                for i = 1:nmin
                    [fm_10bestLogLike{iCom}(i,1)] = bci_fitmodel(fm_10bestParameters{iCom}{i,1},actParamNames,actData,responseLoc,readout);
                end
            case 'fus'
                for i = 1:nmin
                    [fm_10bestLogLike{iCom}(i,1)] = fitModelFus(fm_10bestParameters{iCom}{i,1},actParamNames,actData,responseLoc);
                end
        end
        
        % Find the best parameters combination after fminsearch
        bestSim = find(fm_10bestLogLike{iCom}==min(fm_10bestLogLike{iCom}));
        fm_bestParameters{iCom} = fm_10bestParameters{iCom}{bestSim};
        
        % Produce bci simulations for the best parameters combination
        switch model
            case 'bci'
                [fm_bestlogLike{iCom},all] = bci_fitmodel(fm_bestParameters{iCom},actParamNames,actData,responseLoc,readout);
            case 'fus'
                [fm_bestlogLike{iCom},all] = fitModelFus(fm_bestParameters{iCom},actParamNames,actData,responseLoc);
        end        
        temp = repmat({comLevels(iCom)},size(all));
        [all.comAction] = temp{:};
        bciSimulations = cat(2,bciSimulations,all);
        
        % Sum logLikes over action levels
        logLike_sum = logLike_sum + fm_bestlogLike{iCom};
        
    end
    
    % Computing Bayesian Information Criterion and Akaike Information Criterion
    ndata = size(bciData,1);
    k = numel(actParamNames);
    
    bic = -logLike_sum-0.5*k*log(ndata);
    aic = 2*k + 2*logLike_sum;
    
    %% Saving subject specific bci simulations
    fprintf('\n\nSaving data...\n\n');
    save(fullfile(dataPath,'results',['s' subjID{iSubj} '_' model '_simulations_avert_nosac_sep_best10']),...
        'subjID', 'grid_logLike_com', 'grid_10bestIdx', 'grid_10bestParameters',...
        'grid_10bestlogLike', 'fm_10bestParameters', 'fm_10bestErrorval',...
        'fm_10bestLogLike', 'fm_bestParameters', 'fm_bestlogLike', 'bciSimulations',...
        'logLike_sum', 'ndata', 'k', 'bic', 'aic', '-v7.3');
    
end % End of loop over subjects
