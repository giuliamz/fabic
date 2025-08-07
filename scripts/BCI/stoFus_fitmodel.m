function [logLike,all] = stoFus_fitmodel(parameters,parameterNames,dataVA,responseLoc)
% Simulates responses from the Stochastic Fusion model for a given parameter
% set (see Acerbi et al., 2018, PLoS Comp Bio, cookbook paper).
% 
% DETAILS: 
%   Creating 10000 (e.g.) internal samples and creating a distribution of 
%   potential responses based on these. These distributions are then 
%   compared with the observed distributions during estimation to compute 
%   loglike for this parameter setting.
% USAGE:
%   [logLike,all] = bci_fitmodelFull(parameters,parameterNames,dataVA,responseLoc,readouttype)
% INPUT:
%   parameters (numeric vector): BCI model parameters. I would rather have 
%       used a table instead of using parameters and parameterNames as 
%       separate inputs, but this input has to be compatible with
%       fminsearch. 
%   parameterNames (cell array of strings): BCI model parameter names.
%       Possible values: 'p_common','muP','sigP','kW','sigA','sigV'
%   dataVA (table): behavioural response data with which the 
%       predicted responses of the BCI model will be compared. Given in
%       visual degrees. 
%       Rows = trials, Variables: locA, locV, respA, respV
%   responseLoc is the locations corresponding to the response bottons, e.g.
%       degrees. If a single value, continuous response with a
%       kernel width is used.
%   decisionFun = 1: model averaging, 2: model selection, 3: probability
%       matching
% OUTPUT:
%   logLike:
%   all: 

% This code was contributed by Ulrik Beierholm (Apr 2016), with some code
% left over from Wei Ji Ma and Konrad Koerding at some time 2007-2009.
% Further expanded and adjusted by Uta Noppeney & Mate Aller.
% Finally changed from BCI model to Stochastic Fusion model by Ambra Ferrari.

% Fix the initialising random state - WHY ARE WE USING A FIXED SEED? --> 
% 
% DM: Fixed seed ensures no random noise in different calls of this 
% function. Hence, it helps simple optimizers (e.g., fminsearch) to find 
% the optimal parameter values. However, it slightly reduces the 
% generalizability of the fitted solution, because it depends on the fixed
% seed. BADS (with stochasticity support) also works well with a random 
% seed, and this combination is thus preferred. However, be aware that it 
% may then be impossible to exactly replicate the fitting solutions, unless
% you save the RNG seed when running this function. 

% rand('state',55); 
% randn('state',55);
rng(55);

nIntSamples = 10000;  %% minimum required for model fitting

% Checking dataVA
reqDataVars = {'locA','locV','respA','respV'};
if any(~ismember(reqDataVars,dataVA.Properties.VariableNames))
    error('bci_fitmodel:missingDataVariable',...
        'Missing variable from dataVA in bci_fitmodel.')
end

% Ensure that only data are included where both A and V signals are used
dataVA = dataVA(~isnan(dataVA.locV) & ~isnan(dataVA.locA),:);

% Check responseLoc
if iscolumn(responseLoc)
    responseLoc = responseLoc';
end

% To save computation time, check for what stimulus conditions were
% presented, then simulate each one
conditions = unique(dataVA(:,{'locV','locA'}));
nCond = size(conditions,1);
dataConditions = cell(nCond,1);
for i = 1:nCond
    match = ismember(dataVA(:,{'locV','locA'}),conditions(i,:));
    dataConditions{i} = table2array((dataVA(match,{'respV','respA'})));
end

% Parameters
pFus = parameters(ismember(parameterNames,'pFus'));
sigP = parameters(ismember(parameterNames,'sigP'));
sigA = parameters(ismember(parameterNames,'sigA'));
sigV = parameters(ismember(parameterNames,'sigV'));

% Setting default for muP if it is not specified
if any(ismember(parameterNames,'muP'))
    xP = parameters(ismember(parameterNames,'muP'));
else
    xP = 0;
end

% discrete responses?? i.e. response loc vector e.g. [1 2 3 4 5]
% or continuous responseLoc = 1
if length(responseLoc) == 1
    kernelWidth = parameters(ismember(parameterNames,'kW'));  %% kernel width as a parameter that is fitted to the data
else
    kernelWidth = 1;  % set to 1 for discrete responses
end

% Throw an error if there is a missing parameter
if any(cellfun(@isempty,{pFus,xP,sigP,sigA,sigV,kernelWidth}))
    error('stoFus_fitmodel:missingParameter','Missing parameter in stoFus_fitmodel.')
end

%Throw an error if any of the strictly positive parameters is zero or smaller  
if any([pFus,sigP,sigA,sigV,kernelWidth] <= 0)
    error('fixCrit_fitmodel:nonPositiveParameter','Unexpected non-positive parameter in fixCrit_fitmodel.')
end

% Variances of A and V and prior
varV = sigV^2;
varA = sigA^2;
varP = sigP^2;

% Variances of estimates given common or independent causes
varVA_hat = 1/(1/varV + 1/varA + 1/varP);
varV_hat = 1/(1/varV + 1/varP);
varA_hat = 1/(1/varA + 1/varP);

% % Variances used in computing probability of common or independent causes
% var_common = varV * varA + varV * varP + varA * varP;
% varV_indep = varV + varP;
% varA_indep = varA + varP;

% Initialze variable to collect log-likelihood values for each condition
logLikeCond = NaN(nCond,1);

% Initialize output structure
if nargout > 1
    all = struct;
end

% Simulate responses for each condition
for indCond = 1:nCond
    
    sV = conditions.locV(indCond);%sourcevec(sVind);
    sA = conditions.locA(indCond);%sourcevec(sAind);
    
    % Generation of fake data
    xV = sV + sigV * randn(nIntSamples,1);
    xA = sA + sigA * randn(nIntSamples,1);
    
    % Estimates given common or independent causes
    s_hat_common = (xV/varV + xA/varA + xP/varP) * varVA_hat;
    sV_hat_indep = (xV/varV + xP/varP) * varV_hat;
    sA_hat_indep = (xA/varA + xP/varP) * varA_hat;
    
%     % Probability of common or independent causes
%     quad_common = (xV-xA).^2 * varP + (xV-xP).^2 * varA + (xA-xP).^2 * varV;
%     quadV_indep = (xV-xP).^2;
%     quadA_indep = (xA-xP).^2;
%     
%     % Likelihood of observations (xV,xA) given C, for C=1 and C=2
%     likelihood_common = exp(-quad_common/(2*var_common))/(2*pi*sqrt(var_common));
%     likelihoodV_indep = exp(-quadV_indep/(2*varV_indep))/sqrt(2*pi*varV_indep);
%     likelihoodA_indep = exp(-quadA_indep/(2*varA_indep))/sqrt(2*pi*varA_indep);
%     likelihood_indep =  likelihoodV_indep .* likelihoodA_indep;
%     
%     % Posterior probability of C given observations (xV,xA)
%     post_common = likelihood_common * p_common;
%     post_indep = likelihood_indep * (1-p_common);
%     pC = post_common./(post_common + post_indep);
%     
%     % ----------------------------------------------------------------------
%     % Generate spatial location responses and compute loglike
%     switch decisionFun
%         
%         case 1
%             % Mean of posterior - Model averaging
%             % Overall estimates: weighted averages
%             sV_hat = pC .* s_hat_common + (1-pC) .* sV_hat_indep;
%             sA_hat = pC .* s_hat_common + (1-pC) .* sA_hat_indep;
%             
%         case 2
%             % Model selection instead of model averaging
%             % averaging
%             sV_hat = (pC>0.5).*s_hat_common + (pC<=0.5).*sV_hat_indep;
%             sA_hat = (pC>0.5).*s_hat_common + (pC<=0.5).*sA_hat_indep;
%             
%         case 3
%             % Probability matching
%             thresh = rand(1,1); % to enable probability matching
%             sV_hat = (pC>thresh).*s_hat_common + (pC<=thresh).*sV_hat_indep;
%             sA_hat = (pC>thresh).*s_hat_common + (pC<=thresh).*sA_hat_indep;
%             
%     end
    
    % Apply stochastic fusion --> Ambra Ferrari (28-7-2025)
    sV_hat = (pFus>0.5).*s_hat_common + (pFus<=0.5).*sV_hat_indep;
    sA_hat = (pFus>0.5).*s_hat_common + (pFus<=0.5).*sA_hat_indep;

    % compute predicted responses for discrete and continuous case
    if length(responseLoc) > 1
        % discrete responses
        lengthRespLoc=length(responseLoc);
        
        %find the response location closest to sV_hat and
        %sA_hat, ie with minimum deviation
        [~,tV] = min(abs(repmat(sV_hat,1,lengthRespLoc) - repmat(responseLoc,nIntSamples,1)),[],2);
        [~,tA] = min(abs(repmat(sA_hat,1,lengthRespLoc) - repmat(responseLoc,nIntSamples,1)),[],2);
        sV_pred_resp = responseLoc(tV);
        sA_pred_resp = responseLoc(tA);
    else
        %continuous responses used, no discretisation
        sV_pred_resp = sV_hat;
        sA_pred_resp = sA_hat;
    end
    
    % A, V responses given by participant for particular A,V location
    % combination
    dataV = dataConditions{indCond}(:,1)'; dataV = dataV(~isnan(dataV));
    dataA = dataConditions{indCond}(:,2)'; dataA = dataA(~isnan(dataA));
    
    %  Compute loglike for A and V responses
    if length(responseLoc) > 1
        %discrete case
        %calculate frequencies of predictions, at least 0.00001 (1 out of 100,000)
        freq_predV = max(0.00001,hist(sV_pred_resp, responseLoc)/nIntSamples);
        freq_predA = max(0.00001,hist(sA_pred_resp, responseLoc)/nIntSamples);
        
        %calculate frequencies of actual responses
        freq_dataV = hist(dataV, responseLoc);
        freq_dataA = hist(dataA, responseLoc);
        
        %calculate log-likelihood
        logLikeA = sum(freq_dataA .* log(freq_predA)) ;
        logLikeV = sum(freq_dataV .* log(freq_predV)) ;
        
    else
        %continuous case
        %gaussian kernel distribution for each condition
        %for each condition average over gaussian likelihoods
        
        logLikeA=log(mean(normpdf(repmat(dataA,nIntSamples,1) ,repmat(sA_pred_resp,1,length(dataA)), kernelWidth)));
        logLikeV=log(mean(normpdf(repmat(dataV,nIntSamples,1) ,repmat(sV_pred_resp,1,length(dataV)), kernelWidth)));
    end
    
    %---------------------------------------------------------------------
    % sum loglike across different task responses
%     if isnan(logLikeA) && isnan(logLikeV)
%         logLikeCond(indCond) = NaN;
%     else
%         % If all(isnan(logLikeA)) || all(isnan(logLikeV)), then the sum
%         % is going to be 0 whis means, the likelihood is 1, which is not
%         % true.
%         logLikeCond(indCond) = nansum([nansum(logLikeA) nansum(logLikeV)]);
%     end

    logLikeCond(indCond) = nansum([nansum(logLikeA) nansum(logLikeV)]);
    
    %---------------------------------------------------------------------
    % for a particular parameter setting make plots and save biasses etc.
    if nargout > 1
        
        %store values
        all(indCond).sA_resp = mean(sA_pred_resp);
        all(indCond).sV_resp = mean(sV_pred_resp);
        all(indCond).sV_hat = mean(sV_hat);
        all(indCond).sA_hat = mean(sA_hat);
        all(indCond).s_hat_common = mean(s_hat_common);
        all(indCond).sV_hat_indep = mean(sV_hat_indep);
        all(indCond).sA_hat_indep = mean(sA_hat_indep);
        all(indCond).conditions = conditions(indCond,:);
        all(indCond).logLikeA = logLikeA;
        all(indCond).logLikeV = logLikeV;
        all(indCond).parameters = parameters;
        
        if exist('freq_predV','var')   % if discrete responses
            all(indCond).freq_predV = freq_predV;
            all(indCond).freq_predA = freq_predA;
            all(indCond).freq_dataV = freq_dataV/length(dataV);
            all(indCond).freq_dataA = freq_dataA/length(dataA);            
        end
    end
    
end  % end of loop over conditions

%sum over conditions and turn into negative log likelihood
logLike = -sum(logLikeCond); % neg log like for each attention level
    
end
