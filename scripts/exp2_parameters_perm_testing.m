%% Experiment 2: permutation testing of winning model parameters
% Model: Fixed Criterion Separated

clear;
clc;

%% Exp info

Tail = 'both';
nperm=28;
igroup=1;
group_num=1;

dataPath='<your path>';

% Subjects
subjID = {'70';'71';'72';'73';'75';'76';'77';'78';'79';'80';'81';...
    '83';'84';'85';'86';'87';'88';'89';'93';'94';'95';'98';'99';'101';...
    '102';'103';'104';'106';'107';'111';'112';'113';'114';'115'};
nsubj=length(subjID);

% Model
model = 'fixCrit'; % bci fixCrit

load(fullfile(dataPath,['params_diff_' model '_Com_NCom_exp2.mat']));

%% Effect of action (Com-NCom)

pcommon_mean=mean(pcommon_diff_Com_NCom_exp2);
sigA_mean=mean(sigmaA_diff_Com_NCom_exp2);
sigV_mean=mean(sigmaV_diff_Com_NCom_exp2);

% Gets all the possible permutations (via cartesian product)
for iSub=1:nperm
    sets{iSub} = [-1 1];
end

if nperm == 12
    [a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
    ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:), l(:)];
else
    [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, z,...
        a1, b1, c1, d1, e1] = ndgrid(sets{:});
    ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:), k(:), l(:),...
        m(:), n(:), o(:), p(:), q(:), r(:), s(:), t(:), u(:), v(:), z(:),...
        a1(:), b1(:), c1(:), d1(:), e1(:)];
end

pcommon=pcommon_diff_Com_NCom_exp2(1:nperm);
sigA=sigmaA_diff_Com_NCom_exp2(1:nperm);
sigV=sigmaV_diff_Com_NCom_exp2(1:nperm);

% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp = (ToPermute(iPerm,:))';
    %% main effects and interactions
    pcommonPerms(iPerm,:) = mean(pcommon.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    sigAPerms(iPerm,:) = mean(sigA.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
    sigVPerms(iPerm,:) = mean(sigV.*repmat(tmp,1,group_num,1)); %#ok<*SAGROW>
end

ES(1)=mean(pcommon-mean(pcommonPerms));
ES(2)=mean(sigA-mean(sigAPerms));
ES(3)=mean(sigV-mean(sigVPerms));

CI(1,1)=mean(pcommon-mean(pcommonPerms))-1.96*(std(pcommon-mean(pcommonPerms))/sqrt(nperm));
CI(1,2)=mean(pcommon-mean(pcommonPerms))+1.96*(std(pcommon-mean(pcommonPerms))/sqrt(nperm));
CI(2,1)=mean(sigA-mean(sigAPerms))-1.96*(std(sigA-mean(sigAPerms))/sqrt(nperm));
CI(2,2)=mean(sigA-mean(sigAPerms))+1.96*(std(sigA-mean(sigAPerms))/sqrt(nperm));
CI(3,1)=mean(sigV-mean(sigVPerms))-1.96*(std(sigV-mean(sigVPerms))/sqrt(nperm));
CI(3,2)=mean(sigV-mean(sigVPerms))+1.96*(std(sigV-mean(sigVPerms))/sqrt(nperm));

for i = 1:group_num
    if strcmp(Tail,'left')
        % check the proportion of permutation results that are inferior to
        % the mean of my sample
        P(i,1) = sum(pcommonPerms(:,i)<pcommon_mean(i))/numel(pcommonPerms(:,i));
        P(i,2) = sum(sigAPerms(:,i)<sigA_mean(i))/numel(sigAPerms(:,i));
        P(i,3) = sum(sigVPerms(:,i)<sigV_mean(i))/numel(sigVPerms(:,i));
    elseif strcmp(Tail,'right')
        % same but the other way
        P(i,1) = sum(pcommonPerms(:,i)>pcommon_mean(i))/numel(pcommonPerms(:,i));
        P(i,2) = sum(sigAPerms(:,i)>sigA_mean(i))/numel(sigAPerms(:,i));
        P(i,3) = sum(sigVPerms(:,i)>sigV_mean(i))/numel(sigVPerms(:,i));
    elseif strcmp(Tail,'both')
        % for the 2 tailed just compare to the distribution of absolute value of the distance
        % between the result of each permutation to the mean of all
        % permutation results
        % Then you check the proportion of those distances are superior to
        % the distance between the mean of your sample and the mean of all
        % permutation results
        % P(i) = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
        % Actually not just take the absolute values: the above assumes
        % that your null distribution is symmetric
        P(i,1) = sum( abs(pcommonPerms(:,i)) > abs(pcommon_mean(i) ) )  / numel(pcommonPerms(:,i)) ;
        P(i,2) = sum( abs(sigAPerms(:,i)) > abs(sigA_mean(i) ) )  / numel(sigAPerms(:,i)) ;
        P(i,3) = sum( abs(sigVPerms(:,i)) > abs(sigV_mean(i) ) )  / numel(sigVPerms(:,i)) ;
    end
end

res_list{1}='behav';

for i=1:group_num
    res_list{i,2}=pcommon_mean(i);
    res_list{i,3}=P(i,1);
    res_list{i,4}=sigA_mean(i);
    res_list{i,5}=P(i,2);
    res_list{i,6}=sigV_mean(i);
    res_list{i,7}=P(i,3);
end