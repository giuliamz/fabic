clear;
clc;

% t-test direction (left-tailed, right-tailed, both tails)
Tail = 'both'; % 'left' 'right' 'both'

% settings
%cd('put data path');
nperm=12; % don't go below 10

% load data from script 'WAV_rep_fabic.m': 'fabic_subj_wav_exp1.mat' or 'fabic_subj_wav_exp2.mat'
load('fabic_subj_wav_exp2.mat');
exp='2';
results=table2array(mat);

iroi=1;
roi_num=1;

%% main effects and interactions
% for each effect of interest, we compute differences (e.g. High disparity - Low disparity) for each
% subject and then take the mean across subjects
% each column is a different condition;

% main effect of communication (C - NC)
mainAct(:,iroi)=mean(results(:,1:2:8),2)-mean(results(:,2:2:8),2);
mainAct_mean(iroi,1)=mean(mainAct(:,iroi));
% main effect of report (V - A)
mainRep(:,iroi)=mean(results(:,5:8),2)-mean(results(:,1:4),2);
mainRep_mean(iroi,1)=mean(mainRep(:,iroi));
% min effect of AV disparity (High - Low)
mainDisp(:,iroi)=mean(results(:,[3 4 7 8]),2)-mean(results(:,[1 2 5 6]),2);
mainDisp_mean(iroi,1)=mean(mainDisp(:,iroi));
% action x report (CV - NV - (CA - NA))
intActRep(:,iroi)=mean(results(:,[5 7]),2)-mean(results(:,[6 8]),2)-...
    (mean(results(:,[1 3]),2)-mean(results(:,[2 4]),2));
intActRep_mean(iroi,1)=mean(intActRep(:,iroi));
% action x AV disparity (CH - CL - (NH - NL))
intActDisp(:,iroi)=mean(results(:,[3 7]),2)-mean(results(:,[1 5]),2)-...
    (mean(results(:,[4 8]),2)-mean(results(:,[2 6]),2));
intActDisp_mean(iroi,1)=mean(intActDisp(:,iroi));
% report x AV disparity (VH - VL - (AH - AL))
intRepDisp(:,iroi)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2)-...
    (mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2));
intRepDisp_mean(iroi,1)=mean(intRepDisp(:,iroi));
% action x report x AV disparity ((CVH - CVL - (CAH - CAL)) - ((NVH - NVL) - (NAH - NAL)))
intActRepDisp(:,iroi)=(mean(results(:,7),2)-mean(results(:,5),2)-...
    (mean(results(:,3),2)-mean(results(:,1),2)))-...
    ((mean(results(:,8),2)-mean(results(:,6),2))-...
    (mean(results(:,4),2)-mean(results(:,2),2)));
intActRepDisp_mean(iroi,1)=mean(intActRepDisp(:,iroi));

%% follow-up simple main effects
% here it depends which interactions are significant (p-value);
% for the significant ones, we break them down into simple main effects

% report x disparity interaction

intRepADisp(:,iroi)=mean(results(:,[3 4]),2)-mean(results(:,[1 2]),2); % AH - AL
intRepADisp_mean(iroi,1)=mean(intRepADisp(:,iroi));
intRepVDisp(:,iroi)=mean(results(:,[7 8]),2)-mean(results(:,[5 6]),2); % VH - VL
intRepVDisp_mean(iroi,1)=mean(intRepVDisp(:,iroi));

intDispHRep(:,iroi)=mean(results(:,[7 8]),2)-mean(results(:,[3 4]),2); % VH - AH
intDispHRep_mean(iroi,1)=mean(intDispHRep(:,iroi));
intDispLRep(:,iroi)=mean(results(:,[5 6]),2)-mean(results(:,[1 2]),2); % VL - AL
intDispLRep_mean(iroi,1)=mean(intDispLRep(:,iroi));

% % % action x disparity interaction
%CH - CL
intActCDisp(:,iroi)=mean(results(:,[3 7]),2)-mean(results(:,[1 5]),2); 
intActCDisp_mean(iroi,1)=mean(intActCDisp(:,iroi));
%NH - NL
intActNDisp(:,iroi)=mean(results(:,[4 8]),2)-mean(results(:,[2 6]),2); 
intActNDisp_mean(iroi,1)=mean(intActNDisp(:,iroi));
%NH - CH
intDispHAct(:,iroi)=mean(results(:,[4 8]),2)-mean(results(:,[3 7]),2); 
intDispHAct_mean(iroi,1)=mean(intDispHAct(:,iroi));
%NL - CL
intDispLAct(:,iroi)=mean(results(:,[2 6]),2)-mean(results(:,[1 5]),2); 
intDispLAct_mean(iroi,1)=mean(intDispLAct(:,iroi));

% % % action x response interaction
%CA - CV
intActCRep(:,iroi)=mean(results(:,[1 3]),2)-mean(results(:,[5 7]),2); 
intActCRep_mean(iroi,1)=mean(intActCRep(:,iroi));
%NA - NV
intActNRep(:,iroi)=mean(results(:,[2 4]),2)-mean(results(:,[6 8]),2); 
intActNRep_mean(iroi,1)=mean(intActNRep(:,iroi));
%NA - CA
intRepAAct(:,iroi)=mean(results(:,[2 4]),2)-mean(results(:,[1 3]),2); 
intRepAAct_mean(iroi,1)=mean(intRepAAct(:,iroi));
%NV - CV
intRepVAct(:,iroi)=mean(results(:,[6 8]),2)-mean(results(:,[5 7]),2); 
intRepVAct_mean(iroi,1)=mean(intRepVAct(:,iroi));


%% Create null distribution by permutations
% Gets all the possible permutations (via cartesian product)
for iSub=1:nperm
    sets{iSub} = [-1 1];
end

% with 12 permutations: from a to l (12)
% adjust accordingly
[a, b, c, d, e, f, g, h, i, j, k, l] = ndgrid(sets{:});
ToPermute = [a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:),k(:), l(:)];

mainAct=mainAct(1:nperm);
mainRep=mainRep(1:nperm);
mainDisp=mainDisp(1:nperm);
intActRep=intActRep(1:nperm);
intActDisp=intActDisp(1:nperm);
intRepDisp=intRepDisp(1:nperm);
intActRepDisp=intActRepDisp(1:nperm);

intRepADisp=intRepADisp(1:nperm);
intRepVDisp=intRepVDisp(1:nperm);
intDispHRep=intDispHRep(1:nperm);
intDispLRep=intDispLRep(1:nperm);

intActCDisp=intActCDisp(1:nperm);
intActNDisp=intActNDisp(1:nperm);
intDispHAct=intDispHAct(1:nperm);
intDispLAct=intDispLAct(1:nperm);

intActCRep=intActCRep(1:nperm);
intActNRep=intActNRep(1:nperm);
intRepAAct=intRepAAct(1:nperm);
intRepVAct=intRepVAct(1:nperm);



% Compute the null distributions: one for each permutation
for iPerm = 1:size(ToPermute,1)
    tmp = (ToPermute(iPerm,:))';
    % Each row of Perms contains the mean of a permutation; each column of Perms is a different ROI
    %% main effects and interactions
    PermsAct(iPerm,:) = mean(mainAct.*repmat(tmp,1,roi_num,1)); 
    PermsRep(iPerm,:) = mean(mainRep.*repmat(tmp,1,roi_num,1)); 
    PermsDisp(iPerm,:) = mean(mainDisp.*repmat(tmp,1,roi_num,1)); 
    PermsActRep(iPerm,:) = mean(intActRep.*repmat(tmp,1,roi_num,1));
    PermsActDisp(iPerm,:) = mean(intActDisp.*repmat(tmp,1,roi_num,1)); 
    PermsRepDisp(iPerm,:) = mean(intRepDisp.*repmat(tmp,1,roi_num,1)); 
    PermsActRepDisp(iPerm,:) = mean(intActRepDisp.*repmat(tmp,1,roi_num,1)); 
    % %% follow-up simple main effects
    PermsRepADisp(iPerm,:) = mean(intRepADisp.*repmat(tmp,1,roi_num,1)); 
    PermsRepVDisp(iPerm,:) = mean(intRepVDisp.*repmat(tmp,1,roi_num,1)); 
    PermsDispHRep(iPerm,:) = mean(intDispHRep.*repmat(tmp,1,roi_num,1)); 
    PermsDispLRep(iPerm,:) = mean(intDispLRep.*repmat(tmp,1,roi_num,1)); 
    PermsActCDisp(iPerm,:) = mean(intActCDisp.*repmat(tmp,1,roi_num,1));
    PermsActNDisp(iPerm,:) = mean(intActNDisp.*repmat(tmp,1,roi_num,1)); 
    PermsDispHAct(iPerm,:) = mean(intDispHAct.*repmat(tmp,1,roi_num,1)); 
    PermsDispLAct(iPerm,:) = mean(intDispLAct.*repmat(tmp,1,roi_num,1)); 
    PermsActCRep(iPerm,:) = mean(intActCRep.*repmat(tmp,1,roi_num,1)); 
    PermsActNRep(iPerm,:) = mean(intActNRep.*repmat(tmp,1,roi_num,1)); 
    PermsRepAAct(iPerm,:) = mean(intRepAAct.*repmat(tmp,1,roi_num,1)); 
    PermsRepVAct(iPerm,:) = mean(intRepVAct.*repmat(tmp,1,roi_num,1)); 
end

% calculation of effect size
ES(1)=mean(mainAct-mean(PermsAct));
ES(2)=mean(mainRep-mean(PermsRep));
ES(3)=mean(mainDisp-mean(PermsDisp));
ES(4)=mean(intActRep-mean(PermsActRep));
ES(5)=mean(intActDisp-mean(PermsActDisp));
ES(6)=mean(intRepDisp-mean(PermsRepDisp));
ES(7)=mean(intActRepDisp-mean(PermsActRepDisp));
ES(8)=mean(intRepADisp-mean(PermsRepADisp));
ES(9)=mean(intRepVDisp-mean(PermsRepVDisp));
ES(10)=mean(intDispHRep-mean(PermsDispHRep));
ES(11)=mean(intDispLRep-mean(PermsDispLRep));
ES(12)=mean(intActCDisp-mean(PermsActCDisp));
ES(13)=mean(intActNDisp-mean(PermsActNDisp));
ES(14)=mean(intDispHAct-mean(PermsDispHAct));
ES(15)=mean(intDispLAct-mean(PermsDispLAct));
ES(16)=mean(intActCRep-mean(PermsActCRep));
ES(17)=mean(intActNRep-mean(PermsActNRep));
ES(18)=mean(intRepAAct-mean(PermsRepAAct));
ES(19)=mean(intRepVAct-mean(PermsRepVAct));

% 95% confidence interval
CI(1,1)=mean(mainAct-mean(PermsAct))-1.96*(std(mainAct-mean(PermsAct))/sqrt(size(results,1)));
CI(1,2)=mean(mainAct-mean(PermsAct))+1.96*(std(mainAct-mean(PermsAct))/sqrt(size(results,1)));
CI(2,1)=mean(mainRep-mean(PermsRep))-1.96*(std(mainRep-mean(PermsRep))/sqrt(size(results,1)));
CI(2,2)=mean(mainRep-mean(PermsRep))+1.96*(std(mainRep-mean(PermsRep))/sqrt(size(results,1)));
CI(3,1)=mean(mainDisp-mean(PermsDisp))-1.96*(std(mainDisp-mean(PermsDisp))/sqrt(size(results,1)));
CI(3,2)=mean(mainDisp-mean(PermsDisp))+1.96*(std(mainDisp-mean(PermsDisp))/sqrt(size(results,1)));
CI(4,1)=mean(intActRep-mean(PermsActRep))-1.96*(std(intActRep-mean(PermsActRep))/sqrt(size(results,1)));
CI(4,2)=mean(intActRep-mean(PermsActRep))+1.96*(std(intActRep-mean(PermsActRep))/sqrt(size(results,1)));
CI(5,1)=mean(intActDisp-mean(PermsActDisp))-1.96*(std(intActDisp-mean(PermsActDisp))/sqrt(size(results,1)));
CI(5,2)=mean(intActDisp-mean(PermsActDisp))+1.96*(std(intActDisp-mean(PermsActDisp))/sqrt(size(results,1)));
CI(6,1)=mean(intRepDisp-mean(PermsRepDisp))-1.96*(std(intRepDisp-mean(PermsRepDisp))/sqrt(size(results,1)));
CI(6,2)=mean(intRepDisp-mean(PermsRepDisp))+1.96*(std(intRepDisp-mean(PermsRepDisp))/sqrt(size(results,1)));
CI(7,1)=mean(intActRepDisp-mean(PermsActRepDisp))-1.96*(std(intActRepDisp-mean(PermsActRepDisp))/sqrt(size(results,1)));
CI(7,2)=mean(intActRepDisp-mean(PermsActRepDisp))+1.96*(std(intActRepDisp-mean(PermsActRepDisp))/sqrt(size(results,1)));
CI(8,1)=mean(intRepADisp-mean(PermsRepADisp))-1.96*(std(intRepADisp-mean(PermsRepADisp))/sqrt(size(results,1)));
CI(8,2)=mean(intRepADisp-mean(PermsRepADisp))+1.96*(std(intRepADisp-mean(PermsRepADisp))/sqrt(size(results,1)));
CI(9,1)=mean(intRepVDisp-mean(PermsRepVDisp))-1.96*(std(intRepVDisp-mean(PermsRepVDisp))/sqrt(size(results,1)));
CI(9,2)=mean(intRepVDisp-mean(PermsRepVDisp))+1.96*(std(intRepVDisp-mean(PermsRepVDisp))/sqrt(size(results,1)));
CI(10,1)=mean(intDispHRep-mean(PermsDispHRep))-1.96*(std(intDispHRep-mean(PermsDispHRep))/sqrt(size(results,1)));
CI(10,2)=mean(intDispHRep-mean(PermsDispHRep))+1.96*(std(intDispHRep-mean(PermsDispHRep))/sqrt(size(results,1)));
CI(11,1)=mean(intDispLRep-mean(PermsDispLRep))-1.96*(std(intDispLRep-mean(PermsDispLRep))/sqrt(size(results,1)));
CI(11,2)=mean(intDispLRep-mean(PermsDispLRep))+1.96*(std(intDispLRep-mean(PermsDispLRep))/sqrt(size(results,1)));
CI(12,1)=mean(intActCDisp-mean(PermsActCDisp))-1.96*(std(intActCDisp-mean(PermsActCDisp))/sqrt(size(results,1)));
CI(12,2)=mean(intActCDisp-mean(PermsActCDisp))+1.96*(std(intActCDisp-mean(PermsActCDisp))/sqrt(size(results,1)));
CI(13,1)=mean(intActNDisp-mean(PermsActNDisp))-1.96*(std(intActNDisp-mean(PermsActNDisp))/sqrt(size(results,1)));
CI(13,2)=mean(intActNDisp-mean(PermsActNDisp))+1.96*(std(intActNDisp-mean(PermsActNDisp))/sqrt(size(results,1)));
CI(14,1)=mean(intDispHAct-mean(PermsDispHAct))-1.96*(std(intDispHAct-mean(PermsDispHAct))/sqrt(size(results,1)));
CI(14,2)=mean(intDispHAct-mean(PermsDispHAct))+1.96*(std(intDispHAct-mean(PermsDispHAct))/sqrt(size(results,1)));
CI(15,1)=mean(intDispLAct-mean(PermsDispLAct))-1.96*(std(intDispLAct-mean(PermsDispLAct))/sqrt(size(results,1)));
CI(15,2)=mean(intDispLAct-mean(PermsDispLAct))+1.96*(std(intDispLAct-mean(PermsDispLAct))/sqrt(size(results,1)));
CI(16,1)=mean(intActCRep-mean(PermsActCRep))-1.96*(std(intActCRep-mean(PermsActCRep))/sqrt(size(results,1)));
CI(16,2)=mean(intActCRep-mean(PermsActCRep))+1.96*(std(intActCRep-mean(PermsActCRep))/sqrt(size(results,1)));
CI(17,1)=mean(intActNRep-mean(PermsActNRep))-1.96*(std(intActNRep-mean(PermsActNRep))/sqrt(size(results,1)));
CI(17,2)=mean(intActNRep-mean(PermsActNRep))+1.96*(std(intActNRep-mean(PermsActNRep))/sqrt(size(results,1)));
CI(18,1)=mean(intRepAAct-mean(PermsRepAAct))-1.96*(std(intRepAAct-mean(PermsRepAAct))/sqrt(size(results,1)));
CI(18,2)=mean(intRepAAct-mean(PermsRepAAct))+1.96*(std(intRepAAct-mean(PermsRepAAct))/sqrt(size(results,1)));
CI(19,1)=mean(intRepVAct-mean(PermsRepVAct))-1.96*(std(intRepVAct-mean(PermsRepVAct))/sqrt(size(results,1)));
CI(19,2)=mean(intRepVAct-mean(PermsRepVAct))+1.96*(std(intRepVAct-mean(PermsRepVAct))/sqrt(size(results,1)));


% p-value calculation
for i = 1:roi_num
    if strcmp(Tail,'left')
        % check the proportion of permutation results that are inferior to
        % the mean of my sample
        P(i,1) = sum(PermsAct(:,i)<mainAtt_mean(i))/numel(PermsAct(:,i));
        P(i,2) = sum(PermsRep(:,i)<mainRep_mean(i))/numel(PermsRep(:,i));
        P(i,3) = sum(PermsDisp(:,i)<mainDisp_mean(i))/numel(PermsDisp(:,i));
        P(i,4) = sum(PermsActRep(:,i)<intActRep_mean(i))/numel(PermsActRep(:,i));
        P(i,5) = sum(PermsActDisp(:,i)<intActDisp_mean(i))/numel(PermsActDisp(:,i));
        P(i,6) = sum(PermsRepDisp(:,i)<intRepDisp_mean(i))/numel(PermsRepDisp(:,i));
        P(i,7) = sum(PermsActRepDisp(:,i)<intActRepDisp_mean(i))/numel(PermsActRepDisp(:,i));
    elseif strcmp(Tail,'right')
        % same but the other way
        P(i,1) = sum(PermsAct(:,i)>mainAtt_mean(i))/numel(PermsAct(:,i));
        P(i,2) = sum(PermsRep(:,i)>mainRep_mean(i))/numel(PermsRep(:,i));
        P(i,3) = sum(PermsDisp(:,i)>mainDisp_mean(i))/numel(PermsDisp(:,i));
        P(i,4) = sum(PermsActRep(:,i)>intActRep_mean(i))/numel(PermsActRep(:,i));
        P(i,5) = sum(PermsActDisp(:,i)>intActDisp_mean(i))/numel(PermsActDisp(:,i));
        P(i,6) = sum(PermsRepDisp(:,i)>intRepDisp_mean(i))/numel(PermsRepDisp(:,i));
        P(i,7) = sum(PermsActRepDisp(:,i)>intActRepDisp_mean(i))/numel(PermsActRepDisp(:,i));
    elseif strcmp(Tail,'both')
        % for the 2 tailed we compare to the distribution of absolute value of the distance
        % between the result of each permutation to the mean of all permutation results
        % Then we check the proportion of those distances are superior to 
        % the distance between the mean of the sample and the mean of all
        % permutation results
        % P(i) = sum( abs((Perms(:,i)-mean(Perms(:,i)))) > abs((mean(betas(:,i))-mean(Perms(:,i)))) ) / numel(Perms(:,i)) ;
        % Importantly: the above assumes that your null distribution is symmetric
        P(i,1) = sum( abs(PermsAct(:,i)) > abs(mainAct_mean(i) ) )  / numel(PermsAct(:,i)) ;
        P(i,2) = sum( abs(PermsRep(:,i)) > abs(mainRep_mean(i) ) )  / numel(PermsRep(:,i)) ;
        P(i,3) = sum( abs(PermsDisp(:,i)) > abs(mainDisp_mean(i) ) )  / numel(PermsDisp(:,i)) ;
        P(i,4) = sum( abs(PermsActRep(:,i)) > abs(intActRep_mean(i) ) )  / numel(PermsActRep(:,i)) ;
        P(i,5) = sum( abs(PermsActDisp(:,i)) > abs(intActDisp_mean(i) ) )  / numel(PermsActDisp(:,i)) ;
        P(i,6) = sum( abs(PermsRepDisp(:,i)) > abs(intRepDisp_mean(i) ) )  / numel(PermsRepDisp(:,i)) ;
        P(i,7) = sum( abs(PermsActRepDisp(:,i)) > abs(intActRepDisp_mean(i) ) )  / numel(PermsActRepDisp(:,i)) ;
        P(i,8) = sum( abs(PermsRepADisp(:,i)) > abs(intRepADisp_mean(i) ) )  / numel(PermsRepADisp(:,i)) ;
        P(i,9) = sum( abs(PermsRepVDisp(:,i)) > abs(intRepVDisp_mean(i) ) )  / numel(PermsRepVDisp(:,i)) ;
        P(i,10) = sum( abs(PermsDispHRep(:,i)) > abs(intDispHRep_mean(i) ) )  / numel(PermsDispHRep(:,i)) ;
        P(i,11) = sum( abs(PermsDispLRep(:,i)) > abs(intDispLRep_mean(i) ) )  / numel(PermsDispLRep(:,i)) ;
        P(i,12) = sum( abs(PermsActCDisp(:,i)) > abs(intActCDisp_mean(i) ) )  / numel(PermsActCDisp(:,i)) ;
        P(i,13) = sum( abs(PermsActNDisp(:,i)) > abs(intActNDisp_mean(i) ) )  / numel(PermsActNDisp(:,i)) ;
        P(i,14) = sum( abs(PermsDispHAct(:,i)) > abs(intDispHAct_mean(i) ) )  / numel(PermsDispHAct(:,i)) ;
        P(i,15) = sum( abs(PermsDispLAct(:,i)) > abs(intDispLAct_mean(i) ) )  / numel(PermsDispLAct(:,i)) ;
        P(i,16) = sum( abs(PermsActCRep(:,i)) > abs(intActCRep_mean(i) ) )  / numel(PermsActCRep(:,i)) ;
        P(i,17) = sum( abs(PermsActNRep(:,i)) > abs(intActNRep_mean(i) ) )  / numel(PermsActNRep(:,i)) ;
        P(i,18) = sum( abs(PermsRepAAct(:,i)) > abs(intRepAAct_mean(i) ) )  / numel(PermsRepAAct(:,i)) ;
        P(i,19) = sum( abs(PermsRepVAct(:,i)) > abs(intRepVAct_mean(i) ) )  / numel(PermsRepVAct(:,i)) ;

    end
end

% res_list{1}='behav'; % main effects and interactions
% res_list2{1}='behav'; % simple main effects
% 
% for i=1:roi_num
%     % main effects and interactions
%     res_list{i,2}=mainAct_mean(i);
%     res_list{i,3}=P(i,1);
%     res_list{i,4}=mainRep_mean(i);
%     res_list{i,5}=P(i,2);
%     res_list{i,6}=mainDisp_mean(i);
%     res_list{i,7}=P(i,3);
%     res_list{i,8}=intActRep_mean(i);
%     res_list{i,9}=P(i,4);
%     res_list{i,10}=intActDisp_mean(i);
%     res_list{i,11}=P(i,5);
%     res_list{i,12}=intRepDisp_mean(i);
%     res_list{i,13}=P(i,6);
%     res_list{i,14}=intActRepDisp_mean(i);
%     res_list{i,15}=P(i,7);
%     % simple main effects  
%     % res_list2{i,2}=intRepADisp_mean(i);
%     % res_list2{i,3}=P(i,8);
%     % res_list2{i,4}=intRepVDisp_mean(i);
%     % res_list2{i,5}=P(i,9);
%     % res_list2{i,6}=intDispHRep_mean(i);
%     % res_list2{i,7}=P(i,10);
%     % res_list2{i,8}=intDispLRep_mean(i);
%     % res_list2{i,9}=P(i,11);
% end

%% table with results
ptab{1,1}='PermsAct';
ptab{2,1}=P(1);
ptab{1,2}='PermsRep';
ptab{2,2}=P(2);
ptab{1,3}='PermsDisp';
ptab{2,3}=P(3);
ptab{1,4}='PermsActRep';
ptab{2,4}=P(4);
ptab{1,5}='PermsActDisp';
ptab{2,5}=P(5);
ptab{1,6}='PermsRepDisp';
ptab{2,6}=P(6);
ptab{1,7}='PermsActRepDisp';
ptab{2,7}=P(7);
ptab{1,8}='PermsRepADisp';
ptab{2,8}=P(8);
ptab{1,9}='PermsRepVDisp';
ptab{2,9}=P(9);
ptab{1,10}='PermsDispHRep';
ptab{2,10}=P(10);
ptab{1,11}='PermsDispLRep';
ptab{2,11}=P(11);
ptab{1,12}='PermsActCDisp';
ptab{2,12}=P(12);
ptab{1,13}='PermsActNDisp';
ptab{2,13}=P(13);
ptab{1,14}='PermsDispHAct';
ptab{2,14}=P(14);
ptab{1,15}='PermsDispLAct';
ptab{2,15}=P(15);
ptab{1,16}='PermsActCRep';
ptab{2,16}=P(16);
ptab{1,17}='PermsActNRep';
ptab{2,17}=P(17);
ptab{1,18}='PermsRepAAct';
ptab{2,18}=P(18);
ptab{1,19}='PermsRepVAct';
ptab{2,19}=P(19);

filename=['wav_pv_exp',exp,'.xlsx'];
writecell(ptab,filename)


