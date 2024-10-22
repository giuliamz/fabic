clear;
close all;
clc;

exp = '2'; %exp 1/2
data_path = [pwd, '\data\experiment', exp];

if strcmp(exp,'1')
    %list of participants for Experiment1
    subjID={'05';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';...
       '38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';...
       '51';'52';'53';'54';'55';'56';'57';'58';'59'};
elseif strcmp(exp,'2')
    %list of participants for Experiment2
    subjID={'70';'71';'72';'73';'75';'76';'77';'78';'79';'80';'81';'83';...
        '84';'85';'86';'87';'88';'89';'93';'94';'95';'98';'99';'101';...
        '102';'103';'104';'106';'107';'111';'112';'113';'114';'115'}; 
end

nsubj=length(subjID);

%% compute reported positions in congruent trials pooled over all action 
%% and res modality conditions across all subjects

for s=1:nsubj

    %load log files
    Tablename = [data_path, '\Log_s', subjID{s}, '.mat'];  
    load(Tablename);

    Rtab=table2cell(LogFile);
    S.TrialTot = size(Rtab,1);
    S.order.vis_loc = cell2mat(Rtab(:,1)); %vis stim position (-9°,0°,9°)
    S.order.aud_loc = cell2mat(Rtab(:,2)); %aud stim position (-9°,0°,9°)
    S.order.response = cell2mat(Rtab(:,3)); %response modality: 1=aud, 2=vis
    S.order.comm = cell2mat(Rtab(:,4)); %action intention: 1=com, 2=ncom
    S.response = cell2mat(Rtab(:,5)); %responses (-9°,0°,9°)
    S.rt = cell2mat(Rtab(:,6)); %reaction times

    for i=1:S.TrialTot
    if S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==-9 
        cong_left(i)=S.response(i);
    else
        cong_left(i)=NaN;
    end

    if S.order.aud_loc(i)==0 && S.order.vis_loc(i)==0 
        cong_centre(i)=S.response(i);
    else
        cong_centre(i)=NaN;
    end

    if S.order.aud_loc(i)==9 && S.order.vis_loc(i)==9 
        cong_right(i)=S.response(i);
    else
        cong_right(i)=NaN;
    end


    end

    mean_cong_left(s)=mean(cong_left(~isnan(cong_left)));
    mean_cong_centre(s)=mean(cong_centre(~isnan(cong_centre)));
    mean_cong_right(s)=mean(cong_right(~isnan(cong_right)));

end


mean_cong_left_tot=mean(mean_cong_left); %% E1:-8.589 E2:-8.668 
mean_cong_centre_tot=mean(mean_cong_centre); %% E1:-0.667 E2:-0.639 
mean_cong_right_tot=mean(mean_cong_right); %% E1:8.208 E2:7.849

%% compute wav
mean_cong_left_tot= -8.668; 
mean_cong_centre_tot= -0.639;
mean_cong_right_tot= 7.849;

for s=1:nsubj

    %load log files
    Tablename = [data_path, '\Log_s', subjID{s}, '.mat'];  
    load(Tablename);

    Rtab=table2cell(LogFile);
    Rtab=table2cell(LogFile);
    S.TrialTot = size(Rtab,1);
    S.order.vis_loc = cell2mat(Rtab(:,1)); %vis stim position (-9°,0°,9°)
    S.order.aud_loc = cell2mat(Rtab(:,2)); %aud stim position (-9°,0°,9°)
    S.order.response = cell2mat(Rtab(:,3)); %response modality: 1=aud, 2=vis
    S.order.comm = cell2mat(Rtab(:,4)); %action intention: 1=com, 2=ncom
    S.response = cell2mat(Rtab(:,5)); %responses (-9°,0°,9°)
    S.rt = cell2mat(Rtab(:,6)); %reaction times

    for i=1:S.TrialTot

        %wav pos aud -9 (L) vis 0 (C)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==0 
            WAV_AL_VC_COM_AUD(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_centre_tot-mean_cong_left_tot);
        else
            WAV_AL_VC_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==0 
            WAV_AL_VC_COM_VIS(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_left_tot-mean_cong_centre_tot);
        else
            WAV_AL_VC_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==0 
            WAV_AL_VC_NC_AUD(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_centre_tot-mean_cong_left_tot);
        else
            WAV_AL_VC_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==0 
            WAV_AL_VC_NC_VIS(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_left_tot-mean_cong_centre_tot);
        else
            WAV_AL_VC_NC_VIS(i)=NaN;
        end

        %wav pos aud -9 (L) vid 9 (R)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==9
            WAV_AL_VR_COM_AUD(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_right_tot-mean_cong_left_tot);
        else
            WAV_AL_VR_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==9
            WAV_AL_VR_COM_VIS(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_left_tot-mean_cong_right_tot);
        else
            WAV_AL_VR_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==9 
            WAV_AL_VR_NC_AUD(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_right_tot-mean_cong_left_tot);
        else
            WAV_AL_VR_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==-9 && S.order.vis_loc(i)==9 
            WAV_AL_VR_NC_VIS(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_left_tot-mean_cong_right_tot);
        else
            WAV_AL_VR_NC_VIS(i)=NaN;
        end

        %wav pos aud 0 (C) vid -9 (L)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==-9 
            WAV_AC_VL_COM_AUD(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_left_tot-mean_cong_centre_tot);
        else
            WAV_AC_VL_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==-9 
            WAV_AC_VL_COM_VIS(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_centre_tot-mean_cong_left_tot);
        else
            WAV_AC_VL_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==-9 && ...
                ~isempty(S.response(i)) && any(S.rt(i)>100) && any(S.rt(i)<1500) 
            WAV_AC_VL_NC_AUD(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_left_tot-mean_cong_centre_tot);
        else
            WAV_AC_VL_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==-9 
            WAV_AC_VL_NC_VIS(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_centre_tot-mean_cong_left_tot);
        else
            WAV_AC_VL_NC_VIS(i)=NaN;
        end

        %wav pos aud 0 (C) vid 9 (R)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==9 
            WAV_AC_VR_COM_AUD(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_right_tot-mean_cong_centre_tot);
        else
            WAV_AC_VR_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==9 
            WAV_AC_VR_COM_VIS(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_centre_tot-mean_cong_right_tot);
        else
            WAV_AC_VR_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==9 
            WAV_AC_VR_NC_AUD(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_right_tot-mean_cong_centre_tot);
        else
            WAV_AC_VR_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==0 && S.order.vis_loc(i)==9 
            WAV_AC_VR_NC_VIS(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_centre_tot-mean_cong_right_tot);
        else
            WAV_AC_VR_NC_VIS(i)=NaN;
        end

        %wav pos aud 9 (R) vid -9 (L)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==-9 
            WAV_AR_VL_COM_AUD(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_left_tot-mean_cong_right_tot);
        else
            WAV_AR_VL_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==-9 
            WAV_AR_VL_COM_VIS(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_right_tot-mean_cong_left_tot);
        else
            WAV_AR_VL_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==-9 
            WAV_AR_VL_NC_AUD(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_left_tot-mean_cong_right_tot);
        else
            WAV_AR_VL_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==-9 
            WAV_AR_VL_NC_VIS(i)=(S.response(i)-mean_cong_left_tot)/(mean_cong_right_tot-mean_cong_left_tot);
        else
            WAV_AR_VL_NC_VIS(i)=NaN;
        end

        %wav pos aud 9 (R) vid 0 (C)
        if S.order.comm(i)==1 && S.order.response(i)==1 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==0 
            WAV_AR_VC_COM_AUD(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_centre_tot-mean_cong_right_tot);
        else
            WAV_AR_VC_COM_AUD(i)=NaN;
        end
        if S.order.comm(i)==1 && S.order.response(i)==2 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==0 
            WAV_AR_VC_COM_VIS(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_right_tot-mean_cong_centre_tot);
        else
            WAV_AR_VC_COM_VIS(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==1 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==0 
            WAV_AR_VC_NC_AUD(i)=(S.response(i)-mean_cong_right_tot)/(mean_cong_centre_tot-mean_cong_right_tot);
        else
            WAV_AR_VC_NC_AUD(i)=NaN;
        end
        if S.order.comm(i)==2 && S.order.response(i)==2 && S.order.aud_loc(i)==9 && S.order.vis_loc(i)==0 
            WAV_AR_VC_NC_VIS(i)=(S.response(i)-mean_cong_centre_tot)/(mean_cong_right_tot-mean_cong_centre_tot);
        else
            WAV_AR_VC_NC_VIS(i)=NaN;
        end


    end

%meav wav for each condition for each subject
mean_WAV_AL_VC_COM_AUD(s)=mean(rmmissing(WAV_AL_VC_COM_AUD));
mean_WAV_AL_VC_COM_VIS(s)=mean(rmmissing(WAV_AL_VC_COM_VIS));
mean_WAV_AL_VC_NC_AUD(s)=mean(rmmissing(WAV_AL_VC_NC_AUD));
mean_WAV_AL_VC_NC_VIS(s)=mean(rmmissing(WAV_AL_VC_NC_VIS));

mean_WAV_AL_VR_COM_AUD(s)=mean(rmmissing(WAV_AL_VR_COM_AUD));
mean_WAV_AL_VR_COM_VIS(s)=mean(rmmissing(WAV_AL_VR_COM_VIS));
mean_WAV_AL_VR_NC_AUD(s)=mean(rmmissing(WAV_AL_VR_NC_AUD));
mean_WAV_AL_VR_NC_VIS(s)=mean(rmmissing(WAV_AL_VR_NC_VIS));

mean_WAV_AC_VL_COM_AUD(s)=mean(rmmissing(WAV_AC_VL_COM_AUD));
mean_WAV_AC_VL_COM_VIS(s)=mean(rmmissing(WAV_AC_VL_COM_VIS));
mean_WAV_AC_VL_NC_AUD(s)=mean(rmmissing(WAV_AC_VL_NC_AUD));
mean_WAV_AC_VL_NC_VIS(s)=mean(rmmissing(WAV_AC_VL_NC_VIS));

mean_WAV_AC_VR_COM_AUD(s)=mean(rmmissing(WAV_AC_VR_COM_AUD));
mean_WAV_AC_VR_COM_VIS(s)=mean(rmmissing(WAV_AC_VR_COM_VIS));
mean_WAV_AC_VR_NC_AUD(s)=mean(rmmissing(WAV_AC_VR_NC_AUD));
mean_WAV_AC_VR_NC_VIS(s)=mean(rmmissing(WAV_AC_VR_NC_VIS));

mean_WAV_AR_VL_COM_AUD(s)=mean(rmmissing(WAV_AR_VL_COM_AUD));
mean_WAV_AR_VL_COM_VIS(s)=mean(rmmissing(WAV_AR_VL_COM_VIS));
mean_WAV_AR_VL_NC_AUD(s)=mean(rmmissing(WAV_AR_VL_NC_AUD));
mean_WAV_AR_VL_NC_VIS(s)=mean(rmmissing(WAV_AR_VL_NC_VIS));

mean_WAV_AR_VC_COM_AUD(s)=mean(rmmissing(WAV_AR_VC_COM_AUD));
mean_WAV_AR_VC_COM_VIS(s)=mean(rmmissing(WAV_AR_VC_COM_VIS));
mean_WAV_AR_VC_NC_AUD(s)=mean(rmmissing(WAV_AR_VC_NC_AUD));
mean_WAV_AR_VC_NC_VIS(s)=mean(rmmissing(WAV_AR_VC_NC_VIS));


end

%across subjects mean
meantot_WAV_AL_VC_COM_AUD=mean(mean_WAV_AL_VC_COM_AUD);
meantot_WAV_AL_VC_COM_VIS=mean(mean_WAV_AL_VC_COM_VIS);
meantot_WAV_AL_VC_NC_AUD=mean(mean_WAV_AL_VC_NC_AUD);
meantot_WAV_AL_VC_NC_VIS=mean(mean_WAV_AL_VC_NC_VIS);

meantot_WAV_AL_VR_COM_AUD=mean(mean_WAV_AL_VR_COM_AUD);
meantot_WAV_AL_VR_COM_VIS=mean(mean_WAV_AL_VR_COM_VIS);
meantot_WAV_AL_VR_NC_AUD=mean(mean_WAV_AL_VR_NC_AUD);
meantot_WAV_AL_VR_NC_VIS=mean(mean_WAV_AL_VR_NC_VIS);

meantot_WAV_AC_VL_COM_AUD=mean(mean_WAV_AC_VL_COM_AUD);
meantot_WAV_AC_VL_COM_VIS=mean(mean_WAV_AC_VL_COM_VIS);
meantot_WAV_AC_VL_NC_AUD=mean(mean_WAV_AC_VL_NC_AUD);
meantot_WAV_AC_VL_NC_VIS=mean(mean_WAV_AC_VL_NC_VIS);

meantot_WAV_AC_VR_COM_AUD=mean(mean_WAV_AC_VR_COM_AUD);
meantot_WAV_AC_VR_COM_VIS=mean(mean_WAV_AC_VR_COM_VIS);
meantot_WAV_AC_VR_NC_AUD=mean(mean_WAV_AC_VR_NC_AUD);
meantot_WAV_AC_VR_NC_VIS=mean(mean_WAV_AC_VR_NC_VIS);

meantot_WAV_AR_VL_COM_AUD=mean(mean_WAV_AR_VL_COM_AUD);
meantot_WAV_AR_VL_COM_VIS=mean(mean_WAV_AR_VL_COM_VIS);
meantot_WAV_AR_VL_NC_AUD=mean(mean_WAV_AR_VL_NC_AUD);
meantot_WAV_AR_VL_NC_VIS=mean(mean_WAV_AR_VL_NC_VIS);

meantot_WAV_AR_VC_COM_AUD=mean(mean_WAV_AR_VC_COM_AUD);
meantot_WAV_AR_VC_COM_VIS=mean(mean_WAV_AR_VC_COM_VIS);
meantot_WAV_AR_VC_NC_AUD=mean(mean_WAV_AR_VC_NC_AUD);
meantot_WAV_AR_VC_NC_VIS=mean(mean_WAV_AR_VC_NC_VIS);

% st dev & sem
stdtot_WAV_AL_VC_COM_AUD=std(mean_WAV_AL_VC_COM_AUD);
semtot_WAV_AL_VC_COM_AUD=std(mean_WAV_AL_VC_COM_AUD)/sqrt(size(mean_WAV_AL_VC_COM_AUD,2));
stdtot_WAV_AL_VC_COM_VIS=std(mean_WAV_AL_VC_COM_VIS);
semtot_WAV_AL_VC_COM_VIS=std(mean_WAV_AL_VC_COM_VIS)/sqrt(size(mean_WAV_AL_VC_COM_VIS,2));
stdtot_WAV_AL_VC_NC_AUD=std(mean_WAV_AL_VC_NC_AUD);
semtot_WAV_AL_VC_NC_AUD=std(mean_WAV_AL_VC_NC_AUD)/sqrt(size(mean_WAV_AL_VC_NC_AUD,2));
stdtot_WAV_AL_VC_NC_VIS=mean(mean_WAV_AL_VC_NC_VIS);
semtot_WAV_AL_VC_NC_VIS=mean(mean_WAV_AL_VC_NC_VIS)/sqrt(size(mean_WAV_AL_VC_NC_VIS,2));


stdtot_WAV_AL_VR_COM_AUD=std(mean_WAV_AL_VR_COM_AUD);
semtot_WAV_AL_VR_COM_AUD=std(mean_WAV_AL_VR_COM_AUD)/sqrt(size(mean_WAV_AL_VR_COM_AUD,2));
stdtot_WAV_AL_VR_COM_VIS=mean(mean_WAV_AL_VR_COM_VIS);
semtot_WAV_AL_VR_COM_VIS=mean(mean_WAV_AL_VR_COM_VIS)/sqrt(size(mean_WAV_AL_VR_COM_VIS,2));
stdtot_WAV_AL_VR_NC_AUD=std(mean_WAV_AL_VR_NC_AUD);
semtot_WAV_AL_VR_NC_AUD=std(mean_WAV_AL_VR_NC_AUD)/sqrt(size(mean_WAV_AL_VR_NC_AUD,2));
stdtot_WAV_AL_VR_NC_VIS=mean(mean_WAV_AL_VR_NC_VIS);
semtot_WAV_AL_VR_NC_VIS=mean(mean_WAV_AL_VR_NC_VIS)/sqrt(size(mean_WAV_AL_VR_NC_VIS,2));


stdtot_WAV_AC_VL_COM_AUD=std(mean_WAV_AC_VL_COM_AUD);
semtot_WAV_AC_VL_COM_AUD=std(mean_WAV_AC_VL_COM_AUD)/sqrt(size(mean_WAV_AC_VL_COM_AUD,2));
stdtot_WAV_AC_VL_COM_VIS=mean(mean_WAV_AC_VL_COM_VIS);
semtot_WAV_AC_VL_COM_VIS=mean(mean_WAV_AC_VL_COM_VIS)/sqrt(size(mean_WAV_AC_VL_COM_VIS,2));
stdtot_WAV_AC_VL_NC_AUD=std(mean_WAV_AC_VL_NC_AUD);
semtot_WAV_AC_VL_NC_AUD=std(mean_WAV_AC_VL_NC_AUD)/sqrt(size(mean_WAV_AC_VL_NC_AUD,2));
stdtot_WAV_AC_VL_NC_VIS=mean(mean_WAV_AC_VL_NC_VIS);
semtot_WAV_AC_VL_NC_VIS=mean(mean_WAV_AC_VL_NC_VIS)/sqrt(size(mean_WAV_AC_VL_NC_VIS,2));


stdtot_WAV_AC_VR_COM_AUD=std(mean_WAV_AC_VR_COM_AUD);
semtot_WAV_AC_VR_COM_AUD=std(mean_WAV_AC_VR_COM_AUD)/sqrt(size(mean_WAV_AC_VR_COM_AUD,2));
stdtot_WAV_AC_VR_COM_VIS=mean(mean_WAV_AC_VR_COM_VIS);
semtot_WAV_AC_VR_COM_VIS=mean(mean_WAV_AC_VR_COM_VIS)/sqrt(size(mean_WAV_AC_VR_COM_VIS,2));
stdtot_WAV_AC_VR_NC_AUD=std(mean_WAV_AC_VR_NC_AUD);
semtot_WAV_AC_VR_NC_AUD=std(mean_WAV_AC_VR_NC_AUD)/sqrt(size(mean_WAV_AC_VR_NC_AUD,2));
stdtot_WAV_AC_VR_NC_VIS=mean(mean_WAV_AC_VR_NC_VIS);
semtot_WAV_AC_VR_NC_VIS=mean(mean_WAV_AC_VR_NC_VIS)/sqrt(size(mean_WAV_AC_VR_NC_VIS,2));


stdtot_WAV_AR_VL_COM_AUD=std(mean_WAV_AR_VL_COM_AUD);
semtot_WAV_AR_VL_COM_AUD=std(mean_WAV_AR_VL_COM_AUD)/sqrt(size(mean_WAV_AR_VL_COM_AUD,2));
stdtot_WAV_AR_VL_COM_VIS=mean(mean_WAV_AR_VL_COM_VIS);
semtot_WAV_AR_VL_COM_VIS=mean(mean_WAV_AR_VL_COM_VIS)/sqrt(size(mean_WAV_AR_VL_COM_VIS,2));
stdtot_WAV_AR_VL_NC_AUD=std(mean_WAV_AR_VL_NC_AUD);
semtot_WAV_AR_VL_NC_AUD=std(mean_WAV_AR_VL_NC_AUD)/sqrt(size(mean_WAV_AR_VL_NC_AUD,2));
stdtot_WAV_AR_VL_NC_VIS=mean(mean_WAV_AR_VL_NC_VIS);
semtot_WAV_AR_VL_NC_VIS=mean(mean_WAV_AR_VL_NC_VIS)/sqrt(size(mean_WAV_AR_VL_NC_VIS,2));


stdtot_WAV_AR_VC_COM_AUD=std(mean_WAV_AR_VC_COM_AUD);
semtot_WAV_AR_VC_COM_AUD=std(mean_WAV_AR_VC_COM_AUD)/sqrt(size(mean_WAV_AR_VC_COM_AUD,2));
stdtot_WAV_AR_VC_COM_VIS=mean(mean_WAV_AR_VC_COM_VIS);
semtot_WAV_AR_VC_COM_VIS=mean(mean_WAV_AR_VC_COM_VIS)/sqrt(size(mean_WAV_AR_VC_COM_VIS,2));
stdtot_WAV_AR_VC_NC_AUD=std(mean_WAV_AR_VC_NC_AUD);
semtot_WAV_AR_VC_NC_AUD=std(mean_WAV_AR_VC_NC_AUD)/sqrt(size(mean_WAV_AR_VC_NC_AUD,2));
stdtot_WAV_AR_VC_NC_VIS=mean(mean_WAV_AR_VC_NC_VIS);
semtot_WAV_AR_VC_NC_VIS=mean(mean_WAV_AR_VC_NC_VIS)/sqrt(size(mean_WAV_AR_VC_NC_VIS,2));

%pool tegether high and low disparity conditions
wav_low_com_aud = (meantot_WAV_AC_VL_COM_AUD + meantot_WAV_AC_VR_COM_AUD + meantot_WAV_AL_VC_COM_AUD + meantot_WAV_AR_VC_COM_AUD)/4;
sem_low_com_aud = (semtot_WAV_AC_VL_COM_AUD + semtot_WAV_AC_VR_COM_AUD + semtot_WAV_AL_VC_COM_AUD + semtot_WAV_AR_VC_COM_AUD)/4;
wav_low_nc_aud = (meantot_WAV_AC_VL_NC_AUD + meantot_WAV_AC_VR_NC_AUD + meantot_WAV_AL_VC_NC_AUD + meantot_WAV_AR_VC_NC_AUD)/4;
sem_low_nc_aud = (semtot_WAV_AC_VL_NC_AUD + semtot_WAV_AC_VR_NC_AUD + semtot_WAV_AL_VC_NC_AUD + semtot_WAV_AR_VC_NC_AUD)/4;
wav_high_com_aud = (meantot_WAV_AL_VR_COM_AUD + meantot_WAV_AR_VL_COM_AUD)/2;
sem_high_com_aud = (semtot_WAV_AL_VR_COM_AUD + semtot_WAV_AR_VL_COM_AUD)/2;
wav_high_nc_aud = (meantot_WAV_AL_VR_NC_AUD + meantot_WAV_AR_VL_NC_AUD)/2;
sem_high_nc_aud = (semtot_WAV_AL_VR_NC_AUD + semtot_WAV_AR_VL_NC_AUD)/2;
wav_low_com_vis = 1 - ((meantot_WAV_AC_VL_COM_VIS + meantot_WAV_AC_VR_COM_VIS + meantot_WAV_AL_VC_COM_VIS + meantot_WAV_AR_VC_COM_VIS)/4);
sem_low_com_vis = (semtot_WAV_AC_VL_COM_VIS + semtot_WAV_AC_VR_COM_VIS + semtot_WAV_AL_VC_COM_VIS + semtot_WAV_AR_VC_COM_VIS)/4;
wav_low_nc_vis = 1 - ((meantot_WAV_AC_VL_NC_VIS + meantot_WAV_AC_VR_NC_VIS + meantot_WAV_AL_VC_NC_VIS + meantot_WAV_AR_VC_NC_VIS)/4);
sem_low_nc_vis = (semtot_WAV_AC_VL_NC_VIS + semtot_WAV_AC_VR_NC_VIS + semtot_WAV_AL_VC_NC_VIS + semtot_WAV_AR_VC_NC_VIS)/4;
wav_high_com_vis = 1 - ((meantot_WAV_AL_VR_COM_VIS + meantot_WAV_AR_VL_COM_VIS)/2);
sem_high_com_vis = (semtot_WAV_AL_VR_COM_VIS + semtot_WAV_AR_VL_COM_VIS)/2;
wav_high_nc_vis = 1 - ((meantot_WAV_AL_VR_NC_VIS + meantot_WAV_AR_VL_NC_VIS)/2);
sem_high_nc_vis = (semtot_WAV_AL_VR_NC_VIS + semtot_WAV_AR_VL_NC_VIS)/2;

%table with average and sem for each condition
Res(:,1)={'wav_low_com_aud';'sem_low_com_aud';'wav_low_nc_aud';'sem_low_nc_aud';'wav_high_com_aud';'sem_high_com_aud';...
    'wav_high_nc_aud';'sem_high_nc_aud';'wav_low_com_vis';'sem_low_com_vis';'wav_low_nc_vis';'sem_low_nc_vis';...
    'wav_high_com_vis';'sem_high_com_vis';'wav_high_nc_vis';'sem_high_nc_vis'};
Res(:,2)={wav_low_com_aud;sem_low_com_aud;wav_low_nc_aud;sem_low_nc_aud';wav_high_com_aud;sem_high_com_aud;...
    wav_high_nc_aud;sem_high_nc_aud;wav_low_com_vis;sem_low_com_vis;wav_low_nc_vis;sem_low_nc_vis;...
    wav_high_com_vis;sem_high_com_vis;wav_high_nc_vis;sem_high_nc_vis};
filename=['fabic_wav_exp',exp];
save(filename, 'Res') 

%% statistical analysis
% create n x m mat, n=nsubj, m=ncond (2 resp mod x 2 AV disp x 2 action)

mean_low_com_aud = (mean_WAV_AC_VL_COM_AUD + mean_WAV_AC_VR_COM_AUD + mean_WAV_AL_VC_COM_AUD + mean_WAV_AR_VC_COM_AUD)/4;
mean_high_com_aud = (mean_WAV_AL_VR_COM_AUD + mean_WAV_AR_VL_COM_AUD)/2;
mean_low_nc_aud = (mean_WAV_AC_VL_NC_AUD + mean_WAV_AC_VR_NC_AUD + mean_WAV_AL_VC_NC_AUD + mean_WAV_AR_VC_NC_AUD)/4;
mean_high_nc_aud = (mean_WAV_AL_VR_NC_AUD + mean_WAV_AR_VL_NC_AUD)/2;
mean_low_com_vis = 1-((mean_WAV_AC_VL_COM_VIS + mean_WAV_AC_VR_COM_VIS + mean_WAV_AL_VC_COM_VIS + mean_WAV_AR_VC_COM_VIS)/4);
mean_high_com_vis = 1-((mean_WAV_AL_VR_COM_VIS + mean_WAV_AR_VL_COM_VIS)/2);
mean_low_nc_vis = 1-((mean_WAV_AC_VL_NC_VIS + mean_WAV_AC_VR_NC_VIS + mean_WAV_AL_VC_NC_VIS + mean_WAV_AR_VC_NC_VIS)/4);
mean_high_nc_vis = 1-((mean_WAV_AL_VR_NC_VIS + mean_WAV_AR_VL_NC_VIS)/2);

mat=[mean_low_com_aud',mean_low_nc_aud',mean_high_com_aud',mean_high_nc_aud',mean_low_com_vis',mean_low_nc_vis',mean_high_com_vis',mean_high_nc_vis'];

mat=array2table(mat);

mat.Properties.VariableNames{1} = 'mean_low_com_aud';
mat.Properties.VariableNames{2} = 'mean_low_nc_aud';
mat.Properties.VariableNames{3} = 'mean_high_com_aud';
mat.Properties.VariableNames{4} = 'mean_high_nc_aud';
mat.Properties.VariableNames{5} = 'mean_low_com_vis';
mat.Properties.VariableNames{6} = 'mean_low_nc_vis';
mat.Properties.VariableNames{7} = 'mean_high_com_vis';
mat.Properties.VariableNames{8} = 'mean_high_nc_vis';

filename=['fabic_subj_wav_exp',exp];
save(filename, 'mat')
