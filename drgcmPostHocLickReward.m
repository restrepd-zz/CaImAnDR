%% drgcmPostHocLickReward.m
%
% Needs as an input the ouput file from drgCaImAn_dropc
%
close all
clear all

odor_on_dt=4;
break_dt=2;
RA_dt=2;
no_RAs=1;

figNo=0;

% Read .mat file

[choiceFileName,choiceBatchPathName] = uigetfile({'*_batch_per_trial.mat'},'Select the output from drgcmCaImAn_batch_dropc_rhd');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
load(choiceFileName)

%Plot the percent correct as scored by the operant conditioning computer

%Calculate percent correct
sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

if handles_per_trial.no_trials>sliding_window
    
    %Note: Because this is a reversal I am moving the window for calculation of perCorr to the right by nine points
    for ii=1:handles_per_trial.no_trials-sliding_window+1
        no_Hits=sum([handles_per_trial.trial(ii:ii+sliding_window-1).hit]);
        no_CRs=sum([handles_per_trial.trial(ii:ii+sliding_window-1).CR]);
        perCorr(ii+sliding_window-1)=100*(no_Hits+no_CRs)/sliding_window;
    end
    
    
    perCorr(1:sliding_window)=perCorr(sliding_window+1);
    
    %Note, this is here so that perCorr=0 is included in the 0-10 % bin.
    perCorr(perCorr==0)=0.00001;
    
    %Plot percent correct vs trial
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    hold on
    
    jj_low=find(perCorr<max_percent_low_beh);
    plot(jj_low,perCorr(jj_low),'ob')
    hold on
    jj_high=find(perCorr>min_precent_high_beh);
    plot(jj_high,perCorr(jj_high),'or')
    
    jj_mid=find((perCorr<=min_precent_high_beh)&(perCorr>=max_percent_low_beh));
    plot(jj_mid,perCorr(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot([0 handles_per_trial.no_trials],[50 50],'-k')
    
    
    title(['Percent correct vs. trial number ' choiceFileName(18:end-2)])
    xlabel('Trial number')
    ylabel('Percent correct')
    ylim([0 100])
    
end

%Plot the licks for the original events
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

acq_rate=handles_per_trial.acq_rate;
dt_before=handles_per_trial.dt_before;
delta_odor=handles_per_trial.delta_odor;
delta_odor_on_reinf_on=handles_per_trial.delta_odor_on_reinf_on;
delta_reinf=handles_per_trial.delta_reinf;

dHit_lick_traces=[];
Hitii_lick=0;
dCR_lick_traces=[];
CRii_lick=0;
dFA_lick_traces=[];
FAii_lick=0;
dMiss_lick_traces=[];
Missii_lick=0;

for ii=1:handles_per_trial.no_trials
    dall_lick_traces(ii,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    if handles_per_trial.trial(ii).hit==1
        Hitii_lick=Hitii_lick+1;
        dHit_lick_traces(Hitii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if handles_per_trial.trial(ii).miss==1
        Missii_lick=Missii_lick+1;
        dMiss_lick_traces(Missii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if handles_per_trial.trial(ii).CR==1
        CRii_lick=CRii_lick+1;
        dCR_lick_traces(CRii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if handles_per_trial.trial(ii).FA==1
        FAii_lick=FAii_lick+1;
        dFA_lick_traces(FAii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
end

szalllick=size(dall_lick_traces);
time_licks=([1:szalllick(2)]/(acq_rate/20))-dt_before;
per99=prctile(dall_lick_traces(:),99.9);
per1=prctile(dall_lick_traces(:),1);

mean_Hit_licks=zeros(1,szalllick(2));
mean_Miss_licks=zeros(1,szalllick(2));
mean_FA_licks=zeros(1,szalllick(2));
mean_CR_licks=zeros(1,szalllick(2));

y_shift=0;

%Plot CR lick traces
for ii=1:CRii_lick
    plot(time_licks,dCR_lick_traces(ii,:)+y_shift,'-b')
    mean_CR_licks(1,:)= mean_CR_licks(1,:) + ((dCR_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Plot FA lick traces
for ii=1:FAii_lick
    plot(time_licks,dFA_lick_traces(ii,:)+y_shift,'-m')
    mean_FA_licks(1,:)= mean_FA_licks(1,:) + ((dFA_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Plot Miss lick traces
for ii=1:Missii_lick
    plot(time_licks,dMiss_lick_traces(ii,:)+y_shift,'-c')
    mean_Miss_licks(1,:)= mean_Miss_licks(1,:) + ((dMiss_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%PLot Hit lick traces
for ii=1:Hitii_lick
    plot(time_licks,dHit_lick_traces(ii,:)+y_shift,'-r')
    mean_Hit_licks(1,:)= mean_Hit_licks(1,:) + ((dHit_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Odor on markers
plot([0 0],[0 y_shift],'-k')
plot([mean(delta_odor) mean(delta_odor)],[0 y_shift],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 y_shift],'-r')
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 y_shift],'-r')

title('Licks for events scored by the olfactometer')

%Re-score the licks
new_hits=[];
new_CRs=[];
new_Miss=[];
new_FAs=[];

for trNo=1:handles_per_trial.no_trials
    new_Hits(trNo)=0;
    new_CRs(trNo)=0;
    new_Miss(trNo)=0;
    new_FAs(trNo)=0;
    if (handles_per_trial.trial(trNo).hit==1)||(handles_per_trial.trial(trNo).miss==1)
        RAs_with_licks=0;
        for iiRA=1:no_RAs
            if sum((handles_per_trial.trial(trNo).these_lick_times>break_dt+(iiRA-1)*RA_dt)&(handles_per_trial.trial(trNo).these_lick_times<=break_dt+iiRA*RA_dt))
                RAs_with_licks=RAs_with_licks+1;
            end
        end
        if RAs_with_licks==no_RAs
            new_Hits(trNo)=1;
        else
            new_Miss(trNo)=1;
        end
    else
        RAs_with_licks=0;
        for iiRA=1:no_RAs
            if sum((handles_per_trial.trial(trNo).these_lick_times>break_dt+(iiRA-1)*RA_dt)&(handles_per_trial.trial(trNo).these_lick_times<=break_dt+iiRA*RA_dt))
                RAs_with_licks=RAs_with_licks+1;
            end
        end
        if RAs_with_licks~=no_RAs
            new_CRs(trNo)=1;
        else
            new_FAs(trNo)=1;
        end
    end
end

%Plot re-scored percent correct
perCorrnew=[];

if handles_per_trial.no_trials>sliding_window
    %Note: Because this is a reversal I am moving the window for calculation of perCorrnew to the right by nine points
    for ii=1:handles_per_trial.no_trials-sliding_window+1
        no_Hits=sum([new_Hits(ii:ii+sliding_window-1)]);
        no_CRs=sum([handles_per_trial.trial(ii:ii+sliding_window-1).CR]);
        perCorrnew(ii+sliding_window-1)=100*(no_Hits+no_CRs)/sliding_window;
    end
    
    
    perCorrnew(1:sliding_window)=perCorrnew(sliding_window+1);
    
    %Note, this is here so that perCorrnew=0 is included in the 0-10 % bin.
    perCorrnew(perCorrnew==0)=0.00001;
    
    %Plot percent correct vs trial
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    hold on
    
    jj_low=find(perCorrnew<max_percent_low_beh);
    plot(jj_low,perCorrnew(jj_low),'ob')
    hold on
    jj_high=find(perCorrnew>min_precent_high_beh);
    plot(jj_high,perCorrnew(jj_high),'or')
    
    jj_mid=find((perCorrnew<=min_precent_high_beh)&(perCorrnew>=max_percent_low_beh));
    plot(jj_mid,perCorrnew(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
    hold on
    plot([0 handles_per_trial.no_trials],[50 50],'-k')
    
    
    title(['Re-calculated percent correct vs. trial number ' choiceFileName(18:end-2)])
    xlabel('Trial number')
    ylabel('Percent correct')
    ylim([0 100])
    
end

%Plot re-scored licks
dHit_lick_traces=[];
Hitii_lick=0;
dCR_lick_traces=[];
CRii_lick=0;
dFA_lick_traces=[];
FAii_lick=0;
dMiss_lick_traces=[];
Missii_lick=0;

for ii=1:handles_per_trial.no_trials
    dall_lick_traces(ii,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    if new_Hits(ii)==1
        Hitii_lick=Hitii_lick+1;
        dHit_lick_traces(Hitii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if new_Miss(ii)==1
        Missii_lick=Missii_lick+1;
        dMiss_lick_traces(Missii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if new_CRs(ii)==1
        CRii_lick=CRii_lick+1;
        dCR_lick_traces(CRii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
    if new_FAs(ii)==1
        FAii_lick=FAii_lick+1;
        dFA_lick_traces(FAii_lick,:)=decimate(handles_per_trial.trial(ii).lick_trace,20)';
    end
end

figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
hold on

szalllick=size(dall_lick_traces);
time_licks=([1:szalllick(2)]/(acq_rate/20))-dt_before;
per99=prctile(dall_lick_traces(:),99.9);
per1=prctile(dall_lick_traces(:),1);

mean_Hit_licks=zeros(1,szalllick(2));
mean_Miss_licks=zeros(1,szalllick(2));
mean_FA_licks=zeros(1,szalllick(2));
mean_CR_licks=zeros(1,szalllick(2));

y_shift=0;

%Plot CR lick traces
for ii=1:CRii_lick
    plot(time_licks,dCR_lick_traces(ii,:)+y_shift,'-b')
    mean_CR_licks(1,:)= mean_CR_licks(1,:) + ((dCR_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Plot FA lick traces
for ii=1:FAii_lick
    plot(time_licks,dFA_lick_traces(ii,:)+y_shift,'-m')
    mean_FA_licks(1,:)= mean_FA_licks(1,:) + ((dFA_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Plot Miss lick traces
for ii=1:Missii_lick
    plot(time_licks,dMiss_lick_traces(ii,:)+y_shift,'-c')
    mean_Miss_licks(1,:)= mean_Miss_licks(1,:) + ((dMiss_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Plot Hit lick traces
for ii=1:Hitii_lick
    plot(time_licks,dHit_lick_traces(ii,:)+y_shift,'-r')
    mean_Hit_licks(1,:)= mean_Hit_licks(1,:) + ((dHit_lick_traces(ii,:)-per1)/(per99-per1));
    y_shift=y_shift+1.2*(per99-per1);
end

%Odor on markers
plot([0 0],[0 y_shift],'-k')
plot([mean(delta_odor) mean(delta_odor)],[0 y_shift],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 y_shift],'-r')
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 y_shift],'-r')

title('Licks for events scored post hoc')

pffft=1

