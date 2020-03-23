function drgCaImAnBatchPerSessionEventsPerTriallickvsdFF(choiceBatchPathName,choiceFileName)

%This code does a linear discriminant analysis for spm data
% Needs a choices file such as drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum
% Needs the output files from drgCaImAn_batch_dropc.m
warning('off')

close all
clear all

handles_outs.lick_slopes=[];
handles_outs.dFF_slopes=[];
handles_outs.no_lick_slopes=0;
handles_outs.no_dFF_slopes=0;
handles_outs.slope_triggered_LR=[];

min_trials=20;

tic

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAnBatchPerSessionReversalPerTrial run for ' choiceFileName '\n\n']);



addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;

all_lda_events=[];
all_lda_input_timecourse=[];
no_timepoints=2000000;
all_lda_events_miss_FA=[];

lick_times=[];
no_licks=[];
dLickTraces=[];

for filNum=1:caimanhandles.caimandr_choices.no_files
    
    %Read the file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
    
    
    first_num_odor_trials(filNum)=num_odor_trials+1;
    
    
    for trNo=1:no_odor_trials
        
        %Save epoch
        num_odor_trials=num_odor_trials+1;
        
        %Save lda
        all_lda_events{num_odor_trials}=lda_event{trNo};
        szlit=size(lda_input_timecourse);
        all_lda_input_timecourse(1:length(time_to_eventLDA),1:szlit(2),num_odor_trials)=lda_input_timecourse(:,:,trNo);
        all_lda_no_comp(num_odor_trials)=szlit(2);
        all_lda_fileNo(num_odor_trials)=filNum;
        no_timepoints=min([no_timepoints length(time_to_eventLDA)]);
        
        if epoch_per_trial(trNo)==6
            %Hit
            epochs_per_trial(1,num_odor_trials)=1;
            epochs_per_trial(2:4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=1;
            
            %Was dF/F calculated?
            if sum(which_trial_Hit==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Hit_traces(which_trial_Hit==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save the licks for this trial
                this_Hitii_lick=which_Hitii_lick(find(which_trial_Hit==trNo,1));
                these_Hitii_lick_times=[];
                these_Hitii_lick_times=Hit_lick_times(this_Hitii_lick,1:Hit_no_lick_times(this_Hitii_lick));
                if ~isempty(these_Hitii_lick_times)
                    lick_times(num_odor_trials,1:length(these_Hitii_lick_times))=these_Hitii_lick_times;
                    no_licks(num_odor_trials)=length(these_Hitii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
                dLickTraces(num_odor_trials,:)=dHit_lick_traces(this_Hitii_lick,:);
                
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Hit_traces(which_trial_Hit==trNo,win),2);
                    szhit=size(Hit_traces(which_trial_Hit==trNo,win));
                    %                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Hitii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Hitii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=1;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
                %Calculate the average snip for this trial
                Hitii=handles_out.Hit_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(Hitii).hit_traces);
                num_traces=handles_out.trialNo(Hitii).trace_numHit;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(Hitii).trace_numHit
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(Hitii).hit_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventHit;
                
                
                
            end
        end
        
        if epoch_per_trial(trNo)==7
            %Miss
            epochs_per_trial(2,num_odor_trials)=1;
            epochs_per_trial(1,num_odor_trials)=0;
            epochs_per_trial(3:4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=2;
            
            %Was dF/F calculated?
            if sum(which_trial_Miss==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Miss_traces(which_trial_Miss==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_Missii_lick=which_Missii_lick(find(which_trial_Miss==trNo,1));
                these_Missii_lick_times=[];
                these_Missii_lick_times=Miss_lick_times(this_Missii_lick,1:Miss_no_lick_times(this_Missii_lick));
                if ~isempty(these_Missii_lick_times)
                    lick_times(num_odor_trials,1:length(these_Missii_lick_times))=these_Missii_lick_times;
                    no_licks(num_odor_trials)=length(these_Missii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
                dLickTraces(num_odor_trials,:)=dMiss_lick_traces(this_Missii_lick,:);
                
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Miss_traces(which_trial_Miss==trNo,win),2);
                    %                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Missii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Missii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=2;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
                %Calculate the average snip for this trial
                Missii=handles_out.Miss_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(Missii).miss_traces);
                num_traces=handles_out.trialNo(Missii).trace_numMiss;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(Missii).trace_numMiss
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(Missii).miss_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventMiss;
            end
        end
        
        if epoch_per_trial(trNo)==8
            %FA
            epochs_per_trial(3,num_odor_trials)=1;
            epochs_per_trial(1:2,num_odor_trials)=0;
            epochs_per_trial(4,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=4;
            
            %Was dF/F calculated?
            if sum(which_trial_FA==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(FA_traces(which_trial_FA==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_FAii_lick=which_FAii_lick(find(which_trial_FA==trNo,1));
                these_FAii_lick_times=[];
                these_FAii_lick_times=FA_lick_times(this_FAii_lick,1:FA_no_lick_times(this_FAii_lick));
                if ~isempty(these_FAii_lick_times)
                    lick_times(num_odor_trials,1:length(these_FAii_lick_times))=these_FAii_lick_times;
                    no_licks(num_odor_trials)=length(these_FAii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
                dLickTraces(num_odor_trials,:)=dFA_lick_traces(this_FAii_lick,:);
                
                for winNo=1:szwins(1)
                    %Calculate dFFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(FA_traces(which_trial_FA==trNo,win),2);
                    %                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_FAii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_FAii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=3;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
                %Calculate the average snip for this trial
                FAii=handles_out.FA_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(FAii).FA_traces);
                num_traces=handles_out.trialNo(FAii).trace_numFA;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(FAii).trace_numFA
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(FAii).FA_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventFA;
            end
        end
        
        if epoch_per_trial(trNo)==9
            %CR
            epochs_per_trial(4,num_odor_trials)=1;
            epochs_per_trial(1:3,num_odor_trials)=0;
            all_lda_events_miss_FA(num_odor_trials)=3;
            
            %Was dF/F calculated?
            if sum(which_trial_CR==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(CR_traces(which_trial_CR==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                
                %Save lick times
                this_CRii_lick=which_CRii_lick(find(which_trial_CR==trNo,1));
                these_CRii_lick_times=[];
                these_CRii_lick_times=CR_lick_times(this_CRii_lick,1:CR_no_lick_times(this_CRii_lick));
                if ~isempty(these_CRii_lick_times)
                    lick_times(num_odor_trials,1:length(these_CRii_lick_times))=these_CRii_lick_times;
                    no_licks(num_odor_trials)=length(these_CRii_lick_times);
                else
                    no_licks(num_odor_trials)=0;
                end
                dLickTraces(num_odor_trials,:)=dCR_lick_traces(this_CRii_lick,:);
                
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(CR_traces(which_trial_CR==trNo,win),2);
                    %                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_CRii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_CRii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                
                %Calculate the average snip for this trial
                CRii=handles_out.CR_trial_no(trNo);
                no_time_points=length(handles_out.componentNo(1).trialNo(CRii).CR_traces);
                num_traces=handles_out.trialNo(CRii).trace_numCR;
                these_traces=zeros(num_traces,no_time_points);
                for trace_num=1:handles_out.trialNo(CRii).trace_numCR
                    these_traces(trace_num,:)=handles_out.componentNo(trace_num).trialNo(CRii).CR_traces;
                end
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventCR;
                
                epochs_per_trial_dFF(num_odor_trials_dFF)=4;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
        
    end
    noROIs(filNum)=szhit(1);
end

%Trim the time course for LDA
all_lda_input_timecourse=all_lda_input_timecourse(1:no_timepoints,:,:);
time_to_eventLDA=time_to_eventLDA(1,1:no_timepoints);

%Calculate percent correct
sliding_window=min_trials; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

%Note: Because this is a reversal I am moving the window for calculation of perCorr to the right by nine points
for ii=1:num_odor_trials-sliding_window+1
    no_Hits=sum(epochs_per_trial(1,ii:ii+sliding_window-1)==1);
    no_CRs=sum(epochs_per_trial(4,ii:ii+sliding_window-1)==1);
    perCorr(ii+sliding_window-1)=100*(no_Hits+no_CRs)/sliding_window;
end


perCorr(1:sliding_window)=perCorr(2*sliding_window+1);

%Note, this is here so that perCorr=0 is included in the 0-10 % bin.
perCorr(perCorr==0)=0.00001;




%Plot percent correct vs trial
figNo=1;

try
    close(figNo)
catch
end

hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])

jj_low=find(perCorr<max_percent_low_beh);
plot(jj_low,perCorr(jj_low),'ob')
hold on
jj_high=find(perCorr>min_precent_high_beh);
plot(jj_high,perCorr(jj_high),'or')

jj_mid=find((perCorr<=min_precent_high_beh)&(perCorr>=max_percent_low_beh));
plot(jj_mid,perCorr(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot([0 num_odor_trials],[50 50],'-k')

%Draw the boundaries of each file
for filNum=2:caimanhandles.caimandr_choices.no_files
    %     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
    if isfield(caimanhandles.caimandr_choices,'start_reversal')
        if caimanhandles.caimandr_choices.start_reversal==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_gogo')
        if caimanhandles.caimandr_choices.start_gogo==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Go-go','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_session')
        if sum(caimanhandles.caimandr_choices.start_session==filNum)>0
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
            %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
end

title(['Percent correct vs. trial number ' choiceFileName(18:end-2)])
xlabel('Trial number')
ylabel('Percent correct')
ylim([0 100])


if ~isfield(caimanhandles.caimandr_choices,'start_reversal')
    caimanhandles.caimandr_choices.start_reversal=200;
end

%For reversals plot violin plot of percent
if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
    %Keep track of the percent correct
    %Plot percent correct vs trial
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
    hold on
    
    %Parameters for violin plot
    edges=0:5:100;
    rand_offset=0.8;
     
    %Trials before reversal
    x_val=1;
    handles_out2.pctPerWin(1).pct=perCorr(1:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1);
    bar(x_val,mean(handles_out2.pctPerWin(1).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(1).pct_mean, handles_out2.pctPerWin(1).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(1).pct,edges,x_val,rand_offset,'k',2);
    
    %Trials after reversal
    x_val=2;
    handles_out2.pctPerWin(2).pct=perCorr(first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+15:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+65);
    bar(x_val,mean(handles_out2.pctPerWin(2).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(2).pct_mean, handles_out2.pctPerWin(2).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(2).pct,edges,x_val,rand_offset,'k',2);
    
    %Trials at end
    x_val=3;
    handles_out2.pctPerWin(3).pct=perCorr(end-100:end);
    bar(x_val,mean(handles_out2.pctPerWin(3).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(3).pct_mean, handles_out2.pctPerWin(3).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(3).pct,edges,x_val,rand_offset,'k',2);
    
    %Draw lines between points
    plot([1 2 3],[mean(handles_out2.pctPerWin(1).pct) mean(handles_out2.pctPerWin(2).pct) mean(handles_out2.pctPerWin(3).pct)],'-k')
    
    ylim([0 120])
    xlim([0.3 3.7])
    ylabel('Percent correct')
    xticks([1 2 3 5 6 7])
    xticklabels({'Before','After','End'})
    
    
    
    window_labels{1}='Before';
    window_labels{2}='After';
    window_labels{3}='End';
    
    fprintf(1, ['\n\nranksum p values for percent correct windows\n\n']);
    
    no_pvals=0;
    for ii=1:3
        for jj=ii+1:3
            no_pvals=no_pvals+1;
            if (adtest(handles_out2.pctPerWin(ii).pct)==1)||(adtest(handles_out2.pctPerWin(ii).pct)==1)
                p_vals_corr(no_pvals)=ranksum(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
                fprintf(1, ['p values ranksum for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
            else
                [h p_vals_corr(no_pvals)]=ttest2(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
                fprintf(1, ['p values t test for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
            end
            
        end
    end
    
    pFDRcorr=drsFDRpval(p_vals_corr);
    fprintf(1, ['pFDR for significant difference percent correct  = %d\n\n'],pFDRcorr);
end

%plot number of ROIs
figNo=figNo+1;
try
    close figNo
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
hold on
plot([1:caimanhandles.caimandr_choices.no_files],noROIs,'-ok')
ylim([0 1.2*max(noROIs)])
title('Number of ROIs per file')
xlabel('File no')
ylabel('No ROIs')

%Plot norm dFF for each window
if caimanhandles.caimandr_choices.start_reversal>0
    if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
        tr_reversal=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal);
    else
        tr_reversal=200000;
    end
else
    tr_reversal=0;
end


for winNo=1:szwins(1)
    figNo=figNo+1;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .05+0.3*(winNo-1) .5 .25])
    hold on
    
    for dFF_trNo=1:num_odor_trials_dFF
        
        %Plot odor 1 (S+ forward)
        subplot(2,1,1)
        hold on
        if dFF_trNo<tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
            if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
        
        
        %Plot odor 2 (S- forward)
        subplot(2,1,2)
        hold on
        if dFF_trNo>=tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
            if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
    end
    
    for spno=1:2
        subplot(2,1,spno)
        if isfield(caimanhandles.caimandr_choices,'start_reversal')
            filNum=caimanhandles.caimandr_choices.start_reversal;
            if (filNum>0)&(filNum<caimanhandles.caimandr_choices.no_files)
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Reversal','Color','k','FontSize',18)
            end
        end
        
        if isfield(caimanhandles.caimandr_choices,'start_gogo')
            filNum=caimanhandles.caimandr_choices.start_gogo;
            if (filNum>0)&(filNum<caimanhandles.caimandr_choices.no_files)
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Go-go','Color','k','FontSize',18)
            end
        end
        
        if isfield(caimanhandles.caimandr_choices,'start_session')
            if length(caimanhandles.caimandr_choices.start_session)>=2
                for sessionNo=2:length(caimanhandles.caimandr_choices.start_session)
                    filNum=caimanhandles.caimandr_choices.start_session(sessionNo);
                    plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k')
                    %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
                end
            end
        end
        
        plot([1 num_odor_trials_dFF],[0 0], '-k')
        
        if spno==1
            title(['Odor 1 (S+ forward)'])
        else
            title(['Odor 2 (S- forward)'])
        end
        xlabel('Trial number')
        ylabel('Normalized dF/F')
        ylim([prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))])
    end
    suptitle(['Normalized dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])
end

%Do a linear discriminant analysis
for winNo=1:szwins(1)
    handles_sig.win(winNo).ii_for_sig=0;
end

if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
    %This is a forward run
    total_trial_windows=3;
    supertitle_description{1}=' percent correct <65';
    supertitle_description{2}=' percent correct >=65&<80';
    supertitle_description{3}=' percent correct >=80';
else
    %Forward and reverse
    total_trial_windows=2;
    supertitle_description{1}='dF/F LDA analysis for trials before reversal';
    supertitle_description{2}='dF/F LDA analysis for trials after reversal at end of the session';
end

firstFig=figNo+1;
maxPC1=-2000;
maxPC2=-2000;
minPC1=20000;
minPC2=20000;

trial_window_description{1}='percent correct <65';
trial_window_description{2}='percent correct >=65&<80';
trial_window_description{3}='percent correct >=80';

lick_threshold=20; %This is the threshold to exclude the runs where Ming was adding water manually
t_odor_on=0;
t_odor_off=4;
all_lda_events_CR_zero_f=[];


for no_trial_windows=1:total_trial_windows
     handles_outs.slope_triggered_LR(no_trial_windows).included=0;
    handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
    dFF_trial_mask=[];
    lick_excluded_trials=[];
    jj=0;
    events_miss_FA=[];
    which_trials_in_PCA=[];
    
    if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
        
        fprintf(1, '\n\nPCA processed for dF/F for trials before reversal \n');
        pct_windows=[45 65;65 80;80 100.1];
        
        for ii=1:num_odor_trials_dFF
            if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                lick_excluded_trials(ii)=0;
                if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
                    dFF_trial_mask(ii)=1;
                    jj=jj+1;
                    handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                    events{jj,1}=all_lda_events{ii};
                    events_miss_FA(jj)=all_lda_events_miss_FA(ii);
                    which_trials_in_PCA(jj)=ii;
                else
                    dFF_trial_mask(ii)=0;
                end
            else
                if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
                    lick_excluded_trials(ii)=1;
                else
                    lick_excluded_trials(ii)=0;
                end
                dFF_trial_mask(ii)=0;
            end
        end
    else
        if no_trial_windows==1
            %Forward trials
            fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
            
            for ii=1:num_odor_trials_dFF
                if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                    lick_excluded_trials(ii)=0;
                    if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
                        dFF_trial_mask(ii)=1;
                        jj=jj+1;
                        handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                        events{jj,1}=all_lda_events{ii};
                    else
                        dFF_trial_mask(ii)=0;
                    end
                else
                    dFF_trial_mask(ii)=0
                    if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
                        lick_excluded_trials(ii)=1;
                    else
                        lick_excluded_trials(ii)=0;
                    end
                end
            end
            
        else
            %Trials at end of reversal
            fprintf(1, '\n\nLDA processed for dF/F for trials after reversal \n');
            for ii=1:num_odor_trials_dFF
                if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                    lick_excluded_trials(ii)=0;
                    if (trial_dFF(ii)>=max(trial_dFF)-100)
                        dFF_trial_mask(ii)=1;
                        jj=jj+1;
                        handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                        events{jj,1}=all_lda_events{ii};
                    else
                        dFF_trial_mask(ii)=0;
                    end
                else
                    if (trial_dFF(ii)>=max(trial_dFF)-100)
                        lick_excluded_trials(ii)=1;
                    else
                        lick_excluded_trials(ii)=0;
                    end
                    dFF_trial_mask(ii)=0;
                end
            end
        end
        
    end
    
    
    
    if sum(dFF_trial_mask)>0
        %Do lick analysis
        %dt_lick=0.3;
        
        dt_lick=time_to_eventLDA(2)-time_to_eventLDA(1);
        time_licks=time_to_eventLDA;
        delta_t_gauss=2; %seconds
        no_conv_points_lick=ceil(delta_t_gauss/(time_licks(2)-time_licks(1)));
        no_conv_points_dFF=ceil(delta_t_gauss/(time_to_eventLDA(2)-time_to_eventLDA(1)));
        
        
        %First figure out the lick threshold to exclude trials wheren Ming
        %gave the animal water manually during the odor on
        all_licks_per_dt_per_trial=zeros(num_odor_trials,ceil((dt_after+dt_before)/dt_lick));
        for trial_no=1:num_odor_trials
            for ii_lick=1:no_licks(trial_no)
                all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
            end
        end
        
        
        hit_odor_lick_freq=zeros(1,sum(all_lda_events_miss_FA==1));
        hit_odor_lick_freq(1,:)=sum(all_licks_per_dt_per_trial(all_lda_events_miss_FA==1,(time_licks>=t_odor_on)&(time_licks<=t_odor_off)),2)/(t_odor_off-t_odor_on);
        
        %         lick_threshold=mean(hit_odor_lick_freq)+2*std(hit_odor_lick_freq);
        %         lick_threshold=20;
        %          lick_threshold=200;
        
        %Plot the lick frequency for S+ and S-
        Splick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        Smlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        
        sp_trno=0;
        sm_trno=0;
        
        for trial_no=1:num_odor_trials
            if dFF_trial_mask(trial_no)==1
                %if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                if strcmp(all_lda_events{trial_no},'S+')
                    %S+
                    for ii_lick=1:no_licks(trial_no)
                        Splick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Splick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                    end
                    sp_trno=sp_trno+1;
                else
                    %S-
                    for ii_lick=1:no_licks(trial_no)
                        Smlick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Smlick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                    end
                    sm_trno=sm_trno+1;
                end
                %end
            end
        end
        
        Splick_freq=(Splick_freq/(sp_trno*dt_lick));
        Smlick_freq=(Smlick_freq/(sm_trno*dt_lick));
        
        %Convolve lick_freq using a window of 0.9 sec
        %         no_conv_points=3;
        %         conv_win=ones(1,no_conv_points);
        
        %Convolve lick_freq using a Gaussian window
        conv_win=gausswin(no_conv_points_lick);
        
        Splick_freq=conv(Splick_freq,conv_win,'same')/sum(conv_win);
        Smlick_freq=conv(Smlick_freq,conv_win,'same')/sum(conv_win);
        
        
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        
        hold on
        
        lfreqmax=20;
        
        p1=plot(time_licks,Smlick_freq(1:length(time_licks)),'-b','LineWidth',2);
        p2=plot(time_licks,Splick_freq(1:length(time_licks)),'-r','LineWidth',2);
        
        %Odor on markers
        plot([0 0],[-1 lfreqmax],'-k')
        odorhl=plot([0 mean(delta_odor)],[-0.5 -0.5],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[-1 lfreqmax],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[-1 lfreqmax],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-0.5 -0.5],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-1 lfreqmax],'-r')
        
        switch no_trial_windows
            case 1
                title(['Lick frequency for percent correct <65']);
            case 2
                title(['Lick frequency for percent correct >=65&<80']);
            case 3
                title(['Lick frequency for percent correct >80']);
        end
        
        xlabel('Time (sec)')
        ylabel('Lick frequency (Hz)')
        legend([p1 p2],'S+','S-')
        xlim([-10 20])
        ylim([-1 lfreqmax])
        
        %Now plot the p value for the difference in licks between S+ and S-
        %Get lick p values
        try
            
            no_pvals=0;
            p_val_Sp_Sm=[];
            time_p_lick=[];
            
            for ii=1:length(time_licks)-1
                
                sp_trno=0;
                sm_trno=0;
                this_Sp=[];
                this_Sm=[];
                
                for trial_no=1:num_odor_trials
                    %                     if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                    if dFF_trial_mask(trial_no)==1
                        if strcmp(all_lda_events{trial_no},'S+')
                            %S+
                            sp_trno=sp_trno+1;
                            this_Sp(sp_trno,1)=0;
                            for ii_lick=1:no_licks(trial_no)
                                if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                    this_Sp(sp_trno,1)=1;
                                end
                            end
                        else
                            %S-
                            sm_trno=sm_trno+1;
                            this_Sm(sm_trno,1)=0;
                            for ii_lick=1:no_licks(trial_no)
                                if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                    this_Sm(sm_trno,1)=1;
                                end
                            end
                        end
                    end
                    %                     end
                end
                
                
                no_pvals=no_pvals+1;
                
                if (~isempty(this_Sm))&(~isempty(this_Sp))
                    p_val_Sp_Sm(no_pvals)=ranksum(this_Sm,this_Sp);
                else
                    p_val_Sp_Sm(no_pvals)=1;
                end
                
                %ranksum gives NaN if the values are all the same
                if isnan(p_val_Sp_Sm(no_pvals))
                    p_val_Sp_Sm(no_pvals)=1;
                end
                
                time_p_lick(no_pvals)= (time_licks(ii)+time_licks(ii+1))/2;
                
            end
            
            logpmin=-20;
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            figure(figNo)
            plot(time_p_lick,log10(p_val_Sp_Sm),'-k','LineWidth',2)
            hold on
            plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
            
            %Odor on markers
            plot([0 0],[logpmin 0.5],'-k')
            odorhl=plot([0 mean(delta_odor)],[logpmin+0.5 logpmin+0.5],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],[logpmin 0.5],'-k')
            
            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[logpmin 0.5],'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin+0.5 logpmin+0.5],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin 0.5],'-r')
            
            
            switch no_trial_windows
                case 1
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct <65']);
                case 2
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct >=65&<80']);
                case 3
                    title(['log(p value) for the difference in licks Sp vs Sm for percent correct >80']);
            end
            
            xlabel('Time (sec)')
            ylabel('log10(p value)')
            ylim([logpmin 0.5])
            xlim([-10 20])
            
        catch
        end
        
        %Plot the licks
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        
        per99=prctile(dLickTraces(:),99.9);
        per1=prctile(dLickTraces(:),1);
        szalllick=size(dLickTraces);
        time_licksd=([1:szalllick(2)]/(acq_rate/20))-dt_before;
        y_shift=0;
        
        %Plot Sp lick traces
        for trial_no=1:num_odor_trials
            if strcmp(all_lda_events{trial_no},'S+')
                if dFF_trial_mask(trial_no)==1
                    plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-r')
                    y_shift=y_shift+1.2*(per99-per1);
                else
                    if lick_excluded_trials(trial_no)==1
                        plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-k')
                        y_shift=y_shift+1.2*(per99-per1);
                    end
                end
            end
            
            
            
        end
        
        %Plot Sm lick traces
        for trial_no=1:num_odor_trials
            if strcmp(all_lda_events{trial_no},'S-')
                if dFF_trial_mask(trial_no)==1
                    plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-b')
                    y_shift=y_shift+1.2*(per99-per1);
                else
                    if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))>lick_threshold
                        plot(time_licksd,dLickTraces(trial_no,:)+y_shift,'-k')
                        y_shift=y_shift+1.2*(per99-per1);
                    end
                end
            end
        end
        
        %Odor on markers
        plot([0 0],[0 y_shift],'-k')
        plot([mean(delta_odor) mean(delta_odor)],[0 y_shift],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 y_shift],'-r')
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 y_shift],'-r')
        
        switch no_trial_windows
            case 1
                title(['lick traces per trial for percent correct <65']);
            case 2
                title(['lick traces per trial for percent correct >=65&<80']);
            case 3
                title(['lick traces per trial for percent correct >80']);
        end
        
        xlabel('Time (sec)')
        
        xlim([-10 20])
        
        
        
        
        
        %Plot each trial for dF/F and lick frequency
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .25 .9 .5])
        
        %Plot each trial forthe derivative of dF/F and lick frequency
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        hFig=figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .25 .9 .5])
        
        
        
        maxdFF=-200;
        mindFF=200;
        
        maxdFF_dx=-200;
        mindFF_dx=200;
        
        
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                conv_win=gausswin(no_conv_points_dFF);
                
                this_conv_dFF=[];
                this_conv_dFF=conv(mean_snip_dFF(trNo,:),conv_win,'same')/sum(conv_win);
                
                this_pct95=prctile(this_conv_dFF(1:132),95);
                maxdFF=max([maxdFF this_pct95]);
                
                this_pct5=prctile(this_conv_dFF(1:132),5);
                mindFF=min([mindFF this_pct5]);
                
                
                %Calculate the derivative of dF/F
                this_conv_dFF_dx=gradient(this_conv_dFF);
                
                this_pct95=prctile(this_conv_dFF_dx(1:132),95);
                maxdFF_dx=max([maxdFF_dx this_pct95]);
                
                this_pct5=prctile(this_conv_dFF_dx(1:132),5);
                mindFF_dx=min([mindFF_dx this_pct5]);
                
            end
        end
        
        ymax=maxdFF+0.1*(maxdFF-mindFF);
        ymin=mindFF-0.1*(maxdFF-mindFF);
        
        ymax_dx=maxdFF_dx+0.1*(maxdFF_dx-mindFF_dx);
        ymin_dx=mindFF_dx-0.1*(maxdFF_dx-mindFF_dx);
        
        t_offset=0;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                evNo=all_lda_events_miss_FA(trNo);
                
                %dF/F
                figure(figNo-1)
                subplot(2,1,1)
                hold on
                
                conv_win=gausswin(no_conv_points_dFF);
                
                this_conv_dFF=[];
                this_conv_dFF=conv(mean_snip_dFF(trNo,:),conv_win,'same')/sum(conv_win);
                
                
                switch evNo
                    case 1
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','r','LineWidth',2);
                    case 2
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','c','LineWidth',2);
                    case 3
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','b','LineWidth',2);
                    case 4
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF(1:132)','m','LineWidth',2);
                end
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
                
                
                
                %derivative of dF/F
                figure(figNo)
                subplot(2,1,1)
                hold on
                
                %Calculate the derivative of dF/F
                this_conv_dFF_dx=gradient(this_conv_dFF);
                
                switch evNo
                    case 1
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','r','LineWidth',2);
                    case 2
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','c','LineWidth',2);
                    case 3
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','b','LineWidth',2);
                    case 4
                        plot(time_to_eventLDA(1:132)'+t_offset,this_conv_dFF_dx(1:132)','m','LineWidth',2);
                end
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin_dx ymax_dx],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin_dx ymax_dx],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin_dx ymax_dx],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx ymax_dx],'-r')
                
                
                %                                 if no_trial_windows==3
                %Calculate the slopes of the three different windows
                handles_outs.no_dFF_slopes=handles_outs.no_dFF_slopes+1;
                
                %Before odor window
                t=time_to_eventLDA((time_to_eventLDA>=-mean(delta_odor))&(time_to_eventLDA<=0))';
                dFF=mean_snip_dFF(trNo,(time_to_eventLDA>=-mean(delta_odor))&(time_to_eventLDA<=0))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.dFF_slopes(handles_outs.no_dFF_slopes,1)=lm.Coefficients{2,1};
                handles_outs.dFF_mean(handles_outs.no_dFF_slopes,1)=mean(dFF);
                
                %Odor window
                t=time_to_eventLDA((time_to_eventLDA>=0)&(time_to_eventLDA<=mean(delta_odor)))';
                dFF=mean_snip_dFF(trNo,(time_to_eventLDA>=0)&(time_to_eventLDA<=mean(delta_odor)))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.dFF_slopes(handles_outs.no_dFF_slopes,2)=lm.Coefficients{2,1};
                handles_outs.dFF_mean(handles_outs.no_dFF_slopes,2)=mean(dFF);
                
                %Reinforcement window
                t=time_to_eventLDA((time_to_eventLDA>=mean(delta_odor_on_reinf_on))&(time_to_eventLDA<=mean(delta_odor_on_reinf_on)+3))';
                dFF=mean_snip_dFF(trNo,(time_to_eventLDA>=mean(delta_odor_on_reinf_on))&(time_to_eventLDA<=mean(delta_odor_on_reinf_on)+3))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.dFF_slopes(handles_outs.no_dFF_slopes,3)=lm.Coefficients{2,1};
                handles_outs.dFF_mean(handles_outs.no_dFF_slopes,3)=mean(dFF);
                
                
                handles_outs.dFF_derivatives(handles_outs.no_dFF_slopes,1:132)=this_conv_dFF_dx(1:132);
                handles_outs.conv_dFF(handles_outs.no_dFF_slopes,1:132)=this_conv_dFF(1:132);
                handles_outs.time_to_eventLDA=time_to_eventLDA(1:132);
                %                                 end
                
                handles_outs.dFF_slopes(handles_outs.no_dFF_slopes,3)=lm.Coefficients{2,1};
%                 handles_outs.no_dFF_slopes=handles_outs.no_dFF_slopes+1;
%                 handles_outs.dFF_derivatives(handles_outs.no_dFF_slopes,1:132)=this_conv_dFF_dx(1:132);
%                 
                handles_outs.time_to_eventLDA=time_to_eventLDA(1:132);
                
                t_offset=t_offset+35;
                
            end
            
            
        end
        
        %         xlim([-5 t_offset-25])
        figure(figNo-1)
        xlim([-50 350])
        ylim([ymin ymax])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('dF/F')
        title(['Timecourse for dF/F for ' supertitle_description{no_trial_windows}])
        
        
        figure(figNo)
        xlim([-50 350])
        ylim([ymin_dx ymax_dx])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('Derivative of dF/F 1/sec')
        title(['Timecourse for derivative dF/F for ' supertitle_description{no_trial_windows}])
        
        %lick frequency
        
        
        maxlick=-200;
        minlick=200;
        
        maxlick_dx=-200;
        minlick_dx=200;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                
                this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                end
                
                %Convolve lick_freq using a window of 0.9 sec
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;
                
                this_pct95=prctile(lick_freq,95);
                maxlick=max([maxlick this_pct95]);
                
                this_pct5=prctile(lick_freq,5);
                minlick=min([minlick this_pct5]);
                
                %Calculate the derivative of lick_freq
                lick_freq_dx=gradient(lick_freq);
                
                this_pct95=prctile(lick_freq_dx,95);
                maxlick_dx=max([maxlick_dx this_pct95]);
                
                this_pct5=prctile(lick_freq_dx,5);
                minlick_dx=min([minlick_dx this_pct5]);
                
                
            end
        end
        
        ymax=maxlick+0.1*(maxlick-minlick);
        ymin=minlick-0.1*(maxlick-minlick);
        
        ymax_dx=maxlick_dx+0.1*(maxlick_dx-minlick_dx);
        ymin_dx=minlick_dx-0.1*(maxlick_dx-minlick_dx);
        
        t_offset=0;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                
                evNo=all_lda_events_miss_FA(trNo);
                
                
                this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                end
                
                
                %Convolve lick_freq using a flat window
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;
                
                
                
                figure(figNo-1)
                subplot(2,1,2)
                hold on
                
                switch evNo
                    case 1
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'r','LineWidth',2);
                    case 2
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'c','LineWidth',2);
                    case 3
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'b','LineWidth',2);
                    case 4
                        plot(time_licks(1:132)+t_offset,lick_freq(1:132),'m','LineWidth',2);
                end
                
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[0 ymax],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
                
               
                %Markers for lick segments
                plot([0+t_offset 0+t_offset],[-1.7 -0.3],'-r','LineWidth',3)
                plot([2+t_offset 2+t_offset],[-1.7 -0.3],'-r','LineWidth',3)
                plot([4+t_offset 4+t_offset],[-1.7 -0.3],'-r','LineWidth',3)
                
                 %Plot the lick trace
                time_mask=(time_licksd<=time_licks(132));
                plot(time_licksd(time_mask)+t_offset,-1.5 +(dLickTraces(trNo,time_mask)-per1)/(per99-per1),'-k')
                
                
                %Calculate the derivative of lick_freq
                lick_freq_dx=gradient(lick_freq);
                
                %Plot the derivative
                figure(figNo)
                subplot(2,1,2)
                hold on
                
                switch evNo
                    case 1
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'r','LineWidth',2);
                    case 2
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'c','LineWidth',2);
                    case 3
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'b','LineWidth',2);
                    case 4
                        plot(time_licks(1:132)+t_offset,lick_freq_dx(1:132),'m','LineWidth',2);
                end
                
                
                %Odor on markers
                plot([0+t_offset 0+t_offset],[ymin_dx ymax_dx],'-k')
                odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-k','LineWidth',5);
                plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin_dx ymax_dx],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin_dx ymax_dx],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx + 0.1*(ymax_dx-ymin_dx) ymin_dx + 0.1*(ymax_dx-ymin_dx)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin_dx ymax_dx],'-r')
                
                
                %Calculate the slopes of the three different windows
                handles_outs.no_lick_slopes=handles_outs.no_lick_slopes+1;
                
                %Before odor window
                t=time_licks((time_licks>=-mean(delta_odor))&(time_licks<=0))';
                dFF=lick_freq((time_licks>=-mean(delta_odor))&(time_licks<=0))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.lick_slopes(handles_outs.no_lick_slopes,1)=lm.Coefficients{2,1};
                handles_outs.mean_lick_freq(handles_outs.no_lick_slopes,1)=mean(dFF);
                
                %Odor window
                t=time_licks((time_licks>=0)&(time_licks<=mean(delta_odor)))';
                dFF=lick_freq((time_licks>=0)&(time_licks<=mean(delta_odor)))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.lick_slopes(handles_outs.no_lick_slopes,2)=lm.Coefficients{2,1};
                handles_outs.mean_lick_freq(handles_outs.no_lick_slopes,2)=mean(dFF);
                
                %Reinforcement window
                t=time_licks((time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3))';
                dFF=lick_freq((time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3))';
                tbl= table(t,dFF,'VariableNames',{'t','dFF'});
                lm = fitlm(tbl,'dFF~t');
                handles_outs.lick_slopes(handles_outs.no_lick_slopes,3)=lm.Coefficients{2,1};
                handles_outs.mean_lick_freq(handles_outs.no_lick_slopes,3)=mean(dFF);
                
               
%                 handles_outs.no_lick_slopes=handles_outs.no_lick_slopes+1;
                handles_outs.lick_derivatives(handles_outs.no_lick_slopes,1:132)=lick_freq_dx(1:132);
                handles_outs.conv_lick(handles_outs.no_lick_slopes,1:132)=lick_freq(1:132);
                handles_outs.time_licks=time_licks(1:132);
                
                t_offset=t_offset+35;
                
            end
        end
        
        %         xlim([-5 t_offset-25])
        figure(figNo-1)
        xlim([-50 350])
%         ylim([ymin ymax])
        ylim([-2 10])
        xlabel('Time (sec)')
        ylabel('Lick rate (Hz)')
        title(['Timecourse for lick rate for ' supertitle_description{no_trial_windows}])
        
        figure(figNo)
        xlim([-50 350])
        ylim([ymin_dx ymax_dx])
        %         ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('Derivative of lick rate (Hz/sec)')
        title(['Timecourse for the drivative of lick rate for ' supertitle_description{no_trial_windows}])
        
        %Stop here to save the example plots
        pffft=1;
        
        
        
        
        %Plot the timecourse for dF/F
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        
        hold on
        
        maxdFF=-200;
        mindFF=200;
        
        for evNo=[4 2 3 1]
            %Find out how many trials for this event
            
            these_dFF_trials=logical(dFF_trial_mask)&(all_lda_events_miss_FA==evNo);
            
            
            if sum(these_dFF_trials)>0
                
                
                if sum(these_dFF_trials)>2
                    CI=[];
                    CI = bootci(1000, {@mean, mean_snip_dFF(these_dFF_trials,:)})';
                    CI(:,1)=mean(mean_snip_dFF(these_dFF_trials,:))'-CI(:,1);
                    CI(:,2)=CI(:,2)-mean(mean_snip_dFF(these_dFF_trials,:))';
                    
                    maxdFF=max([maxdFF,max(CI(:,1)+mean(mean_snip_dFF(these_dFF_trials,:))')]);
                    mindFF=min([mindFF,min(mean(mean_snip_dFF(these_dFF_trials,:))'-CI(:,2))]);
                    
                    switch evNo
                        case 1
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'r');
                        case 2
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'c');
                        case 3
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'b');
                        case 4
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', CI(1:length(time_to_eventLDA),:), 'm');
                    end
                else
                    switch evNo
                        case 1
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))','r');
                        case 2
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', 'c');
                        case 3
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))', 'b');
                        case 4
                            plot(time_to_eventLDA',mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA)))',  'm');
                    end
                end
            end
            
        end
        
        ymax=maxdFF+0.1*(maxdFF-mindFF);
        ymin=mindFF-0.1*(maxdFF-mindFF);
        
        %Odor on markers
        plot([0 0],[ymin ymax],'-k')
        odorhl=plot([0 mean(delta_odor)],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[ymin ymax],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[ymin ymax],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[ymin ymax],'-r')
        
        xlim([-10 20])
        %         ylim([ymin ymax])
        ylim([-0.4 1])
        xlabel('Time (sec)')
        ylabel('dF/F')
        title(['Timecourse for dF/F for ' supertitle_description{no_trial_windows}])
        
        %Calculate and plot the mean lick rate and dFF
        %Note: This overwrites the other percent correct
        handles_outs.window=[];
        for evNo=1:6
            handles_outs.no_lick_freq_choice_win(evNo)=0;
            for winNo=1:length(caimanhandles.caimandr_choices.wins)
                handles_outs.window(winNo).no_lick_freq_choice_win(evNo)=0;
            end
        end
        
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                actual_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    if floor((lick_times(trNo,ii_lick)+dt_before)/dt_lick)>0
                        actual_lick_freq( floor((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=actual_lick_freq( floor((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                    end
                end
                
                
                %Convolve lick_freq using a flat window
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(actual_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;
                
                
                evNo=all_lda_events_miss_FA(trNo);
                handles_outs.no_lick_freq_choice_win(evNo)=handles_outs.no_lick_freq_choice_win(evNo)+1;
                
                for winNo=1:length(caimanhandles.caimandr_choices.wins)
                    if winNo==2
                        %On ourpose I make a time_mask covering the two 2 sec time seegments for decision making
                        this_time_mask=(time_licks>=0)&(time_licks<=4);
                    else
                        this_time_mask=(time_licks>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_licks<=caimanhandles.caimandr_choices.wins(winNo,2));
                    end
                    this_lick_freq=actual_lick_freq(this_time_mask)/dt';
                    handles_outs.window(winNo).event(evNo).lick_freq_choice_win(handles_outs.no_lick_freq_choice_win(evNo))=mean(this_lick_freq);
                    this_dFF=mean_snip_dFF(trNo,(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_eventLDA<=caimanhandles.caimandr_choices.wins(winNo,2)))';
                    handles_outs.window(winNo).event(evNo).dFF_mean_choice_win(handles_outs.no_lick_freq_choice_win(evNo))=mean(this_dFF);
                    
                    %Calculate mean lick frequency for the two 4 sec periods
                    %                 this_lick_freq4=lick_freq((time_licks>=0)&(time_licks<=4))';
                    handles_outs.window(winNo).no_lick_freq_choice_win(evNo)=handles_outs.window(winNo).no_lick_freq_choice_win(evNo)+1;
                    handles_outs.window(winNo).event(evNo).dFF_choice_win(handles_outs.window(winNo).no_lick_freq_choice_win(evNo))=mean(this_dFF);

                    if evNo==3
                        time_mask=(time_licks>=0)&(time_licks<=4);  %Time mask covering the two 2 sec time seegments for decision making
                        this_actual_lick_freq=mean(actual_lick_freq(time_mask))/dt;
                        if this_actual_lick_freq>0
                            %There are licks in the two 2 sec segments
                            all_lda_events_CR_zero_f(trNo)=1;
                            handles_outs.window(winNo).no_lick_freq_choice_win(5)=handles_outs.window(winNo).no_lick_freq_choice_win(5)+1;
                            handles_outs.window(winNo).event(5).lick_freq_choice_win(handles_outs.window(winNo).no_lick_freq_choice_win(5))=mean(this_lick_freq);
                            handles_outs.window(winNo).event(5).dFF_choice_win(handles_outs.window(winNo).no_lick_freq_choice_win(5))=mean(this_dFF);
                        else
                            %No licks
                            all_lda_events_CR_zero_f(trNo)=2;
                            handles_outs.window(winNo).no_lick_freq_choice_win(6)=handles_outs.window(winNo).no_lick_freq_choice_win(6)+1;
                            handles_outs.window(winNo).event(6).lick_freq_choice_win(handles_outs.window(winNo).no_lick_freq_choice_win(6))=mean(this_lick_freq);
                            handles_outs.window(winNo).event(6).dFF_choice_win(handles_outs.window(winNo).no_lick_freq_choice_win(6))=mean(this_dFF);
                        end
                    else
                        all_lda_events_CR_zero_f(trNo)=-1;
                    end
                end
            end
        end
        
        if no_trial_windows==3
            pffft=1;
        end
        
        %Plot the lick frequency
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        ii_bar=0;
        
        for winNo=1:length(caimanhandles.caimandr_choices.wins)
            
            for evNo=1:4
                
                if handles_outs.no_lick_freq_choice_win(evNo)>0
                    these_lfs=zeros(length(handles_outs.window(winNo).event(evNo).lick_freq_choice_win),1);
                    these_lfs(:,1)=handles_outs.window(winNo).event(evNo).lick_freq_choice_win;
                    this_mean_lf=mean(these_lfs);
                    switch evNo
                        case 1
                            bar(ii_bar,this_mean_lf, 'r');
                        case 2
                            bar(ii_bar,this_mean_lf, 'c');
                        case 3
                            bar(ii_bar,this_mean_lf, 'b');
                        case 4
                            bar(ii_bar,this_mean_lf, 'm');
                    end
                    
                    if handles_outs.window(winNo).no_lick_freq_choice_win(evNo)>2
                        CI=[];
                        CI = bootci(1000, {@mean, these_lfs})';
                        plot([ii_bar ii_bar],CI,'-k','LineWidth',2) 
                    end
                    ii_bar=ii_bar+1;
                end
               
            end
             ii_bar=ii_bar+1;
        end
        title('Lick frequency for the different events')
        

        %Plot the lick frequency for CR zero licks vs CR
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        ii_bar=0;
        
        for winNo=1:length(caimanhandles.caimandr_choices.wins)
            
            for evNo_CR=5:6
                
                if handles_outs.window(winNo).no_lick_freq_choice_win(evNo_CR)>0
                    these_lfs=zeros(length(handles_outs.window(winNo).event(evNo_CR).lick_freq_choice_win),1);
                    these_lfs(:,1)=handles_outs.window(winNo).event(evNo_CR).lick_freq_choice_win;
                    this_mean_lf=mean(these_lfs);
                    switch evNo_CR
                        case 5
                            bar(ii_bar,this_mean_lf, 'r');
                        case 6
                            bar(ii_bar,this_mean_lf, 'b');
                    end
                    
                    if handles_outs.window(winNo).no_lick_freq_choice_win(evNo_CR)>2
                        CI=[];
                        CI = bootci(1000, {@mean, these_lfs})';
                        plot([ii_bar ii_bar],CI,'-k','LineWidth',2) 
                    end
                    ii_bar=ii_bar+1;
                end
               
            end
             ii_bar=ii_bar+1;
        end
        title('Lick frequency for CR no licks=blue, CR licks=red')
        
        %Plot dFF for CR zero licks vs CR licks
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        hold on
        
        ii_bar=0;
        
        for winNo=1:length(caimanhandles.caimandr_choices.wins)
            
            for evNo_CR=5:6
                
                if handles_outs.window(winNo).no_lick_freq_choice_win(evNo_CR)>0
                    these_dffs=zeros(length(handles_outs.window(winNo).event(evNo_CR).dFF_choice_win),1);
                    these_dffs(:,1)=handles_outs.window(winNo).event(evNo_CR).dFF_choice_win;
                    this_mean_dFF=mean(these_dffs);
                    switch evNo_CR
                           case 5
                            bar(ii_bar,this_mean_dFF, 'r');
                        case 6
                            bar(ii_bar,this_mean_dFF, 'b');
                    end
                    
                    if handles_outs.window(winNo).no_lick_freq_choice_win(evNo_CR)>2
                        CI=[];
                        CI = bootci(1000, {@mean, these_dffs})';
                        plot([ii_bar ii_bar],CI,'-k','LineWidth',2) 
                    end
                    ii_bar=ii_bar+1;
                end
               
            end
             ii_bar=ii_bar+1;
        end
        title('dFF for CR no licks=blue, CR licks=red')
        
        
        %Now do the analysis for the licks aligned to areas with large changes in
        %dF/F
        dFF_slope_thr=0.03;
        lick_freq_dFFslope_triggered=[];
        lick_freq_dx_dFFslope_triggered=[];
        dFF_dx_dFFslope_triggered=[];
        slope_triggered_events_miss_FA=[];
        dFF_dFFslope_triggered=[];
        dFFslope_triggered_t=[];
        ii_dFF_dx=0;
        no_bins=21;
        dt=time_to_eventLDA(2)-time_to_eventLDA(1);
        t_delta=-dt*((no_bins-1)/2):dt:dt*((no_bins-1)/2);
        trial_no=0;
        
        for trNo=1:length(dFF_trial_mask)
            
            if dFF_trial_mask(trNo)
                
                evNo=all_lda_events_miss_FA(trNo);
                
                %dF/F
                conv_win=gausswin(no_conv_points_dFF);
                
                this_conv_dFF=[];
                this_conv_dFF=conv(mean_snip_dFF(trNo,:),conv_win,'same')/sum(conv_win);
                this_conv_dFF=this_conv_dFF(1:132);
                
                %Calculate the derivative of dF/F
                this_conv_dFF_dx=gradient(this_conv_dFF);
                
                this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                for ii_lick=1:no_licks(trNo)
                    this_lick_freq(ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trNo,ii_lick)+dt_before)/dt_lick))+1;
                end
                
                
                %Convolve lick_freq using a flat window
                %                 no_conv_points=4;
                %                 conv_win=ones(1,no_conv_points);
                
                %Convolve lick_freq using a Gaussian window
                conv_win=gausswin(no_conv_points_lick);
                
                lick_freq=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
                lick_freq=lick_freq/dt_lick;
                lick_freq=lick_freq(1:132);
                
                %Calculate the derivative of lick_freq
                lick_freq_dx=gradient(lick_freq);
                
                trial_no=trial_no+1;
              
                
                %First calculate the PSTH for the licks asssociated with large slopes
                ii=find(time_to_eventLDA>=0,1,'first');
                
%               these_licks_high_dFF_dx=zeros(1,no_bins);
                this_dFFslope_triggered_dFF=zeros(1,no_bins);
                this_dFFslope_triggered_lick_f=zeros(1,no_bins);
                this_dFFslope_triggered_lick_f_dx=zeros(1,no_bins);
                this_dFFslope_triggered_dFF_dx=zeros(1,no_bins);
                
                
                %I amlooking for the first large slope change in dFF within the odor period
                if ~isempty(find(this_conv_dFF_dx(ii:end)>=dFF_slope_thr,1,'first'))
                    
                    %Do lick accounting
                    this_ii=find(this_conv_dFF_dx(ii:end)>=dFF_slope_thr,1,'first');
                    
                    if (this_ii+ii-1-((no_bins-1)/2)>=1)
                        if ((this_ii+ii-1+((no_bins-1)/2)<=length(this_conv_dFF_dx)))&(time_to_eventLDA(this_ii+ii-1)<=mean(delta_odor_on_reinf_on))
                            ii_dFF_dx= ii_dFF_dx+1;
                            
                            this_time=time_to_eventLDA(this_ii+ii-1);
                            
                            dFFslope_triggered_t(ii_dFF_dx)=this_time;
                            
                            this_dFFslope_triggered_dFF(1,:)=this_conv_dFF(this_ii+ii-1-((no_bins-1)/2):this_ii+ii-1+((no_bins-1)/2));
                            this_dFFslope_triggered_lick_f(1,:)=lick_freq(this_ii+ii-1-((no_bins-1)/2):this_ii+ii-1+((no_bins-1)/2));
                            
                            this_dFFslope_triggered_dFF_dx(1,:)=this_conv_dFF_dx(this_ii+ii-1-((no_bins-1)/2):this_ii+ii-1+((no_bins-1)/2));
                            this_dFFslope_triggered_lick_f_dx(1,:)=lick_freq_dx(this_ii+ii-1-((no_bins-1)/2):this_ii+ii-1+((no_bins-1)/2));
                            
                            lick_freq_dFFslope_triggered(ii_dFF_dx,:)=this_dFFslope_triggered_lick_f;
                            dFF_dFFslope_triggered(ii_dFF_dx,:)=this_dFFslope_triggered_dFF;
                            lick_freq_dx_dFFslope_triggered(ii_dFF_dx,:)=this_dFFslope_triggered_lick_f_dx;
                            dFF_dx_dFFslope_triggered(ii_dFF_dx,:)=this_dFFslope_triggered_dFF_dx;
                            slope_triggered_events_miss_FA(ii_dFF_dx)=all_lda_events_miss_FA(trNo);

                            %This is here to browse the choice of slope-triggered events
                            view_event=1;
                            if view_event==1
                                figNo=figNo+1;
                                try
                                    close(figNo)
                                catch
                                end
                                figure(figNo)
                                
                                hold on
                                
                                %dFF
                                subplot(4,1,1)
                                hold on
                                
                                switch evNo
                                    case 1
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF(1:132)','r','LineWidth',2);
                                    case 2
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF(1:132)','c','LineWidth',2);
                                    case 3
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF(1:132)','b','LineWidth',2);
                                    case 4
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF(1:132)','m','LineWidth',2);
                                end
                                plot([time_to_eventLDA(this_ii+ii-1) time_to_eventLDA(this_ii+ii-1)],[-0.75 2.5],'-k')
                                xlim([time_to_eventLDA(1) time_to_eventLDA(132)])
                                ylim([-0.75 2.5])
                                title('dF/F')
                                
                                %dFF derivative
                                subplot(4,1,2)
                                hold on
                                
                                switch evNo
                                    case 1
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF_dx(1:132)','r','LineWidth',2);
                                    case 2
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF_dx(1:132)','c','LineWidth',2);
                                    case 3
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF_dx(1:132)','b','LineWidth',2);
                                    case 4
                                        plot(time_to_eventLDA(1:132)',this_conv_dFF_dx(1:132)','m','LineWidth',2);
                                end
                                plot([time_to_eventLDA(this_ii+ii-1) time_to_eventLDA(this_ii+ii-1)],[-0.1 0.1],'-k')
                                xlim([time_to_eventLDA(1) time_to_eventLDA(132)])
                                ylim([-0.1 0.1])
                                title('Derivative of dF/F')
                                
                                %lick frequency
                                subplot(4,1,3)
                                hold on
                                
                                switch evNo
                                    case 1
                                        plot(time_to_eventLDA(1:132)',lick_freq(1:132)','r','LineWidth',2);
                                    case 2
                                        plot(time_to_eventLDA(1:132)',lick_freq(1:132)','c','LineWidth',2);
                                    case 3
                                        plot(time_to_eventLDA(1:132)',lick_freq(1:132)','b','LineWidth',2);
                                    case 4
                                        plot(time_to_eventLDA(1:132)',lick_freq(1:132)','m','LineWidth',2);
                                end
                                plot([time_to_eventLDA(this_ii+ii-1) time_to_eventLDA(this_ii+ii-1)],[-1 15],'-k')
                                xlim([time_to_eventLDA(1) time_to_eventLDA(132)])
                                ylim([-1 15])
                                title('Lick frequency')
                                
                                %dFF derivative
                                subplot(4,1,4)
                                hold on
                                
                                switch evNo
                                    case 1
                                        plot(time_to_eventLDA(1:132)',lick_freq_dx(1:132)','r','LineWidth',2);
                                    case 2
                                        plot(time_to_eventLDA(1:132)',lick_freq_dx(1:132)','c','LineWidth',2);
                                    case 3
                                        plot(time_to_eventLDA(1:132)',lick_freq_dx(1:132)','b','LineWidth',2);
                                    case 4
                                        plot(time_to_eventLDA(1:132)',lick_freq_dx(1:132)','m','LineWidth',2);
                                end
                                plot([time_to_eventLDA(this_ii+ii-1) time_to_eventLDA(this_ii+ii-1)],[-2 2],'-k')
                                xlim([time_to_eventLDA(1) time_to_eventLDA(132)])
                                ylim([-2 2])
                                title('Derivative of lick frequency')
                                
                                suptitle(['Trial no ', num2str(trNo)])
                                
                                figNo=figNo-1;
                                
                                pffft=1;
                            end
                            
                            
                            
                        end
                        
                        
                    end
                    
                end
                
                
            end
        end
                  
        handles_outs.slope_triggered_LR(no_trial_windows).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).lick_freq_dx_dFFslope_triggered=lick_freq_dx_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).dFF_dFFslope_triggered=dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).dFF_dx_dFFslope_triggered=dFF_dx_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).slope_triggered_events_miss_FA=slope_triggered_events_miss_FA;
        
        handles_outs.slope_triggered_LR(no_trial_windows).dFFslope_triggered_t=dFFslope_triggered_t;
        

        %Now plot the slope-triggered lick frequency
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        
        hold on
        
        t_delta=-dt*((no_bins-1)/2):dt:dt*((no_bins-1)/2);
        
        handles_outs.slope_triggered_LR(no_trial_windows).included=1;
        handles_outs.t_delta=t_delta;
        
        %S+
        these_ii=(slope_triggered_events_miss_FA==2)|(slope_triggered_events_miss_FA==1);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(5).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(5).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
        
        %S-
        these_ii=(slope_triggered_events_miss_FA==3)|(slope_triggered_events_miss_FA==4);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(6).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(6).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
       
        %Miss
        these_ii=(slope_triggered_events_miss_FA==2);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(2).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(2).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, lick_freq_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_lick_freq_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_lick_freq_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_lick_freq_dFFslope_triggered',  CI_lickf_slope_triggered', 'c');
        else
            plot(t_delta',mean_lick_freq_dFFslope_triggered','c')
        end
        
        %FA
        these_ii=(slope_triggered_events_miss_FA==4);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(4).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(4).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, lick_freq_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_lick_freq_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_lick_freq_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_lick_freq_dFFslope_triggered',  CI_lickf_slope_triggered', 'm');
        else
            plot(t_delta',mean_lick_freq_dFFslope_triggered','m')
        end
        
        %CR
        these_ii=(slope_triggered_events_miss_FA==3);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(3).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(3).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, lick_freq_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_lick_freq_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_lick_freq_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_lick_freq_dFFslope_triggered',  CI_lickf_slope_triggered', 'b');
        else
            plot(t_delta',mean_lick_freq_dFFslope_triggered','b')
        end
        
        %Hit
        these_ii=(slope_triggered_events_miss_FA==1);
        mean_lick_freq_dFFslope_triggered=mean(lick_freq_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(1).mean_lick_freq_dFFslope_triggered=mean_lick_freq_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(1).lick_freq_dFFslope_triggered=lick_freq_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, lick_freq_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_lick_freq_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_lick_freq_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_lick_freq_dFFslope_triggered',  CI_lickf_slope_triggered', 'r');
        else
            plot(t_delta',mean_lick_freq_dFFslope_triggered','r')
        end
        
            
        xlabel('Time (sec)')
        ylabel('Lick frequency')
        title(['dF/F slope-triggered-lick frequency ' supertitle_description{no_trial_windows}])
        
        
        %Now plot the slope-triggered dFF
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        
        hold on
        
        t_delta=-dt*((no_bins-1)/2):dt:dt*((no_bins-1)/2);
        
        handles_outs.slope_triggered_LR(no_trial_windows).included=1;
        
        %S+
        these_ii=(slope_triggered_events_miss_FA==2)|(slope_triggered_events_miss_FA==1);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(5).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(5).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        %S-
        these_ii=(slope_triggered_events_miss_FA==3)|(slope_triggered_events_miss_FA==4);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(6).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(6).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        %Miss
        these_ii=(slope_triggered_events_miss_FA==2);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(2).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(2).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, dFF_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_dFF_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_dFF_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_dFF_dFFslope_triggered',  CI_lickf_slope_triggered', 'c');
        else
            plot(t_delta',mean_dFF_dFFslope_triggered','c')
        end
        
        
     
        
        %FA
        these_ii=(slope_triggered_events_miss_FA==4);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(4).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(4).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, dFF_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_dFF_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_dFF_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_dFF_dFFslope_triggered',  CI_lickf_slope_triggered', 'm');
        else
            plot(t_delta',mean_dFF_dFFslope_triggered','m')
        end
        
          %CR
        these_ii=(slope_triggered_events_miss_FA==3);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(3).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(3).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, dFF_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_dFF_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_dFF_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_dFF_dFFslope_triggered',  CI_lickf_slope_triggered', 'b');
        else
            plot(t_delta',mean_dFF_dFFslope_triggered','b')
        end
        
        %Hit
        these_ii=(slope_triggered_events_miss_FA==1);
        mean_dFF_dFFslope_triggered=mean(dFF_dFFslope_triggered(these_ii,:),1);
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(1).mean_dFF_dFFslope_triggered=mean_dFF_dFFslope_triggered;
        handles_outs.slope_triggered_LR(no_trial_windows).epoch(1).dFF_dFFslope_triggered=dFF_dFFslope_triggered(these_ii,:);
        
        if sum(these_ii)>2
            CI_lickf_slope_triggered=bootci(1000, @mean, dFF_dFFslope_triggered(these_ii,:),1);
            CI_lickf_slope_triggered(1,:)=mean_dFF_dFFslope_triggered-CI_lickf_slope_triggered(1,:);
            CI_lickf_slope_triggered(2,:)=CI_lickf_slope_triggered(2,:)-mean_dFF_dFFslope_triggered;
            
            [hlCR, hpCR] = boundedline(t_delta',mean_dFF_dFFslope_triggered',  CI_lickf_slope_triggered', 'r');
        else
            plot(t_delta',mean_dFF_dFFslope_triggered','r')
        end
        
            
    
            
        xlabel('Time (sec)')
        ylabel('dF/F')
        title(['dF/F ' supertitle_description{no_trial_windows}])
        
       
        
    end
    pffft=1;
end



%
% for fNo=firstFig:firstFig+2*szwins(1)-1
%     figure(fNo)
%     xlim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
%     ylim([minPC2-0.1*(maxPC2-minPC2) maxPC2+0.1*(maxPC2-minPC2)])
% end

%Plot the slopes
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo)
set(hFig, 'units','normalized','position',[.15 .25 .5 .25])

maxy=-200;
miny=200;

maxx=-200;
minx=200;

for ii=1:3
    maxx=max([maxx max(handles_outs.lick_slopes(:,ii))]);
    minx=min([minx min(handles_outs.lick_slopes(:,ii))]);

    maxy=max([maxy max(handles_outs.dFF_slopes(:,ii))]);
    miny=min([miny min(handles_outs.dFF_slopes(:,ii))]);
end


hold on

subplot(1,3,1)
hold on
plot( handles_outs.lick_slopes(:,1), handles_outs.dFF_slopes(:,1),'ob')
[rho,pval] = corr(handles_outs.lick_slopes(:,1),handles_outs.dFF_slopes(:,1));
fprintf(1, ['rho and p value for dFF slopes vs lick slopes for pre-odor window  = %d, %d\n\n'],rho,pval);
xlabel('slope for lick frequency Hz/sec')
ylabel('slope for dF/F')
handles_outs.slope_rho(1)=rho;
handles_outs.slope_pval(1)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Pre-odor')

%Least squares fit
x=handles_outs.lick_slopes(:,1);
y=handles_outs.dFF_slopes(:,1);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)


subplot(1,3,2)
hold on
plot( handles_outs.lick_slopes(:,2), handles_outs.dFF_slopes(:,2),'om')
[rho,pval] = corr(handles_outs.lick_slopes(:,2),handles_outs.dFF_slopes(:,2));
fprintf(1, ['rho and p value for dFF slopes vs lick slopes for odor window  = %d, %d\n\n'],rho,pval);
xlabel('slope for lick frequency Hz/sec')
ylabel('slope for dF/F')
handles_outs.slope_rho(2)=rho;
handles_outs.slope_pval(2)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Odor')

%Least squares fit
x=handles_outs.lick_slopes(:,2);
y=handles_outs.dFF_slopes(:,2);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)


subplot(1,3,3)
hold on
plot( handles_outs.lick_slopes(:,3), handles_outs.dFF_slopes(:,3),'or')
[rho,pval] = corr(handles_outs.lick_slopes(:,3),handles_outs.dFF_slopes(:,3));
fprintf(1, ['rho and p value for dFF slopes vs lick slopes for reinforcement window  = %d, %d\n\n'],rho,pval);
xlabel('slope for lick frequency Hz/sec')
ylabel('slope for dF/F')
handles_outs.slope_rho(3)=rho;
handles_outs.slope_pval(3)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Reinforcement')

%Least squares fit
x=handles_outs.lick_slopes(:,3);
y=handles_outs.dFF_slopes(:,3);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)



%Plot mean values for dFF and lick rate
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo)
set(hFig, 'units','normalized','position',[.15 .25 .5 .25])

maxy=-200;
miny=200;

maxx=-200;
minx=200;

for ii=1:3
    maxx=max([maxx max(handles_outs.mean_lick_freq(:,ii))]);
    minx=min([minx min(handles_outs.mean_lick_freq(:,ii))]);

    maxy=max([maxy max(handles_outs.dFF_mean(:,ii))]);
    miny=min([miny min(handles_outs.dFF_mean(:,ii))]);
end


hold on

subplot(1,3,1)
hold on
plot( handles_outs.mean_lick_freq(:,1), handles_outs.dFF_mean(:,1),'ob')
[rho,pval] = corr(handles_outs.mean_lick_freq(:,1),handles_outs.dFF_mean(:,1));
fprintf(1, ['rho and p value for mean dFF vs mean lick rate for pre-odor window  = %d, %d\n\n'],rho,pval);
xlabel('lick frequency Hz/sec')
ylabel('mean dF/F')
handles_outs.mean_rho(1)=rho;
handles_outs.mean_pval(1)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Pre-odor')

%Least squares fit
x=handles_outs.mean_lick_freq(:,1);
y=handles_outs.dFF_mean(:,1);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)


subplot(1,3,2)
hold on
plot( handles_outs.mean_lick_freq(:,2), handles_outs.dFF_mean(:,2),'om')
[rho,pval] = corr(handles_outs.mean_lick_freq(:,2),handles_outs.dFF_mean(:,2));
fprintf(1, ['rho and p value for mean dFF  vs lick rate for odor window  = %d, %d\n\n'],rho,pval);
xlabel('lick frequency Hz/sec')
ylabel('mean dF/F')
handles_outs.mean_rho(2)=rho;
handles_outs.mean_pval(2)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Odor')

%Least squares fit
x=handles_outs.mean_lick_freq(:,2);
y=handles_outs.dFF_mean(:,2);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)


subplot(1,3,3)
hold on
plot( handles_outs.mean_lick_freq(:,3), handles_outs.dFF_mean(:,3),'or')
[rho,pval] = corr(handles_outs.mean_lick_freq(:,3),handles_outs.dFF_mean(:,3));
fprintf(1, ['rho and p value for mean dFF vs lick rates for reinforcement window  = %d, %d\n\n'],rho,pval);
xlabel('lick frequency Hz/sec')
ylabel('mean dF/F')
handles_outs.mean_rho(3)=rho;
handles_outs.mean_pval(3)=pval;
xlim([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)])
ylim([miny-0.1*(maxy-miny) maxy+0.1*(maxy-miny)])
title('Reinforcement')

%Least squares fit
x=handles_outs.mean_lick_freq(:,3);
y=handles_outs.dFF_mean(:,3);
tbl= table(x,y,'VariableNames',{'x','y'});
lm = fitlm(tbl,'y~x');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};
plot([minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)], this_slope*[minx-0.1*(maxx-minx) maxx+0.1*(maxx-minx)]+this_intercept,'-k','LineWIdth',3)



%Plot the derivative plot
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.1 .25 .75 .25])

% %For the entire timecourse
% subplot(1,4,1)
% hold on
% 
lick_derivatives=[];
dFF_derivatives=[];
for trNo=1:handles_outs.no_lick_slopes
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,:)];
    dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,:)];
    plot(handles_outs.lick_derivatives(trNo,:),handles_outs.dFF_derivatives(trNo,:),'.b')
end

ymax=max(dFF_derivatives)+0.1*(max(dFF_derivatives)-min(dFF_derivatives));
ymin=min(dFF_derivatives)-0.1*(max(dFF_derivatives)-min(dFF_derivatives));

xmax=max(lick_derivatives)+0.1*(max(lick_derivatives)-min(lick_derivatives));
xmin=min(lick_derivatives)-0.1*(max(lick_derivatives)-min(lick_derivatives));
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative  = %d, %d\n\n'],rho,pval);
% % handles_outs.rho(1)=rho;
% % handles_outs.pval(1)=pval;
% xlabel('lick frequency derivativeHz/sec')
% ylabel('dF/F derivative 1/sec')
% title('All times')
% xlim([xmin xmax])
% ylim([ymin ymax])

%Now do the three windows
time_licks=[];
time_licks=handles_outs.time_licks;

%Before odor window
subplot(1,3,1)
hold on

time_mask=(time_licks>=-mean(delta_odor))&(time_licks<=0);
lick_derivatives=[];
dFF_derivatives=[];
for trNo=1:handles_outs.no_lick_slopes
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
    dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
    plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
end


tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(lick_derivatives',dFF_derivatives');
fprintf(1, ['rho and p value for dFF derivative vs lick derivative before odor = %d, %d\n\n'],rho,pval);
handles_outs.DtdFF_rho(1)=rho;
handles_outs.DtdFF_pval(1)=pval;
xlabel('lick frequency derivativeHz/sec')
ylabel('dF/F derivative 1/sec')
title('Before odor')
xlim([xmin xmax])
ylim([ymin ymax])

%Odor window
subplot(1,3,2)
hold on


time_mask=(time_licks>=0)&(time_licks<=mean(delta_odor));
lick_derivatives=[];
dFF_derivatives=[];
for trNo=1:handles_outs.no_lick_slopes
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
    dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
    plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
end


tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(lick_derivatives',dFF_derivatives');
fprintf(1, ['rho and p value for dFF derivative vs lick derivative during odor = %d, %d\n\n'],rho,pval);
handles_outs.DtdFF_rho(2)=rho;
handles_outs.DtdFF_pval(2)=pval;
xlabel('lick frequency derivativeHz/sec')
ylabel('dF/F derivative 1/sec')
title('During odor')
xlim([xmin xmax])
ylim([ymin ymax])

%Reinforcement window
subplot(1,3,3)
hold on

time_mask=(time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3);
lick_derivatives=[];
dFF_derivatives=[];
for trNo=1:handles_outs.no_lick_slopes
    lick_derivatives=[lick_derivatives handles_outs.lick_derivatives(trNo,time_mask)];
    dFF_derivatives=[dFF_derivatives handles_outs.dFF_derivatives(trNo,time_mask)];
    plot(handles_outs.lick_derivatives(trNo,time_mask),handles_outs.dFF_derivatives(trNo,time_mask),'.b')
end


tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(lick_derivatives',dFF_derivatives');
fprintf(1, ['rho and p value for dFF derivative vs lick derivative during reinforcement = %d, %d\n\n'],rho,pval);
handles_outs.DtdFF_rho(3)=rho;
handles_outs.DtdFF_pval(3)=pval;
xlabel('lick frequency derivativeHz/sec')
ylabel('dF/F derivative 1/sec')
title('Reinforcement')
xlim([xmin xmax])
ylim([ymin ymax])


%Plot the dFF vs lick freq plot
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.1 .25 .75 .25])

% %For the entire timecourse
% subplot(1,4,1)
% hold on
% 
conv_lick=[];
conv_dFF=[];
for trNo=1:handles_outs.no_lick_slopes
    conv_lick=[conv_lick handles_outs.conv_lick(trNo,:)];
    conv_dFF=[conv_dFF handles_outs.conv_dFF(trNo,:)];
    plot(handles_outs.conv_lick(trNo,:),handles_outs.conv_dFF(trNo,:),'.b')
end

ymax=max(conv_dFF)+0.1*(max(conv_dFF)-min(conv_dFF));
ymin=min(conv_dFF)-0.1*(max(conv_dFF)-min(conv_dFF));

xmax=max(conv_lick)+0.1*(max(conv_lick)-min(conv_lick));
xmin=min(conv_lick)-0.1*(max(conv_lick)-min(conv_lick));
% 
% tbl= table(lick_derivatives',dFF_derivatives','VariableNames',{'lick_derivatives','dFF_derivatives'});
% lm = fitlm(tbl,'dFF_derivatives~lick_derivatives');
% this_slope=lm.Coefficients{2,1};
% this_intercept=lm.Coefficients{1,1};
% 
% plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)
% 
% [rho,pval] = corr(lick_derivatives',dFF_derivatives');
% fprintf(1, ['rho and p value for dFF derivative vs lick derivative  = %d, %d\n\n'],rho,pval);
% % handles_outs.rho(1)=rho;
% % handles_outs.pval(1)=pval;
% xlabel('lick frequency derivativeHz/sec')
% ylabel('dF/F derivative 1/sec')
% title('All times')
% xlim([xmin xmax])
% ylim([ymin ymax])

%Now do the three windows
time_licks=[];
time_licks=handles_outs.time_licks;

%Before odor window
subplot(1,3,1)
hold on

time_mask=(time_licks>=-mean(delta_odor))&(time_licks<=0);
conv_lick=[];
conv_dFF=[];
for trNo=1:handles_outs.no_lick_slopes
    conv_lick=[conv_lick handles_outs.conv_lick(trNo,time_mask)];
    conv_dFF=[conv_dFF handles_outs.conv_dFF(trNo,time_mask)];
    plot(handles_outs.conv_lick(trNo,time_mask),handles_outs.conv_dFF(trNo,time_mask),'.b')
end


tbl= table(conv_lick',conv_dFF','VariableNames',{'conv_lick','conv_dFF'});
lm = fitlm(tbl,'conv_dFF~conv_lick');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(conv_lick',conv_dFF');
fprintf(1, ['rho and p value for dFF vs lick frequency before odor = %d, %d\n\n'],rho,pval);
handles_outs.dFF_rho(1)=rho;
handles_outs.dFF_pval(1)=pval;
xlabel('lick frequency Hz')
ylabel('dF/F')
title('Before odor')
xlim([xmin xmax])
ylim([ymin ymax])

%Odor window
subplot(1,3,2)
hold on


time_mask=(time_licks>=0)&(time_licks<=mean(delta_odor));
conv_lick=[];
conv_dFF=[];
for trNo=1:handles_outs.no_lick_slopes
    conv_lick=[conv_lick handles_outs.conv_lick(trNo,time_mask)];
    conv_dFF=[conv_dFF handles_outs.conv_dFF(trNo,time_mask)];
    plot(handles_outs.conv_lick(trNo,time_mask),handles_outs.conv_dFF(trNo,time_mask),'.b')
end


tbl= table(conv_lick',conv_dFF','VariableNames',{'conv_lick','conv_dFF'});
lm = fitlm(tbl,'conv_dFF~conv_lick');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(conv_lick',conv_dFF');
fprintf(1, ['rho and p value for dFF  vs lick frequency during odor = %d, %d\n\n'],rho,pval);
handles_outs.dFF_rho(2)=rho;
handles_outs.dFF_pval(2)=pval;
xlabel('lick frequency Hz')
ylabel('dF/F')
title('During odor')
xlim([xmin xmax])
ylim([ymin ymax])

%Reinforcement window
subplot(1,3,3)
hold on

time_mask=(time_licks>=mean(delta_odor_on_reinf_on))&(time_licks<=mean(delta_odor_on_reinf_on)+3);
conv_lick=[];
conv_dFF=[];
for trNo=1:handles_outs.no_lick_slopes
    conv_lick=[conv_lick handles_outs.conv_lick(trNo,time_mask)];
    conv_dFF=[conv_dFF handles_outs.conv_dFF(trNo,time_mask)];
    plot(handles_outs.conv_lick(trNo,time_mask),handles_outs.conv_dFF(trNo,time_mask),'.b')
end


tbl= table(conv_lick',conv_dFF','VariableNames',{'conv_lick','conv_dFF'});
lm = fitlm(tbl,'conv_dFF~conv_lick');
this_slope=lm.Coefficients{2,1};
this_intercept=lm.Coefficients{1,1};

plot([xmin xmax],this_slope*[xmin xmax]+this_intercept,'-r','LineWidth',2)

[rho,pval] = corr(conv_lick',conv_dFF');
fprintf(1, ['rho and p value for dFF vs lick frequency during reinforcement = %d, %d\n\n'],rho,pval);
handles_outs.dFF_rho(3)=rho;
handles_outs.dFF_pval(3)=pval;
xlabel('lick frequency Hz')
ylabel('dF/F')
title('Reinforcement')
xlim([xmin xmax])
ylim([ymin ymax])



save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) '_slopes.mat'],'handles_outs')


pffft=1;
