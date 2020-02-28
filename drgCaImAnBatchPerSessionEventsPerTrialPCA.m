function drgCaImAnBatchPerSessionEventsPerTrialPCA(choiceBatchPathName,choiceFileName)

%This code does a linear discriminant analysis for spm data
% Needs a choices file such as drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum
% Needs the output files from drgCaImAn_batch_dropc.m
warning('off')

close all
clear all
 
min_trials=20;

tic
 
if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAnBatchPerSessionReversalPerTrialPCA run for ' choiceFileName '\n\n']);

 

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

%Override if the user chooses a different min_trials
if isfield(caimanhandles.caimandr_choices,'min_trials')
    min_trials=caimanhandles.caimandr_choices.min_trials;
end

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
                num_comps_per_trial(num_odor_trials_dFF)=size(these_traces,1);
                mean_snip_dFF_per_comp(num_odor_trials_dFF,1:size(these_traces,1),1:no_time_points)=these_traces;
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
                num_comps_per_trial(num_odor_trials_dFF)=size(these_traces,1);
                mean_snip_dFF_per_comp(num_odor_trials_dFF,1:size(these_traces,1),1:no_time_points)=these_traces;
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
                num_comps_per_trial(num_odor_trials_dFF)=size(these_traces,1);
                mean_snip_dFF_per_comp(num_odor_trials_dFF,1:size(these_traces,1),1:no_time_points)=these_traces;
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
                mean_snip_dFF_per_comp(num_odor_trials_dFF,1:size(these_traces,1),1:no_time_points)=these_traces;
                num_comps_per_trial(num_odor_trials_dFF)=size(these_traces,1);
                mean_snip_dFF(num_odor_trials_dFF,1:no_time_points)=mean(these_traces,1);
                CI_snip_dFF(num_odor_trials_dFF,1:2,1:no_time_points)=bootci(1000, @mean, these_traces);
                time(num_odor_trials_dFF).time_to_event=handles_out.time_to_eventCR;
                
                epochs_per_trial_dFF(num_odor_trials_dFF)=4;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
  
    end
%     noROIs(filNum)=szhit(1);
noROIs(filNum)=szlit(2);
    
end

%Trim the time course for LDA
all_lda_input_timecourse=all_lda_input_timecourse(1:no_timepoints,:,:);
time_to_eventLDA=time_to_eventLDA(1,1:no_timepoints);

%Calculate percent correct
sliding_window=20; %Trials for determination of behavioral performance
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

% %For reversals plot violin plot of percent
% if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
%     %Keep track of the percent correct
%     %Plot percent correct vs trial
%     figNo=figNo+1;
%     try
%         close(figNo)
%     catch
%     end
%     
%     hFig1 = figure(figNo);
%     set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
%     hold on
%     
%     %Parameters for violin plot
%     edges=0:5:100;
%     rand_offset=0.8;
%     
%     %Trials before reversal
%     x_val=1;
%     handles_out2.pctPerWin(1).pct=perCorr(1:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1);
%     bar(x_val,mean(handles_out2.pctPerWin(1).pct),'FaceColor',[0.7 0.7 0.7])
%     [handles_out2.pctPerWin(1).pct_mean, handles_out2.pctPerWin(1).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(1).pct,edges,x_val,rand_offset,'k',2);
%     
%     %Trials after reversal
%     x_val=2;
%     handles_out2.pctPerWin(2).pct=perCorr(first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+15:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+65);
%     bar(x_val,mean(handles_out2.pctPerWin(2).pct),'FaceColor',[0.7 0.7 0.7])
%     [handles_out2.pctPerWin(2).pct_mean, handles_out2.pctPerWin(2).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(2).pct,edges,x_val,rand_offset,'k',2);
%     
%     %Trials at end
%     x_val=3;
%     handles_out2.pctPerWin(3).pct=perCorr(end-100:end);
%     bar(x_val,mean(handles_out2.pctPerWin(3).pct),'FaceColor',[0.7 0.7 0.7])
%     [handles_out2.pctPerWin(3).pct_mean, handles_out2.pctPerWin(3).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(3).pct,edges,x_val,rand_offset,'k',2);
%     
%     %Draw lines between points
%     plot([1 2 3],[mean(handles_out2.pctPerWin(1).pct) mean(handles_out2.pctPerWin(2).pct) mean(handles_out2.pctPerWin(3).pct)],'-k')
%     
%     ylim([0 120])
%     xlim([0.3 3.7])
%     ylabel('Percent correct')
%     xticks([1 2 3 5 6 7])
%     xticklabels({'Before','After','End'})
%     
%     
%     
%     window_labels{1}='Before';
%     window_labels{2}='After';
%     window_labels{3}='End';
%     
%     fprintf(1, ['\n\nranksum p values for percent correct windows\n\n']);
%     
%     no_pvals=0;
%     for ii=1:3
%         for jj=ii+1:3
%             no_pvals=no_pvals+1;
%             if (adtest(handles_out2.pctPerWin(ii).pct)==1)||(adtest(handles_out2.pctPerWin(ii).pct)==1)
%                 p_vals_corr(no_pvals)=ranksum(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
%                 fprintf(1, ['p values ranksum for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
%             else
%                 [h p_vals_corr(no_pvals)]=ttest2(handles_out2.pctPerWin(ii).pct,handles_out2.pctPerWin(jj).pct);
%                 fprintf(1, ['p values t test for ' window_labels{ii} ' vs. ' window_labels{jj} ' =%d\n'],p_vals_corr(no_pvals));
%             end
%             
%         end
%     end
%     
%     pFDRcorr=drsFDRpval(p_vals_corr);
%     fprintf(1, ['pFDR for significant difference percent correct  = %d\n\n'],pFDRcorr);
% end

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
    supertitle_description{1}='dF/F LDA analysis for trials with percent correct <65';
    supertitle_description{2}='dF/F LDA analysis for trials with percent correct >=65&<80';
    supertitle_description{3}='dF/F LDA analysis for trials with percent correct >=80';
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

for no_trial_windows=1:total_trial_windows
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
                        events_miss_FA(jj)=all_lda_events_miss_FA(ii);
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
        dt_lick=0.3;
        
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
        no_conv_points=3;
        conv_win=ones(1,no_conv_points);
        Splick_freq=conv(Splick_freq,conv_win,'same')/no_conv_points;
        Smlick_freq=conv(Smlick_freq,conv_win,'same')/no_conv_points;
        
        
        
        figNo=figNo+1
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        
        hold on
        
        lfreqmax=20;
        
        p1=plot(time_licks,Smlick_freq,'-b','LineWidth',2);
        p2=plot(time_licks,Splick_freq,'-r','LineWidth',2);
        
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
            
            figNo=figNo+1
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
        figNo=figNo+1
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
        
        
        
        fprintf(1,'For S+ vs S-\n')
        
        
        for time_point=1:length(time_to_eventLDA)
            handles_par(no_trial_windows).lda(time_point).test_out_per_timepoint=[];
            handles_par(no_trial_windows).lda(time_point).shuffled_out_per_timepoint=[];
            handles_par(no_trial_windows).lda(time_point).discriminant_correct=[];
            handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled=[];
            handles_par(no_trial_windows).lda(time_point).auROC=[];
            handles_par(no_trial_windows).lda(time_point).trials_processed_mask=[];
        end
        
        keep_lick=ones(1,length(dFF_trial_mask));
        for trial_no=1:num_odor_trials
            if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))>=lick_threshold
                keep_lick(trial_no)=0;
            end
        end
        vetted_mask=logical(dFF_trial_mask);
        Nall=sum(dFF_trial_mask&keep_lick);
        
        no_timepoints=length(time_to_eventLDA);
        for time_point=1:no_timepoints
            
            %dFF per trial per component
            measurements=zeros(Nall,max(all_lda_no_comp(vetted_mask)));
            this_all_lda_input_timecourse=zeros(Nall,max(all_lda_no_comp(vetted_mask)))';
            this_all_lda_input_timecourse(:,:)=all_lda_input_timecourse(time_point,1:max(all_lda_no_comp(vetted_mask)),vetted_mask);
            measurements(:,:)=this_all_lda_input_timecourse';
            
            which_file=[];
            which_file=all_lda_fileNo(vetted_mask);
            no_comps=[];
            no_comps=all_lda_no_comp(vetted_mask);
            
            
            
            %For each of the files do the PCA
            %         for ii=1:Nall
            par_t_out(time_point).trial_window(no_trial_windows).no_files_processed=0;
            for fileNo=unique(which_file)
                
                %If there are enough trials process the LDA
                N=sum(fileNo==which_file);
                
                
                if N>=min_trials
                    par_t_out(time_point).trial_window(no_trial_windows).no_files_processed=par_t_out(time_point).trial_window(no_trial_windows).no_files_processed+1;
                    par_t_out(time_point).trial_window(no_trial_windows).file_processed(par_t_out(time_point).trial_window(no_trial_windows).no_files_processed).principal_components=[];
                    %Partition the data into training and test set
                    
                    %Extract the measurements for this file
                    these_measurements=[];
                    these_measurements=measurements(which_file==fileNo,1:no_comps(find(which_file==fileNo,1,'first')));
                    
                    %What odor are these trials?
                    try
                        noEvs=0;
                        for jj=1:Nall
                            if (which_file(jj)==fileNo)
                                noEvs=noEvs+1;
                                if strcmp(events{jj},'S+')
                                    par_t_out(time_point).trial_window(no_trial_windows).file_processed(par_t_out(time_point).trial_window(no_trial_windows).no_files_processed).event(noEvs)=1;
                                else
                                    par_t_out(time_point).trial_window(no_trial_windows).file_processed(par_t_out(time_point).trial_window(no_trial_windows).no_files_processed).event(noEvs)=0;
                                end
                                par_t_out(time_point).trial_window(no_trial_windows).file_processed(par_t_out(time_point).trial_window(no_trial_windows).no_files_processed).event_miss_FA(noEvs)=events_miss_FA(jj);
                            end
                        end
                    catch
                    end
                    
                    [coeff,par_t_out(time_point).trial_window(no_trial_windows).file_processed(par_t_out(time_point).trial_window(no_trial_windows).no_files_processed).principal_components,latent]=pca(these_measurements);
                    
                end
            end
            
        end
        
        %Show the result of the PCA for each window
        if par_t_out(time_point).trial_window(no_trial_windows).no_files_processed>0
            for winNo=1:szwins(1)
                
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                
                hold on
                
                pca1min=-10;
                pca1max=15;
                
                
                win=(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_eventLDA<=caimanhandles.caimandr_choices.wins(winNo,2));
                
                for time_point=1:no_timepoints
                    
                    if win(time_point)==1
                        
                        %Plot these principal components
                        no_correct=0;
                        no_wrong=0;
                        for fileNo=1:par_t_out(time_point).trial_window(no_trial_windows).no_files_processed
                            these_pcs=[];
                            these_pcs=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).principal_components;
                            
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==1;
                            no_correct=no_correct+sum(these_events);
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'or')
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==2;
                            no_wrong=no_wrong+sum(these_events);
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'oc')
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==3;
                            no_correct=no_correct+sum(these_events);
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'ob')
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==4;
                            no_wrong=no_wrong+sum(these_events);
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'om')
                            
                            maxPC1=max([maxPC1 max(these_pcs(:,1))]);
                            maxPC2=max([maxPC2 max(these_pcs(:,2))]);
                            minPC1=min([minPC1 min(these_pcs(:,1))]);
                            minPC2=min([minPC2 min(these_pcs(:,2))]);
                        end
                        
                    end
                end
                
                
                xlabel('PC1')
                ylabel('PC2')
                xlim([-10 15])
                ylim([-6 8])
                title(['PCA for timepoints from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' secs for ' supertitle_description{no_trial_windows}])
                
                if (winNo==1)&(no_trial_windows==3)
                    pffft=1;
                end
            end
            
            %Show the result of the PCA for each window
            for winNo=1:szwins(1)
                
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                figure(figNo)
                
                
                hold on
                
                win=(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_eventLDA<=caimanhandles.caimandr_choices.wins(winNo,2));
                
                for time_point=1:no_timepoints
                    
                    if win(time_point)==1
                        
                        %Plot these principal components
                        for fileNo=1:par_t_out(time_point).trial_window(no_trial_windows).no_files_processed
                            these_pcs=[];
                            these_pcs=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).principal_components;
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event==1;
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'or','MarkerFaceColor','r')
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event==0;
                            plot(these_pcs(these_events,1),these_pcs(these_events,2),'ob','MarkerFaceColor','b')
                            maxPC1=max([maxPC1 max(these_pcs(:,1))]);
                            maxPC2=max([maxPC2 max(these_pcs(:,2))]);
                            minPC1=min([minPC1 min(these_pcs(:,1))]);
                            minPC2=min([minPC2 min(these_pcs(:,2))]);
                        end
                        
                    end
                end
                
                
                
                xlabel('PC1')
                ylabel('PC2')
                xlim([-10 15])
                ylim([-6 8])
                title(['PCA for timepoints from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' secs for ' supertitle_description{no_trial_windows}])
                
            end
            fprintf(1, ['percent correct = %d for trial window with ' trial_window_description{no_trial_windows} '\n'],100*no_correct/(no_correct+no_wrong));
            
            
            %Plot the timecourse for PCA1
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            figure(figNo)
            
            
            hold on
            
            
            for evNo=[4 2 3 1]
                %Find out how many trials for this event
                no_all_trials=0;
                for fileNo=1:par_t_out(1).trial_window(no_trial_windows).no_files_processed
                    these_events=[];
                    these_events=par_t_out(1).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==evNo;
                    no_all_trials=no_all_trials+sum(these_events);
                end
                
                if no_all_trials>0
                    PC1=zeros(no_timepoints,no_all_trials);
                    for time_point=1:no_timepoints
                        no_trials_ii=0;
                        for fileNo=1:par_t_out(time_point).trial_window(no_trial_windows).no_files_processed
                            these_events=[];
                            these_events=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).event_miss_FA==evNo;
                            these_pcs=[];
                            these_pcs=par_t_out(time_point).trial_window(no_trial_windows).file_processed(fileNo).principal_components;
                            PC1(time_point,no_trials_ii+1:no_trials_ii+sum(these_events))=these_pcs(these_events,1)';
                            no_trials_ii=no_trials_ii+sum(these_events);
                        end
                    end
                    
                    if no_all_trials>2
                        CI=[];
                        CI = bootci(1000, {@mean, PC1'})';
                        CI(:,1)=mean(PC1,2)-CI(:,1);
                        CI(:,2)=CI(:,2)-mean(PC1,2);
                        
                        switch evNo
                            case 1
                                [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(PC1,2), CI, 'r');
                            case 2
                                [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(PC1,2), CI, 'c');
                            case 3
                                [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(PC1,2), CI, 'b');
                            case 4
                                [hlCR, hpCR] = boundedline(time_to_eventLDA',mean(PC1,2), CI, 'm');
                        end
                    else
                        switch evNo
                            case 1
                                plot(time_to_eventLDA',mean(PC1,2),'r');
                            case 2
                                plot(time_to_eventLDA',mean(PC1,2), 'c');
                            case 3
                                plot(time_to_eventLDA',mean(PC1,2), 'b');
                            case 4
                                plot(time_to_eventLDA',mean(PC1,2),  'm');
                        end
                    end
                end
                
            end
            
            %Odor on markers
            plot([0 0],[pca1min pca1max],'-k')
            odorhl=plot([0 mean(delta_odor)],[pca1min + 0.1*(pca1max-pca1min) pca1min + 0.1*(pca1max-pca1min)],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],[pca1min pca1max],'-k')
            
            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pca1min pca1max],'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pca1min + 0.1*(pca1max-pca1min) pca1min + 0.1*(pca1max-pca1min)],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pca1min pca1max],'-r')
            
            
            xlabel('Time (sec)')
            ylabel('PC1')
            ylim([-10 15])
            title(['Timecourse for PCA1 for ' supertitle_description{no_trial_windows}])
            
        end
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
            
            these_dFF_trials=which_trials_in_PCA(events_miss_FA==evNo);
            
            
            if length(these_dFF_trials)>0
                
                this_mean=mean(mean_snip_dFF(these_dFF_trials,:),1)';
                this_mean=this_mean(1:length(time_to_eventLDA'),1);
                if length(these_dFF_trials)>2
                    CI=[];
                    CI = bootci(1000, {@mean, mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA'))})';
                    CI(:,1)=this_mean-CI(:,1);
                    CI(:,2)=CI(:,2)-this_mean;
                    
                    maxdFF=max([maxdFF,max(CI(:,1)+this_mean)]);
                    mindFF=min([mindFF,min(this_mean-CI(:,2))]);
                    
                    switch evNo
                        case 1
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'r');
                        case 2
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'c');
                        case 3
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'b');
                        case 4
                            [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'm');
                    end
                else
                    
                    switch evNo
                        case 1
                            plot(time_to_eventLDA',this_mean,'r');
                        case 2
                            plot(time_to_eventLDA',this_mean, 'c');
                        case 3
                            plot(time_to_eventLDA',this_mean, 'b');
                        case 4
                            plot(time_to_eventLDA',this_mean,  'm');
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
        
        %Plot the timecourse for dF/F for Hits and lick-excluded trials
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        figure(figNo)
        
        
        hold on
        
        maxdFF=-200;
        mindFF=200;
        
        evNo=1
        %Find out how many trials for this event
        
        these_dFF_trials=which_trials_in_PCA(events_miss_FA==evNo);
        
        this_mean=mean(mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA')))';
        
        if length(these_dFF_trials)>0
            
            
            if length(these_dFF_trials)>2
                CI=[];
                CI = bootci(1000, {@mean, mean_snip_dFF(these_dFF_trials,1:length(time_to_eventLDA'))})';
                CI(:,1)=this_mean-CI(:,1);
                CI(:,2)=CI(:,2)-this_mean;
                
                maxdFF=max([maxdFF,max(CI(:,1)+this_mean)]);
                mindFF=min([mindFF,min(this_mean-CI(:,2))]);
                
                switch evNo
                    case 1
                        [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'r');
                    case 2
                        [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'c');
                    case 3
                        [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'b');
                    case 4
                        [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'm');
                end
            else
                switch evNo
                    case 1
                        plot(time_to_eventLDA',this_mean,'r');
                    case 2
                        plot(time_to_eventLDA',this_mean, 'c');
                    case 3
                        plot(time_to_eventLDA',this_mean, 'b');
                    case 4
                        plot(time_to_eventLDA',this_mean,  'm');
                end
            end
        end
        
        
        
        %lick-excluded trials
        this_mean=mean(mean_snip_dFF(logical(lick_excluded_trials),1:length(time_to_eventLDA')),1)';
        
        if sum(lick_excluded_trials)>0
            
            
            if sum(lick_excluded_trials)>2
                CI=[];
                CI = bootci(1000, {@mean, mean_snip_dFF(logical(lick_excluded_trials),1:length(time_to_eventLDA'))})';
                CI(:,1)=this_mean-CI(:,1);
                CI(:,2)=CI(:,2)-this_mean;
                
                maxdFF=max([maxdFF,max(CI(:,1)+this_mean)]);
                mindFF=min([mindFF,min(this_mean-CI(:,2))]);
                
                
                [hlCR, hpCR] = boundedline(time_to_eventLDA',this_mean, CI, 'k');
                
            else
                
                plot(time_to_eventLDA',this_mean,'k');
                
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
        
        pffft=1;
        
        
    end
    pffft=1
end
 
% for fNo=firstFig:firstFig+2*szwins(1)-1
%     figure(fNo)
%     xlim([minPC1-0.1*(maxPC1-minPC1) maxPC1+0.1*(maxPC1-minPC1)])
%     ylim([minPC2-0.1*(maxPC2-minPC2) maxPC2+0.1*(maxPC2-minPC2)])
% end

pffft=1
