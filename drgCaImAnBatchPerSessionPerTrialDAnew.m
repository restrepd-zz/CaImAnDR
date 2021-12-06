function drgCaImAnBatchPerSessionPerTrialDAnew(choiceBatchPathName,choiceFileName)

% This code does a linear discriminant analysis for spm data for Figure 3 of
% Ma et al 2020
% Needs a choices file such as drgCaImAnChoicesDiego20180910_mmPVG04_Cerebellum
% Needs the output files from drgCaImAn_batch_dropc.m
warning('off')

close all
clear all
 
%Note: not all classifiers are enabled
classifier_names{1}='LDA';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='NN';
classifier_names{5}='Tree';


min_trials=20;
% min_trials=19;
 
ii_for_sig=0;
ii_for_sig_spike=0;

%Note that I am not using lick threshold here. I am setting it to 100
lick_threshold=100; %This is the threshold to exclude the runs where Ming was adding water manually
t_odor_on=0;
t_odor_off=4;

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

%Override if the user chooses a different min_trials
if isfield(caimanhandles.caimandr_choices,'min_trials')
    min_trials=caimanhandles.caimandr_choices.min_trials;
end

this_cost=caimanhandles.caimandr_choices.this_cost;
MLalgo=caimanhandles.caimandr_choices.MLalgo;
p_thr=caimanhandles.caimandr_choices.p_thr

gcp


%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;

all_lda_events=[]; 
all_lda_input_timecourse=[];
no_timepoints=2000000;

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

%For reversals plot violin plot of percent
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
    sgtitle(['Normalized dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])
end

%Do a linear discriminant analysis


if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
    %This is a forward run
    total_trial_windows=3;
    supertitle_description{1}='dF/F decoding analysis for trials with percent correct <65';
    supertitle_description{2}='dF/F decoding analysis for trials with percent correct >=65&<80';
    supertitle_description{3}='dF/F decoding analysis for trials with percent correct >=80';
else
    %Forward and reverse
    total_trial_windows=2;
    supertitle_description{1}='dF/F LDA analysis for trials before reversal';
    supertitle_description{2}='dF/F LDA analysis for trials after reversal at end of the session';
end

%Now calculate spikes
all_lda_spike_input_timecourse=zeros(size(all_lda_input_timecourse,1),size(all_lda_input_timecourse,2),size(all_lda_input_timecourse,3));
for ii_file=1:max(all_lda_fileNo)
    first_trNo=find(all_lda_fileNo==ii_file,1,'first');
    this_file_all_lda_input_timecourse=zeros(length(time_to_eventLDA),all_lda_no_comp(first_trNo),sum(all_lda_fileNo==ii_file));  
    this_file_all_lda_input_timecourse(:,:,:)=all_lda_input_timecourse(:,1:all_lda_no_comp(first_trNo),all_lda_fileNo==ii_file);
    
    %Threshold out the spikes
    std_traces=std(this_file_all_lda_input_timecourse(:));
    thr=5*std_traces;
    
    this_file_all_lda_spike_input_timecourse=zeros(length(time_to_eventLDA),all_lda_no_comp(first_trNo),sum(all_lda_fileNo==ii_file));
    for ii_comp=1:all_lda_no_comp(first_trNo)
        these_dFF_traces=zeros(length(time_to_eventLDA),sum(all_lda_fileNo==ii_file));
        these_dFF_spikes=zeros(length(time_to_eventLDA),sum(all_lda_fileNo==ii_file));
        these_dFF_traces(:,:)=this_file_all_lda_input_timecourse(:,ii_comp,:);
        baseline_dF=prctile(these_dFF_traces(:),5);
        these_dFF_spikes=these_dFF_traces>baseline_dF+thr;
        this_file_all_lda_spike_input_timecourse(:,ii_comp,:)=these_dFF_spikes;
    end
    
    all_lda_spike_input_timecourse(:,1:all_lda_no_comp(first_trNo),all_lda_fileNo==ii_file)=this_file_all_lda_spike_input_timecourse;
end
    


for no_trial_windows=1:total_trial_windows
    handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
    dFF_trial_mask=[];
    lick_excluded_trials=[];
    jj=0;
    events=[];
    
    if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
        
        fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
        pct_windows=[45 65;65 80;80 100.1];
        
        for ii=1:num_odor_trials_dFF
            if sum((lick_times(ii,1:no_licks(ii))>=t_odor_on)&(lick_times(ii,1:no_licks(ii))<=t_odor_off))<lick_threshold
                lick_excluded_trials(ii)=0;
                if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
                    dFF_trial_mask(ii)=1;
                    jj=jj+1;
                    handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
                    events{jj,1}=all_lda_events{ii};
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
    
    Nall=sum(dFF_trial_mask);
    
    if sum(dFF_trial_mask)>0
        
        %Do lick analysis
        %Plot the lick frequency for S+ and S-
        dt_lick=0.3;
        Splick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        Smlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        sp_trno=0;
        sm_trno=0;
        
        for trial_no=1:num_odor_trials
            if dFF_trial_mask(trial_no)==1
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
            end
        end
        
        Splick_freq=(Splick_freq/(sp_trno*dt_lick));
        Smlick_freq=(Smlick_freq/(sm_trno*dt_lick));
        
        %Convolve lick_freq using a window of 0.9 sec
        no_conv_points=3;
        conv_win=ones(1,no_conv_points);
        Splick_freq=conv(Splick_freq,conv_win,'same')/no_conv_points;
        Smlick_freq=conv(Smlick_freq,conv_win,'same')/no_conv_points;
        
        figNo=figNo+1;
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
        
        
        fprintf(1,'For S+ vs S-\n')
        
        
        for time_point=1:length(time_to_eventLDA)
            handles_par(no_trial_windows).lda(time_point).test_out_per_timepoint=[];
            handles_par(no_trial_windows).lda(time_point).shuffled_out_per_timepoint=[];
            handles_par(no_trial_windows).lda(time_point).discriminant_correct=[];
            handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled=[];
            handles_par(no_trial_windows).lda(time_point).auROC=[];
            handles_par(no_trial_windows).lda(time_point).trials_processed_mask=[];
        end
        
        no_timepoints=length(time_to_eventLDA);
        
        %Masks for S+ and S-
        which_file=[];
        which_file=all_lda_fileNo(logical(dFF_trial_mask));
        
        these_files=[min(which_file):max(which_file)];
        Sp_mask=zeros(max(which_file),Nall);
        Sm_mask=zeros(max(which_file),Nall);
        allTr_mask=zeros(max(which_file),Nall);

        for ww=1:length(events)
            if strcmp(events{ww},'S+')
                Sp_mask(which_file(ww),ww)=1;
            else
                Sm_mask(which_file(ww),ww)=1;
            end
            allTr_mask(which_file(ww),ww)=1;
        end
        
        no_comps=[];
        no_comps=all_lda_no_comp(logical(dFF_trial_mask));
        
       
        
        %Make en exclussion mask based on differences between S+ and S-
        comps_included=zeros(max(which_file),max(all_lda_no_comp(logical(dFF_trial_mask))));
        for ii_file=these_files
            all_t_this_all_lda_input_timecourse=zeros(no_timepoints,max(all_lda_no_comp(logical(dFF_trial_mask))),Nall);
            all_t_this_all_lda_input_timecourse(:,:,:)=all_lda_input_timecourse(:,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
            for ii_comp=1:size(all_t_this_all_lda_input_timecourse,2)
                
                sp_dFF=zeros(sum(time_to_eventLDA>10),sum(Sp_mask(ii_file,:)));
                this_Sp_mask=zeros(1,Nall);
                this_Sp_mask(1,:)=Sp_mask(ii_file,:);
                sp_dFF(:,:)=all_t_this_all_lda_input_timecourse(time_to_eventLDA>10,ii_comp,logical(this_Sp_mask));
                
                sm_dFF=zeros(sum(time_to_eventLDA>10),sum(Sm_mask(ii_file,:)));
                this_Sm_mask=zeros(1,Nall);
                this_Sm_mask(1,:)=Sm_mask(ii_file,:);
                sm_dFF(:,:)=all_t_this_all_lda_input_timecourse(time_to_eventLDA>10,ii_comp,logical(this_Sm_mask));
                if (sum(this_Sp_mask>0))&(sum(this_Sm_mask>0))
                    [h,p]=ttest2(sm_dFF(:),sp_dFF(:));
                    
                    if p<p_thr
                        comps_included(ii_file,ii_comp)=1;
                    end
                end
            end

            ii=find(this_Sp_mask==1,1,'first');
            fprintf(1,'For file No %d the number of comps included was %d out of %d\n', ii_file,sum(comps_included(ii_file,:)),no_comps(ii))
        end
        
        
        for ii_file=these_files
            %Get the splus and sminus traces for each file
            %Note that we do not align the components across sessions
            %Also, calculate the baseline for each component and the threshold per file
            
            this_Sp_mask=zeros(1,Nall);
            this_Sp_mask(1,:)=Sp_mask(ii_file,:);
            
            this_Sm_mask=zeros(1,Nall);
            this_Sm_mask(1,:)=Sm_mask(ii_file,:);
            
            these_comps_included=zeros(1,size(comps_included,2));
            these_comps_included(1,:)=comps_included(ii_file,:);
            
            szSp=size(splus_traces);
            szSm=size(sminus_traces);
            %time_to_event=([1:szSm(2)]*dt-dt_before);
            time_to_eventSm=([1:szSm(2)]*dt-dt_before);
            time_to_eventSp=([1:szSp(2)]*dt-dt_before);
            
            if ((sum(this_Sp_mask==1)>4)&(sum(this_Sm_mask==1)>4)&(sum(these_comps_included==1)>2))
                
                these_sminus_traces=zeros(sum(these_comps_included==1)*sum(this_Sm_mask==1),length(time_to_eventSm));
                these_splus_traces=zeros(sum(these_comps_included==1)*sum(this_Sp_mask==1),length(time_to_eventSp));
                
                
                tic
                for time_point=1:no_timepoints
                    measurements=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))));
                    this_all_lda_input_timecourse=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))))';
                    this_all_lda_input_timecourse(:,:)=all_lda_input_timecourse(time_point,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
                    measurements(:,:)=this_all_lda_input_timecourse';
                    
                    %Enter S+
                    ii_sp=0;
                    for jj=1:length(this_Sp_mask)
                        if this_Sp_mask(jj)==1
                            these_splus_traces(ii_sp+1:ii_sp+sum(these_comps_included==1),time_point)=measurements(jj,logical(these_comps_included));
                            ii_sp=ii_sp+sum(these_comps_included==1);
                        end
                    end
                    
                    %Enter S-
                    ii_sm=0;
                    for jj=1:length(this_Sm_mask)
                        if this_Sm_mask(jj)==1
                            these_sminus_traces(ii_sm+1:ii_sm+sum(these_comps_included==1),time_point)=measurements(jj,logical(these_comps_included));
                            ii_sm=ii_sm+sum(these_comps_included==1);
                        end
                    end
                end
                toc
                
                tic
                %Show the average calcium dynamics for the cells chosen
                %S+, S-, all snips
                CIsm = bootci(1000, @mean, these_sminus_traces);
                meansm=mean(these_sminus_traces,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                CIsp = bootci(1000, @mean, these_splus_traces);
                meansp=mean(these_splus_traces,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;
                
                
                %First plot the average Splus and Sminus
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                hFig = figure(figNo);
                
                hold on
                
                
                
                pct1=prctile([mean(these_sminus_traces,1)'; mean(these_splus_traces(:,1:szSp(2)),1)'],1);
                pct99=prctile([mean(these_sminus_traces,1)'; mean(these_splus_traces(:,1:szSp(2)),1)'],99);
                
                
                
                [hlsm, hpsm] = boundedline(time_to_eventSm',mean(these_sminus_traces,1)', CIsm', 'b');
                [hlsp, hpsp] = boundedline(time_to_eventSp',mean(these_splus_traces,1)', CIsp', 'r');
                
                %Odor on markers
                plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
                odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
                plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
                
                switch no_trial_windows
                    case 1
                        title(['Ca changes for percent <65%% and file ' num2str(ii_file)])
                    case 2
                        title(['Ca changes for percent >65%% and <80%% and file ' num2str(ii_file)])
                    case 3
                        title(['Ca changes for percent >80%% and file ' num2str(ii_file)])
                end
                
                legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
                xlabel('Time (sec)')
                ylabel('dF/F')
                ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
                xlim([-10 19.8])
                toc
                
                %Calculate p values
                p_dFF=[];
                for ii=1:min([length(time_to_eventSm) length(time_to_eventSp)])
                    [h p_dFF(ii)]=ttest2(these_sminus_traces(:,ii),these_splus_traces(:,ii));
                end
                
                if length(time_to_eventSm)<length(time_to_eventSp)
                    time_p=time_to_eventSm;
                else
                    time_p=time_to_eventSp; 
                end
                logpmin=-20;
                
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                figure(figNo)
                plot(time_p,log10(p_dFF),'-k','LineWidth',2)
                hold on
                plot([time_p(1) time_p(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
                
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
                        title(['log(p value) for dFF Sp vs Sm for percent <65, file ' num2str(ii_file)]);
                    case 2
                        title(['log(p value) for dFF Sp vs Sm for percent  >=65&<80, file ' num2str(ii_file)]);
                    case 3
                        title(['log(p value) for dFF Sp vs Sm for percent  >80, file ' num2str(ii_file)]);
                end
                
                xlabel('Time (sec)')
                ylabel('log10(p value)')
                ylim([logpmin 0.5])
                xlim([-10 20])
            end
        end
        
        %Now do discriminant analysis
        for time_point=1:no_timepoints
            
            scores=[];
            trials_processed=[];
            correct_predict=[];
            correct_predict_shuffled=[];
            lda_no_trials=[];
            
            
            if time_point==1
                maxN=0;
                for ii=1:Nall
                    if sum(which_file(ii)==which_file)>maxN
                        maxN=sum(which_file(ii)==which_file);
                    end
                end
                
                
                switch no_trial_windows
                    case 1
                        fprintf(1,'\n\nFor percent <65%% max number of trials per file= %d\n\n',maxN)
                    case 2
                        fprintf(1,'\n\nFor percent >65%% and <80%% max number of trials per file= %d\n\n',maxN)
                    case 3
                        fprintf(1,'\n\nFor percent >80%% max number of trials per file= %d\n\n',maxN)
                end
                
              
            end
            
            measurements=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))));
            this_all_lda_input_timecourse=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))))';
            this_all_lda_input_timecourse(:,:)=all_lda_input_timecourse(time_point,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
            measurements(:,:)=this_all_lda_input_timecourse';
            
            parfor ii=1:Nall
                
                %If there are enough trials process the LDA
                N=sum(which_file(ii)==which_file);
                these_comps_included=zeros(1,size(comps_included,2));
                these_comps_included(1,:)=comps_included(which_file(ii),:);
                
                if (N>=min_trials)&(sum(these_comps_included)>2)
                    trials_processed(ii)=1;
                    %Partition the data into training and test sets.
                    
                    %Create input and target vectors leaving one trial out
                    %For per_input each column has the dF/F for one trial
                    %each row is a single time point for dF/F for one of the cells
                    %For per_target the top row is 1 if the odor is S+ and 0 if it is
                    %S-, and row 2 has 1 for S-
                    idxTrn=ones(Nall,1)&(which_file==which_file(ii))';
                    idxTrn(ii)=0;
                    idxTest=zeros(Nall,1);
                    idxTest(ii)=1;
                    
                    %Store the training data in a table.
                    tblTrn=[];
%                     tblTrn = array2table(measurements(logical(idxTrn),1:no_comps(ii)));
                    tblTrn = array2table(measurements(logical(idxTrn),logical(these_comps_included)));
                    these_events=[];
                    noEvs=0;
                    all_these_events=[];
                    noallEvs=0;
                    for jj=1:Nall
                        if (which_file(jj)==which_file(ii))&(jj~=ii)
                            noEvs=noEvs+1;
                            these_events{noEvs}=events{jj};
                        end
                        if (which_file(jj)==which_file(ii))
                            noallEvs=noallEvs+1;
                            all_these_events{noallEvs}=events{jj};
                            if jj==ii
                                ii_event=noallEvs;
                            end
                        end
                    end
                    tblTrn.Y = these_events';
                     
                    %Train a discriminant analysis model using the training set and default options.
                    %By default this is a regularized linear discriminant analysis (LDA)
%                     Mdl = fitcdiscr(tblTrn,'Y');
                    
                    switch MLalgo
                        case 1
                            Mdl = fitcdiscr(tblTrn,'Y','Cost',this_cost);
                        case 2
                            Mdl = fitcsvm(tblTrn,'Y','Cost',this_cost);
                        case 3
                            Mdl = fitcnb(tblTrn,'Y','Cost',this_cost);
                        case 4
                            Mdl = fitcnet(tblTrn,'Y');
                        case 5
                            Mdl = fitctree(tblTrn,'Y','Cost',this_cost);
                    end
                    
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %                     [label,score] = predict(Mdl,measurements(logical(idxTest),1:no_comps(ii)));
                    [label,score] = predict(Mdl,measurements(logical(idxTest),logical(these_comps_included)));
                    
                    %label is the predicted label, and score is the predicted class
                    %posterior probability
                    scores(ii)=score(1);
                    lda_no_trials(ii)=N;
                    correct_predict(ii)=strcmp(events{ii},label);
                    
                    ii_shuffled=randperm(N);
                    correct_predict_shuffled(ii)=strcmp(all_these_events{ii_shuffled(ii_event)},label);
                    
                else
                    trials_processed(ii)=0;
                end
            end
            
            if sum(trials_processed)>=1
                handles_par(no_trial_windows).lda(time_point).trials_processed_mask=logical(trials_processed);
                
                events_processed=[];
                jj=0;
                for ii=1:Nall
                    if trials_processed(ii)==1
                        jj=jj+1;
                        events_processed{jj}=events{ii};
                    end
                end
                
                trials_processed=logical(trials_processed);
                scores_processed=[];
                scores_processed=scores(trials_processed);
                
                %Calculate auROC
                [X,Y,T,AUC] = perfcurve(events_processed,scores_processed','S+');
                
                handles_par(no_trial_windows).lda(time_point).auROC=AUC-0.5;
                handles_par(no_trial_windows).lda(time_point).discriminant_correct=100*sum(correct_predict(trials_processed))/sum(trials_processed);
                handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled=100*sum(correct_predict_shuffled(trials_processed))/sum(trials_processed);
                handles_par(no_trial_windows).lda(time_point).N=lda_no_trials(trials_processed);
                fprintf(1, 'LDA percent correct classification %d (for timepoint %d out of %d)\n',100*sum(correct_predict(trials_processed))/sum(trials_processed),time_point,no_timepoints);
                
            else
                fprintf(1, 'LDA percent correct classification not calculated because there were not enough trials in this file for timepoint %d out of %d\n',time_point,no_timepoints);
            end
        end
        
        %Plot the timecourse
        if sum(trials_processed)>=1
            
            
            discriminant_correct=zeros(1,no_timepoints);
            discriminant_correct_shuffled=zeros(1,no_timepoints);
            auROC=zeros(1,no_timepoints);
            
            for time_point=1:no_timepoints
                discriminant_correct(time_point)= handles_par(no_trial_windows).lda(time_point).discriminant_correct;
                discriminant_correct_shuffled(time_point)= handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled;
                auROC(time_point)=handles_par(no_trial_windows).lda(time_point).auROC;
            end
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            figure(figNo)
            
            %         subplot(1,2,1)
            hold on
            
            
            per95=prctile(discriminant_correct_shuffled(1,:),95);
            per5=prctile(discriminant_correct_shuffled(1,:),5);
            CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
            [hlCR, hpCR] = boundedline([time_to_eventLDA(1) time_to_eventLDA(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
            
            plot(time_to_eventLDA',discriminant_correct(1,:),'-k')
            
            %Odor on markers
            plot([0 0],[0 110],'-k')
            odorhl=plot([0 mean(delta_odor)],[32 32],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],[0 110],'-k')
            
            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 110],'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[32 32],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 110],'-r')
            
            
            ylim([00 110])
            
            xlabel('Time (sec)')
            ylabel('Percent correct')
            title(supertitle_description{no_trial_windows})
            
            
            %Save data for the tests of significance for the LDA
            for winNo=1:szwins(1)
                
                %discriminant correct within windows
                ii_for_sig=ii_for_sig+1;
                
                win=(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_eventLDA<=caimanhandles.caimandr_choices.wins(winNo,2));
                handles_sig.win(ii_for_sig).discriminant_correct=discriminant_correct(win);
                if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                    switch no_trial_windows
                        case 1
                            handles_sig.win(ii_for_sig).description=['<65 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                        case 2
                            handles_sig.win(ii_for_sig).description=['>=65&<80 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                        case 3
                            handles_sig.win(ii_for_sig).description=['>=80 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    end
                else
                    if no_trial_windows==1
                        handles_sig.win(ii_for_sig).description=['before reversal from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    else
                        handles_sig.win(ii_for_sig).description=['after reversal from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    end
                end
            end
            
            %Save the shuffled data
            %Note that we use all the time points for the shuffled data on
            %purpose!
            ii_for_sig=ii_for_sig+1;
            handles_sig.win(ii_for_sig).discriminant_correct=discriminant_correct_shuffled;
            handles_sig.win(ii_for_sig).win=winNo;
            handles_sig.win(ii_for_sig).pct=no_trial_windows;
            if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                switch no_trial_windows
                    case 1
                        handles_sig.win(ii_for_sig).description=['<65 shuffled'];
                    case 2
                        handles_sig.win(ii_for_sig).description=['>=65&<80 shuffled'];
                    case 3
                        handles_sig.win(ii_for_sig).description=['>=80 shuffled'];
                end
            else
                if no_trial_windows==1
                    handles_sig.win(ii_for_sig).description=['before reversal (shuffled)'];
                else
                    handles_sig.win(ii_for_sig).description=['after reversal (shuffled)'];
                end
            end
            
            
            handles_sig.ii_for_sig=ii_for_sig;
            
            
        end
        
        %Now do spikes
        %all_lda_spike_input_timecourse
        for ii_file=these_files
            %Get the splus and sminus traces for each file
            %Note that we do not align the components across sessions
            %Also, calculate the baseline for each component and the threshold per file
            
            this_Sp_mask=zeros(1,Nall);
            this_Sp_mask(1,:)=Sp_mask(ii_file,:);
            
            this_Sm_mask=zeros(1,Nall);
            this_Sm_mask(1,:)=Sm_mask(ii_file,:);
            
            these_comps_included=zeros(1,size(comps_included,2));
            these_comps_included(1,:)=comps_included(ii_file,:);
            
            szSp=size(splus_traces);
            szSm=size(sminus_traces);
            %time_to_event=([1:szSm(2)]*dt-dt_before);
            time_to_eventSm=([1:szSm(2)]*dt-dt_before);
            time_to_eventSp=([1:szSp(2)]*dt-dt_before);
            
            if ((sum(this_Sp_mask==1)>4)&(sum(this_Sm_mask==1)>4)&(sum(these_comps_included==1)>2))
                
                these_sminus_traces=zeros(sum(these_comps_included==1)*sum(this_Sm_mask==1),length(time_to_eventSm));
                these_splus_traces=zeros(sum(these_comps_included==1)*sum(this_Sp_mask==1),length(time_to_eventSp));
                
                
                tic
                for time_point=1:no_timepoints
                    measurements=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))));
                    this_all_lda_spike_input_timecourse=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))))';
                    this_all_lda_spike_input_timecourse(:,:)=all_lda_spike_input_timecourse(time_point,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
                    measurements(:,:)=this_all_lda_spike_input_timecourse';
                    
                    %Enter S+
                    ii_sp=0;
                    for jj=1:length(this_Sp_mask)
                        if this_Sp_mask(jj)==1
                            these_splus_traces(ii_sp+1:ii_sp+sum(these_comps_included==1),time_point)=measurements(jj,logical(these_comps_included));
                            ii_sp=ii_sp+sum(these_comps_included==1);
                        end
                    end
                    
                    %Enter S-
                    ii_sm=0;
                    for jj=1:length(this_Sm_mask)
                        if this_Sm_mask(jj)==1
                            these_sminus_traces(ii_sm+1:ii_sm+sum(these_comps_included==1),time_point)=measurements(jj,logical(these_comps_included));
                            ii_sm=ii_sm+sum(these_comps_included==1);
                        end
                    end
                end
                toc
                
                tic
                %Show the average calcium dynamics for the cells chosen
                %S+, S-, all snips
                CIsm = bootci(1000, @mean, these_sminus_traces);
                meansm=mean(these_sminus_traces,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                CIsp = bootci(1000, @mean, these_splus_traces);
                meansp=mean(these_splus_traces,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;
                
                
                %First plot the average Splus and Sminus
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                hFig = figure(figNo);
                
                hold on
                
                
                
                pct1=prctile([mean(these_sminus_traces,1)'; mean(these_splus_traces(:,1:szSp(2)),1)'],1);
                pct99=prctile([mean(these_sminus_traces,1)'; mean(these_splus_traces(:,1:szSp(2)),1)'],99);
                
                
                
                [hlsm, hpsm] = boundedline(time_to_eventSm',mean(these_sminus_traces,1)', CIsm', 'b');
                [hlsp, hpsp] = boundedline(time_to_eventSp',mean(these_splus_traces,1)', CIsp', 'r');
                
                %Odor on markers
                plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
                odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
                plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
                
                %Reinforcement markers
                plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
                reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
                plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
                
                switch no_trial_windows
                    case 1
                        title(['Ca spikes for percent <65%% and file ' num2str(ii_file)])
                    case 2
                        title(['Ca spikes for percent >65%% and <80%% and file ' num2str(ii_file)])
                    case 3
                        title(['Ca spikes for percent >80%% and file ' num2str(ii_file)])
                end
                
                legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
                xlabel('Time (sec)')
                ylabel('dF/F')
                ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
                xlim([-10 19.8])
                toc
                
                %Calculate p values
                p_dFF=[];
                for ii=1:min([length(time_to_eventSm) length(time_to_eventSp)])
                    [h p_dFF(ii)]=ttest2(these_sminus_traces(:,ii),these_splus_traces(:,ii));
                end
                
                if length(time_to_eventSm)<length(time_to_eventSp)
                    time_p=time_to_eventSm;
                else
                    time_p=time_to_eventSp; 
                end
                logpmin=-20;
                
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                figure(figNo)
                plot(time_p,log10(p_dFF),'-k','LineWidth',2)
                hold on
                plot([time_p(1) time_p(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
                
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
                        title(['log(p value) for dFF spikes Sp vs Sm for percent <65, file ' num2str(ii_file)]);
                    case 2
                        title(['log(p value) for dFF spikes Sp vs Sm for percent  >=65&<80, file ' num2str(ii_file)]);
                    case 3
                        title(['log(p value) for dFF spikes Sp vs Sm for percent  >80, file ' num2str(ii_file)]);
                end
                
                xlabel('Time (sec)')
                ylabel('log10(p value)')
                ylim([logpmin 0.5])
                xlim([-10 20])
            end
        end
        
        %Now do discriminant analysis
        for time_point=1:no_timepoints
            
            scores=[];
            trials_processed=[];
            correct_predict=[];
            correct_predict_shuffled=[];
            lda_no_trials=[];
            
            
            if time_point==1
                maxN=0;
                for ii=1:Nall
                    if sum(which_file(ii)==which_file)>maxN
                        maxN=sum(which_file(ii)==which_file);
                    end
                end
                
                
                switch no_trial_windows
                    case 1
                        fprintf(1,'\n\nFor percent <65%% max number of trials per file= %d\n\n',maxN)
                    case 2
                        fprintf(1,'\n\nFor percent >65%% and <80%% max number of trials per file= %d\n\n',maxN)
                    case 3
                        fprintf(1,'\n\nFor percent >80%% max number of trials per file= %d\n\n',maxN)
                end
                
              
            end
            
            measurements=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))));
            this_all_lda_spike_input_timecourse=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))))';
            this_all_lda_spike_input_timecourse(:,:)=all_lda_spike_input_timecourse(time_point,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
            measurements(:,:)=this_all_lda_spike_input_timecourse';
            
            parfor ii=1:Nall
                
                %If there are enough trials process the LDA
                N=sum(which_file(ii)==which_file);
                these_comps_included=zeros(1,size(comps_included,2));
                these_comps_included(1,:)=comps_included(which_file(ii),:);
                
                if (N>=min_trials)&(sum(these_comps_included)>2)
                    trials_processed(ii)=1;
                    %Partition the data into training and test sets.
                    
                    %Create input and target vectors leaving one trial out
                    %For per_input each column has the dF/F for one trial
                    %each row is a single time point for dF/F for one of the cells
                    %For per_target the top row is 1 if the odor is S+ and 0 if it is
                    %S-, and row 2 has 1 for S-
                    idxTrn=ones(Nall,1)&(which_file==which_file(ii))';
                    idxTrn(ii)=0;
                    idxTest=zeros(Nall,1);
                    idxTest(ii)=1;
                    
                    %Store the training data in a table.
                    tblTrn=[];
%                     tblTrn = array2table(measurements(logical(idxTrn),1:no_comps(ii)));
                    tblTrn = array2table(measurements(logical(idxTrn),logical(these_comps_included)));
                    these_events=[];
                    noEvs=0;
                    all_these_events=[];
                    noallEvs=0;
                    for jj=1:Nall
                        if (which_file(jj)==which_file(ii))&(jj~=ii)
                            noEvs=noEvs+1;
                            these_events{noEvs}=events{jj};
                        end
                        if (which_file(jj)==which_file(ii))
                            noallEvs=noallEvs+1;
                            all_these_events{noallEvs}=events{jj};
                            if jj==ii
                                ii_event=noallEvs;
                            end
                        end
                    end
                    tblTrn.Y = these_events';
                     
                    %Train a discriminant analysis model using the training set and default options.
                    %By default this is a regularized linear discriminant analysis (LDA)
%                     Mdl = fitcdiscr(tblTrn,'Y');
                    
                    switch MLalgo
                        case 1
                            Mdl = fitcdiscr(tblTrn,'Y','Cost',this_cost);
                        case 2
                            Mdl = fitcsvm(tblTrn,'Y','Cost',this_cost);
                        case 3
                            Mdl = fitcnb(tblTrn,'Y','Cost',this_cost);
                        case 4
                            Mdl = fitcnet(tblTrn,'Y');
                        case 5
                            Mdl = fitctree(tblTrn,'Y','Cost',this_cost);
                    end
                    
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %                     [label,score] = predict(Mdl,measurements(logical(idxTest),1:no_comps(ii)));
                    [label,score] = predict(Mdl,measurements(logical(idxTest),logical(these_comps_included)));
                    
                    %label is the predicted label, and score is the predicted class
                    %posterior probability
                    scores(ii)=score(1);
                    lda_no_trials(ii)=N;
                    correct_predict(ii)=strcmp(events{ii},label);
                    
                    ii_shuffled=randperm(N);
                    correct_predict_shuffled(ii)=strcmp(all_these_events{ii_shuffled(ii_event)},label);
                    
                else
                    trials_processed(ii)=0;
                end
            end
            
            if sum(trials_processed)>=1
                handles_par(no_trial_windows).lda(time_point).trials_processed_mask=logical(trials_processed);
                
                events_processed=[];
                jj=0;
                for ii=1:Nall
                    if trials_processed(ii)==1
                        jj=jj+1;
                        events_processed{jj}=events{ii};
                    end
                end
                
                trials_processed=logical(trials_processed);
                scores_processed=[];
                scores_processed=scores(trials_processed);
                
                %Calculate auROC
                [X,Y,T,AUC] = perfcurve(events_processed,scores_processed','S+');
                
                handles_par(no_trial_windows).lda(time_point).auROC=AUC-0.5;
                handles_par(no_trial_windows).lda(time_point).discriminant_correct=100*sum(correct_predict(trials_processed))/sum(trials_processed);
                handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled=100*sum(correct_predict_shuffled(trials_processed))/sum(trials_processed);
                handles_par(no_trial_windows).lda(time_point).N=lda_no_trials(trials_processed);
                fprintf(1, 'Decoding accuracy for dFF spikes %d (for timepoint %d out of %d)\n',100*sum(correct_predict(trials_processed))/sum(trials_processed),time_point,no_timepoints);
                
            else
                fprintf(1, 'Decoding accuracy for dFF spikes not calculated because there were not enough trials in this file for timepoint %d out of %d\n',time_point,no_timepoints);
            end
        end
        
        %Plot the timecourse
        if sum(trials_processed)>=1
            
            
            discriminant_correct=zeros(1,no_timepoints);
            discriminant_correct_shuffled=zeros(1,no_timepoints);
            auROC=zeros(1,no_timepoints);
            
            for time_point=1:no_timepoints
                discriminant_correct(time_point)= handles_par(no_trial_windows).lda(time_point).discriminant_correct;
                discriminant_correct_shuffled(time_point)= handles_par(no_trial_windows).lda(time_point).discriminant_correct_shuffled;
                auROC(time_point)=handles_par(no_trial_windows).lda(time_point).auROC;
            end
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            figure(figNo)
            
            %         subplot(1,2,1)
            hold on
            
            
            per95=prctile(discriminant_correct_shuffled(1,:),95);
            per5=prctile(discriminant_correct_shuffled(1,:),5);
            CIsh=[mean(discriminant_correct_shuffled(1,:))-per5 per95-mean(discriminant_correct_shuffled(1,:))]';
            [hlCR, hpCR] = boundedline([time_to_eventLDA(1) time_to_eventLDA(end)],[mean(discriminant_correct_shuffled(1,:)) mean(discriminant_correct_shuffled(1,:))], CIsh', 'r');
            
            plot(time_to_eventLDA',discriminant_correct(1,:),'-k')
            
            %Odor on markers
            plot([0 0],[0 110],'-k')
            odorhl=plot([0 mean(delta_odor)],[32 32],'-k','LineWidth',5);
            plot([mean(delta_odor) mean(delta_odor)],[0 110],'-k')
            
            %Reinforcement markers
            plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 110],'-r')
            reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[32 32],'-r','LineWidth',5);
            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 110],'-r')
            
            
            ylim([00 110])
            
            xlabel('Time (sec)')
            ylabel('Percent correct')
            title(['Spike' supertitle_description{no_trial_windows}])
            
            
            %Save data for the tests of significance for the LDA
            for winNo=1:szwins(1)
                
                %discriminant correct within windows
                ii_for_sig_spike=ii_for_sig_spike+1;
                
                win=(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_eventLDA<=caimanhandles.caimandr_choices.wins(winNo,2));
                handles_sig_spike.win(ii_for_sig_spike).discriminant_correct=discriminant_correct(win);
                if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                    switch no_trial_windows
                        case 1
                            handles_sig_spike.win(ii_for_sig_spike).description=['<65 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                        case 2
                            handles_sig_spike.win(ii_for_sig_spike).description=['>=65&<80 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                        case 3
                            handles_sig_spike.win(ii_for_sig_spike).description=['>=80 from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    end
                else
                    if no_trial_windows==1
                        handles_sig_spike.win(ii_for_sig_spike).description=['before reversal from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    else
                        handles_sig_spike.win(ii_for_sig_spike).description=['after reversal from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec'];
                    end
                end
            end
            
            %Save the shuffled data
            %Note that we use all the time points for the shuffled data on
            %purpose!
            ii_for_sig_spike=ii_for_sig_spike+1;
            handles_sig_spike.win(ii_for_sig_spike).discriminant_correct=discriminant_correct_shuffled;
            handles_sig_spike.win(ii_for_sig_spike).win=winNo;
            handles_sig_spike.win(ii_for_sig_spike).pct=no_trial_windows;
            if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                switch no_trial_windows
                    case 1
                        handles_sig_spike.win(ii_for_sig_spike).description=['<65 shuffled'];
                    case 2
                        handles_sig_spike.win(ii_for_sig_spike).description=['>=65&<80 shuffled'];
                    case 3
                        handles_sig_spike.win(ii_for_sig_spike).description=['>=80 shuffled'];
                end
            else
                if no_trial_windows==1
                    handles_sig_spike.win(ii_for_sig_spike).description=['before reversal (shuffled)'];
                else
                    handles_sig_spike.win(ii_for_sig_spike).description=['after reversal (shuffled)'];
                end
            end
            
            
            handles_sig_spike.ii_for_sig_spike=ii_for_sig_spike;
            
            
        end
        
    end %end for if sum(dFF_trial_mask)>0
end %end for for no_trial_windows =

%Do test of significance for the LDA

fprintf(1, 'Tests of significance for difference in percent correct LDA\n')
fprintf(1, [choiceFileName '\n']);
fprintf(1, ['Note: For shuffled trials we use all time points\n\n'])
p_vals_LDA=[];
no_LDA_pvals=0;
for ii=1:ii_for_sig
    for jj=ii+1:ii_for_sig
        
        no_LDA_pvals=no_LDA_pvals+1;
        if (length(handles_sig.win(ii).discriminant_correct)<4)||(length(handles_sig.win(jj).discriminant_correct)<4)
            %adtest does not work with n<4. In that case go the
            %safe way ranksum
            p_vals_LDA(no_LDA_pvals)=ranksum(handles_sig.win(ii).discriminant_correct,handles_sig.win(jj).discriminant_correct);
            fprintf(1, ['p values ranksum for ' handles_sig.win(ii).description ' vs. ' handles_sig.win(jj).description ' =%d\n'],p_vals_LDA(no_LDA_pvals));
        else
            if (adtest(handles_sig.win(ii).discriminant_correct)==1)||(adtest(handles_sig.win(jj).discriminant_correct)==1)
                p_vals_LDA(no_LDA_pvals)=ranksum(handles_sig.win(ii).discriminant_correct,handles_sig.win(jj).discriminant_correct);
                fprintf(1, ['p values ranksum for ' handles_sig.win(ii).description ' vs. ' handles_sig.win(jj).description ' =%d\n'],p_vals_LDA(no_LDA_pvals));
            else
                [h p_vals_LDA(no_LDA_pvals)]=ttest2(handles_sig.win(ii).discriminant_correct,handles_sig.win(jj).discriminant_correct);
                fprintf(1, ['p values t test for ' handles_sig.win(ii).description ' vs. ' handles_sig.win(jj).description ' =%d\n'],p_vals_LDA(no_LDA_pvals));
            end
        end
        
        
    end
end
  
pFDRLDA=drsFDRpval(p_vals_LDA);
fprintf(1, ['\npFDR for significant difference percent correct  = %d\n\n'],pFDRLDA);

save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) '_lda.mat'],'handles_par','handles_sig','handles_sig_spike')
  
pffft=1
