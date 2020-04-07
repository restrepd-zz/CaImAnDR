function drgCaImAnBatchPerSessionPerTrialGoGoDeltaReward

% This function calculates the per trial dFF timecourse for go-no go sessions
%for the go-go expeeriment with changes in reward
%
% The input is a series of CalmAn_batch_pre_per.mat files with the CaImAn
% data for dFF for each ROI, the data on epochs and licks. 
% Each training session includes several of these files
% The name and location of these files and some choice parameters are 
% entered in a drgCaImAnChoices file
% caimanhandles.caimandr_choices.start_reversal is the file number for the 
% start of a reversal
% Processing of the data is different if there is a reversal
%
%
% Needs a choices file such as drgCaImAnChoices_20180515_mmPVG02_Cerebellum.m
% Needs the CalmAn_batch_pre_per.mat files output files from drgCaImAn_batch_dropc.m

%
warning('off')

close all
clear all

dry_time=4;
wet_time=7.5;

these_colors{1}='b';
these_colors{2}='r';
these_colors{3}='m';
these_colors{4}='k';
these_colors{5}='g';
these_colors{6}='c';
these_colors{7}='y';

these_lines{1}='-b';
these_lines{2}='-r';
these_lines{3}='-m';
these_lines{4}='-k';
these_lines{5}='-g';
these_lines{6}='-c';
these_lines{7}='-y';
these_lines{8}='--k';




tic
 
[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSessionPerTrial run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;
files_per_trial=[];

all_lda_events=[]; 
all_lda_input_timecourse=[];

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
        files_per_trial(num_odor_trials)=filNum;
        
        %Save lda
%         all_lda_events{num_odor_trials}=lda_event{trNo};
%         szlit=size(lda_input_timecourse);
%         all_lda_input_timecourse(1:length(time_to_event),1:szlit(2),num_odor_trials)=lda_input_timecourse(:,:,trNo);
%         all_lda_no_comp=szlit(2);
        
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

%Calculate percent correct
sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

%Note: I am moving the window for calculation of perCorr to the right by nine points
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

%Now calculate the per ROI dFFs
for filNum=1:caimanhandles.caimandr_choices.no_files
        
    %Read the file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
     
     
    first_num_odor_trials(filNum)=num_odor_trials+1;
    
    for trNo=1:no_odor_trials
        
    end
    
end

%If no start_reversal file we enter a large file number so that the data
%are all processed for forward go-no go runs
if ~isfield(caimanhandles.caimandr_choices,'start_reversal')
    caimanhandles.caimandr_choices.start_reversal=2000;
end

%For reversals plot violin plot of percent
if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
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
    [handles_out2.pctPerWin(1).pct_mean, handles_out2.pctPerWin(1).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(1).pct,edges,x_val,rand_offset,'k','k',2);
    
    %Trials after reversal
    x_val=2;
    handles_out2.pctPerWin(2).pct=perCorr(first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+15:first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+65);
    bar(x_val,mean(handles_out2.pctPerWin(2).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(2).pct_mean, handles_out2.pctPerWin(2).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(2).pct,edges,x_val,rand_offset,'k','k',2);
    
    %Trials at end
    x_val=3;
    handles_out2.pctPerWin(3).pct=perCorr(end-100:end);
    bar(x_val,mean(handles_out2.pctPerWin(3).pct),'FaceColor',[0.7 0.7 0.7])
    [handles_out2.pctPerWin(3).pct_mean, handles_out2.pctPerWin(3).pct_CI]=drgViolinPoint(handles_out2.pctPerWin(3).pct,edges,x_val,rand_offset,'k','k',2);
    
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

glm_dFF=[];
glm_ii=0;
for winNo=2:szwins(1)
    figNo=figNo+1;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .05+0.3*(winNo-1) .5 .25])
    hold on
    
    last_file=1;
    
    for dFF_trNo=1:num_odor_trials_dFF
        
        
        if files_per_trial(dFF_trNo)~=last_file
            subplot(2,1,1)
            hold on
            plot([trial_dFF(dFF_trNo)-0.5 trial_dFF(dFF_trNo)-0.5], [-0.5 2],'-k')
            subplot(2,1,2)
            hold on
            plot([trial_dFF(dFF_trNo)-0.5 trial_dFF(dFF_trNo)-0.5], [-0.5 2],'-k')
            last_file=files_per_trial(dFF_trNo);
        end
        
        
        %Plot odor 1 (S+ forward)
        subplot(2,1,1)
        hold on
        %         if dFF_trNo<tr_reversal

        
        if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
            
            glm_dFF.data(glm_ii+1)=mean_win_dFF(winNo,dFF_trNo);
            glm_dFF.winNo(glm_ii+1)=winNo;
            glm_dFF.spm(glm_ii+1)=1;
            glm_dFF.reward(glm_ii+1)=caimanhandles.caimandr_choices.ul_reward(files_per_trial(dFF_trNo));
            glm_dFF.odor_pair(glm_ii+1)=caimanhandles.caimandr_choices.odor_pair(files_per_trial(dFF_trNo));
            glm_ii=glm_ii+1;
            
            %Confidence interval
            this_CI=zeros(1,2);
            this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
            plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
            
            
            switch caimanhandles.caimandr_choices.odor_pair(files_per_trial(dFF_trNo))
                case 1
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                case 2
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                case 3
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
            end
        end
        
        
        
        %             %If Hit or Miss
        %             if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
        %                 %Confidence interval
        %                 this_CI=zeros(1,2);
        %                 this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
        %                 plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==1
        %                     %Hit
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==4
        %                     %CR
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==2
        %                     %Miss
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==3
        %                     %FA
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
        %                 end
        %             end
        %         else
        %             %If CR or FA
        %              if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
        %                 %Confidence interval
        %                 this_CI=zeros(1,2);
        %                 this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
        %                 plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==1
        %                     %Hit
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==4
        %                     %CR
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==2
        %                     %Miss
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==3
        %                     %FA
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
        %                 end
        %             end
        %
        %         end
        
        
        
        
        %Plot odor 2 (S- forward)
        subplot(2,1,2)
        hold on
        %         if dFF_trNo>=tr_reversal
        %If Hit or Miss
        if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
            
            glm_dFF.data(glm_ii+1)=mean_win_dFF(winNo,dFF_trNo);
            glm_dFF.winNo(glm_ii+1)=winNo;
            glm_dFF.spm(glm_ii+1)=2;
            glm_dFF.reward(glm_ii+1)=caimanhandles.caimandr_choices.ul_reward(files_per_trial(dFF_trNo));
            glm_dFF.odor_pair(glm_ii+1)=caimanhandles.caimandr_choices.odor_pair(files_per_trial(dFF_trNo));
            glm_ii=glm_ii+1;
            
            %Confidence interval
            this_CI=zeros(1,2);
            this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
            plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
            
             
            switch caimanhandles.caimandr_choices.odor_pair(files_per_trial(dFF_trNo))
                case 1
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                case 2
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                case 3
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
            end
            
            
            %             %Confidence interval
            %             this_CI=zeros(1,2);
            %             this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
            %             plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
            %
            %             if epochs_per_trial_dFF(dFF_trNo)==1
            %                 %Hit
            %                 plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
            %             end
            %
            %             if epochs_per_trial_dFF(dFF_trNo)==4
            %                 %CR
            %                 plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
            %             end
            %
            %             if epochs_per_trial_dFF(dFF_trNo)==2
            %                 %Miss
            %                 plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
            %             end
            %
            %             if epochs_per_trial_dFF(dFF_trNo)==3
            %                 %FA
            %                 plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
            %             end
        end
        %         else
        %             %If CR or FA
        %              if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
        %                 %Confidence interval
        %                 this_CI=zeros(1,2);
        %                 this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
        %                 plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==1
        %                     %Hit
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==4
        %                     %CR
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==2
        %                     %Miss
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
        %                 end
        %
        %                 if epochs_per_trial_dFF(dFF_trNo)==3
        %                     %FA
        %                     plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
        %                 end
        %             end
        %
        %         end
    end
    
    for spno=1:2
        subplot(2,1,spno)
%         if isfield(caimanhandles.caimandr_choices,'start_reversal')
%             filNum=caimanhandles.caimandr_choices.start_reversal;
%             if (filNum>0)&(filNum<caimanhandles.caimandr_choices.no_files)
%                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
%                 text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Reversal','Color','k','FontSize',18)
%             end
%         end
%         
%         if isfield(caimanhandles.caimandr_choices,'start_gogo')
%             filNum=caimanhandles.caimandr_choices.start_gogo;
%             if filNum>0
%                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
%                 text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Go-go','Color','k','FontSize',18)
%             end
%         end
%         
%         if isfield(caimanhandles.caimandr_choices,'start_session')
%             if length(caimanhandles.caimandr_choices.start_session)>=2
%                 for sessionNo=2:length(caimanhandles.caimandr_choices.start_session)
%                     filNum=caimanhandles.caimandr_choices.start_session(sessionNo);
%                     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k')
%                     %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
%                 end
%             end
%         end
        
        plot([1 num_odor_trials_dFF],[0 0], '-k')
        
        if spno==1
            title(['Odor 1 (S+ forward)'])
        else
            title(['Odor 2 (S- forward)'])
        end
        xlabel('Trial number')
        ylabel('dF/F')
        ylim([prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))])
    end
    suptitle(['dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])

end


%Now plot all the timecourses per reward

for ii_odor_pair=1:max(caimanhandles.caimandr_choices.odor_pair)
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig1 = figure(figNo);
    set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
    hold on
    
    these_ul=caimanhandles.caimandr_choices.ul_reward(caimanhandles.caimandr_choices.odor_pair==ii_odor_pair);
    
    for ii_ul=1:length(these_ul)
        this_file=(caimanhandles.caimandr_choices.ul_reward==these_ul(ii_ul))&(caimanhandles.caimandr_choices.odor_pair==ii_odor_pair);
        this_file_ii=find(this_file==1,1,'first');
        
        these_mean_dFFs=zeros(sum(files_per_trial==this_file_ii),size(mean_snip_dFF,2));
        these_mean_dFFs(:,:)=mean_snip_dFF(files_per_trial==this_file_ii,:);
        this_mean_dFF=zeros(1,size(mean_snip_dFF,2));
        this_mean_dFF(1,:)=mean(these_mean_dFFs,1);
        these_CIs=bootci(1000, @mean, these_mean_dFFs);
        these_CIs(1,:)=this_mean_dFF- these_CIs(1,:);
        these_CIs(2,:)=(these_CIs(2,:)-this_mean_dFF);
        this_time=time(1).time_to_event(1:size(mean_snip_dFF,2))';
        plot(this_time,this_mean_dFF', these_lines{ii_ul})
        %             boundedline(this_time,this_mean_dFF', these_CIs', these_colors{ii_ul});
    end
    
    %Odor on markers
    lowdFF=-0.5;
    highdFF=3;
    plot([0 0],[lowdFF highdFF],'-k')
    odorhl=plot([0 mean(delta_odor)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[lowdFF highdFF],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF highdFF],'-r')
    
    
    ylim([lowdFF highdFF])
    xlim([-10 19.8])
    xlabel('sec')
    ylabel('dF/F')
    
    title(['dF/F timecourse ' caimanhandles.caimandr_choices.odor_pair_names{ii_odor_pair}])
    
end


glm_dFF.data(glm_ii+1)=mean_win_dFF(winNo,dFF_trNo);
glm_dFF.winNo(glm_ii+1)=winNo;
glm_dFF.spm(glm_ii+1)=2;
glm_dFF.reward(glm_ii+1)=caimanhandles.caimandr_choices.ul_reward(files_per_trial(dFF_trNo));
glm_dFF.odor_pair(glm_ii+1)=caimanhandles.caimandr_choices.odor_pair(files_per_trial(dFF_trNo));
glm_ii=glm_ii+1;

%Perform the glm for LDA percent correct
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_dFF.data',glm_dFF.winNo',glm_dFF.spm',glm_dFF.reward',glm_dFF.odor_pair',...
    'VariableNames',{'dFF','window','spm','reward','odor_pair'});
mdl = fitglm(tbl,'dFF~window+spm+reward+odor_pair+window*spm*reward*odor_pair'...
    ,'CategoricalVars',[2,3,5])

%Note: For 20200219mmPVG6f06_Cerebellum the glm yields no difference for spm and odor_pair. Because of that below I
%average both odors in each condition

%Now plot all the dFF timecourses merging all odor_pairs
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .3 .3])
hold on

unique_ul=unique(caimanhandles.caimandr_choices.ul_reward);
glm_mean_dFF=[];
glm_ii=0;
winNo=3;
all_dFFs=[];
all_uls=[];
all_CIs=zeros(2,length(unique_ul));
all_mean_dFFs=zeros(1,length(unique_ul));

glm_mean_wd_dFF=[];
glm_wd_ii=0;
all_dFFs_dry=[];
all_dFFs_wet=[];
all_CIs_dry=zeros(2,length(unique_ul));
all_mean_dFFs_dry=zeros(1,length(unique_ul));
all_CIs_wet=zeros(2,length(unique_ul));
all_mean_dFFs_wet=zeros(1,length(unique_ul));

this_time=time(1).time_to_event(1:size(mean_snip_dFF,2))';

ii_dry=find(this_time>=dry_time,1,'first');
ii_wet=find(this_time>=wet_time,1,'first');

for ii_ul=1:length(unique_ul)
    these_files=find(caimanhandles.caimandr_choices.ul_reward==unique_ul(ii_ul));
    this_no_trials=0;
    for jj=1:length(these_files)
        this_no_trials=this_no_trials+sum(files_per_trial==these_files(jj));
    end
    
    these_mean_dFFs=zeros(this_no_trials,size(mean_snip_dFF,2));
    this_ii=0;
    for jj=1:length(these_files)
        these_mean_dFFs(this_ii+1:sum(files_per_trial==these_files(jj))+this_ii,:)=mean_snip_dFF(files_per_trial==these_files(jj),:);
        this_ii=sum(files_per_trial==these_files(jj))+this_ii;
    end
    this_mean_dFF=zeros(1,size(mean_snip_dFF,2));
    this_mean_dFF(1,:)=mean(these_mean_dFFs,1);
    these_CIs=bootci(1000, @mean, these_mean_dFFs);
   
    %Find max and save the mean and CI
    [maxdFF,max_ii]=max(this_mean_dFF);
    all_mean_dFFs(ii_ul)=this_mean_dFF(max_ii);
    all_CIs(:,ii_ul)=these_CIs(:,max_ii);
    
    glm_mean_dFF.data(glm_ii+1:glm_ii+this_ii)=these_mean_dFFs(:,max_ii);
    glm_mean_dFF.ul(glm_ii+1:glm_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_ii=glm_ii+this_ii;
    
    all_dFFs=[all_dFFs these_mean_dFFs(:,max_ii)'];
    all_uls=[all_uls unique_ul(ii_ul)*ones(this_ii,1)'];
    
    %Save dry and wet
    all_mean_dFFs_dry(ii_ul)=this_mean_dFF(ii_dry);
    all_CIs_dry(:,ii_ul)=these_CIs(:,ii_dry);
    all_mean_dFFs_wet(ii_ul)=this_mean_dFF(ii_wet);
    all_CIs_wet(:,ii_ul)=these_CIs(:,ii_wet);
    
    glm_mean_wd_dFF.data(glm_wd_ii+1:glm_wd_ii+this_ii)=these_mean_dFFs(:,ii_dry);
    glm_mean_wd_dFF.ul(glm_wd_ii+1:glm_wd_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_mean_wd_dFF.dry_wet(glm_wd_ii+1:glm_wd_ii+this_ii)=zeros(1,this_ii);
    glm_wd_ii=glm_wd_ii+this_ii;
    
    glm_mean_wd_dFF.data(glm_wd_ii+1:glm_wd_ii+this_ii)=these_mean_dFFs(:,ii_wet);
    glm_mean_wd_dFF.ul(glm_wd_ii+1:glm_wd_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_mean_wd_dFF.dry_wet(glm_wd_ii+1:glm_wd_ii+this_ii)=ones(1,this_ii);
    glm_wd_ii=glm_wd_ii+this_ii;
    
    all_dFFs_dry=[all_dFFs_dry these_mean_dFFs(:,ii_dry)'];
    all_dFFs_wet=[all_dFFs_wet these_mean_dFFs(:,ii_wet)'];
    
    
    %Now plot the timecourse
    these_CIs(1,:)=this_mean_dFF- these_CIs(1,:);
    these_CIs(2,:)=(these_CIs(2,:)-this_mean_dFF);

    plot(this_time,this_mean_dFF', these_lines{ii_ul})
    
end

%Odor on markers
lowdFF=-0.5;
highdFF=3;
plot([0 0],[lowdFF highdFF],'-k')
odorhl=plot([0 mean(delta_odor)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[lowdFF highdFF],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF highdFF],'-r')


ylim([lowdFF highdFF])
xlim([-10 19.8])
ylim([-0.5 2])
xlabel('sec')
ylabel('dF/F')

title(['dF/F timecourse merging all odor pairs'])

%Plot dFF vs ul
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(unique_ul,all_mean_dFFs,'ok')

for ii=1:length(unique_ul)
    plot([unique_ul(ii) unique_ul(ii)],all_CIs(:,ii),'-k')
end

f=fit(unique_ul',all_mean_dFFs','poly1');
plot(f,unique_ul',all_mean_dFFs')

ylim([1.1 1.8])
xlim([-0.5 35])

xlabel('Reward volume (ul)')
ylabel('dF/F peak')


%Perform the glm for dFF peak vs ul
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_mean_dFF.data',glm_mean_dFF.ul',...
    'VariableNames',{'dFF','ul'});
mdl = fitglm(tbl,'dFF~ul')


%Plot dFF vs ul for wet vs dry
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(unique_ul,all_mean_dFFs_dry,'or')
plot(unique_ul,all_mean_dFFs_wet,'ob')

for ii=1:length(unique_ul)
    plot([unique_ul(ii) unique_ul(ii)],all_CIs_dry(:,ii),'-r')
    plot([unique_ul(ii) unique_ul(ii)],all_CIs_wet(:,ii),'-b')
end

f=fit(unique_ul',all_mean_dFFs_dry','poly1');
plot(f,'r')

f=fit(unique_ul',all_mean_dFFs_wet','poly1');
plot(f,'b')

ylim([0.9 1.5])
xlim([-0.5 35])

xlabel('Reward volume (ul)')
ylabel('dF/F peak')
title('dF/F peak vs. reward volume for wet vs dry')


%Perform the glm for dFF for ul and wet vs dry
fprintf(1, ['\n\nglm for dF/F ul and dry vs wet\n'])
tbl = table(glm_mean_wd_dFF.data',glm_mean_wd_dFF.ul',glm_mean_wd_dFF.dry_wet',...
    'VariableNames',{'dFF','ul','dry_wet'});
mdl = fitglm(tbl,'dFF~ul+dry_wet+ul*dry_wet','CategoricalVars',[3])
 

%Do lick analysis
%dt_lick=0.3;

dt_lick=time_to_eventLDA(2)-time_to_eventLDA(1);
time_licks=([1:ceil((dt_after+dt_before)/dt_lick)]*dt_lick)-dt_before;
no_lick_timepoints=length(time_licks);
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



%Plot the lick frequency for S+ and S-
lick_freq_per_trial=zeros(num_odor_trials,ceil((dt_after+dt_before)/dt_lick));


for trial_no=1:num_odor_trials
    
    this_lick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
    for ii_lick=1:no_licks(trial_no)
        this_lick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=this_lick_freq( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
    end
    
    
    %Convolve lick_freq using a Gaussian window
    conv_win=gausswin(no_conv_points_lick);
    
    this_lick_freq_conv=conv(this_lick_freq,conv_win,'same')/sum(conv_win);
    this_lick_freq_conv=this_lick_freq_conv/dt_lick;
    
    
    lick_freq_per_trial(trial_no,:)=this_lick_freq_conv;
    
end



%Now plot all the lick timecourses merging all odor_pairs
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

unique_ul=unique(caimanhandles.caimandr_choices.ul_reward);
glm_mean_licks=[];
glm_ii=0;
winNo=3;
all_lickfs=[];
all_uls_licks=[];
all_CIlick=zeros(2,length(unique_ul));
all_mean_lickfs=zeros(1,length(unique_ul));

glm_mean_wd_licks=[];
glm_wd_ii=0;
all_lickfs_dry=[];
all_lickfs_wet=[];
all_CIlick_dry=zeros(2,length(unique_ul));
all_mean_lickfs_dry=zeros(1,length(unique_ul));
all_CIlick_wet=zeros(2,length(unique_ul));
all_mean_lickfs_wet=zeros(1,length(unique_ul));

this_time=time(1).time_to_event(1:no_time_points)';

ii_dry=find(this_time>=dry_time,1,'first');
ii_wet=find(this_time>=wet_time,1,'first');

for ii_ul=1:length(unique_ul)
    these_files=find(caimanhandles.caimandr_choices.ul_reward==unique_ul(ii_ul));
    this_no_trials=0;
    for jj=1:length(these_files)
        this_no_trials=this_no_trials+sum(files_per_trial==these_files(jj));
    end
    
    these_mean_lickfs=zeros(this_no_trials,no_lick_timepoints);
    this_ii=0;
    for jj=1:length(these_files)
        these_mean_lickfs(this_ii+1:sum(files_per_trial==these_files(jj))+this_ii,:)=lick_freq_per_trial(files_per_trial==these_files(jj),:);
        this_ii=sum(files_per_trial==these_files(jj))+this_ii;
    end
    this_mean_lickf=zeros(1,no_lick_timepoints);
    this_mean_lickf(1,:)=mean(these_mean_lickfs,1);
    these_CIlicks=bootci(1000, @mean, these_mean_lickfs);
   
    %Find max and save the mean and CI
    [maxlickf,max_ii]=max(this_mean_lickf);
    all_mean_lickfs(ii_ul)=this_mean_lickf(max_ii);
    all_CIlick(:,ii_ul)=these_CIlicks(:,max_ii);
    
    glm_mean_licks.data(glm_ii+1:glm_ii+this_ii)=these_mean_lickfs(:,max_ii);
    glm_mean_licks.ul(glm_ii+1:glm_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_ii=glm_ii+this_ii;
    
    all_lickfs=[all_lickfs these_mean_lickfs(:,max_ii)'];
    all_uls_licks=[all_uls_licks unique_ul(ii_ul)*ones(this_ii,1)'];
    
    
    %Save dry and wet
    all_mean_lickfs_dry(ii_ul)=this_mean_lickf(ii_dry);
    all_CIlick_dry(:,ii_ul)=these_CIlicks(:,ii_dry);
    all_mean_lickfs_wet(ii_ul)=this_mean_lickf(ii_wet);
    all_CIlick_wet(:,ii_ul)=these_CIlicks(:,ii_wet);
    
    glm_mean_wd_licks.data(glm_wd_ii+1:glm_wd_ii+this_ii)=these_mean_lickfs(:,ii_dry);
    glm_mean_wd_licks.ul(glm_wd_ii+1:glm_wd_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_mean_wd_licks.dry_wet(glm_wd_ii+1:glm_wd_ii+this_ii)=zeros(1,this_ii);
    glm_wd_ii=glm_wd_ii+this_ii;
    
    glm_mean_wd_licks.data(glm_wd_ii+1:glm_wd_ii+this_ii)=these_mean_lickfs(:,ii_wet);
    glm_mean_wd_licks.ul(glm_wd_ii+1:glm_wd_ii+this_ii)=unique_ul(ii_ul)*ones(1,this_ii);
    glm_mean_wd_licks.dry_wet(glm_wd_ii+1:glm_wd_ii+this_ii)=ones(1,this_ii);
    glm_wd_ii=glm_wd_ii+this_ii;
    
    all_lickfs_dry=[all_lickfs_dry these_mean_lickfs(:,ii_dry)'];
    all_lickfs_wet=[all_lickfs_wet these_mean_lickfs(:,ii_wet)'];
    
    
    these_CIlicks(1,:)=this_mean_lickf- these_CIlicks(1,:);
    these_CIlicks(2,:)=(these_CIlicks(2,:)-this_mean_lickf);

    
    plot(time_licks,this_mean_lickf', these_lines{ii_ul})
    
end

%Odor on markers
lowdFF=0;
highdFF=20;
plot([0 0],[lowdFF highdFF],'-k')
odorhl=plot([0 mean(delta_odor)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-k','LineWidth',5);
plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')

%Reinforcement markers
plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[lowdFF highdFF],'-r')
reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF+0.05*(highdFF-lowdFF) lowdFF+0.05*(highdFF-lowdFF)],'-r','LineWidth',5);
plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[lowdFF highdFF],'-r')


ylim([lowdFF highdFF])
xlim([-10 19.8])
xlabel('sec')
ylabel('Lick frequency')

title(['Lick frequency timecourse merging all odor pairs'])

%Plot dFF vs ul
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(unique_ul,all_mean_lickfs,'ok')

for ii=1:length(unique_ul)
    plot([unique_ul(ii) unique_ul(ii)],all_CIlick(:,ii),'-k')
end

f=fit(unique_ul',all_mean_lickfs','poly2')
plot(f,unique_ul',all_mean_lickfs')

ylim([13 16])

xlabel('Reward volume (ul)')
ylabel('Lick frequency')


%Perform the glm for lick peak vs ul
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_mean_licks.data',glm_mean_licks.ul',...
    'VariableNames',{'lick_frequency','ul'});
mdl = fitglm(tbl,'lick_frequency~ul')


%Plot lickf vs ul for wet vs dry
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(unique_ul,all_mean_lickfs_dry,'or')
plot(unique_ul,all_mean_lickfs_wet,'ob')

for ii=1:length(unique_ul)
    plot([unique_ul(ii) unique_ul(ii)],all_CIlick_dry(:,ii),'-r')
    plot([unique_ul(ii) unique_ul(ii)],all_CIlick_wet(:,ii),'-b')
end

f=fit(unique_ul',all_mean_lickfs_dry','poly1');
plot(f,'r')

f=fit(unique_ul',all_mean_lickfs_wet','poly1');
plot(f,'b')

ylim([7 14])
xlim([-0.5 35])

xlabel('Reward volume (ul)')
ylabel('Lick frequency peak')


%Perform the glm for dFF for ul and wet vs dry
fprintf(1, ['\n\nglm for lick frequency ul and dry vs wet\n'])
tbl = table(glm_mean_wd_licks.data',glm_mean_wd_licks.ul',glm_mean_wd_licks.dry_wet',...
    'VariableNames',{'lick_f','ul','dry_wet'});
mdl = fitglm(tbl,'lick_f~ul+dry_wet+ul*dry_wet','CategoricalVars',[3])


%Now plot the correlation of peak lick frequency with peak dF/F 
%on a trial by trial basis
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(all_dFFs, all_lickfs,'ob')
f=fit(all_dFFs',all_lickfs','poly1')
plot(f,all_dFFs',all_lickfs')
[rho, pval]=corr(all_dFFs', all_lickfs');
fprintf(1, ['\n\nrho %d and p value %d for peak lick frequency vs. peak dF/F\n'],rho,pval)
xlabel('dF/F')
ylabel('Lick frequency (Hz)')
title('Peak dF/F and lick frequency')

%Now plot the correlation of dry lick frequency with dF/F 
%on a trial by trial basis
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(all_dFFs_dry, all_lickfs_dry,'ob')
f=fit(all_dFFs_dry',all_lickfs_dry','poly1')
plot(f,all_dFFs_dry',all_lickfs_dry')
[rho, pval]=corr(all_dFFs_dry', all_lickfs_dry');
fprintf(1, ['\n\nrho %d and p value %d for dry lick frequency vs. dF/F\n'],rho,pval)
xlabel('dF/F')
ylabel('Lick frequency (Hz)')
title('Dry dF/F and lick frequency')

%Now plot the correlation of wet lick frequency with dF/F 
%on a trial by trial basis
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on

plot(all_dFFs_wet, all_lickfs_wet,'ob')
f=fit(all_dFFs_wet',all_lickfs_wet','poly1')
plot(f,all_dFFs_wet',all_lickfs_wet')
[rho, pval]=corr(all_dFFs_wet', all_lickfs_wet');
fprintf(1, ['\n\nrho %d and p value %d for wet lick frequency vs. dF/F\n'],rho,pval)
xlabel('dF/F')
ylabel('Lick frequency (Hz)')
title('Wet dF/F and lick frequency')

 
%Now plot the correlation of mean peak lick frequency with mean peak dF/F 
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on


for ii_ul=1:length(unique_ul)
   plot(all_mean_dFFs(ii_ul),all_mean_lickfs(ii_ul), ['o' these_colors{ii_ul}])
   plot(all_CIs(:,ii_ul), [all_mean_lickfs(ii_ul) all_mean_lickfs(ii_ul)],'-k')
   plot([all_mean_dFFs(ii_ul) all_mean_dFFs(ii_ul)], all_CIlick(:,ii_ul),'-k')
end


f=fit(all_mean_dFFs',all_mean_lickfs','poly1')
plot(f,all_mean_dFFs',all_mean_lickfs')
[rho, pval]=corr(all_mean_dFFs', all_mean_lickfs');
fprintf(1, ['\n\nrho %d and p value %d for mean lick frequency vs. mean dF/F\n'],rho,pval)
xlabel('Mean dF/F')
ylabel('Mean lick frequency (Hz)')


%Now plot the correlation of mean dry lick frequency with mean peak dF/F 
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on


for ii_ul=1:length(unique_ul)
   plot(all_mean_dFFs_dry(ii_ul),all_mean_lickfs_dry(ii_ul), ['o' these_colors{ii_ul}])
   plot(all_CIs_dry(:,ii_ul), [all_mean_lickfs_dry(ii_ul) all_mean_lickfs_dry(ii_ul)],'-k')
   plot([all_mean_dFFs_dry(ii_ul) all_mean_dFFs_dry(ii_ul)], all_CIlick_dry(:,ii_ul),'-k')
end


f=fit(all_mean_dFFs_dry',all_mean_lickfs_dry','poly1')
plot(f,'k')
[rho, pval]=corr(all_mean_dFFs_dry', all_mean_lickfs_dry');
fprintf(1, ['\n\nrho %d and p value %d for mean lick frequency vs. mean dF/F dry\n'],rho,pval)
xlabel('Mean dF/F')
ylabel('Mean lick frequency (Hz)')
title('dry')


%Now plot the correlation of mean wet lick frequency with mean peak dF/F 
figNo=figNo+1;
try
    close(figNo)
catch
end
hFig1 = figure(figNo);
set(hFig1, 'units','normalized','position',[.25 .25 .25 .25])
hold on


for ii_ul=1:length(unique_ul)
   plot(all_mean_dFFs_wet(ii_ul),all_mean_lickfs_wet(ii_ul), ['o' these_colors{ii_ul}])
   plot(all_CIs_wet(:,ii_ul), [all_mean_lickfs_wet(ii_ul) all_mean_lickfs_wet(ii_ul)],'-k')
   plot([all_mean_dFFs_wet(ii_ul) all_mean_dFFs_wet(ii_ul)], all_CIlick_wet(:,ii_ul),'-k')
end


f=fit(all_mean_dFFs_wet',all_mean_lickfs_wet','poly1')
plot(f,'k')
[rho, pval]=corr(all_mean_dFFs_wet', all_mean_lickfs_wet');
fprintf(1, ['\n\nrho %d and p value %d for mean lick frequency vs. mean dF/F wet\n'],rho,pval)
xlabel('Mean dF/F')
ylabel('Mean lick frequency (Hz)')
title('wet')

% save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName],'handles_out2')
 
pffft=1
