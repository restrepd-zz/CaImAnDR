function drgCaImAnBatchPerSessionPerTrial

% This function calculates the per trial dFF timecourse for go-no go sessions
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
            if filNum>0
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


%Now show violin plot for the different time windows for dFF calculated
%within each time window

%Parameters for violin plot
prct1=prctile(mean_win_dFF(:),1);
prct99=prctile(mean_win_dFF(:),99);
edges=prct1:(prct99-prct1)/20:prct99;
rand_offset=0.8;

handles_out2.mean_win_dFF=mean_win_dFF;
handles_out2.prct1=prct1;
handles_out2.edges=edges;
handles_out2.rand_offset=rand_offset;

if caimanhandles.caimandr_choices.start_reversal<caimanhandles.caimandr_choices.no_files
    for winNo=1:szwins(1)
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.25 .2+0.1*(winNo-1) .5 .5])
        
        
        
        %Before reversal
        dFF_trial_mask=[];
        for ii=1:num_odor_trials_dFF
            if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
                dFF_trial_mask(ii)=1;
            else
                dFF_trial_mask(ii)=0;
            end
        end
        
        %Odor 1 subplot
        subplot(2,1,1)
        hold on
        
        %Hit
        x_val=1;
        sample_no=1;
        decription_of_sample{sample_no}='Odor 1 Hit before';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==1));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==1))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==1));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 0])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Miss
        x_val=2;
        sample_no=2;
        decription_of_sample{sample_no}='Odor 1 Miss before';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==2));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==2))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==2));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 1 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Odor 2 subplot
        subplot(2,1,2)
        hold on
        
        %CR
        x_val=1;
        sample_no=3;
        decription_of_sample{sample_no}='Odor 2 CR before';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==4));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==4))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==4));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %FA
        x_val=2;
        sample_no=4;
        decription_of_sample{sample_no}='Odor 2 FA before';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==3));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==3))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==3));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %After reversal
        dFF_trial_mask=[];
        for ii=1:num_odor_trials_dFF
            if (trial_dFF(ii)>=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+15)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)+65)
                dFF_trial_mask(ii)=1;
            else
                dFF_trial_mask(ii)=0;
            end
        end
        
        %Odor 2 subplot
        subplot(2,1,2)
        hold on
        
        %Hit
        x_val=4;
        sample_no=5;
        decription_of_sample{sample_no}='Odor 2 Hit after';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==1));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==1))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==1));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 0])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Miss
        x_val=5;
        sample_no=6;
        decription_of_sample{sample_no}='Odor 2 Miss after';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==2));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==2))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==2));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 1 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Odor 1 subplot
        subplot(2,1,1)
        hold on
        
        %CR
        x_val=4;
        sample_no=7;
        decription_of_sample{sample_no}='Odor 1 CR after';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==4));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==4))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==4));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %FA
        x_val=5;
        sample_no=8;
        decription_of_sample{sample_no}='Odor 1 FA after';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==3));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==3))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==3));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %At end
        dFF_trial_mask=[];
        for ii=1:num_odor_trials_dFF
            if (trial_dFF(ii)>=num_odor_trials-100)&(trial_dFF(ii)<=num_odor_trials)
                dFF_trial_mask(ii)=1;
            else
                dFF_trial_mask(ii)=0;
            end
        end
        
        %Odor 2 subplot
        subplot(2,1,2)
        hold on
        
        %Hit
        x_val=7;
        sample_no=9;
        decription_of_sample{sample_no}='Odor 2 Hit at end';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==1));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==1))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==1));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 0])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Miss
        x_val=8;
        sample_no=10;
        decription_of_sample{sample_no}='Odor 2 Miss at end';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==2));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==2))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==2));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 1 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %Odor 1 subplot
        subplot(2,1,1)
        hold on
        
        %CR
        x_val=7;
        sample_no=11;
        decription_of_sample{sample_no}='Odor 1 CR at end';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==4));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==4))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==4));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[0 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        %FA
        x_val=8;
        sample_no=12;
        decription_of_sample{sample_no}='Odor 1 FA at end';
        no_trials_per_sample(sample_no)=sum(dFF_trial_mask&(epochs_per_trial_dFF==3));
        if sum(dFF_trial_mask&(epochs_per_trial_dFF==3))>=2
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=1;
            handles_out2.dFFPerWin(winNo,sample_no).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==3));
            bar(x_val,mean(handles_out2.dFFPerWin(winNo,sample_no).dFF),'FaceColor',[1 0 1])
            [handles_out2.dFFPerWin(winNo,sample_no).dFF_mean, handles_out2.dFFPerWin(winNo,sample_no).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,sample_no).dFF,edges,x_val,rand_offset,'k','k',2);
        else
            handles_out2.dFFPerWin(winNo,sample_no).dFFexists=0;
        end
        
        subplot(2,1,1)
        ylim([prct1 prct99])
        xlim([0 9])
        ylabel('dF/F')
        xticks([2.5 8.5 14.5])
        xticklabels({'Before','After','End'})
        
        subplot(2,1,2)
        ylim([prct1 prct99])
        xlim([0 9])
        ylabel('dF/F')
        xticks([1.5 4.5 7.5])
        xticklabels({'Before','After','End'})
        
        suptitle(['dF/F before and after reversal for window ' num2str(winNo)])
        
        %Do test of significance
        fprintf(1, ['Tests of significance for difference in dFF for window No %d\n\n'],winNo)
        p_vals_dFF=[];
        no_dFF_pvals=0;
        for ii=1:12
            for jj=ii+1:12
                if (handles_out2.dFFPerWin(winNo,ii).dFFexists==1)&(handles_out2.dFFPerWin(winNo,jj).dFFexists==1)
                    no_dFF_pvals=no_dFF_pvals+1;
                    if (length(handles_out2.dFFPerWin(winNo,jj).dFF)<4)||(length(handles_out2.dFFPerWin(winNo,ii).dFF)<4)
                        %adtest does not work with n<4. In that case go the
                        %safe way ranksum
                        p_vals_dFF(no_dFF_pvals)=ranksum(handles_out2.dFFPerWin(winNo,ii).dFF,handles_out2.dFFPerWin(winNo,jj).dFF);
                        fprintf(1, ['p values ranksum for ' decription_of_sample{ii} ' vs. ' decription_of_sample{jj} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                    else
                        if (adtest(handles_out2.dFFPerWin(winNo,ii).dFF)==1)||(adtest(handles_out2.dFFPerWin(winNo,jj).dFF)==1)
                            p_vals_dFF(no_dFF_pvals)=ranksum(handles_out2.dFFPerWin(winNo,ii).dFF,handles_out2.dFFPerWin(winNo,jj).dFF);
                            fprintf(1, ['p values ranksum for ' decription_of_sample{ii} ' vs. ' decription_of_sample{jj} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        else
                            [h p_vals_dFF(no_dFF_pvals)]=ttest2(handles_out2.dFFPerWin(winNo,ii).dFF,handles_out2.dFFPerWin(winNo,jj).dFF);
                            fprintf(1, ['p values t test for ' decription_of_sample{ii} ' vs. ' decription_of_sample{jj} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        end
                    end
                    
                end
            end
        end
        
        pFDRdFF=drsFDRpval(p_vals_dFF);
        fprintf(1, ['pFDR for significant difference percent correct  = %d\n\n'],pFDRdFF);
        
        
    end
else
    %If there are no reversals do a violin plot of dFF vs percent
    for winNo=1:szwins(1)
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.25 .2+0.1*(winNo-1) .5 .5])
         
         
        pct_windows=[45 65;65 80;80 100.1];
        
        for pct_win_ii=1:3
            
            dFF_trial_mask=((perCorr>=pct_windows( pct_win_ii,1))&(perCorr<pct_windows( pct_win_ii,2)));
            x_val=pct_win_ii*2;
            
            %Hit
            subplot(2,2,1)
            hold on
            epochNo=1;
            handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFFno=sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
            if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>=2
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFFexists=1;
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo,dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                [handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean, handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,edges,x_val,rand_offset,'r','k',2);
            else
                if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>0
                    handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                    plot(2*pct_win_ii,handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
                end
            end
            
            %Miss
            subplot(2,2,2)
            hold on
            epochNo=2;
            if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>=2
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFFexists=1;
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                [handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean, handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,edges,x_val,rand_offset,'c','k',2);
            else
                if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>0
                    handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                    plot(2*pct_win_ii,handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
                end
            end
            
            
            %CR
            subplot(2,2,3)
            hold on
            epochNo=4;
            if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>=2
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFFexists=1;
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                [handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean, handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,edges,x_val,rand_offset,'b','k',2);
            else
                if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>0
                    handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                    plot(2*pct_win_ii,handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
                end
            end
            
            %FA
            subplot(2,2,4)
            hold on
            epochNo=3;
            if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>=2
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFFexists=1;
                handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                [handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean, handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_CI]=drgViolinPoint(handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,edges,x_val,rand_offset,'m','k',2);
            else
                if sum(dFF_trial_mask&(epochs_per_trial_dFF==epochNo))>0
                    handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF=mean_win_dFF(winNo, dFF_trial_mask&(epochs_per_trial_dFF==epochNo));
                    plot(2*pct_win_ii,handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF,'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
                end
            end
        end
        
        for epochNo=1:4
             
            mean_dFF=[];
            x_val=[];
            ii=0;
            for pct_win_ii=1:3
                if ~isempty(handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean)
                    ii=ii+1;
                    x_val(ii)=2*pct_win_ii;
                    mean_dFF(ii)=handles_out2.dFFPerWin(winNo,epochNo,pct_win_ii).dFF_mean;
                end
            end
            switch epochNo
                case 1
                    subplot(2,2,1)
                    if ii>1
                        plot(x_val,mean_dFF,'-r')
                    end
                    ylim([prct1 prct99])
                    xlim([1 7])
                    ylabel('dF/F')
                    xticks([2 4 6])
                    xticklabels({'<65','>=65&<80','>=80'})
                    xlabel('Percent correct')
                case 2
                    subplot(2,2,2)
                    if ii>1
                        plot(x_val,mean_dFF,'-c')
                    end
                    ylim([prct1 prct99])
                    xlim([1 7])
                    ylabel('dF/F')
                    xticks([2 4 6])
                    xticklabels({'<65','>=65&<80','>=80'})
                    xlabel('Percent correct')
                case 3
                    subplot(2,2,4)
                    if ii>1
                        plot(x_val,mean_dFF,'-m')
                    end
                    ylim([prct1 prct99])
                    xlim([1 7])
                    ylabel('dF/F')
                    xticks([2 4 6])
                    xticklabels({'<65','>=65&<80','>=80'})
                    xlabel('Percent correct')
                case 4
                    subplot(2,2,3)
                    if ii>1
                        plot(x_val,mean_dFF,'-b')
                    end
                    ylim([prct1 prct99])
                    xlim([1 7])
                    ylabel('dF/F')
                    xticks([2 4 6])
                    xticklabels({'<65','>=65&<80','>=80'})
                    xlabel('Percent correct')
            end
            
        end
        
        
        suptitle(['dF/F vs percent correct for window ' num2str(winNo) ' ' choiceFileName])
        
        per_description{1}='<65';
        per_description{2}='>=65&<80';
        per_description{3}='>=80';
        
        epoch_description{1}='Hit';
        epoch_description{2}='Miss';
        epoch_description{3}='CR';
        epoch_description{4}='FA';
        
        %Do test of significance
        fprintf(1, ['Tests of significance for difference in dFF for window No %d\n\n'],winNo)
        
        p_vals_dFF=[];
        no_dFF_pvals=0;
        for epoch1=1:4
            for epoch2=epoch1+1:4
                for per_ii1=1:3
                    
                    no_dFF_pvals=no_dFF_pvals+1;
                    if (~isempty(handles_out2.dFFPerWin(winNo,epoch2,per_ii1).dFF))&(~isempty(handles_out2.dFFPerWin(winNo,epoch1,per_ii1).dFF))
                        [p_vals_dFF(no_dFF_pvals), r_or_t]=drg_ranksum_or_ttest(handles_out2.dFFPerWin(winNo,epoch1,per_ii1).dFF,handles_out2.dFFPerWin(winNo,epoch2,per_ii1).dFF);
                        if r_or_t==0
                            fprintf(1, ['p values ranksum for ' epoch_description{epoch1} ' ' per_description{per_ii1} ' vs. ' epoch_description{epoch2} ' ' per_description{per_ii1} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        else
                            fprintf(1, ['p values t test for ' epoch_description{epoch1} ' ' per_description{per_ii1} ' vs. ' epoch_description{epoch2} ' ' per_description{per_ii1} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        end
                    end
                    
                end
                
            end
        end
        
        for epoch1=1:4
            
            for per_ii1=1:3
                for per_ii2=per_ii1+1:3
                    no_dFF_pvals=no_dFF_pvals+1;
                    if (~isempty(handles_out2.dFFPerWin(winNo,epoch1,per_ii2).dFF))&(~isempty(handles_out2.dFFPerWin(winNo,epoch1,per_ii1).dFF))
                        [p_vals_dFF(no_dFF_pvals), r_or_t]=drg_ranksum_or_ttest(handles_out2.dFFPerWin(winNo,epoch1,per_ii1).dFF,handles_out2.dFFPerWin(winNo,epoch1,per_ii2).dFF);
                        if r_or_t==0
                            fprintf(1, ['p values ranksum for ' epoch_description{epoch1} ' ' per_description{per_ii1} ' vs. ' epoch_description{epoch1} ' ' per_description{per_ii2} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        else
                            fprintf(1, ['p values t test for ' epoch_description{epoch1} ' ' per_description{per_ii1} ' vs. ' epoch_description{epoch1} ' ' per_description{per_ii2} ' =%d\n'],p_vals_dFF(no_dFF_pvals));
                        end
                        
                    end
                end
            end
            
            
        end
        
        pFDRdFF=drsFDRpval(p_vals_dFF);
        fprintf(1, ['pFDR for significant difference percent correct  = %d\n\n'],pFDRdFF);
        
    end
end
 
%Plot bounded lines for S+ and S- snips for the different percent windows
if caimanhandles.caimandr_choices.start_reversal>caimanhandles.caimandr_choices.no_files
    
    %Calculate the means and CIs
        timepoints=200000;
        for dFF_trNo=1:num_odor_trials_dFF
            timepoints=min([timepoints length(time(dFF_trNo).time_to_event)]);
        end
        
        pct_windows=[45 65;65 80;80 100.1];
        
        all_mean_snip_dFFs=zeros(2,3,timepoints);
        all_CI_snip_dFFs=zeros(2,3,2,timepoints);
            
        for pct_win_ii=1:3
            dFF_percent_mask=((perCorr>=pct_windows( pct_win_ii,1))&(perCorr<pct_windows( pct_win_ii,2)));

            %S+
            dFF_epoch_mask=(epochs_per_trial_dFF==1)|(epochs_per_trial_dFF==2);
            if sum(dFF_percent_mask&dFF_epoch_mask)>2
                all_mean_snip_dFFs(1,pct_win_ii,:)=mean(mean_snip_dFF(dFF_percent_mask&dFF_epoch_mask,1:timepoints),1);
                all_CI_snip_dFFs(1,pct_win_ii,:,:)=bootci(1000, @mean, mean_snip_dFF(dFF_percent_mask&dFF_epoch_mask,1:timepoints));
            end
             
            %S-
            dFF_epoch_mask=(epochs_per_trial_dFF==3)|(epochs_per_trial_dFF==4);
            if sum(dFF_percent_mask&dFF_epoch_mask)>2
                all_mean_snip_dFFs(2,pct_win_ii,:)=mean(mean_snip_dFF(dFF_percent_mask&dFF_epoch_mask,1:timepoints),1);
                all_CI_snip_dFFs(2,pct_win_ii,:,:)=bootci(1000, @mean, mean_snip_dFF(dFF_percent_mask&dFF_epoch_mask,1:timepoints));
            end
        end
        
        this_time=time(1).time_to_event(1:timepoints)';
        
        highdFF=prctile(all_CI_snip_dFFs(:),99)+0.1*(prctile(all_CI_snip_dFFs(:),99)-prctile(all_CI_snip_dFFs(:),1));
        lowdFF=prctile(all_CI_snip_dFFs(:),1)-0.1*(prctile(all_CI_snip_dFFs(:),99)-prctile(all_CI_snip_dFFs(:),1));
        
        for pct_win_ii=1:3
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
                
            hFig = figure(figNo);
            set(hFig, 'units','normalized','position',[.25 .2+0.1*(winNo-1) .5 .5])
            hold on
            
            %S-
            this_mean_snip_dFF=zeros(timepoints,1);
            this_mean_snip_dFF(:,1)=all_mean_snip_dFFs(2,pct_win_ii,:);
            this_CI_snip_dFF=zeros(2,timepoints);
            this_CI_snip_dFF(:,:)=all_CI_snip_dFFs(2,pct_win_ii,:,:);
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=(this_CI_snip_dFF(:,2)-this_mean_snip_dFF);
            boundedline(this_time,this_mean_snip_dFF, this_CI_snip_dFF, 'b');
            this_mean_snip_dFFsp=this_mean_snip_dFF;
            
            %S+
            this_mean_snip_dFF=zeros(timepoints,1);
            this_mean_snip_dFF(:,1)=all_mean_snip_dFFs(1,pct_win_ii,:);
            this_CI_snip_dFF=zeros(2,timepoints);
            this_CI_snip_dFF(:,:)=all_CI_snip_dFFs(1,pct_win_ii,:,:);
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=(this_CI_snip_dFF(:,2)-this_mean_snip_dFF);
            boundedline(this_time,this_mean_snip_dFF, this_CI_snip_dFF, 'r');
            this_mean_snip_dFFsm=this_mean_snip_dFF;
            
            %Odor on markers
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
            title(['Trials with percent correct ' per_description{pct_win_ii}])
            
            if pct_win_ii==3
                save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) 'gcamp.mat'],'this_time','this_mean_snip_dFFsm','this_mean_snip_dFFsp')
            end
                
        end
    
end
 
%Now plot snips for all trials
    
%Initailize the figures
figNoOd1=figNo+1;
subPOd1=0;
try
    close (figNoOd1)
catch
end

hFigOd1 = figure(figNoOd1);
set(hFigOd1, 'units','normalized','position',[.05 .1 0.4 0.7])


figNoOd2=figNo+2;
subPOd2=0;
try
    close (figNoOd2)
catch
end

hFigOd2 = figure(figNoOd2);
set(hFigOd2, 'units','normalized','position',[.55 .1 0.4 0.7])

lowdFF=prctile(mean_snip_dFF(:),0.01)-0.1*(prctile(mean_snip_dFF(:),99)-prctile(mean_snip_dFF(:),1));
highdFF=prctile(mean_snip_dFF(:),99.99)+0.1*(prctile(mean_snip_dFF(:),99)-prctile(mean_snip_dFF(:),1));

for dFF_trNo=1:num_odor_trials_dFF
    
    if dFF_trNo==80
        pffft=1;
    end
    
    %Plot odor 1 (S+ forward)
    if dFF_trNo<tr_reversal
        
        %Forward
        if epochs_per_trial_dFF(dFF_trNo)==1
             
            %Hit, this is odor 1
             
            subPOd1=subPOd1+1;
            figure(figNoOd1)
            
            if subPOd1>12
                suptitle('Odor 1 S+ in forward')
                subPOd1=1;
                figNoOd1=figNoOd1+2;
                
                try
                    close (figNoOd1)
                catch
                end
                
                hFigOd1 = figure(figNoOd1);
                set(hFigOd1, 'units','normalized','position',[.05 .1 0.4 0.7])
            end
            
            subplot(4,3,subPOd1)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'r');
            
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['fwd ' num2str(dFF_trNo)])
            
            %Used to generate the examples
            if dFF_trNo==52
%                 pffft=1;
%                 figure(51)
%                 hold on
%                 
%                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
%                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=this_CI_snip_dFF';
%                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
%                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
%                 
%                 
%                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'r');
%                 plot([0 0],[lowdFF highdFF],'-k')
%                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
%                 
%                 
%                 ylim([lowdFF highdFF])
%                 xlim([-10 19.8])
%                 xlabel('sec')
%                 ylabel('dF/F')
%                 title(['fwd ' num2str(dFF_trNo)])
            end
            
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==4
            %CR this is Odor 2
            
            subPOd2=subPOd2+1;
            figure(figNoOd2)
            
            if subPOd2>12
                suptitle('Odor 2 S- in forward')
                subPOd2=1;
                figNoOd2=figNoOd2+2;
                
                try
                    close (figNoOd2)
                catch
                end
                
                hFigOd2 = figure(figNoOd2);
                set(hFigOd2, 'units','normalized','position',[.5 .1 0.4 0.7])
            end
            
            subplot(4,3,subPOd2)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'b');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['fwd ' num2str(dFF_trNo)])
            
            if dFF_trNo==51
%                 figure(50)
%                 hold on
%                 
%                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
%                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=this_CI_snip_dFF';
%                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
%                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
%                 
%                 
%                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'b'); %Odor on markers
%                 plot([0 0],[lowdFF highdFF],'-k')
%                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
%                 
%                 
%                 
%                 ylim([lowdFF highdFF])
%                 xlim([-10 19.8])
%                 xlabel('sec')
%                 ylabel('dF/F')
%                 title(['fwd ' num2str(dFF_trNo)])
                
            end
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==2
            %Miss, this is Odor 1
            
            subPOd1=subPOd1+1;
            figure(figNoOd1)
            
            if subPOd1>12
                suptitle('Odor 1 S+ in forward')
                subPOd1=1;
                figNoOd1=figNoOd1+2;
                
                try
                    close (figNoOd1)
                catch
                end
                
                hFigOd1 = figure(figNoOd1);
                set(hFigOd1, 'units','normalized','position',[.05 .1 0.4 0.7])
            end
            
            subplot(4,3,subPOd1)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'c');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['fwd ' num2str(dFF_trNo)])
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==3
            %FA this is Odor 2
            
            subPOd2=subPOd2+1;
            figure(figNoOd2)
            
            if subPOd2>12
                suptitle('Odor 2 S- in forward')
                subPOd2=1;
                figNoOd2=figNoOd2+2;
                
                try
                    close (figNoOd2)
                catch
                end
                
                hFigOd2 = figure(figNoOd2);
                set(hFigOd2, 'units','normalized','position',[.5 .1 0.4 0.7])
            end
           
            subplot(4,3,subPOd2)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'm');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['fwd ' num2str(dFF_trNo)])
        end
        
    else
        
        %Reverse
        
        if epochs_per_trial_dFF(dFF_trNo)==1
            
            %Hit, this is odor 2
            
            subPOd2=subPOd2+1;
             figure(figNoOd2)
             
            if (subPOd2>12)
                suptitle('Odor 2 S+ after reversal')
                subPOd2=1;
                figNoOd2=figNoOd2+2;
                
                try
                    close (figNoOd2)
                catch
                end
                
                hFigOd2 = figure(figNoOd2);
                set(hFigOd2, 'units','normalized','position',[.5 .1 0.4 0.7])
            end
            
            subplot(4,3,subPOd2)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'r');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['\color{red}rev ' num2str(dFF_trNo)])
            
            if dFF_trNo==438
                pffft=1;
                %                 figure(55)
                %                 hold on
                %
                %                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
                %                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
                %                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
                %                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
                %                 this_CI_snip_dFF=this_CI_snip_dFF';
                %                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
                %                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
                %
                %
                %                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'r');
                %                 %Odor on markers
                %                 plot([0 0],[lowdFF highdFF],'-k')
                %                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
                %
                %                 ylim([lowdFF highdFF])
                %                 xlim([-10 19.8])
                %                 xlabel('sec')
                %                 ylabel('dF/F')
                %                 title(['\color{red}rev ' num2str(dFF_trNo)])
            end
            
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==4
            
            %CR this is Odor 1
            
            subPOd1=subPOd1+1;
            figure(figNoOd1)
            
            if (subPOd1>12)
                suptitle('Odor 1 S- after reversal')
                subPOd1=1;
                figNoOd1=figNoOd1+2;
                
                try
                    close (figNoOd1)
                catch
                end
                
                hFigOd1 = figure(figNoOd1);
                set(hFigOd1, 'units','normalized','position',[.05 .1 0.4 0.7])
            end
            
            subplot(4,3,subPOd1)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'b');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['\color{red}rev ' num2str(dFF_trNo)])
            
            if dFF_trNo==439
                pffft=1;
%                 figure(57)
%                 hold on
%                 
%                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
%                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=this_CI_snip_dFF';
%                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
%                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
%                 
%                 
%                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'b');
%                 
%                 %Odor on markers
%                 plot([0 0],[lowdFF highdFF],'-k')
%                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
%                 
%                 ylim([lowdFF highdFF])
%                 xlim([-10 19.8])
%                 xlabel('sec')
%                 ylabel('dF/F')
%                 title(['\color{red}rev ' num2str(dFF_trNo)])
            end
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==2
            %Miss, this is Odor 2
            
            subPOd2=subPOd2+1;
            figure(figNoOd2)
            
            if (subPOd2>12)
                suptitle('Odor 2 S+ after reversal')
                subPOd2=1;
                figNoOd2=figNoOd2+2;
                
                try
                    close (figNoOd2)
                catch
                end
                
                hFigOd2 = figure(figNoOd2);
                set(hFigOd2, 'units','normalized','position',[.5 .1 0.4 0.7])
            end
            
            
            subplot(4,3,subPOd2)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'c');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['\color{red}rev' num2str(dFF_trNo)])
            if dFF_trNo==91
                pffft=1;
%                 figure(53)
%                 hold on
%                 
%                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
%                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=this_CI_snip_dFF';
%                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
%                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
%                 
%                 
%                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'c');
%                 
%                 %Odor on markers
%                 plot([0 0],[lowdFF highdFF],'-k')
%                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
%                 
%                 ylim([lowdFF highdFF])
%                 xlim([-10 19.8])
%                 xlabel('sec')
%                 ylabel('dF/F')
%                 title(['\color{red}rev' num2str(dFF_trNo)])
            end
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==3
            %FA this is Odor 1
         
            subPOd1=subPOd1+1;
            figure(figNoOd1)
            
            if (subPOd1>12)
                suptitle('Odor 1 S- after reversal')
                subPOd1=1;
                figNoOd1=figNoOd1+2;
                
                try
                    close (figNoOd1)
                catch
                end
                
                hFigOd1 = figure(figNoOd1);
                set(hFigOd1, 'units','normalized','position',[.05 .1 0.4 0.7])
            end
            
            
            subplot(4,3,subPOd1)
            hold on
            
            this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
            this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
            this_CI_snip_dFF=this_CI_snip_dFF';
            this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
            this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
            
            
            boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'm');
            
            ylim([lowdFF highdFF])
            xlim([-10 19.8])
            xlabel('sec')
            ylabel('dF/F')
            title(['\color{red}rev' num2str(dFF_trNo)])
            
            if dFF_trNo==92
                pffft=1;
%                 figure(54)
%                 hold on
%                 
%                 this_mean_snip_dFF=zeros(length(time(dFF_trNo).time_to_event),1);
%                 this_mean_snip_dFF(:,1)=mean_snip_dFF(dFF_trNo,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=zeros(2,length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF(:,:)=CI_snip_dFF(dFF_trNo,:,1:length(time(dFF_trNo).time_to_event));
%                 this_CI_snip_dFF=this_CI_snip_dFF';
%                 this_CI_snip_dFF(:,1)=this_mean_snip_dFF-this_CI_snip_dFF(:,1);
%                 this_CI_snip_dFF(:,2)=this_CI_snip_dFF(:,2)-this_mean_snip_dFF;
%                 
%                 
%                 boundedline(time(dFF_trNo).time_to_event',this_mean_snip_dFF, this_CI_snip_dFF, 'm');
%                 
%                 %Odor on markers
%                 plot([0 0],[lowdFF highdFF],'-k')
%                 plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
%                 
%                 ylim([lowdFF highdFF])
%                 xlim([-10 19.8])
%                 xlabel('sec')
%                 ylabel('dF/F')
%                 title(['\color{red}rev' num2str(dFF_trNo)])
            end
        end
        
    end
    
   
    
    %Odor on markers
    plot([0 0],[lowdFF highdFF],'-k')
    plot([mean(delta_odor) mean(delta_odor)],[lowdFF highdFF],'-k')
    
    
    
end

hFigOd1 = figure(figNoOd1);
suptitle('Odor 1 S- after reversal')

figure(figNoOd2)
suptitle('Odor 2 S+ after reversal')

%Show all responses for each odor in a pseudocolor for two ranges of trials

save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName],'handles_out2')
 
pffft=1
