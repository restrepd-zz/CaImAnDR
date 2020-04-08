function drgCaImAnBatchMultiSessionLicksv3(choiceBatchPathName,choiceFileName)

%This code does a lick analysis for sessions where we record the behavior and licks
%for optogenetic modulation of behavior/licks by turning MLIs with
%optogenetics
% Needs a choices file such as
% drgCaImAn_multichoices_PVhM4Di_all_files_08192019.m
% Needs the output files from drgCaImAn_batch_dropc_no_microscope.m
%
% Use display_choice=1 for
% drgCaImAn_multichoices_PVhM4Di_all_files_05052019.m 
%
%Used display_choice=2 for drgCaImAn_multichoices_PVhM4Di_all_files_08192019
warning('off')

close all
clear all

display_choice=2;

lick_threshold=80; %This is the threshold to exclude the runs where Ming was adding water manually, this was set at 20
t_odor_on=0;
t_odor_off=4;

min_trials=20; %Number of trials for the sliding window for percent correct

tic

if nargin==0
    [choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_multichoices*.m'},'Select the .m file with all the choices for analysis');
end

fprintf(1, ['\ndrgCaImAnBatchMultiSessionLicks run for ' choiceFileName '\n\n']);



addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

mice_included=unique(caimanhandles.caimandr_choices.mouse);
for ii_mouse=1:length(mice_included)
    if mice_included(ii_mouse)~=0
        for filNum=1:caimanhandles.caimandr_choices.no_files
            if caimanhandles.caimandr_choices.mouse(filNum)~=0
                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_smtrials=0;
                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_sptrials=0;
                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).splicks=zeros(179,1000);
                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).smlicks=zeros(179,1000);
            end
        end
        for hM4DNo=1:4
           caimanhandles.per_mouse(mice_included(ii_mouse)).hM4D(hM4DNo).no_sptrials=0;
           caimanhandles.per_mouse(mice_included(ii_mouse)).hM4D(hM4DNo).no_smtrials=0;
        end
    end
end



% szwins=size(caimanhandles.caimandr_choices.wins);

%Read the files and calculate the dF/F in each window
per_file_lick_freq_sp=[];
per_file_lick_freq_sm=[];
per_file_lick_log10_p_val=[];
no_files_included=0;
per_file_hM4D=[];
per_file_mouse=[];
per_file_trial_window=[];
per_file_decision_time=[];
per_file_decision_time_reinf=[];
per_file_fractional_lick_time_on_sp=[];
per_file_fractional_lick_time_on_sm=[];

for filNum=1:caimanhandles.caimandr_choices.no_files
    
    num_odor_trials=0;
    epochs_per_trial=[];
    num_odor_trials_dFF=0;
    
    all_lda_events=[];
    all_lda_events_miss_FA=[];
    all_opto_on_per_trial=[];
    
    lick_times=[];
    no_licks=[];
    dLickTraces=[];
    
    handles_outs.lick_slopes=[];
    handles_outs.dFF_slopes=[];
    handles_outs.no_lick_slopes=0;
    handles_outs.no_dFF_slopes=0;
    handles_outs.no_vels=0;
    handles_outs.velocities=[];
    handles_outs.no_accels=0;
    handles_outs.acceleration=[];
    handles_outs.epochs_per_trial=[];
    handles_outs.perCorr=[];
    
    %Variables for optical flow
    skip_ii=19;
    baseline_ii=200;
    
    
    
    %Read the CalmAn_batch file
    if isfile([caimanhandles.caimandr_choices.pathName{filNum} caimanhandles.caimandr_choices.fileName{filNum}])
        load([caimanhandles.caimandr_choices.pathName{filNum} caimanhandles.caimandr_choices.fileName{filNum}])
        
        if no_odor_trials>20
            
            first_num_odor_trials(filNum)=num_odor_trials+1;
            
            
            for trNo=1:no_odor_trials
                
                %Save epoch
                num_odor_trials=num_odor_trials+1;
                all_lda_events{num_odor_trials}=lda_event{trNo};
                all_opto_on_per_trial(num_odor_trials)=opto_on_per_trial(trNo);
                
                
                if epoch_per_trial(trNo)==6
                    %Hit
                    epochs_per_trial(1,num_odor_trials)=1;
                    epochs_per_trial(2:4,num_odor_trials)=0;
                    all_lda_events_miss_FA(num_odor_trials)=1;
                    
                    
                    if sum(which_trial_Hit==trNo)>0
                        
                        
                        %Save the licks for this trial
                        this_Hitii_lick=find(which_trial_Hit==trNo,1);
                        these_Hitii_lick_times=[];
                        these_Hitii_lick_times=Hit_lick_times(this_Hitii_lick,1:Hit_no_lick_times(this_Hitii_lick));
                        if ~isempty(these_Hitii_lick_times)
                            lick_times(num_odor_trials,1:length(these_Hitii_lick_times))=these_Hitii_lick_times;
                            no_licks(num_odor_trials)=length(these_Hitii_lick_times);
                        else
                            no_licks(num_odor_trials)=0;
                        end
                        dLickTraces(num_odor_trials,:)=dHit_lick_traces(this_Hitii_lick,:);
                        
                        epochs_per_trial_dFF(num_odor_trials)=1;
                        
                    end
                end
                
                if epoch_per_trial(trNo)==7
                    %Miss
                    epochs_per_trial(2,num_odor_trials)=1;
                    epochs_per_trial(1,num_odor_trials)=0;
                    epochs_per_trial(3:4,num_odor_trials)=0;
                    all_lda_events_miss_FA(num_odor_trials)=2;
                    
                    
                    if sum(which_trial_Miss==trNo)>0
                        
                        %Save lick times
                        this_Missii_lick=find(which_trial_Miss==trNo,1);
                        these_Missii_lick_times=[];
                        these_Missii_lick_times=Miss_lick_times(this_Missii_lick,1:Miss_no_lick_times(this_Missii_lick));
                        if ~isempty(these_Missii_lick_times)
                            lick_times(num_odor_trials,1:length(these_Missii_lick_times))=these_Missii_lick_times;
                            no_licks(num_odor_trials)=length(these_Missii_lick_times);
                        else
                            no_licks(num_odor_trials)=0;
                        end
                        dLickTraces(num_odor_trials,:)=dMiss_lick_traces(this_Missii_lick,:);
                        
                        epochs_per_trial_dFF(num_odor_trials)=2;
                        
                    end
                end
                
                if epoch_per_trial(trNo)==8
                    %FA
                    epochs_per_trial(3,num_odor_trials)=1;
                    epochs_per_trial(1:2,num_odor_trials)=0;
                    epochs_per_trial(4,num_odor_trials)=0;
                    all_lda_events_miss_FA(num_odor_trials)=4;
                    
                    if sum(which_trial_FA==trNo)>0
                        
                        %Save lick times
                        this_FAii_lick=find(which_trial_FA==trNo,1);
                        these_FAii_lick_times=[];
                        these_FAii_lick_times=FA_lick_times(this_FAii_lick,1:FA_no_lick_times(this_FAii_lick));
                        if ~isempty(these_FAii_lick_times)
                            lick_times(num_odor_trials,1:length(these_FAii_lick_times))=these_FAii_lick_times;
                            no_licks(num_odor_trials)=length(these_FAii_lick_times);
                        else
                            no_licks(num_odor_trials)=0;
                        end
                        dLickTraces(num_odor_trials,:)=dFA_lick_traces(this_FAii_lick,:);
                        
                        epochs_per_trial_dFF(num_odor_trials)=3;
                        
                    end
                end
                
                if epoch_per_trial(trNo)==9
                    %CR
                    epochs_per_trial(4,num_odor_trials)=1;
                    epochs_per_trial(1:3,num_odor_trials)=0;
                    all_lda_events_miss_FA(num_odor_trials)=3;
                    
                    if sum(which_trial_CR==trNo)>0
                        
                        %Save lick times
                        this_CRii_lick=find(which_trial_CR==trNo,1);
                        these_CRii_lick_times=[];
                        these_CRii_lick_times=CR_lick_times(this_CRii_lick,1:CR_no_lick_times(this_CRii_lick));
                        if ~isempty(these_CRii_lick_times)
                            lick_times(num_odor_trials,1:length(these_CRii_lick_times))=these_CRii_lick_times;
                            no_licks(num_odor_trials)=length(these_CRii_lick_times);
                        else
                            no_licks(num_odor_trials)=0;
                        end
                        dLickTraces(num_odor_trials,:)=dCR_lick_traces(this_CRii_lick,:);
                        
                        
                        epochs_per_trial_dFF(num_odor_trials)=4;
                        
                    end
                end
                
            end
            
            
            
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
            
            
            perCorr(1:sliding_window)=perCorr(sliding_window+1);
            
            %Note, this is here so that perCorr=0 is included in the 0-10 % bin.
            perCorr(perCorr==0)=0.00001;
            
            
            %Plot percent correct vs trial
            close all
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
            
            %     %Draw the boundaries of each file
            %     for filNum=2:caimanhandles.caimandr_choices.no_files
            %         %     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
            %         if isfield(caimanhandles.caimandr_choices,'start_reversal')
            %             if caimanhandles.caimandr_choices.start_reversal==filNum
            %                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            %                 text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
            %             end
            %         end
            %         if isfield(caimanhandles.caimandr_choices,'start_gogo')
            %             if caimanhandles.caimandr_choices.start_gogo==filNum
            %                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            %                 text(first_num_odor_trials(filNum)+2,80,'Go-go','Color','k','FontSize',18)
            %             end
            %         end
            %         if isfield(caimanhandles.caimandr_choices,'start_session')
            %             if sum(caimanhandles.caimandr_choices.start_session==filNum)>0
            %                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
            %                 %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
            %             end
            %         end
            %     end
            
            title(['Percent correct vs. trial number ' choiceFileName(18:end-2)])
            xlabel('Trial number')
            ylabel('Percent correct')
            ylim([0 100])
            
            %Save the percent correct
            caimanhandles.per_file(filNum).file_processed=1;
            caimanhandles.per_file(filNum).perCorr=perCorr;
            caimanhandles.per_file(filNum).jj_low=jj_low;
            caimanhandles.per_file(filNum).jj_high=jj_high;
            caimanhandles.per_file(filNum).jj_mid=jj_mid;
            caimanhandles.per_file(filNum).no_odor_trials=no_odor_trials;
            
            %Now plot the intertrial interval to see whether laser on changes the
            %subsequent ITI
            ITI=per_trial_odor_on_time(2:end)-per_trial_odor_on_time(1:end-1);
            ITItrials=[1:no_odor_trials-1];
            
            figNo=figNo+1;
            try
                close(figNo)
            catch
            end
            
            hFig1 = figure(figNo);
            set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])
            
            %Plot ITI for trials with no previous laser
            plot(ITItrials(opto_on_per_trial(1:end-1)==0),ITI(opto_on_per_trial(1:end-1)==0),'ob')
            hold on
            plot(ITItrials(opto_on_per_trial(1:end-1)==1),ITI(opto_on_per_trial(1:end-1)==1),'or')
            title('ITI')
            xlabel('Trial No')
            legend('No laser','Laser')
            ylabel('ITI (sec)')
            
            %Save ITI
            caimanhandles.per_file(filNum).ITI=ITI;
            
            
            if ~isfield(caimanhandles.caimandr_choices,'start_reversal')
                caimanhandles.caimandr_choices.start_reversal=200;
            end
            
            
            
            
            %Now
            % for winNo=1:szwins(1)
            %     handles_sig.win(winNo).ii_for_sig=0;
            % end
            
            if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                %This is a forward run
                total_trial_windows=3;
                supertitle_description{1}=' percent correct <65';
                supertitle_description{2}=' percent correct >=65&<80';
                supertitle_description{3}=' percent correct >=75';
            else
                %Forward and reverse
                total_trial_windows=2;
                supertitle_description{1}=' trials before reversal';
                supertitle_description{2}=' trials after reversal at end of the session';
            end
            
            
            
            trial_window_description{1}='percent correct <65';
            trial_window_description{2}='percent correct >=65&<80';
            trial_window_description{3}='percent correct >=75';
            
            
            
            for no_trial_windows=1:total_trial_windows
                
                caimanhandles.per_file(filNum).trial_window_processed(no_trial_windows)=0;
                
                dFF_trial_mask=[];
                lick_excluded_trials=[];
                jj=0;
                events_miss_FA=[];
                which_trials_in_PCA=[];
                
                if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
                    
                    %             fprintf(1, '\n\nPCA processed for dF/F for trials before reversal \n');
%                     pct_windows=[45 65;65 80;80 100.1];
                    pct_windows=[45 65;65 75;75 100.1];
                    
                    for ii=1:num_odor_trials
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
                        
                        for ii=1:num_odor_trials
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
                    
                    
                    %First plot the lick frequency
                    dt_lick=1/6; %6 Hz
                    time_licks=([1:ceil((dt_after+dt_before)/dt_lick)]*dt_lick)-dt_before;
                    delta_t_gauss=2; %seconds
                    no_conv_points_lick=ceil(delta_t_gauss/(time_licks(2)-time_licks(1)));
                    no_conv_points_dFF=ceil(delta_t_gauss/(1/6));
                    
                    
                    
                    %Figure out the lick threshold to exclude trials wheren Ming
                    %gave the animal water manually during the odor on
                    all_licks_per_dt_per_trial=zeros(num_odor_trials,ceil((dt_after+dt_before)/dt_lick));
                    for trial_no=1:num_odor_trials
                        for ii_lick=1:no_licks(trial_no)
                            all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=all_licks_per_dt_per_trial(trial_no, ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                        end
                    end
                    
                    %
                    %         hit_odor_lick_freq=zeros(1,sum(all_lda_events_miss_FA==1));
                    %         hit_odor_lick_freq(1,:)=sum(all_licks_per_dt_per_trial(all_lda_events_miss_FA==1,(time_licks>=t_odor_on)&(time_licks<=t_odor_off)),2)/(t_odor_off-t_odor_on);
                    %
                    %         lick_threshold=mean(hit_odor_lick_freq)+2*std(hit_odor_lick_freq);
                    %         lick_threshold=20;
                    %          lick_threshold=200;
                    
                    %Calculate the lick frequency for S+ and S-
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
                    
                    Splick_freq_no_conv=Splick_freq;
                    Smlick_freq_no_conv=Smlick_freq;
                    
                    %Convolve lick_freq using a window of 0.9 sec
                    %         no_conv_points=3;
                    %         conv_win=ones(1,no_conv_points);
                    
                    %Convolve lick_freq using a Gaussian window
                    conv_win=gausswin(no_conv_points_lick);
                    
                    Splick_freq=conv(Splick_freq,conv_win,'same')/sum(conv_win);
                    Smlick_freq=conv(Smlick_freq,conv_win,'same')/sum(conv_win);
                    
                    
                    %Plot the lick frequency
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
                    
                    caimanhandles.per_file(filNum).trial_window_processed(no_trial_windows)=1;
                    caimanhandles.per_file(filNum).trial_window(no_trial_windows).time_licks=time_licks;
                    caimanhandles.per_file(filNum).trial_window(no_trial_windows).Smlick_freq=Smlick_freq;
                    caimanhandles.per_file(filNum).trial_window(no_trial_windows).Splick_freq=Splick_freq;
                    caimanhandles.delta_odor_on_reinf_on=mean(delta_odor_on_reinf_on);
                    caimanhandles.delta_reinf=mean(delta_reinf);
                    caimanhandles.delta_odor=mean(delta_odor);
                    
                  
                    
                    %If there were a subset of trials with light on for optogenetics
                    %show the frequency with light on and off
                    if (sum(all_opto_on_per_trial==1)>0)&(sum(all_opto_on_per_trial==0)>0)
                        
                        %Calculate the lick frequency for S+ and S- with light on or
                        %off
                        Splick_freq_on=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                        Splick_freq_off=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                        Smlick_freq_on=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                        Smlick_freq_off=zeros(1,ceil((dt_after+dt_before)/dt_lick));
                        sp_trno_on=0;
                        sp_trno_off=0;
                        sm_trno_on=0;
                        sm_trno_off=0;
                        
                        for trial_no=1:num_odor_trials
                            if dFF_trial_mask(trial_no)==1
                                if all_opto_on_per_trial(trial_no)==1
                                    if strcmp(all_lda_events{trial_no},'S+')
                                        %S+
                                        for ii_lick=1:no_licks(trial_no)
                                            Splick_freq_on( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Splick_freq_on( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                                        end
                                        sp_trno_on=sp_trno_on+1;
                                    else
                                        %S-
                                        for ii_lick=1:no_licks(trial_no)
                                            Smlick_freq_on( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Smlick_freq_on( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                                        end
                                        sm_trno_on=sm_trno_on+1;
                                    end
                                else
                                    if strcmp(all_lda_events{trial_no},'S+')
                                        %S+
                                        for ii_lick=1:no_licks(trial_no)
                                            Splick_freq_off( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Splick_freq_off( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                                        end
                                        sp_trno_off=sp_trno_off+1;
                                    else
                                        %S-
                                        for ii_lick=1:no_licks(trial_no)
                                            Smlick_freq_off( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))=Smlick_freq_off( ceil((lick_times(trial_no,ii_lick)+dt_before)/dt_lick))+1;
                                        end
                                        sm_trno_off=sm_trno_off+1;
                                    end
                                end
                            end
                        end
                        
                        if sp_trno_on>0
                            Splick_freq_on=(Splick_freq_on/(sp_trno_on*dt_lick));
                        end
                        if sp_trno_off>0
                            Splick_freq_off=(Splick_freq_off/(sp_trno_off*dt_lick));
                        end
                        if sm_trno_on>0
                            Smlick_freq_on=(Smlick_freq_on/(sm_trno_on*dt_lick));
                        end
                        if sm_trno_off>0
                            Smlick_freq_off=(Smlick_freq_off/(sm_trno_off*dt_lick));
                        end
                        
                        %Convolve lick_freq using a Gaussian window
                        conv_win=gausswin(no_conv_points_lick);
                        
                        Splick_freq_on=conv(Splick_freq_on,conv_win,'same')/sum(conv_win);
                        Smlick_freq_on=conv(Smlick_freq_on,conv_win,'same')/sum(conv_win);
                        Splick_freq_off=conv(Splick_freq_off,conv_win,'same')/sum(conv_win);
                        Smlick_freq_off=conv(Smlick_freq_off,conv_win,'same')/sum(conv_win);
                        
                        
                        %Plot the lick frequency for light on and off
                        figNo=figNo+1;
                        try
                            close(figNo)
                        catch
                        end
                        
                        figure(figNo)
                        
                        hold on
                        
                        lfreqmax=20;
                        
                        p1=plot(time_licks,Smlick_freq_off,'-','LineWidth',2,'Color',[0.7 0.7 1]);
                        p2=plot(time_licks,Smlick_freq_on,'-b','LineWidth',2);
                        p3=plot(time_licks,Splick_freq_off,'-','LineWidth',2,'Color',[1 0.7 0.7]);
                        p4=plot(time_licks,Splick_freq_on,'-r','LineWidth',2);
                        
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
                        legend([p1 p2 p3 p4],'S- light off','S- light on','S+ light off','S+ light on')
                        xlim([-10 20])
                        ylim([-1 lfreqmax])
                        
                    end
                    
                    %Now plot the p value for the difference in licks between S+ and S-
                    %Get lick p values
                    try
                        trials_used=10;
                        no_pvals=0;
                        p_val_Sp_Sm=[];
                        time_p_lick=[];
                        fractional_lick_time_on_sp=zeros(1,length(time_licks));
                        fractional_lick_time_on_sm=zeros(1,length(time_licks));
                        
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
                            
                            if ii==1
                                no_sptrial_start=caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_sptrials+1;
                                no_smtrial_start=caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_smtrials+1;
                                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_sptrials=...
                                    caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_sptrials+sp_trno;
                                caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_smtrials=...
                                    caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_smtrials+sm_trno;
                                no_sptrial_end=caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_sptrials;
                                no_smtrial_end=caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).no_smtrials;
                            end
                            
                            caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).splicks(ii,no_sptrial_start:no_sptrial_end)=this_Sp;
                            caimanhandles.per_mouse(caimanhandles.caimandr_choices.mouse(filNum)).hM4D(caimanhandles.caimandr_choices.hM4D(filNum)).smlicks(ii,no_smtrial_start:no_smtrial_end)=this_Sm;
                            
                            no_pvals=no_pvals+1;
                            p_val_Sp_Sm(no_pvals)=1;
                            
                            if (~isempty(this_Sm))&(~isempty(this_Sp))
                                p_val_done=0;
                                no_these_p_vals=0;
                                ii_s=1;
                                these_pvals=[];
                                while p_val_done==0
                                    if ii_s+trials_used-1>length(this_Sp)
                                        p_val_done=1;
                                    else
                                        if ii_s+trials_used-1>length(this_Sm)
                                            p_val_done=1;
                                        else
                                            no_these_p_vals=no_these_p_vals+1;
                                            these_pvals(no_these_p_vals)=ranksum(this_Sm(ii_s:ii_s+trials_used-1),this_Sp(ii_s:ii_s+trials_used-1));
                                            ii_s=ii_s+trials_used;
                                        end
                                    end
                                end
                                 p_val_Sp_Sm(no_pvals)=mean(these_pvals);
                             else
                                p_val_Sp_Sm(no_pvals)=1;
                            end
                               
                            
%                             if (~isempty(this_Sm))&(~isempty(this_Sp))
%                                 p_val_Sp_Sm(no_pvals)=ranksum(this_Sm,this_Sp);
%                             else
%                                 p_val_Sp_Sm(no_pvals)=1;
%                             end
                            
                            %ranksum gives NaN if the values are all the same
                            if isnan(p_val_Sp_Sm(no_pvals))
                                p_val_Sp_Sm(no_pvals)=1;
                            end
                            
                            time_p_lick(no_pvals)= (time_licks(ii)+time_licks(ii+1))/2;
                            
                            fractional_lick_time_on_sp(ii)=mean(this_Sp);
                            fractional_lick_time_on_sm(ii)=mean(this_Sm);
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
                        
                        
                        %Save the logp data
                        caimanhandles.per_file(filNum).trial_window(no_trial_windows).log10p_val_Sp_Sm=log10(p_val_Sp_Sm);
                        caimanhandles.per_file(filNum).trial_window(no_trial_windows).time_p_lick=time_p_lick;
                        
                        %Find the decision time for odor detection
                        ii_zero=find(time_p_lick>0,1,'first');
                        decision_ii=find(p_val_Sp_Sm(time_p_lick>0)<0.05,1,'first');
                        if ~isempty(decision_ii)
                            decision_time=time_p_lick(decision_ii-1+ii_zero);
                        else
                            decision_time=time_p_lick(end);
                        end
                        %Note, if the decision time is longer than 10 seconds it is assigned 10 +
                        %gaussian noise (so that they do not overlap)
                        if decision_time>10
                            decision_time=10+randn;
                        end
                        caimanhandles.per_file(filNum).trial_window(no_trial_windows).tdecision_time=decision_time;
                        
                        %Find the decision time for water detection (this will only work for <65%)
                        dt_reinf_on=mean(delta_odor_on_reinf_on);
                        ii_zero_reinf=find(time_p_lick>dt_reinf_on,1,'first');
                        decision_ii_reinf=find(p_val_Sp_Sm(time_p_lick>dt_reinf_on)<0.05,1,'first');
                        if ~isempty(decision_ii_reinf)
                            decision_time_reinf=time_p_lick(decision_ii_reinf-1+ii_zero_reinf);
                        else
                            decision_time_reinf=time_p_lick(end);
                        end
                        decision_time_reinf=decision_time_reinf-dt_reinf_on;
                         if decision_time_reinf>10
                            decision_time_reinf=10+randn;
                        end
                        
                        %These if statements are here to exclude NaNs and zeros
                        if (sum(isnan(Splick_freq))==0)&(sum(isnan(Smlick_freq))==0)&(length(this_Sp)>10)&(length(this_Sm)>10)  
                            no_files_included=no_files_included+1;
                            per_file_lick_log10_p_val(no_files_included,1:length(log10(p_val_Sp_Sm)))=log10(p_val_Sp_Sm);
                            per_file_lick_freq_sp(no_files_included,1:length(Splick_freq_no_conv))=Splick_freq_no_conv;
                            per_file_lick_freq_sm(no_files_included,1:length(Smlick_freq_no_conv))=Smlick_freq_no_conv;
                            per_file_hM4D(no_files_included)=caimanhandles.caimandr_choices.hM4D(filNum);
                            per_file_mouse(no_files_included)=caimanhandles.caimandr_choices.mouse(filNum);
                            per_file_noSp(no_files_included)=length(this_Sp);
                            per_file_noSm(no_files_included)=length(this_Sm);
                            per_file_trial_window(no_files_included)=no_trial_windows; 
                            per_file_decision_time(no_files_included)=decision_time;  
                            per_file_decision_time_reinf(no_files_included)=decision_time_reinf;
                            per_file_fractional_lick_time_on_sp(no_files_included,1:length(fractional_lick_time_on_sp))=fractional_lick_time_on_sp;
                            per_file_fractional_lick_time_on_sm(no_files_included,1:length(fractional_lick_time_on_sm))=fractional_lick_time_on_sm;
                        end
                    catch
                        pffft=1;
                    end
                    
                    
                    %If there were a subset of trials with light on for optogenetics
                    %show the p values for light on vs off
                    if (sum(all_opto_on_per_trial==1)>0)&(sum(all_opto_on_per_trial==0)>0)
                        
                        %S+ light on vs off
                        figNo=figNo+1;
                        try
                            close(figNo)
                        catch
                        end
                        
                        figure(figNo)
                        
                        try
                            
                            no_pvals=0;
                            p_vals=[];
                            time_p_lick=[];
                            
                            for ii=1:length(time_licks)-1
                                
                                sp_trno_on=0;
                                sp_trno_off=0;
                                this_Sp_on=[];
                                this_Sp_off=[];
                                
                                for trial_no=1:num_odor_trials
                                    %                     if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                                    if dFF_trial_mask(trial_no)==1
                                        if strcmp(all_lda_events{trial_no},'S+')
                                            if all_opto_on_per_trial(trial_no)==1
                                                %S+ light on
                                                sp_trno_on=sp_trno_on+1;
                                                this_Sp_on(sp_trno_on,1)=0;
                                                for ii_lick=1:no_licks(trial_no)
                                                    if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                                        this_Sp_on(sp_trno_on,1)=1;
                                                    end
                                                end
                                            else
                                                %S+ light off
                                                sp_trno_off=sp_trno_off+1;
                                                this_Sp_off(sp_trno_off,1)=0;
                                                for ii_lick=1:no_licks(trial_no)
                                                    if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                                        this_Sp_off(sp_trno_off,1)=1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    %                     end
                                end
                                
                                
                                no_pvals=no_pvals+1;
                                
                                if (~isempty(this_Sp_on))&(~isempty(this_Sp_off))
                                    p_vals(no_pvals)=ranksum(this_Sp_on,this_Sp_off);
                                else
                                    p_vals(no_pvals)=1;
                                end
                                
                                %ranksum gives NaN if the values are all the same
                                if isnan(p_vals(no_pvals))
                                    p_vals(no_pvals)=1;
                                end
                                
                                time_p_lick(no_pvals)= (time_licks(ii)+time_licks(ii+1))/2;
                                
                            end
                            
                            logpmin=-20;
                            
                            subplot(2,1,1)
                            hold on
                            plot(time_p_lick,log10(p_vals),'-k','LineWidth',2)
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
                                    title(['log(p value) for lick difference S+ light on vs light off for %correct <65']);
                                case 2
                                    title(['log(p value) for lick difference S+ light on vs light off for %correct >=65&<80']);
                                case 3
                                    title(['log(p value) for lick difference S+ light on vs light off for %correct >80']);
                            end
                            
                            xlabel('Time (sec)')
                            ylabel('log10(p value)')
                            ylim([logpmin 0.5])
                            xlim([-10 20])
                            
                        catch
                        end
                        
                        %S- light on vs off
                        try
                            
                            no_pvals=0;
                            p_vals=[];
                            time_p_lick=[];
                            
                            for ii=1:length(time_licks)-1
                                
                                sm_trno_on=0;
                                sm_trno_off=0;
                                this_Sm_on=[];
                                this_Sm_off=[];
                                
                                for trial_no=1:num_odor_trials
                                    %                     if sum((lick_times(trial_no,1:no_licks(trial_no))>=t_odor_on)&(lick_times(trial_no,1:no_licks(trial_no))<=t_odor_off))<lick_threshold
                                    if dFF_trial_mask(trial_no)==1
                                        if strcmp(all_lda_events{trial_no},'S-')
                                            if all_opto_on_per_trial(trial_no)==1
                                                %S+ light on
                                                sm_trno_on=sm_trno_on+1;
                                                this_Sm_on(sm_trno_on,1)=0;
                                                for ii_lick=1:no_licks(trial_no)
                                                    if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                                        this_Sm_on(sm_trno_on,1)=1;
                                                    end
                                                end
                                            else
                                                %S+ light off
                                                sm_trno_off=sm_trno_off+1;
                                                this_Sm_off(sm_trno_off,1)=0;
                                                for ii_lick=1:no_licks(trial_no)
                                                    if (lick_times(trial_no,ii_lick)>=time_licks(ii))&(lick_times(trial_no,ii_lick)<time_licks(ii+1))
                                                        this_Sm_off(sm_trno_off,1)=1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    %                     end
                                end
                                
                                
                                no_pvals=no_pvals+1;
                                
                                if (~isempty(this_Sm_on))&(~isempty(this_Sm_off))
                                    p_vals(no_pvals)=ranksum(this_Sm_on,this_Sm_off);
                                else
                                    p_vals(no_pvals)=1;
                                end
                                
                                %ranksum gives NaN if the values are all the same
                                if isnan(p_vals(no_pvals))
                                    p_vals(no_pvals)=1;
                                end
                                
                                time_p_lick(no_pvals)= (time_licks(ii)+time_licks(ii+1))/2;
                                
                            end
                            
                            logpmin=-20;
                            
                            subplot(2,1,2)
                            hold on
                            plot(time_p_lick,log10(p_vals),'-k','LineWidth',2)
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
                                    title(['log(p value) for lick difference S- light on vs light off for %correct <65']);
                                case 2
                                    title(['log(p value) for lick difference S- light on vs light off for %correct >=65&<80']);
                                case 3
                                    title(['log(p value) for lick difference S- light on vs light off for %correct >80']);
                            end
                            
                            xlabel('Time (sec)')
                            ylabel('log10(p value)')
                            ylim([logpmin 0.5])
                            xlim([-10 20])
                            
                        catch
                        end
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
                    
                    if (no_files_included>0)
                        if (per_file_hM4D(no_files_included)==1)&(per_file_trial_window(no_files_included)==3)
                            pffft=1;
                        end
                    end
                    
                    %Plot lick frequency and the derivative of lick frequency per trial
                    %Plot the licks
                    figNo=figNo+1;
                    try
                        close(figNo)
                    catch
                    end
                    
                    figure(figNo)
                    hold on
                    
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
                            
                            
                            subplot(2,1,1)
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
                            plot([0+t_offset 0+t_offset],[ymin ymax],'-k')
                            odorhl=plot([0+t_offset mean(delta_odor)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-k','LineWidth',5);
                            plot([mean(delta_odor)+t_offset mean(delta_odor)+t_offset],[ymin ymax],'-k')
                            
                            %Reinforcement markers
                            plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+t_offset],[ymin ymax],'-r')
                            reinfhl=plot([mean(delta_odor_on_reinf_on)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin + 0.1*(ymax-ymin) ymin + 0.1*(ymax-ymin)],'-r','LineWidth',5);
                            plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset mean(delta_odor_on_reinf_on)+mean(delta_reinf)+t_offset],[ymin ymax],'-r')
                            
                            %Calculate the derivative of lick_freq
                            lick_freq_dx=gradient(lick_freq);
                            
                            %Plot the derivative
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
                            
                            
                            handles_outs.no_lick_slopes=handles_outs.no_lick_slopes+1;
                            handles_outs.lick_freq(handles_outs.no_lick_slopes,1:132)=lick_freq(1:132);
                            handles_outs.lick_derivatives(handles_outs.no_lick_slopes,1:132)=lick_freq_dx(1:132);
                            handles_outs.time_licks=time_licks(1:132);
                            
                            t_offset=t_offset+35;
                            
                        end
                    end
                    
                    
                    subplot(2,1,1)
                    xlim([-50 350])
                    ylim([ymin ymax])
                    xlabel('Time (sec)')
                    ylabel('LR')
                    title(['Timecourse for lick rate for ' supertitle_description{no_trial_windows}])
                    
                    subplot(2,1,2)
                    xlim([-50 350])
                    ylim([ymin_dx ymax_dx])
                    xlabel('Time (sec)')
                    ylabel('dLR/dt')
                    title(['Timecourse for the drivative of lick rate for ' supertitle_description{no_trial_windows}])
                    
                    
                end
                
            end
            fprintf(1, ['\nProcessed file number %d ' caimanhandles.caimandr_choices.fileName{filNum} '\n'],filNum);
            switch caimanhandles.caimandr_choices.hM4D(filNum)
                case 1
                    fprintf(1, ['\nhM4Di with CNO\n\n'],filNum);
                case 2
                    fprintf(1, ['\nhM4Di no CNO\n\n'],filNum);
                case 3
                    fprintf(1, ['\control with CNO\n\n'],filNum);
                case 4
                    fprintf(1, ['\control no CNO\n\n'],filNum);
            end
            pffft=1;
        else
            caimanhandles.per_file(filNum).file_processed=0;
            fprintf(1, ['\nDid not pocess file number %d ' caimanhandles.caimandr_choices.fileName{filNum} ' number of trials <20\n'],filNum);
        end
        
    else
        caimanhandles.per_file(filNum).file_processed=0;
        fprintf(1, ['\nFile number %d ' caimanhandles.caimandr_choices.fileName{filNum} ' does not exist\n'],filNum);
    end
end
   
figNo
% save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) '_licks.mat'],'handles_outs')
 
%Now plot the figures

%Plot decision time

if display_choice==1
    %Find decision time per mouse per trial window
    for mouseNo=1:caimanhandles.caimandr_choices.no_mice
        these_decision_times=[];
        these_no_decision_time=zeros(1,total_trial_windows);
        for filNum=1:caimanhandles.caimandr_choices.no_files
            for trial_window_no=1:total_trial_windows
                if ~isempty(caimanhandles.per_file(filNum).trial_window_processed)
                    if (caimanhandles.per_file(filNum).trial_window_processed(trial_window_no)==1)&(caimanhandles.caimandr_choices.mouse_no(filNum)==mouseNo)
                        these_no_decision_time(trial_window_no)=these_no_decision_time(trial_window_no)+1;
                        these_decision_times(trial_window_no,these_no_decision_time(trial_window_no))=caimanhandles.per_file(filNum).trial_window(trial_window_no).tdecision_time;
                    end
                end
            end
        end
        %For now sort these and calculate the average of the three smallest
        for trial_window_no=1:total_trial_windows
            if these_no_decision_time(trial_window_no)<3
                caimanhandles.per_mouse(mouseNo).decision_time(trial_window_no)=mean(these_decision_times(trial_window_no,1:these_no_decision_time(trial_window_no)));
            else
                sortedts=sort(these_decision_times(trial_window_no,1:these_no_decision_time(trial_window_no)));
                caimanhandles.per_mouse(mouseNo).decision_time(trial_window_no)=mean(sortedts(1:3));
            end
            if caimanhandles.caimandr_choices.hM4D(mouseNo)==1
                plot(trial_window_no,caimanhandles.per_mouse(mouseNo).decision_time(trial_window_no),'or')
            else
                plot(trial_window_no,caimanhandles.per_mouse(mouseNo).decision_time(trial_window_no),'ob')
            end
        end
    end
    
    %Plot the decision time bar graph
    
    title('Decision times')
    ylabel('Decision time (sec)')
    
    for mouseNo=1:caimanhandles.caimandr_choices.no_mice
        these_decision_times=[];
        these_no_decision_timecoursesSm=[];
        these_no_decision_timecoursesSp=[];
        for filNum=1:caimanhandles.caimandr_choices.no_files
            for trial_window_no=1:total_trial_windows
                if ~isempty(caimanhandles.per_file(filNum).trial_window_processed)
                    if (caimanhandles.per_file(filNum).trial_window_processed(trial_window_no)==1)&(caimanhandles.caimandr_choices.mouse_no(filNum)==mouseNo)
                        these_no_decision_time(trial_window_no)=these_no_decision_time(trial_window_no)+1;
                        these_decision_timecoursesSm(trial_window_no,these_no_decision_time(trial_window_no),:)=caimanhandles.per_file(filNum).trial_window(trial_window_no).Smlick_freq;
                        these_decision_timecoursesSp(trial_window_no,these_no_decision_time(trial_window_no),:)=caimanhandles.per_file(filNum).trial_window(trial_window_no).Splick_freq;
                    end
                end
            end
        end
        %alculate the average of the time courses
        caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSm=zeros(total_trial_windows,size(these_decision_timecoursesSm,3));
        caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSp=zeros(total_trial_windows,size(these_decision_timecoursesSp,3));
        time_licks=caimanhandles.per_file(filNum).trial_window(trial_window_no).time_licks;
        for trial_window_no=1:total_trial_windows
            caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSm(trial_window_no,:)=mean(these_decision_timecoursesSm(trial_window_no,1:these_no_decision_time(trial_window_no),:),2);
            caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSp(trial_window_no,:)=mean(these_decision_timecoursesSp(trial_window_no,1:these_no_decision_time(trial_window_no),:),2);
        end
    end
    
    %PLot lick freq
    for trial_window_no=1:total_trial_windows
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        for mouseNo=1:caimanhandles.caimandr_choices.no_mice
            if caimanhandles.caimandr_choices.hM4D(mouseNo)==1
                phM4DSm=plot(time_licks,caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSm(trial_window_no,:),'-b');
            else
                pcontSm=plot(time_licks,caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSm(trial_window_no,:),'-c');
            end
            if caimanhandles.caimandr_choices.hM4D(mouseNo)==1
                phM4DSp=plot(time_licks,caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSp(trial_window_no,:),'-r');
            else
                pcontSp=plot(time_licks,caimanhandles.per_mouse(mouseNo).mean_decision_timecourseSp(trial_window_no,:),'-m');
            end
        end
        legend([phM4DSp phM4DSm pcontSp pcontSm],{'hM4Di Sp','hM4Di Sm','Control Sp','Control Sm'})
        title(['Lick frequency ' supertitle_description{trial_window_no}])
        ylabel('Frequency (Hz)')
        xlabel('Time (sec)')
    end
    
end

show_CI=1;

if display_choice==2
    
    %Plot lick frequency
    %     for trial_window_no=1:total_trial_windows
    trial_window_no=3;
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .4 .6 .25])
    
    subplot(1,2,1)
    hold on
    
    %Plot the bounded lines
    if show_CI==1
        for hM4DNo=2:-1:1
            
            
            mean_freq_Sp=mean(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            no_trials=size(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            these_per_file_lick_freq_sp=[];
            these_per_file_lick_freq_sp=per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
            CIfreq = bootci(1000, {@mean, these_per_file_lick_freq_sp})';
            CIfreq(:,1)=mean_freq_Sp'-CIfreq(:,1);
            CIfreq(:,2)=CIfreq(:,2)-mean_freq_Sp';
            
            if hM4DNo==1
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sp', CIfreq, 'r');
            else
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sp', CIfreq, 'm');
            end
            
            mean_freq_Sm=mean(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            no_trials=size(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            these_per_file_lick_freq_sm=[];
            these_per_file_lick_freq_sm=per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
            CIfreq = bootci(1000, {@mean, these_per_file_lick_freq_sm})';
            CIfreq(:,1)=mean_freq_Sm'-CIfreq(:,1);
            CIfreq(:,2)=CIfreq(:,2)-mean_freq_Sm';
            
            if hM4DNo==1
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sm', CIfreq, 'b');
            else
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sm', CIfreq, 'c');
            end
        end
    end
    
    %Plot the thick lines and do the glm
    glm_lick_f=[];
    glm_lick_f_ii=0;
    
    glm_sp_lick_f=[];
    glm_sp_lick_f_ii=0;
    
    glm_sm_lick_f=[];
    glm_sm_lick_f_ii=0;
    
    %Plot the thick lines and do glm
    for hM4DNo=2:-1:1
        
        mean_freq_Sp=mean(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        no_trials=size(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        these_per_file_lick_freq_sp=[];
        these_per_file_lick_freq_sp=per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
        for ii=1:no_trials
            glm_lick_f.data(glm_lick_f_ii+1:glm_lick_f_ii+180)=these_per_file_lick_freq_sp(ii,:);
            if hM4DNo==1
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
            else
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=0;
            end
            glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=2;
            glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=2;
            glm_lick_f.time(glm_lick_f_ii+1:glm_lick_f_ii+180)=time_licks;
            glm_lick_f_ii=glm_lick_f_ii+180;
            
            glm_sp_lick_f.data(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=these_per_file_lick_freq_sp(ii,:);
            if hM4DNo==1
                glm_sp_lick_f.CNO(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=1;
            else
                glm_sp_lick_f.CNO(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=0;
            end
            glm_sp_lick_f.genotype(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=2;
            glm_sp_lick_f.time(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=time_licks;
            glm_sp_lick_f_ii=glm_sp_lick_f_ii+180;
        end
        if hM4DNo==1
            p2=plot(time_licks,mean_freq_Sp,'-r','LineWidth',2);
        else
            p4=plot(time_licks,mean_freq_Sp,'-m','LineWidth',2);
        end
        
        mean_freq_Sm=mean(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        no_trials=size(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        these_per_file_lick_freq_sm=[];
        these_per_file_lick_freq_sm=per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
        for ii=1:no_trials
            glm_lick_f.data(glm_lick_f_ii+1:glm_lick_f_ii+180)=these_per_file_lick_freq_sm(ii,:);
            if hM4DNo==1
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=2;
            else
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=0;
                glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=2;
            end
            glm_lick_f.time(glm_lick_f_ii+1:glm_lick_f_ii+180)=time_licks;
            glm_lick_f_ii=glm_lick_f_ii+180;
            
            glm_sm_lick_f.data(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=these_per_file_lick_freq_sm(ii,:);
            if hM4DNo==1
                glm_sm_lick_f.CNO(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=1;
            else
                glm_sm_lick_f.CNO(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=0;
            end
            glm_sm_lick_f.genotype(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=2;
            glm_sm_lick_f.time(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=time_licks;
            glm_sm_lick_f_ii=glm_sm_lick_f_ii+180;
        end
        if hM4DNo==1
            p1=plot(time_licks,mean_freq_Sm,'-b','LineWidth',2);
        else
            p3=plot(time_licks,mean_freq_Sm,'-c','LineWidth',2);
        end
        
    end
    
    %Odor on markers
    plot([0 0],[-2 16],'-k')
    odorhl=plot([0 mean(delta_odor)],[-1 -1],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[-2 16],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[-2 16],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-1 -1],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-2 16],'-r')
    
    
    legend([p1 p2 p3 p4],{'CNO Sm','CNO Sp','No CNO Sm','No CNO Sp'})
    title('hM4Di')
    ylim([-2 16])
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    
    subplot(1,2,2)
    hold on
    
    %Plot the bounded lines
    if show_CI==1
        for hM4DNo=4:-1:3
            
            mean_freq_Sp=mean(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            no_trials=size(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            these_per_file_lick_freq_sp=[];
            these_per_file_lick_freq_sp=per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
            CIfreq = bootci(1000, {@mean, these_per_file_lick_freq_sp})';
            CIfreq(:,1)=mean_freq_Sp'-CIfreq(:,1);
            CIfreq(:,2)=CIfreq(:,2)-mean_freq_Sp';
            
            if hM4DNo==3
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sp', CIfreq, 'r');
            else
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sp', CIfreq, 'm');
            end
            
            mean_freq_Sm=mean(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            no_trials=size(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
            these_per_file_lick_freq_sm=[];
            these_per_file_lick_freq_sm=per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
            CIfreq = bootci(1000, {@mean, these_per_file_lick_freq_sm})';
            CIfreq(:,1)=mean_freq_Sm'-CIfreq(:,1);
            CIfreq(:,2)=CIfreq(:,2)-mean_freq_Sm';
            
            if hM4DNo==3
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sm', CIfreq, 'b');
            else
                [hlCR, hpCR] = boundedline(time_licks,mean_freq_Sm', CIfreq, 'c');
            end
        end
    end
    
    %Plot the thick lines and do glm
    for hM4DNo=4:-1:3
        
        mean_freq_Sp=mean(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        no_trials=size(per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        these_per_file_lick_freq_sp=[];
        these_per_file_lick_freq_sp=per_file_lick_freq_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
        for ii=1:no_trials
            glm_lick_f.data(glm_lick_f_ii+1:glm_lick_f_ii+180)=these_per_file_lick_freq_sp(ii,:);
            if hM4DNo==3
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
            else
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=0;
            end
            glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=2;
            glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
            glm_lick_f.time(glm_lick_f_ii+1:glm_lick_f_ii+180)=time_licks;
            glm_lick_f_ii=glm_lick_f_ii+180;
            
            glm_sp_lick_f.data(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=these_per_file_lick_freq_sp(ii,:);
            if hM4DNo==3
                glm_sp_lick_f.CNO(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=1;
            else
                glm_sp_lick_f.CNO(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=0;
            end
            glm_sp_lick_f.genotype(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=1;
            glm_sp_lick_f.time(glm_sp_lick_f_ii+1:glm_sp_lick_f_ii+180)=time_licks;
            glm_sp_lick_f_ii=glm_sp_lick_f_ii+180;
        end
        if hM4DNo==3
            p2=plot(time_licks,mean_freq_Sp,'-r','LineWidth',2);
        else
            p4=plot(time_licks,mean_freq_Sp,'-m','LineWidth',2);
        end
        
        mean_freq_Sm=mean(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        no_trials=size(per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
        these_per_file_lick_freq_sm=[];
        these_per_file_lick_freq_sm=per_file_lick_freq_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
        for ii=1:no_trials
            glm_lick_f.data(glm_lick_f_ii+1:glm_lick_f_ii+180)=these_per_file_lick_freq_sm(ii,:);
            if hM4DNo==3
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
            else
                glm_lick_f.CNO(glm_lick_f_ii+1:glm_lick_f_ii+180)=0;
                glm_lick_f.spm(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
                glm_lick_f.genotype(glm_lick_f_ii+1:glm_lick_f_ii+180)=1;
            end
            glm_lick_f.time(glm_lick_f_ii+1:glm_lick_f_ii+180)=time_licks;
            glm_lick_f_ii=glm_lick_f_ii+180;
            
            glm_sm_lick_f.data(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=these_per_file_lick_freq_sm(ii,:);
            if hM4DNo==3
                glm_sm_lick_f.CNO(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=1;
            else
                glm_sm_lick_f.CNO(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=0;
            end
            glm_sm_lick_f.genotype(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=1;
            glm_sm_lick_f.time(glm_sm_lick_f_ii+1:glm_sm_lick_f_ii+180)=time_licks;
            glm_sm_lick_f_ii=glm_sm_lick_f_ii+180;
        end
        if hM4DNo==3
            p1=plot(time_licks,mean_freq_Sm,'-b','LineWidth',2);
        else
            p3=plot(time_licks,mean_freq_Sm,'-c','LineWidth',2);
        end
        
        
    end
    
    %Odor on markers
    plot([0 0],[-2 16],'-k')
    odorhl=plot([0 mean(delta_odor)],[-1 -1],'-k','LineWidth',5);
    plot([mean(delta_odor) mean(delta_odor)],[-2 16],'-k')
    
    %Reinforcement markers
    plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[-2 16],'-r')
    reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-1 -1],'-r','LineWidth',5);
    plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[-2 16],'-r')
    
    legend([p1 p2 p3 p4],{'CNO Sm','CNO Sp','No CNO Sm','No CNO Sp'})
    ylim([-2 16])
    title('Control')
    ylabel('Frequency (Hz)')
    xlabel('Time (sec)')
    
    suptitle(['Lick frequency ' supertitle_description{trial_window_no}])
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick frequency timecourse for both spms ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_lick_f.data',glm_lick_f.spm',glm_lick_f.CNO',glm_lick_f.genotype',glm_lick_f.time',...
        'VariableNames',{'lick_frequency','spm','CNO','genotype','time'});
    mdl = fitglm(tbl,'lick_frequency~spm+CNO+genotype+time+spm*CNO*genotype*time'...
        ,'CategoricalVars',[2,3,4])
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick frequency timecourse for S+ ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_sp_lick_f.data',glm_sp_lick_f.CNO',glm_sp_lick_f.genotype',glm_sp_lick_f.time',...
        'VariableNames',{'lick_frequency','CNO','genotype','time'});
    mdl = fitglm(tbl,'lick_frequency~CNO+genotype+time+CNO*genotype*time'...
        ,'CategoricalVars',[2,3])
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick frequency timecourse for S- ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_sm_lick_f.data',glm_sm_lick_f.CNO',glm_sm_lick_f.genotype',glm_sm_lick_f.time',...
        'VariableNames',{'lick_frequency','CNO','genotype','time'});
    mdl = fitglm(tbl,'lick_frequency~CNO+genotype+time+CNO*genotype*time'...
        ,'CategoricalVars',[2,3])
    
    
    %     %Plot lick proabability
    %     for trial_window_no=1:total_trial_windows
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         hFig=figure(figNo);
    %         set(hFig, 'units','normalized','position',[.1 .4 .6 .25])
    %
    %         subplot(1,2,1)
    %         hold on
    %
    %         %Plot the bounded lines
    %         if show_CI==1
    %             for hM4DNo=1:2
    %                 mean_prob_Sm=mean(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_prob_sm=[];
    %                 these_per_file_lick_prob_sm=per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIprob = bootci(1000, {@mean, these_per_file_lick_prob_sm})';
    %                 CIprob(:,1)=mean_prob_Sm'-CIprob(:,1);
    %                 CIprob(:,2)=CIprob(:,2)-mean_prob_Sm';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sm', CIprob, 'b');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sm', CIprob, 'c');
    %                 end
    %
    %                 mean_prob_Sp=mean(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_prob_sp=[];
    %                 these_per_file_lick_prob_sp=per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIprob = bootci(1000, {@mean, these_per_file_lick_prob_sp})';
    %                 CIprob(:,1)=mean_prob_Sp'-CIprob(:,1);
    %                 CIprob(:,2)=CIprob(:,2)-mean_prob_Sp';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sp', CIprob, 'r');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sp', CIprob, 'm');
    %                 end
    %             end
    %         end
    %
    %         %Plot the thick lines and do the glm
    %         glm_lick_prob=[];
    %         glm_lick_prob_ii=0;
    %         for hM4DNo=1:2
    %             mean_prob_Sm=mean(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_prob_sm=[];
    %             these_per_file_lick_prob_sm=per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_prob.data(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=these_per_file_lick_prob_sm(ii,:);
    %                 if hM4DNo==1
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                 else
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=0;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                 end
    %                 glm_lick_prob.time(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=time_licks;
    %                 glm_lick_prob_ii=glm_lick_prob_ii+180;
    %             end
    %             if hM4DNo==1
    %                 p1=plot(time_licks,mean_prob_Sm,'-b','LineWidth',2);
    %             else
    %                 p3=plot(time_licks,mean_prob_Sm,'-c','LineWidth',2);
    %             end
    %
    %             mean_prob_Sp=mean(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_prob_sp=[];
    %             these_per_file_lick_prob_sp=per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_prob.data(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=these_per_file_lick_prob_sp(ii,:);
    %                 if hM4DNo==1
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                 else
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=0;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                 end
    %                 glm_lick_prob.time(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=time_licks;
    %                 glm_lick_prob_ii=glm_lick_prob_ii+180;
    %             end
    %             if hM4DNo==1
    %                 p2=plot(time_licks,mean_prob_Sp,'-r','LineWidth',2);
    %             else
    %                 p4=plot(time_licks,mean_prob_Sp,'-m','LineWidth',2);
    %             end
    %         end
    %         legend([p1 p2 p3 p4],{'CNO Sm','CNO Sp','No CNO Sm','No CNO Sp'})
    %         title('hM4Di')
    %         ylabel('Frequency (Hz)')
    %         xlabel('Time (sec)')
    %
    %         subplot(1,2,2)
    %         hold on
    %
    %         %Plot the bounded lines
    %         if show_CI==1
    %             for hM4DNo=3:4
    %                 mean_prob_Sm=mean(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_prob_sm=[];
    %                 these_per_file_lick_prob_sm=per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIprob = bootci(1000, {@mean, these_per_file_lick_prob_sm})';
    %                 CIprob(:,1)=mean_prob_Sm'-CIprob(:,1);
    %                 CIprob(:,2)=CIprob(:,2)-mean_prob_Sm';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sm', CIprob, 'b');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sm', CIprob, 'c');
    %                 end
    %
    %                 mean_prob_Sp=mean(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_prob_sp=[];
    %                 these_per_file_lick_prob_sp=per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIprob = bootci(1000, {@mean, these_per_file_lick_prob_sp})';
    %                 CIprob(:,1)=mean_prob_Sp'-CIprob(:,1);
    %                 CIprob(:,2)=CIprob(:,2)-mean_prob_Sp';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sp', CIprob, 'r');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_licks,mean_prob_Sp', CIprob, 'm');
    %                 end
    %             end
    %         end
    %
    %         %Plot the thick lines and calcualte the glm
    %         for hM4DNo=3:4
    %             mean_prob_Sm=mean(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_prob_sm=[];
    %             these_per_file_lick_prob_sm=per_file_fractional_lick_time_on_sm((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_prob.data(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=these_per_file_lick_prob_sm(ii,:);
    %                 if hM4DNo==3
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                 else
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=0;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                 end
    %                 glm_lick_prob.time(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=time_licks;
    %                 glm_lick_prob_ii=glm_lick_prob_ii+180;
    %             end
    %             if hM4DNo==3
    %                 p1=plot(time_licks,mean_prob_Sm,'-b','LineWidth',2);
    %             else
    %                 p3=plot(time_licks,mean_prob_Sm,'-c','LineWidth',2);
    %             end
    %             mean_prob_Sp=mean(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_prob_sp=[];
    %             these_per_file_lick_prob_sp=per_file_fractional_lick_time_on_sp((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_prob.data(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=these_per_file_lick_prob_sp(ii,:);
    %                 if hM4DNo==3
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=1;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                 else
    %                     glm_lick_prob.CNO(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=0;
    %                     glm_lick_prob.spm(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                     glm_lick_prob.genotype(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=2;
    %                 end
    %                 glm_lick_prob.time(glm_lick_prob_ii+1:glm_lick_prob_ii+180)=time_licks;
    %                 glm_lick_prob_ii=glm_lick_prob_ii+180;
    %             end
    %             if hM4DNo==3
    %                 p2=plot(time_licks,mean_prob_Sp,'-r','LineWidth',2);
    %             else
    %                 p4=plot(time_licks,mean_prob_Sp,'-m','LineWidth',2);
    %             end
    %         end
    %         legend([p1 p2 p3 p4],{'CNO Sm','CNO Sp','No CNO Sm','No CNO Sp'})
    %         title('Control')
    %         ylabel('Probability')
    %         xlabel('Time (sec)')
    %
    %         suptitle(['Lick probability ' supertitle_description{trial_window_no}])
    %
    %         %Perform the glm
    %         fprintf(1, ['\n\nglm for lick probability timecourse for ' trial_window_description{trial_window_no} '\n'])
    %         tbl = table(glm_lick_prob.data',glm_lick_prob.spm',glm_lick_prob.CNO',glm_lick_prob.genotype',glm_lick_prob.time',...
    %             'VariableNames',{'lick_probability','spm','CNO','genotype','time'});
    %         mdl = fitglm(tbl,'lick_probability~spm+CNO+genotype+time+spm*CNO*genotype*time'...
    %             ,'CategoricalVars',[2,3,4])
    %     end
    
    %     %Plot p value timecourses for lick decision
    %     for trial_window_no=1:total_trial_windows
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         hFig=figure(figNo);
    %         set(hFig, 'units','normalized','position',[.1 .4 .6 .25])
    %
    %         subplot(1,2,1)
    %         hold on
    %
    %         %Plot the bounded lines
    %         if show_CI==1
    %             for hM4DNo=1:2
    %                 mean_log10_pval=mean(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_log10_p_val=[];
    %                 these_per_file_lick_log10_p_val=per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIpval = bootci(1000, {@mean, these_per_file_lick_log10_p_val})';
    %                 CIpval(:,1)=mean_log10_pval'-CIpval(:,1);
    %                 CIpval(:,2)=CIpval(:,2)-mean_log10_pval';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_p_lick,mean_log10_pval', CIpval, 'r');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_p_lick,mean_log10_pval', CIpval, 'b');
    %                 end
    %
    %             end
    %         end
    %
    %         %Plot the thick lines and glm
    %         glm_lick_p=[];
    %         glm_lick_p_ii=0;
    %
    %         for hM4DNo=1:2
    %             mean_log10_pval=mean(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_log10_p_val=[];
    %             these_per_file_lick_log10_p_val=per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_p.data(glm_lick_p_ii+1:glm_lick_p_ii+179)=these_per_file_lick_log10_p_val(ii,:);
    %                 if hM4DNo==1
    %                     glm_lick_p.CNO(glm_lick_p_ii+1:glm_lick_p_ii+179)=1;
    %                     glm_lick_p.genotype(glm_lick_p_ii+1:glm_lick_p_ii+179)=1;
    %                 else
    %                     glm_lick_p.CNO(glm_lick_p_ii+1:glm_lick_p_ii+179)=0;
    %                     glm_lick_p.genotype(glm_lick_p_ii+1:glm_lick_p_ii+179)=1;
    %                 end
    %                 glm_lick_p.time(glm_lick_p_ii+1:glm_lick_p_ii+179)=time_p_lick;
    %                 glm_lick_p_ii=glm_lick_p_ii+179;
    %             end
    %             if hM4DNo==1
    %                 p1=plot(time_p_lick,mean_log10_pval,'-r','LineWidth',2);
    %             else
    %                 p2=plot(time_p_lick,mean_log10_pval,'-b','LineWidth',2);
    %             end
    %         end
    %
    %         plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
    %
    %         %Odor on markers
    %         plot([0 0],[logpmin 0.5],'-k')
    %         odorhl=plot([0 mean(delta_odor)],[logpmin+0.5 logpmin+0.5],'-k','LineWidth',5);
    %         plot([mean(delta_odor) mean(delta_odor)],[logpmin 0.5],'-k')
    %
    %         %Reinforcement markers
    %         plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[logpmin 0.5],'-r')
    %         reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin+0.5 logpmin+0.5],'-r','LineWidth',5);
    %         plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin 0.5],'-r')
    %
    %
    %         legend([p1 p2],{'CNO','No CNO'})
    %         title(['h4DMi'])
    %         ylim([-20 1])
    %         ylabel('Log10(pvalue)')
    %         xlabel('Time (sec)')
    %
    %         subplot(1,2,2)
    %         hold on
    %
    %         %Plot the bounded lines
    %         if show_CI==1
    %             for hM4DNo=3:4
    %                 mean_log10_pval=mean(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 no_trials=size(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %                 these_per_file_lick_log10_p_val=[];
    %                 these_per_file_lick_log10_p_val=per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %                 CIpval = bootci(1000, {@mean, these_per_file_lick_log10_p_val})';
    %                 CIpval(:,1)=mean_log10_pval'-CIpval(:,1);
    %                 CIpval(:,2)=CIpval(:,2)-mean_log10_pval';
    %
    %                 if hM4DNo==1
    %                     [hlCR, hpCR] = boundedline(time_p_lick,mean_log10_pval', CIpval, 'r');
    %                 else
    %                     [hlCR, hpCR] = boundedline(time_p_lick,mean_log10_pval', CIpval, 'b');
    %                 end
    %
    %             end
    %         end
    %
    %         for hM4DNo=3:4
    %             mean_log10_pval=mean(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             no_trials=size(per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:),1);
    %             these_per_file_lick_log10_p_val=[];
    %             these_per_file_lick_log10_p_val=per_file_lick_log10_p_val((per_file_hM4D==hM4DNo)&(per_file_trial_window==trial_window_no),:);
    %             for ii=1:no_trials
    %                 glm_lick_p.data(glm_lick_p_ii+1:glm_lick_p_ii+179)=these_per_file_lick_log10_p_val(ii,:);
    %                 if hM4DNo==3
    %                     glm_lick_p.CNO(glm_lick_p_ii+1:glm_lick_p_ii+179)=1;
    %                     glm_lick_p.genotype(glm_lick_p_ii+1:glm_lick_p_ii+179)=2;
    %                 else
    %                     glm_lick_p.CNO(glm_lick_p_ii+1:glm_lick_p_ii+179)=0;
    %                     glm_lick_p.genotype(glm_lick_p_ii+1:glm_lick_p_ii+179)=2;
    %                 end
    %                 glm_lick_p.time(glm_lick_p_ii+1:glm_lick_p_ii+179)=time_p_lick;
    %                 glm_lick_p_ii=glm_lick_p_ii+179;
    %             end
    %             if hM4DNo==3
    %                 p1=plot(time_p_lick,mean_log10_pval,'-r','LineWidth',2);
    %             else
    %                 p2=plot(time_p_lick,mean_log10_pval,'-b','LineWidth',2);
    %             end
    %         end
    %
    %         plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)],'-r','LineWidth',1)
    %
    %         %Odor on markers
    %         plot([0 0],[logpmin 0.5],'-k')
    %         odorhl=plot([0 mean(delta_odor)],[logpmin+0.5 logpmin+0.5],'-k','LineWidth',5);
    %         plot([mean(delta_odor) mean(delta_odor)],[logpmin 0.5],'-k')
    %
    %         %Reinforcement markers
    %         plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[logpmin 0.5],'-r')
    %         reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin+0.5 logpmin+0.5],'-r','LineWidth',5);
    %         plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[logpmin 0.5],'-r')
    %
    %
    %         legend([p1 p2],{'CNO','No CNO'})
    %         ylim([-20 1])
    %         title(['Control'])
    %         ylabel('Frequency (Hz)')
    %         xlabel('Time (sec)')
    %
    %         suptitle(['Log10(p value) decision making for ' supertitle_description{trial_window_no}])
    %
    %         %Perform the glm
    %         fprintf(1, ['\n\nglm for decision making p value timecourse for ' trial_window_description{trial_window_no} '\n'])
    %         tbl = table(glm_lick_p.data',glm_lick_p.CNO',glm_lick_p.genotype',glm_lick_p.time',...
    %             'VariableNames',{'lick_p_value','CNO','genotype','time'});
    %         mdl = fitglm(tbl,'lick_p_value~CNO+genotype+time+CNO*genotype*time'...
    %             ,'CategoricalVars',[2,3])
    %     end
    
    %     %Plot the number of trials per session
    %     for trial_window_no=1:total_trial_windows
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         hFig=figure(figNo);
    %         set(hFig, 'units','normalized','position',[.1 .4 .6 .25])
    %
    %         %Plot cumultaive histo for S+ trials
    %         subplot(1,2,1)
    %         hold on
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSp((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'r')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSp((per_file_hM4D==2)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'b')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSp((per_file_hM4D==3)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'m')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSp((per_file_hM4D==4)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'c')
    %
    %
    %         legend('h4DMi CNO','h4DMi no CNO','control CNO','control no CNO')
    %         title(['S+'])
    %         ylim([0 1.1])
    %         ylabel('Probability')
    %         xlabel('Number of trials')
    %
    %         %Plot cumultaive histo for S- trials
    %         subplot(1,2,2)
    %         hold on
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSm((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'r')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSm((per_file_hM4D==2)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'b')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSm((per_file_hM4D==3)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'m')
    %
    %         [f_trNo,x_trNo] = drg_ecdf(per_file_noSm((per_file_hM4D==4)&(per_file_trial_window==trial_window_no)));
    %         plot(x_trNo,f_trNo,'c')
    %
    %
    %         legend('h4DMi CNO','h4DMi no CNO','control CNO','control no CNO')
    %         title(['S-'])
    %         ylim([0 1.1])
    %         ylabel('Probability')
    %         xlabel('Number of trials')
    %
    %         suptitle(['Number of trials per session' trial_window_description{trial_window_no}])
    %
    %     end
    
    %     %Now plot decision times
    %     for trial_window_no=1:total_trial_windows
    %
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         figure(figNo)
    %         hold on
    %
    %         glm_dt_ii=0;
    %         glm_dt=[];
    %
    %         %Control, no CNO
    %         these_dts=[];
    %         these_dts=per_file_decision_time((per_file_hM4D==4)&(per_file_trial_window==trial_window_no));
    %         bar(1,mean(these_dts),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    %         CIdt = bootci(1000, {@mean, these_dts})';
    %         plot([1 1],CIdt,'-k','LineWidth',2)
    %         plot(ones(1,length(these_dts)),these_dts,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dt.data(glm_dt_ii+1:glm_dt_ii+length(these_dts))=these_dts;
    %         glm_dt.genotype(glm_dt_ii+1:glm_dt_ii+length(these_dts))=1;
    %         glm_dt.CNO(glm_dt_ii+1:glm_dt_ii+length(these_dts))=0;
    %         glm_dt_ii=glm_dt_ii+length(these_dts);
    %
    %         %Control, CNO
    %         these_dts=[];
    %         these_dts=per_file_decision_time((per_file_hM4D==3)&(per_file_trial_window==trial_window_no));
    %         bar(2,mean(these_dts),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    %         CIdt = bootci(1000, {@mean, these_dts})';
    %         plot([2 2],CIdt,'-k','LineWidth',2)
    %         plot(2*ones(1,length(these_dts)),these_dts,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dt.data(glm_dt_ii+1:glm_dt_ii+length(these_dts))=these_dts;
    %         glm_dt.genotype(glm_dt_ii+1:glm_dt_ii+length(these_dts))=1;
    %         glm_dt.CNO(glm_dt_ii+1:glm_dt_ii+length(these_dts))=1;
    %         glm_dt_ii=glm_dt_ii+length(these_dts);
    %
    %         %hM4D, no CNO
    %         these_dts=[];
    %         these_dts=per_file_decision_time((per_file_hM4D==2)&(per_file_trial_window==trial_window_no));
    %         bar(4,mean(these_dts),'FaceColor','b','EdgeColor','b')
    %         CIdt = bootci(1000, {@mean, these_dts})';
    %         plot([4 4],CIdt,'-k','LineWidth',2)
    %         plot(4*ones(1,length(these_dts)),these_dts,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dt.data(glm_dt_ii+1:glm_dt_ii+length(these_dts))=these_dts;
    %         glm_dt.genotype(glm_dt_ii+1:glm_dt_ii+length(these_dts))=2;
    %         glm_dt.CNO(glm_dt_ii+1:glm_dt_ii+length(these_dts))=0;
    %         glm_dt_ii=glm_dt_ii+length(these_dts);
    %
    %         %hM4D, CNO
    %         these_dts=[];
    %         these_dts=per_file_decision_time((per_file_hM4D==1)&(per_file_trial_window==trial_window_no));
    %         bar(5,mean(these_dts),'FaceColor','r','EdgeColor','r')
    %         CIdt = bootci(1000, {@mean, these_dts})';
    %         plot([5 5],CIdt,'-k','LineWidth',2)
    %         plot(5*ones(1,length(these_dts)),these_dts,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dt.data(glm_dt_ii+1:glm_dt_ii+length(these_dts))=these_dts;
    %         glm_dt.genotype(glm_dt_ii+1:glm_dt_ii+length(these_dts))=2;
    %         glm_dt.CNO(glm_dt_ii+1:glm_dt_ii+length(these_dts))=1;
    %         glm_dt_ii=glm_dt_ii+length(these_dts);
    %
    %         xticks([1 2 4 5])
    %         ylim([0 15])
    %         xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    %         title(['Decision times for' trial_window_description{trial_window_no}])
    %         ylabel('Decision time (sec)')
    %
    %         %Perform the glm
    %         fprintf(1, ['\n\nglm for decision times for ' trial_window_description{trial_window_no} '\n'])
    %         tbl = table(glm_dt.data',glm_dt.genotype',glm_dt.CNO',...
    %             'VariableNames',{'decision_times','genotype','CNO'});
    %         mdl = fitglm(tbl,'decision_times~genotype+CNO+genotype*CNO'...
    %             ,'CategoricalVars',[2,3])
    %
    %         fprintf(1, ['Number of sessions and mice:\n'])
    %         fprintf(1, ['hM4d CNO: %d, %d\n'],sum((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)),...
    %             length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    %         fprintf(1, ['hM4d no CNO: %d, %d\n'],sum((per_file_hM4D==2)&(per_file_trial_window==trial_window_no)),...
    %             length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    %         fprintf(1, ['control CNO: %d, %d\n'],sum((per_file_hM4D==3)&(per_file_trial_window==trial_window_no)),...
    %             length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    %         fprintf(1, ['control no CNO: %d, %d\n'],sum((per_file_hM4D==4)&(per_file_trial_window==trial_window_no)),...
    %             length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    %
    %     end
    %
    %     %Now plot decision times for water reinforcement
    %     %This is only done for <65% correct
    %     for trial_window_no=1:1
    %
    %         figNo=figNo+1;
    %         try
    %             close(figNo)
    %         catch
    %         end
    %
    %         figure(figNo)
    %         hold on
    %
    %         glm_dtr_ii=0;
    %         glm_dtr=[];
    %
    %         %Control, no CNO
    %         these_dtrs=[];
    %         these_dtrs=per_file_decision_time_reinf((per_file_hM4D==4)&(per_file_trial_window==trial_window_no));
    %         bar(1,mean(these_dtrs),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    %         CIdt = bootci(1000, {@mean, these_dtrs})';
    %         plot([1 1],CIdt,'-k','LineWidth',2)
    %         plot(1*ones(1,length(these_dtrs)),these_dtrs,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dtr.data(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=these_dtrs;
    %         glm_dtr.genotype(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=1;
    %         glm_dtr.CNO(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=0;
    %         glm_dtr_ii=glm_dtr_ii+length(these_dtrs);
    %
    %         %Control, CNO
    %         these_dtrs=[];
    %         these_dtrs=per_file_decision_time_reinf((per_file_hM4D==3)&(per_file_trial_window==trial_window_no));
    %         bar(2,mean(these_dtrs),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    %         CIdt = bootci(1000, {@mean, these_dtrs})';
    %         plot([2 2],CIdt,'-k','LineWidth',2)
    %         plot(2*ones(1,length(these_dtrs)),these_dtrs,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dtr.data(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=these_dtrs;
    %         glm_dtr.genotype(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=1;
    %         glm_dtr.CNO(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=1;
    %         glm_dtr_ii=glm_dtr_ii+length(these_dtrs);
    %
    %         %hM4D, no CNO
    %         these_dtrs=[];
    %         these_dtrs=per_file_decision_time_reinf((per_file_hM4D==2)&(per_file_trial_window==trial_window_no));
    %         bar(4,mean(these_dtrs),'FaceColor','b','EdgeColor','b')
    %         CIdt = bootci(1000, {@mean, these_dtrs})';
    %         plot([4 4],CIdt,'-k','LineWidth',2)
    %         plot(4*ones(1,length(these_dtrs)),these_dtrs,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dtr.data(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=these_dtrs;
    %         glm_dtr.genotype(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=2;
    %         glm_dtr.CNO(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=0;
    %         glm_dtr_ii=glm_dtr_ii+length(these_dtrs);
    %
    %         %hM4D, CNO
    %         these_dtrs=[];
    %         these_dtrs=per_file_decision_time_reinf((per_file_hM4D==1)&(per_file_trial_window==trial_window_no));
    %         bar(5,mean(these_dtrs),'FaceColor','r','EdgeColor','r')
    %         CIdt = bootci(1000, {@mean, these_dtrs})';
    %         plot([5 5],CIdt,'-k','LineWidth',2)
    %         plot(5*ones(1,length(these_dtrs)),these_dtrs,'o','Color',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
    %
    %         glm_dtr.data(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=these_dtrs;
    %         glm_dtr.genotype(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=2;
    %         glm_dtr.CNO(glm_dtr_ii+1:glm_dtr_ii+length(these_dtrs))=1;
    %         glm_dtr_ii=glm_dtr_ii+length(these_dtrs);
    %
    %         xticks([1 2 4 5])
    %         xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    %         ylim([0 15])
    %         title(['Decision times for water detection for ' trial_window_description{trial_window_no}])
    %         ylabel('Decision time (sec)')
    %
    %         %Perform the glm
    %         fprintf(1, ['\n\nglm for decision times for water detection for ' trial_window_description{trial_window_no} '\n'])
    %         tbl = table(glm_dtr.data',glm_dtr.genotype',glm_dtr.CNO',...
    %             'VariableNames',{'decision_time_reinf','genotype','CNO'});
    %         mdl = fitglm(tbl,'decision_time_reinf~genotype+CNO+genotype*CNO'...
    %             ,'CategoricalVars',[2,3])
    %
    %     end
    %
    
    %Now plot bar graphs for fraction of the time that the animl is licking
    trial_window_no=3;
    
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    figure(figNo)
    
    glm_flt_ii=0;
    glm_flt=[];
    
    %S+ trials
    subplot(1,2,1)
    hold on
    
    %Control, no CNO
    flt_timecourses=per_file_fractional_lick_time_on_sp((per_file_hM4D==4)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(1,mean_prob,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    plot([1 1],CIprob,'-k','LineWidth',2)
    plot(ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %Control, CNO
    flt_timecourses=per_file_fractional_lick_time_on_sp((per_file_hM4D==3)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(2,mean_prob,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    plot([2 2],CIprob,'-k','LineWidth',2)
    plot(2*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %hM4D, no CNO
    flt_timecourses=per_file_fractional_lick_time_on_sp((per_file_hM4D==2)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(4,mean_prob,'b')
    plot([4 4],CIprob,'-k','LineWidth',2)
    plot(4*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %hM4D, CNO
    flt_timecourses=per_file_fractional_lick_time_on_sp((per_file_hM4D==1)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(5,mean_prob,'r')
    plot([5 5],CIprob,'-k','LineWidth',2)
    plot(5*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    xticks([1 2 4 5])
    ylim([0 1.1])
    xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    title('S+')
    ylabel('Lick on time')
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick on time for S+ for ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_flt.data',glm_flt.genotype',glm_flt.CNO',...
        'VariableNames',{'lick_on_time','genotype','CNO'});
    mdl = fitglm(tbl,'lick_on_time~genotype+CNO+genotype*CNO'...
        ,'CategoricalVars',[2,3])
    
    
    
    glm_flt_ii=0;
    glm_flt=[];
    
    %S- trials
    subplot(1,2,2)
    hold on
    
    %Control, no CNO
    flt_timecourses=per_file_fractional_lick_time_on_sm((per_file_hM4D==4)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';;
    bar(1,mean_prob,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    plot([1 1],CIprob,'-k','LineWidth',2)
    plot(ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %Control, CNO
    flt_timecourses=per_file_fractional_lick_time_on_sm((per_file_hM4D==3)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(2,mean_prob,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    plot([2 2],CIprob,'-k','LineWidth',2)
    plot(2*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %hM4D, no CNO
    flt_timecourses=per_file_fractional_lick_time_on_sm((per_file_hM4D==2)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(4,mean_prob,'b')
    plot([4 4],CIprob,'-k','LineWidth',2)
    plot(4*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    %hM4D, CNO
    flt_timecourses=per_file_fractional_lick_time_on_sm((per_file_hM4D==1)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=0.8)&(time_licks<=1.8)),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(5,mean_prob,'r')
    plot([5 5],CIprob,'-k','LineWidth',2)
    plot(5*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    xticks([1 2 4 5])
    ylim([0 1.1])
    xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    title('S-')
    ylabel('Lick on time')
    
    suptitle(['Lick on time for' trial_window_description{trial_window_no}])
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick on time for ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_flt.data',glm_flt.genotype',glm_flt.CNO',...
        'VariableNames',{'lick_on_time','genotype','CNO'});
    mdl = fitglm(tbl,'lick_on_time~genotype+CNO+genotype*CNO'...
        ,'CategoricalVars',[2,3])
    
    
    
    fprintf(1, ['Number of sessions and mice:\n'])
    fprintf(1, ['hM4d CNO: %d, %d\n'],sum((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['hM4d no CNO: %d, %d\n'],sum((per_file_hM4D==2)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['control CNO: %d, %d\n'],sum((per_file_hM4D==3)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['control no CNO: %d, %d\n'],sum((per_file_hM4D==4)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    
    
    
    %Now plot bar graphs for lick frequency
    %0.8-1.8 sec
    t_start=[0.8 5.2];
    t_end=[1.8 6.2];
    t_win_label{1}='odor';
    t_win_label{2}='reward';
    
    
    trial_window_no=3
    ii_t_win=1;
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    figure(figNo)
    
    glm_flt_ii=0;
    glm_flt=[];
    rst_ii=0;
    rst_data=[];
    
    %S+ trials
    subplot(1,2,1)
    hold on
    
    %Control, no CNO
    flt_timecourses=per_file_lick_freq_sp((per_file_hM4D==4)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(1,mean_prob,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    plot([1 1],CIprob,'-k','LineWidth',2)
    plot(ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='Control, no CNO';
    
    %Control, CNO
    flt_timecourses=per_file_lick_freq_sp((per_file_hM4D==3)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(2,mean_prob,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    plot([2 2],CIprob,'-k','LineWidth',2)
    plot(2*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='Control, CNO';
    
    %hM4D, no CNO
    flt_timecourses=per_file_lick_freq_sp((per_file_hM4D==2)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(4,mean_prob,'b')
    plot([4 4],CIprob,'-k','LineWidth',2)
    plot(4*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='hM4Di, no CNO';
    
    %hM4D, CNO
    flt_timecourses=per_file_lick_freq_sp((per_file_hM4D==1)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(5,mean_prob,'r')
    plot([5 5],CIprob,'-k','LineWidth',2)
    plot(5*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='HM4Di, CNO';
    
    xticks([1 2 4 5])
    ylim([0 13])
    xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    title('S+')
    ylabel('Frequency (Hz)')
    
    %Perform the glm
    fprintf(1, ['\n\nglm for frequency for S+ for ' t_win_label{ii_t_win} ' ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_flt.data',glm_flt.genotype',glm_flt.CNO',...
        'VariableNames',{'frequency','genotype','CNO'});
    mdl = fitglm(tbl,'frequency~genotype+CNO+genotype*CNO'...
        ,'CategoricalVars',[2,3])
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for frequency S+ for ' t_win_label{ii_t_win} ' ' trial_window_description{trial_window_no} '\n'])
    [output_data] = drgMutiRanksumorTtest(rst_data);
    
    glm_flt_ii=0;
    glm_flt=[];
    rst_ii=0;
    rst_data=[];
    
    %S- trials
    subplot(1,2,2)
    hold on
    
    %Control, no CNO
    flt_timecourses=per_file_lick_freq_sm((per_file_hM4D==4)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';;
    bar(1,mean_prob,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
    plot([1 1],CIprob,'-k','LineWidth',2)
    plot(ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='Control, no CNO';
    
    %Control, CNO
    flt_timecourses=per_file_lick_freq_sm((per_file_hM4D==3)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(2,mean_prob,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
    plot([2 2],CIprob,'-k','LineWidth',2)
    plot(2*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='Control, CNO';
    
    %hM4D, no CNO
    flt_timecourses=per_file_lick_freq_sm((per_file_hM4D==2)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(4,mean_prob,'b')
    plot([4 4],CIprob,'-k','LineWidth',2)
    plot(4*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='hM4Di, no CNO';
    
    %hM4D, CNO
    flt_timecourses=per_file_lick_freq_sm((per_file_hM4D==1)&(per_file_trial_window==trial_window_no),:);
    flts=mean(flt_timecourses(:,(time_licks>=t_start(ii_t_win))&(time_licks<=t_end(ii_t_win))),2);
    mean_prob=mean(flts);
    CIprob = bootci(1000, {@mean, flts})';
    bar(5,mean_prob,'r')
    plot([5 5],CIprob,'-k','LineWidth',2)
    plot(5*ones(1,length(flts)),flts,'ok')
    
    glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
    glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
    glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
    glm_flt_ii=glm_flt_ii+length(flts);
    
    rst_ii=rst_ii+1;
    rst_data(rst_ii).data=flts;
    rst_data(rst_ii).description='HM4Di, CNO';
    
    xticks([1 2 4 5])
    ylim([0 13])
    xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
    title('S-')
    ylabel('Frequency (Hz)')
    
    suptitle(['Lick frequency for ' t_win_label{ii_t_win} ' ' trial_window_description{trial_window_no}])
    
    %Perform the glm
    fprintf(1, ['\n\nglm for lick frequency for S- for ' t_win_label{ii_t_win} ' ' trial_window_description{trial_window_no} '\n'])
    tbl = table(glm_flt.data',glm_flt.genotype',glm_flt.CNO',...
        'VariableNames',{'frequency','genotype','CNO'});
    mdl = fitglm(tbl,'frequency~genotype+CNO+genotype*CNO'...
        ,'CategoricalVars',[2,3])
    
    %Do the ranksum/t-test
    fprintf(1, ['\n\nRanksum or t-test p values for frequency S- for ' t_win_label{ii_t_win} ' ' trial_window_description{trial_window_no} '\n'])
    [output_data] = drgMutiRanksumorTtest(rst_data);
    
    
    
    fprintf(1, ['Number of sessions and mice:\n'])
    fprintf(1, ['hM4d CNO: %d, %d\n'],sum((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['hM4d no CNO: %d, %d\n'],sum((per_file_hM4D==2)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['control CNO: %d, %d\n'],sum((per_file_hM4D==3)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    fprintf(1, ['control no CNO: %d, %d\n'],sum((per_file_hM4D==4)&(per_file_trial_window==trial_window_no)),...
        length(unique(per_file_mouse((per_file_hM4D==1)&(per_file_trial_window==trial_window_no)))))
    
    pffft=1;
    
end


pffft=1;
