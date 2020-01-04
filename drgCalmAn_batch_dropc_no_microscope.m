%% drgCaImAn_batch_dropc_no_microscope.m
% This is used to process files acquired in a session to study the changes
% in behavior and licks when the MLIs are activated with optogenetics
% Needs as an input the ouput files from the dropcspm_hf.m (.mat) and from
% INTAN (.rhd)
%
close all
clear all

%Choices



% 1 dropc_nsampler_piriform

% 2 dropcspm_hf before 2/23/2018

% 3 dropcspm_hf after 2/24/2018
%
%  handles.dropcData.epochEvent
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR

dropc_program=3;

% Read choices file

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_dropc_choices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles_choice=' choiceFileName(1:end-2) ';'])
handles_choice.choiceFileName=choiceFileName;
handles_choice.choiceBatchPathName=choiceBatchPathName;




for fileNo=handles_choice.first_file:handles_choice.no_files
    try
        fprintf(1, ['\nProcessing  ' num2str(fileNo) '\n']);
        
        if iscell(handles_choice.PathName)
            cd(handles_choice.PathName{fileNo})
        else
            cd(handles_choice.PathName)
        end
        
        
        close all
        
        
        fprintf(1, ['\ndrgCaImAn_batch_dropc_no_microscope run for ' handles_choice.spmFileName{fileNo} '\n\n']);
        
        %Read the dropc file
        handles=[];
        load(handles_choice.spmFileName{fileNo})
        
        %Read the rhd file
        adc_in=[];
        digital_in=[];
        acq_rate=[];
        [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(handles_choice.rhdFileName{fileNo},3);
        
        %Get the events
        odor_on_times=[];
        ootii=0;
        y_shift=100;
        switch dropc_program
            case 1
                for event=2:handles.dropcData.eventIndex
                    plot([handles.dropcData.eventTime(event) handles.dropcData.eventTime(event)], [0 (no_traces+2)*y_shift],...
                        these_lines{handles.dropcData.odorNo(event)},'LineWidth',1)
                end
            case 2
                for event=1:handles.dropcData.allTrialIndex
                    plot([handles.dropcData.allTrialTime(event)-2.5 handles.dropcData.allTrialTime(event)-2.5], [0 (no_traces+2)*y_shift],...
                        these_lines{handles.dropcData.odorType(event)},'LineWidth',1)
                end
            case 3
                %For S+ and S- plot odor on and reinforcement
                for epoch=1:handles.dropcData.epochIndex
                    %Epoch 2 is odor on, 3 is odor off
                    plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
                    if plot_epoch
                        
                        if (handles.dropcData.epochEvent(epoch)==2)
                            ootii=ootii+1;
                            odor_on_times(ootii)=handles.dropcData.epochTime(epoch);
                        end
                    end
                    
                    
                end
        end
        
        
        %Align the rhd times with the olfactometer
        
        %Find the FV, odor on and odor off events in digital_in recorded by INTAN
        ii=1;
        at_end=0;
        odor_on_times_rhd=[];
        FV_times_rhd=[];
        odor_off_times_rhd=[];
        iioon=0;
        iiFV=0;
        iiooff=0;
        opto_on=bitand(digital_in,64);
        digital_in=bitand(digital_in,2+4+8+16);
        while at_end==0
            ii_FV=find(digital_in(ii:end)==6,1,'first');
            if isempty(ii_FV)
                at_end=1;
            else
                %FV
                ii=ii+ii_FV-1;
                iiFV=iiFV+1;
                FV_times_rhd(iiFV)=ii/acq_rate;
                
                %Note: There was a problem with transmission of bit 16 in some of Ming's
                %recordings. As a result, the value following FV(6) was 2
                
                if isempty(find(digital_in(ii:end)==18,1,'first'))
                    %Data missing bit 16?
                    
                    %Odor on
                    ii_odor_on=find(digital_in(ii:end)==2,1,'first');
                    %Odor off
                    ii_odor_off=find(digital_in(ii:end)~=2,1,'first');
                    
                    if (~isempty(ii_odor_on))&(~isempty(ii_odor_off))
                        
                        %Odor on
                        ii=ii+ii_odor_on-1;
                        iioon=iioon+1;
                        odor_on_times_rhd(iioon)=ii/acq_rate;
                        
                        %Odor off
                        
                        ii=ii+ii_odor_off-1;
                        iiooff=iiooff+1;
                        odor_off_times_rhd(iiooff)=ii/acq_rate;
                        
                        ii=ii+1;
                        if ii>=length(digital_in)
                            at_end=1;
                        end
                    else
                        at_end=1;
                    end
                    
                else
                    %Odor on
                    ii_odor_on=find(digital_in(ii:end)==18,1,'first');
                    %Odor off
                    ii_odor_off=find(digital_in(ii:end)<18,1,'first');
                    
                    if (~isempty(ii_odor_on))&(~isempty(ii_odor_off))
                        
                        %Odor on
                        ii=ii+ii_odor_on-1;
                        iioon=iioon+1;
                        odor_on_times_rhd(iioon)=ii/acq_rate;
                        
                        %Odor off
                        
                        ii=ii+ii_odor_off-1;
                        iiooff=iiooff+1;
                        odor_off_times_rhd(iiooff)=ii/acq_rate;
                        
                        ii=ii+1;
                        if ii>=length(digital_in)
                            at_end=1;
                        end
                    else
                        at_end=1;
                    end
                    
                end
            end
        end
        
        %Do opto on
        at_end=0;
        ii_OpOn=0;
        ii=1;
        OpOn_times_rhd=[];
        while at_end==0
            ii_Op_On=find(opto_on(ii:end)==64,1,'first');
            if isempty(ii_Op_On)
                at_end=1;
            else
                %Opto on
                ii=ii+ii_Op_On-1;
                ii_OpOn=ii_OpOn+1;
                OpOn_times_rhd(ii_OpOn)=ii/acq_rate;
                
                %Opto on off
                ii_Op_off=find(opto_on(ii:end)==0,1,'first');
                
                if (~isempty(ii_Op_off))
                    
                    
                    ii=ii+ii_Op_off-1;
                    
                    ii=ii+1;
                    if ii>=length(opto_on)
                        at_end=1;
                    end
                else
                    at_end=1;
                end
            end
        end
        
        %     %Figures for troubleshootimg
        %     %Are the times aligned?
        %     figure(1)
        %     hold on
        %     plot(odor_on_times,ones(1,length(odor_on_times)),'or')
        %     plot(odor_on_times_rhd,1.2*ones(1,length(odor_on_times_rhd)),'ob')
        %
        %     %Is the INTAN reading the digital input from the olfactometer
        %     figure(2)
        %     plot(digital_in(200000:600000))
        
        %Find the alignment of the rhd vs the olfactometer times
        if length(odor_on_times)<length(odor_on_times_rhd)
            sum_delta=[];
            for ii=0:length(odor_on_times_rhd)-length(odor_on_times)
                sum_delta(ii+1)=abs(sum(odor_on_times_rhd(1+ii:ii+length(odor_on_times))-odor_on_times));
            end
            [min_del min_jj]=min(sum_delta);
            odor_on_times_rhd=odor_on_times_rhd(min_jj:min_jj+length(odor_on_times)-1);
        end
        
        
        
        delta_t_rhd=mean(odor_on_times-odor_on_times_rhd);
        
        %     %Find the trials when opto was on
        %     opto_on_per_trial=zeros(1,length(odor_on_times));
        %     for ii=1:length(OpOn_times_rhd)
        %         if sum(abs(OpOn_times_rhd(ii)-odor_on_times_rhd)<3)>0
        %             this_trial=find(abs(OpOn_times_rhd(ii)-odor_on_times_rhd)<3,1,'first');
        %             opto_on_per_trial(1,this_trial)=1;
        %         end
        %     end
        
        
        %Plot the licks recorded by the INTAN (adc_in)
        figure(1)
        time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
        OpOn_times_rhd=OpOn_times_rhd+delta_t_rhd;
        pct998=prctile(adc_in,99.8);
        pct1=prctile(adc_in,1);
        norm_fact=0.8*y_shift/(pct998-pct1);
        
        plot(time_rhd(time_rhd>0),adc_in(time_rhd>0)*norm_fact)
        
        xlabel('time (s)')
        ylabel('Licks')
        title(handles_choice.rhdFileName{fileNo})
        
        dt_before=10;
        dt_after=20;
        
        
        %Plot the responses aligned with the onset of the epochs
        switch dropc_program
            case 3
                dt_odor_onset=0.1085;  %This is the time from FV off to odor entering the nose cone
                
                timesSD=5;
                timesSDodorOn=2.5;
                
                response_points=1*dt_after;
                odor_response_points=1;
                dt_trace=2000000;
                dt_rhd_trace=1000000000000;
                
                
                %Determine the odor on and reinforcement times
                iido=0;
                delta_odor=[];
                iidr=0;
                delta_reinf=[];
                iidro=0;
                delta_odor_on_reinf_on=[];
                
                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                        
                        if (handles.dropcData.epochEvent(epoch)==4)
                            %This is a reinforcement on, find the odor on and
                            %the reinforcement off
                            
                            %Find reinforcement off
                            next_epoch=epoch+1;
                            while handles.dropcData.epochEvent(next_epoch)~=5
                                next_epoch=next_epoch+1;
                            end
                            iidr=iidr+1;
                            delta_reinf(iidr)=handles.dropcData.epochTime(next_epoch)-handles.dropcData.epochTime(epoch);
                            
                            %Find odor on
                            next_epoch=epoch-1;
                            while handles.dropcData.epochEvent(next_epoch)~=2
                                next_epoch=next_epoch-1;
                            end
                            iidro=iidro+1;
                            delta_odor_on_reinf_on(iidro)=handles.dropcData.epochTime(epoch)-handles.dropcData.epochTime(next_epoch);
                        end
                        
                        if (handles.dropcData.epochEvent(epoch)==2)
                            %This is an odor on event, find the next odor off
                            next_epoch=epoch+1;
                            while handles.dropcData.epochEvent(next_epoch)~=3
                                next_epoch=next_epoch+1;
                            end
                            iido=iido+1;
                            delta_odor(iido)=handles.dropcData.epochTime(next_epoch)-handles.dropcData.epochTime(epoch);
                        end
                    end
                end
                
                
                %Extract all licks and determine the threshold
                all_lick_traces=[];
                allii_lick=0;
                
                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                        if (handles.dropcData.epochEvent(epoch)==2)
                            %Odor on
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            allii_lick=allii_lick+1;
                            all_lick_traces(allii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                        end
                    end
                end
                threshold_lick=prctile(all_lick_traces(:),1)+((prctile(all_lick_traces(:),99)-prctile(all_lick_traces(:),1))/2);
                
                %Initialize variables
                
                odor_traces=[];
                od_ii=0;
                
                fvii=0;
                fvii_lick=0;
                no_fv_traces=0;
                fv_traces=[];
                
                sp_odor_response=[];
                splus_traces=[];
                spii=0;
                spii_lick=0;
                smii=0;
                smii_lick=0;
                no_sp_traces=0;
                no_sm_traces=0;
                splus_traces=[];
                sminus_traces=[];
                sm_odor_response=[];
                which_trace_sp=[];
                which_trace_sm=[];
                splus_lick_traces=[];
                sminus_lick_traces=[];
                
                
                %lda input
                lda_input_timecourse=[];
                lda_event=[];
                
                
                %Diefferent epochs
                Hitii_lick=0;
                no_Hit_traces=0;
                Hit_traces=[];
                which_trace_Hit=[];
                which_trial_Hit=[];
                which_Hitii_lick=[];
                Hit_lick_times=[];
                Hit_no_lick_times=0;
                no_Hit_trials=0;
                
                Missii_lick=0;
                no_Miss_traces=0;
                Miss_traces=[];
                which_trace_Miss=[];
                which_trial_Miss=[];
                which_Missii_lick=[];
                Miss_lick_times=[];
                Miss_no_lick_times=0;
                
                FAii_lick=0;
                no_FA_traces=0;
                FA_traces=[];
                which_trace_FA=[];
                which_trial_FA=[];
                which_FAii_lick=[];
                FA_lick_times=[];
                FA_no_lick_times=0;
                
                CRii_lick=0;
                no_CR_traces=0;
                CR_traces=[];
                which_trace_CR=[];
                which_trial_CR=[];
                which_CRii_lick=[];
                CR_lick_times=[];
                CR_no_lick_times=0;
                
                no_odor_trials=0;
                no_spm_odor_trials=0;
                epoch_per_trial=[];
                epoch_time=[];
                opto_on_per_trial=[];
                
                %Find Hits, CRs, etc
                for epoch=1:handles.dropcData.epochIndex
                    if ((handles.dropcData.epochTime(epoch)-dt_before)>0)
                        %This event is processed for Ca
                        
                        %Final valve epoch
                        if (handles.dropcData.epochEvent(epoch)==1)
                            
                            
                            
                            %Get the licks
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after);
                            fvii_lick=fvii_lick+1;
                            fv_lick_traces(fvii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                        end
                        handles_out.no_fv_trials=fvii_lick;
                        
                        %Now do S+ and S-
                        if (handles.dropcData.epochEvent(epoch)==2)
                            
                            
                            if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                                
                                %S plus
                                
                                rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                                spii_lick=spii_lick+1;
                                splus_lick_traces(spii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                                dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                                
                                handles_out.no_sp_trials=spii_lick;
                            else
                                
                                %S minus
                                
                                rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                                smii_lick=smii_lick+1;
                                sminus_lick_traces(smii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                                dt_rhd_trace=min([dt_rhd_trace (1/acq_rate)*sum(rhd_mask)]);
                                
                                handles_out.no_sm_trials=smii_lick;
                            end
                            %end
                        end
                        
                        %Find Hit, CR, FA and Miss
                        
                        %Hit
                        if (handles.dropcData.epochEvent(epoch)==6)
                            no_odor_trials=no_odor_trials+1;
                            epoch_per_trial(no_odor_trials)=6;
                            per_trial_odor_on_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S+';
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            %Is this an opto on trial
                            if sum((OpOn_times_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(OpOn_times_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset))>0
                                opto_on_per_trial(no_odor_trials)=1;
                            else
                                opto_on_per_trial(no_odor_trials)=0;
                            end
                            Hitii_lick=Hitii_lick+1;
                            which_trial_Hit(Hitii_lick)=no_odor_trials;
                            Hit_lick_traces(Hitii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            
                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            Hit_lick_times(Hitii_lick,1:length(these_lick_times))=these_lick_times;
                            Hit_no_lick_times(Hitii_lick)=length(these_lick_times);
                            
                            
                            handles_out.no_Hit_trials=Hitii_lick;
                            handles_out.Hit_trial_no(no_odor_trials)=Hitii_lick;
                        end
                        
                        %Miss
                        if (handles.dropcData.epochEvent(epoch)==7)
                            no_odor_trials=no_odor_trials+1;
                            epoch_per_trial(no_odor_trials)=7;
                            per_trial_odor_on_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S+';
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            %Is this an opto on trial
                            if sum((OpOn_times_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(OpOn_times_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset))>0
                                opto_on_per_trial(no_odor_trials)=1;
                            else
                                opto_on_per_trial(no_odor_trials)=0;
                            end
                            Missii_lick=Missii_lick+1;
                            which_trial_Miss(Missii_lick)=no_odor_trials;
                            Miss_lick_traces(Missii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            
                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            Miss_lick_times(Missii_lick,1:length(these_lick_times))=these_lick_times;
                            Miss_no_lick_times(Missii_lick)=length(these_lick_times);
                            
                            handles_out.no_Miss_trials=Missii_lick;
                            handles_out.Miss_trial_no(no_odor_trials)=Missii_lick;
                        end
                        
                        %FA
                        if (handles.dropcData.epochEvent(epoch)==8)
                            no_odor_trials=no_odor_trials+1;
                            epoch_per_trial(no_odor_trials)=8;
                            per_trial_odor_on_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S-';
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            %Is this an opto on trial
                            if sum((OpOn_times_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(OpOn_times_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset))>0
                                opto_on_per_trial(no_odor_trials)=1;
                            else
                                opto_on_per_trial(no_odor_trials)=0;
                            end
                            FAii_lick=FAii_lick+1;
                            which_trial_FA(FAii_lick)=no_odor_trials;
                            FA_lick_traces(FAii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            
                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            FA_lick_times(FAii_lick,1:length(these_lick_times))=these_lick_times;
                            FA_no_lick_times(FAii_lick)=length(these_lick_times);
                            
                            handles_out.no_FA_trials=FAii_lick;
                            handles_out.FA_trial_no(no_odor_trials)=FAii_lick;
                        end
                        
                        %CR
                        if (handles.dropcData.epochEvent(epoch)==9)
                            no_odor_trials=no_odor_trials+1;
                            epoch_per_trial(no_odor_trials)=9;
                            per_trial_odor_on_time(no_odor_trials)=handles.dropcData.epochTime(epoch);
                            lda_event{no_odor_trials}='S-';
                            rhd_mask=(time_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                &(time_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset);
                            %Is this an opto on trial
                            if sum((OpOn_times_rhd>=handles.dropcData.epochTime(epoch)-dt_before+dt_odor_onset)...
                                    &(OpOn_times_rhd<=handles.dropcData.epochTime(epoch)+dt_after+dt_odor_onset))>0
                                opto_on_per_trial(no_odor_trials)=1;
                            else
                                opto_on_per_trial(no_odor_trials)=0;
                            end
                            CRii_lick=CRii_lick+1;
                            which_trial_CR(CRii_lick)=no_odor_trials;
                            CR_lick_traces(CRii_lick,1:sum(rhd_mask))=adc_in(rhd_mask);
                            
                            these_lick_times=[];
                            these_lick_times=drgGetLicksCaImAn(acq_rate,threshold_lick,adc_in(rhd_mask))-dt_before;
                            CR_lick_times(CRii_lick,1:length(these_lick_times))=these_lick_times;
                            CR_no_lick_times(CRii_lick)=length(these_lick_times);
                            
                            handles_out.no_CR_trials=CRii_lick;
                            handles_out.CR_trial_no(no_odor_trials)=CRii_lick;
                            
                            
                        end
                        
                    end
                end
        end
        
        
        
        %Calculate lick frequency
        dt_lick=0.3;
        lick_t_start=0;
        lick_t_end=3;
        
        %Hit
        Hitlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        Hitlick_per_trial_timecourse=zeros(Hitii_lick,ceil((dt_after+dt_before)/dt_lick));
        
        for ii=1:Hitii_lick
            for ii_lick=1:Hit_no_lick_times(ii)
                Hitlick_freq( ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Hitlick_freq( ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
                Hitlick_per_trial_timecourse(ii,ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Hitlick_per_trial_timecourse(ii,ceil((Hit_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
            these_times=zeros(1,Hit_no_lick_times(ii));
            these_times(1,1:Hit_no_lick_times(ii))=Hit_lick_times(ii,1:Hit_no_lick_times(ii));
            Hitlick_per_trial(ii)=sum((these_times>=lick_t_start)&(these_times<=lick_t_end))/(lick_t_end-lick_t_start);
        end
        
        Hitlick_freq=(Hitlick_freq/(Hitii_lick*dt_lick));
        Hitlick_per_trial_timecourse=Hitlick_per_trial_timecourse/dt_lick;
        
        %Miss
        Misslick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        
        
        for ii=1:Missii_lick
            for ii_lick=1:Miss_no_lick_times(ii)
                Misslick_freq( ceil((Miss_lick_times(ii,ii_lick)+dt_before)/dt_lick))=Misslick_freq( ceil((Miss_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end
        
        Misslick_freq=(Misslick_freq/(Missii_lick*dt_lick));
        
        
        %CR
        CRlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        
        
        for ii=1:CRii_lick
            for ii_lick=1:CR_no_lick_times(ii)
                CRlick_freq( ceil((CR_lick_times(ii,ii_lick)+dt_before)/dt_lick))=CRlick_freq( ceil((CR_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end
        
        CRlick_freq=(CRlick_freq/(CRii_lick*dt_lick));
        
        %FA
        FAlick_freq=zeros(1,ceil((dt_after+dt_before)/dt_lick));
        
        
        for ii=1:FAii_lick
            for ii_lick=1:FA_no_lick_times(ii)
                FAlick_freq( ceil((FA_lick_times(ii,ii_lick)+dt_before)/dt_lick))=FAlick_freq( ceil((FA_lick_times(ii,ii_lick)+dt_before)/dt_lick))+1;
            end
        end
        
        FAlick_freq=(FAlick_freq/(FAii_lick*dt_lick));
        
        time_licks=([1:length(Hitlick_freq)]*dt_lick)-dt_before+dt_lick/2;
        
        
        
        for ii=1:smii_lick
            dsminus_lick_traces(ii,:)=decimate(sminus_lick_traces(ii,:)',20)';
        end
        pct99=max(dsminus_lick_traces(:));
        pct1=prctile(dsminus_lick_traces(:),1);
        dsminus_lick_traces=(dsminus_lick_traces/(pct99-pct1))+4;
        CIsmlick = bootci(1000, @mean, dsminus_lick_traces);
        meansmlick=mean(dsminus_lick_traces,1);
        CIsmlick(1,:)=meansmlick-CIsmlick(1,:);
        CIsmlick(2,:)=CIsmlick(2,:)-meansmlick;
        
        for ii=1:spii_lick
            dsplus_lick_traces(ii,:)=decimate(splus_lick_traces(ii,:)',20)';
        end
        pct99=prctile(dsplus_lick_traces(:),99);
        pct1=prctile(dsplus_lick_traces(:),1);
        dsplus_lick_traces=(dsplus_lick_traces/(pct99-pct1))+5.5;
        CIsplick = bootci(1000, @mean, dsplus_lick_traces);
        meansplick=mean(dsplus_lick_traces,1);
        CIsplick(1,:)=meansplick-CIsplick(1,:);
        CIsplick(2,:)=CIsplick(2,:)-meansplick;
        
        
        %Decimate the licks
        
        %CR lick traces
        dCR_lick_traces=[];
        for ii=1:CRii_lick
            dCR_lick_traces(ii,:)=decimate(CR_lick_traces(ii,:)',20)';
        end
        
        %FA lick traces
        dFA_lick_traces=[];
        for ii=1:FAii_lick
            dFA_lick_traces(ii,:)=decimate(FA_lick_traces(ii,:)',20)';
        end
        
        %Miss lick traces
        dMiss_lick_traces=[];
        for ii=1:Missii_lick
            dMiss_lick_traces(ii,:)=decimate(Miss_lick_traces(ii,:)',20)';
        end
        
        %Hit lick traces
        dHit_lick_traces=[];
        for ii=1:Hitii_lick
            dHit_lick_traces(ii,:)=decimate(Hit_lick_traces(ii,:)',20)';
        end
        
        
        %Save the calculated data
        rhdname=handles_choice.rhdFileName{fileNo};
        save_name=[rhdname(1:end-4) '_batch_licks.mat']
        
        
        save(save_name, 'Hitii_lick','no_Hit_traces','Hit_traces','which_trace_Hit','which_trial_Hit',...
            'Missii_lick','no_Miss_traces','Miss_traces','which_trace_Miss','which_trial_Miss',...
            'FAii_lick','no_FA_traces','FA_traces','which_trace_FA','which_trial_FA',...
            'CRii_lick','no_CR_traces','CR_traces','which_trace_CR','which_trial_CR',...
            'no_odor_trials','epoch_per_trial','epoch_time',...
            'dt_before','delta_odor','delta_odor_on_reinf_on','delta_reinf',...
            'dt_after','dt_odor_onset','splus_traces','sminus_traces','sp_odor_response',...
            'sm_odor_response','opto_on_per_trial',...
            'which_CRii_lick','which_FAii_lick','which_Hitii_lick','which_Missii_lick',...
            'CR_lick_times','FA_lick_times','Hit_lick_times','Miss_lick_times',...
            'CR_no_lick_times','FA_no_lick_times','Hit_no_lick_times','Miss_no_lick_times',...
            'time_licks','FAlick_freq','CRlick_freq','Hitlick_freq','Misslick_freq',...
            'handles','odor_traces','all_lick_traces','acq_rate',...
            'y_shift','time_rhd','adc_in','handles_out',...
            'lda_input_timecourse','lda_event','dHit_lick_traces'...
            ,'dCR_lick_traces','dMiss_lick_traces','dFA_lick_traces','per_trial_odor_on_time')
        
        
        
        %Plot the licks
        figure(2)
        hold on
        
        for ii=1:allii_lick
            dall_lick_traces(ii,:)=decimate(all_lick_traces(ii,:)',20)';
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
            dCR_lick_traces(ii,:)=decimate(CR_lick_traces(ii,:)',20)';
            plot(time_licks,dCR_lick_traces(ii,:)+y_shift,'-b')
            mean_CR_licks(1,:)= mean_CR_licks(1,:) + ((dCR_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end
        
        %Plot FA lick traces
        for ii=1:FAii_lick
            dFA_lick_traces(ii,:)=decimate(FA_lick_traces(ii,:)',20)';
            plot(time_licks,dFA_lick_traces(ii,:)+y_shift,'-m')
            mean_FA_licks(1,:)= mean_FA_licks(1,:) + ((dFA_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end
        
        %Plot Miss lick traces
        for ii=1:Missii_lick
            dMiss_lick_traces(ii,:)=decimate(Miss_lick_traces(ii,:)',20)';
            plot(time_licks,dMiss_lick_traces(ii,:)+y_shift,'-c')
            mean_Miss_licks(1,:)= mean_Miss_licks(1,:) + ((dMiss_lick_traces(ii,:)-per1)/(per99-per1));
            y_shift=y_shift+1.2*(per99-per1);
        end
        
        %PLot Hit lick traces
        for ii=1:Hitii_lick
            dHit_lick_traces(ii,:)=decimate(Hit_lick_traces(ii,:)',20)';
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
        
        
        %Plot lick frequency
        figure(3)
        hold on
        time_licks=([1:length(Hitlick_freq)]*dt_lick)-dt_before;
        plot(time_licks,Hitlick_freq,'-r')
        plot(time_licks,Misslick_freq,'-c')
        plot(time_licks,CRlick_freq,'-b')
        plot(time_licks,FAlick_freq,'-m')
        
        title('Lick frequency, r=Hit, c=Miss, b=CR, m=FA')
        xlabel('Time (sec)')
        ylabel('Lick frequency (Hz)')
        
        
        %Get lick p values
        try
            dt_lick_pval=0.1;
            no_pvals=0;
            p_val_Hit_CR=[];
            for ii=1:dt_lick_pval*acq_rate:sum(rhd_mask)-dt_lick_pval*acq_rate
                %Hit vs CR
                this_CR=zeros(CRii_lick,1);
                for jj=1:CRii_lick
                    if sum(CR_lick_traces(jj,ii:ii+dt_lick_pval*acq_rate)>threshold_lick)>=1
                        this_CR(jj,1)=1;
                    else
                        this_CR(jj,1)=0;
                    end
                end
                
                this_Hit=zeros(Hitii_lick,1);
                for jj=1:Hitii_lick
                    if sum(Hit_lick_traces(jj,ii:ii+dt_lick_pval*acq_rate)>threshold_lick)>=1
                        this_Hit(jj,1)=1;
                    else
                        this_Hit(jj,1)=0;
                    end
                end
                
                no_pvals=no_pvals+1;
                
                if (~isempty(this_CR))&(~isempty(this_CR))
                    p_val_Hit_CR(no_pvals)=ranksum(this_CR,this_Hit);
                else
                    p_val_Hit_CR(no_pvals)=1;
                end
                
                %ranksum gives NaN if the values are all the same
                if isnan(p_val_Hit_CR(no_pvals))
                    p_val_Hit_CR(no_pvals)=1;
                end
                
            end
            
            time_p_lick=[-dt_before+(dt_lick_pval/2):dt_lick_pval:dt_after-(dt_lick_pval)];
            
            
            figure(4)
            plot(time_p_lick,log10(p_val_Hit_CR))
            hold on
            plot([time_p_lick(1) time_p_lick(end)],[log10(0.05) log10(0.05)])
            title('log(p value) for the difference in licks Hit vs CR')
            xlabel('Time (sec)')
            ylabel('log10(p value)')
            
        catch
        end
        
    catch
        fprintf(1, ['\ndrgCaImAn_batch_dropc_no_microscope was not run for ' handles_choice.spmFileName{fileNo} '\n\n']);
    end
    
end

pfft=1;


















