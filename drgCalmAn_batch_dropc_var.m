%% drgCaImAn_batch_dropc.m
%
% Needs as an input the ouput file from drgCaImAn_dropc
%
close all
clear all

%Choices
do_warp=0;     %1=warped components from a reference file

%Other default variables
plot_raw=1; %Used plot_raw=1 for Ming's data
ref_win=[-5 -1.5];


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
 
cd(handles_choice.PathName)
 
reverse_fileNo=3;

for fileNo=handles_choice.first_file:handles_choice.no_files
    
    fprintf(1, ['\nProcessing  ' num2str(fileNo) '\n']);
    
    Yr=[];
    A_or=[];
    C_or=[];
    b2=[];
    f2=[];
    Cn=[];
    options=[];
    if do_warp==1
        [Yr,A_or,C_or,b2,f2,Cn,options] = drgCaImAn_get_warped_components(handles_choice.CaImAn_scriptFileName{fileNo},handles_choice.CaImAn_referenceFileName);
    else
        load(handles_choice.CaImAn_scriptFileName{fileNo})
    end
    
    fnameca=handles_choice.CaImAn_scriptFileName{fileNo};
    
    close all
    
    dt=1/options.fr;
    raw=[];
    inferred=[];
    try
        [raw,inferred]=drgGetCAtraces(Yr,A_or,C_or,b2,f2,Cn,options);
    catch
        pffft=1
    end
    
    % Should we use raw or inferred?
    if plot_raw==1
        traces=raw;
    else
        traces=inferred;
    end
    
    sz_traces=size(traces);
    no_traces=sz_traces(1);
    no_images=sz_traces(2);
    
    fprintf(1, ['\ndrgCaImAn_dropc run for ' fnameca '\n\n']);
     
    %Read the dropc file
    handles=[];
    load(handles_choice.spmFileName{fileNo})
     
    %Read the rhd file
    adc_in=[];
    digital_in=[];
    acq_rate=[];
    [adc_in,digital_in,acq_rate]=drg_read_Intan_RHD2000_file(handles_choice.rhdFileName{fileNo},3);
    
    
    % Plot the traces
    these_lines{1}='-b';
    these_lines{2}='-r';
    these_lines{3}='-m';
    these_lines{8}='-g';
    these_lines{5}='-y';
    these_lines{6}='-k';
    these_lines{7}='-c';
    these_lines{4}='-k';
    
    try
        close 1
    catch
    end
    
    hFig1 = figure(1);
    set(hFig1, 'units','normalized','position',[.05 .1 .85 .8])
    
    
    hold on
    
    % Determine the y spacing of the traces
    y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));
    
    %Plot the event lines
    odor_on_times=[];
    ootii=0;
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
                    if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                        plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                            '-r','LineWidth',1)
                    else
                        plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces+2)*y_shift],...
                            '-b','LineWidth',1)
                    end
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
    
    %Plot the licks recorded by the INTAN (adc_in)
    time_rhd=([1:length(digital_in)]/acq_rate)+delta_t_rhd;
    pct998=prctile(adc_in,99.8);
    pct1=prctile(adc_in,1);
    norm_fact=0.8*y_shift/(pct998-pct1);
    
    plot(time_rhd(time_rhd>0),adc_in(time_rhd>0)*norm_fact)
    
    %Plot the traces
    time=[1:no_images]*dt;
    for trNo=1:no_traces
        % for trNo=1:20
        plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
    end
    
    
    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    
    
    xlabel('time (s)')
    ylabel('deltaF/F')
    title(fnameca(1:end-4))
    
    out_traces=traces;
    %Exclude changes during odor on
    for noOdOn=1:ootii
        this_ii=floor(odor_on_times(noOdOn)/dt);
        from_ii=this_ii-floor(5/dt);
        if from_ii<1
            from_ii=1;
        end
        to_ii=this_ii+floor(15/dt);
        if to_ii>length(time)
            to_ii=length(time);
        end
        for trNo=1:no_traces
            out_traces(trNo,from_ii:to_ii)=out_traces(trNo,from_ii);
        end
    end
    
    %Calculate variance
    var_window=5;
    k=floor(var_window/dt);
    
    dFFvar_out=movvar(out_traces,k);
    dFFvar=movvar(traces,k);
    %plot(dFFvar)
    
    freq=0.1:0.1:2;
    for trNo=1:no_traces
        [S,f,t,P]=spectrogram(detrend(double(traces(trNo,:))),20,1,freq,1/dt);
        if trNo==1
            allP=ones(no_traces,length(f),length(t));
            allP(1,:,:)=P;
        else
            allP(trNo,:,:)=P;
        end
    end
    logP=ones(length(f),length(t));
    logP(:,:)=mean(10*log10(allP),1);
    
        for trNo=1:no_traces
        [S,f,t,P_out]=spectrogram(detrend(double(out_traces(trNo,:))),20,1,freq,1/dt);
        if trNo==1
            allP_out=ones(no_traces,length(f),length(t));
            allP_out(1,:,:)=P_out;
        else
            allP_out(trNo,:,:)=P_out;
        end
    end
    logP_out=ones(length(f),length(t));
    logP_out(:,:)=mean(10*log10(allP_out),1);
    
    if fileNo==handles_choice.first_file
        cum_dFFvar=mean(dFFvar,1);
        cum_dFFvar_out=mean(dFFvar_out,1);
        all_times=time;
        allodor_on_times=odor_on_times;
        all_logP=logP;
        all_logP_times=t;
        all_logP_out=logP_out;
        
    else
        this_dFFvar=mean(dFFvar,1);
        cum_dFFvar=[cum_dFFvar this_dFFvar];
        this_dFFvar_out=mean(dFFvar_out,1);
        cum_dFFvar_out=[cum_dFFvar_out this_dFFvar_out];
        last_time=all_times(end);
        if fileNo==reverse_fileNo
            reverse_time=last_time+dt;
        end
        all_times=[all_times time+last_time+dt];
        allodor_on_times=[allodor_on_times odor_on_times+last_time+dt];
        all_logP=[all_logP logP];
        all_logP_out=[all_logP_out logP_out];
        last_logP_time=all_logP_times(end);
         if fileNo==reverse_fileNo
            reverse_logP_time=last_logP_time+1;
        end
        all_logP_times=[all_logP_times t+last_logP_time+1];
    end
end

figure(2)
plot(all_times,cum_dFFvar,'-k')
hold on
plot([reverse_time reverse_time],[min(cum_dFFvar) max(cum_dFFvar)],'-r')
ylabel('Variance')
xlabel('Time (sec)')
title(['Variance for ' choiceFileName])

hFig3=figure(3)
set(hFig3, 'units','normalized','position',[.07 .1 .75 .3])
plot(all_times,cum_dFFvar_out,'-k')
hold on
plot([reverse_time reverse_time],[min(cum_dFFvar) max(cum_dFFvar)],'-r')
ylabel('Variance')
xlabel('Time (sec)')
title(['Variance excluding odor events for ' choiceFileName])



log_P_ref=zeros(length(f),length(all_logP_times));
reft=(all_logP_times<reverse_logP_time);
log_P_ref(:,:)=repmat(mean(all_logP(:,reft),2)',length(all_logP_times),1)';

logPmax=prctile(all_logP(:)-log_P_ref(:),99);
logPmin=prctile(all_logP(:)-log_P_ref(:),1);
if (logPmax-logPmin)>10
    logPmin=logPmax-10;
end

hFig4 = figure(4);
set(hFig4, 'units','normalized','position',[.07 .1 .75 .3])
drg_pcolor(repmat(all_logP_times,length(freq),1)',repmat(freq,length(all_logP_times),1),(all_logP-log_P_ref)')
hold on
plot([reverse_logP_time reverse_logP_time],[0.1 2],'-k','LineWidth',3)

colormap jet
shading interp
caxis([logPmin logPmax]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title('Power (dB) timecourse ')




log_P_ref_out=zeros(length(f),length(all_logP_times));
reft=(all_logP_times<reverse_logP_time);
log_P_ref_out(:,:)=repmat(mean(all_logP_out(:,reft),2)',length(all_logP_times),1)';

logPmax_out=prctile(all_logP_out(:)-log_P_ref_out(:),99);
logPmin_out=prctile(all_logP_out(:)-log_P_ref_out(:),1);
if (logPmax_out-logPmin_out)>10
    logPmin_out=logPmax_out-10;
end

hFig5 = figure(5);
set(hFig5, 'units','normalized','position',[.07 .1 .75 .3])
drg_pcolor(repmat(all_logP_times,length(freq),1)',repmat(freq,length(all_logP_times),1),(all_logP_out-log_P_ref_out)')
hold on
plot([reverse_logP_time reverse_logP_time],[0.1 2],'-k','LineWidth',3)

colormap jet
shading interp
caxis([logPmin_out logPmax_out]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)');
title('Power (dB) timecourse excluding odor events')