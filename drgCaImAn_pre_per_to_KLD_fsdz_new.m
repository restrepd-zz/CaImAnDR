function handles_out2=drgCaImAn_pre_per_to_KLD_fsdz_new(pre_perBatchPathName, pre_perFileName, p_threshold, figNo, show_figures,no_sp_sm_trials_to_use,first_sp_sm_trial_no)
%
% this function performs a Kullback–Leibler divergence calculation on the
% data
%

warning('off')


min_no_trials=10;


handles_out2.pre_perBatchPathName=pre_perBatchPathName;
handles_out2.pre_perFileName=pre_perFileName;
% handles_out2.p_threshold=p_threshold;
% handles_out2.MLalgo=MLalgo;
handles_out2.decoding_processed=1;

% if ~iscell(pre_perFileName)
load([pre_perBatchPathName pre_perFileName])

handles_out2.trialNo_sp=zeros(1,handles_out.no_sp_trials);
ii_sp=0;
handles_out2.trialNo_sm=zeros(1,handles_out.no_sm_trials);
ii_sm=0;

for ii_trials=1:no_odor_trials
    if (epoch_per_trial(ii_trials)==6)||(epoch_per_trial(ii_trials)==7)
        %This is an S+ trial
        ii_sp=ii_sp+1;
        handles_out2.trialNo_sp(ii_sp)=ii_trials;
    else
        %This is an S- trial
        ii_sm=ii_sm+1;
        handles_out2.trialNo_sm(ii_sm)=ii_trials;
    end
end

handles_out2.no_odor_trials=no_odor_trials;

if (handles_out.no_sp_trials>=min_no_trials)&(handles_out.no_sm_trials>=min_no_trials)
    
    szSp=size(splus_traces);
    szSm=size(sminus_traces);
    %time_to_event=([1:szSm(2)]*dt-dt_before);
    time_to_eventSm=([1:szSm(2)]*dt-dt_before);
    time_to_eventSp=([1:szSp(2)]*dt-dt_before);
    handles_out2.time_to_eventSm=time_to_eventSm;
    handles_out2.time_to_eventSp=time_to_eventSp;
    
    if show_figures==1
        
        %S+, S-, all snips
        CIsm = bootci(1000, @mean, sminus_traces);
        meansm=mean(sminus_traces,1);
        CIsm(1,:)=meansm-CIsm(1,:);
        CIsm(2,:)=CIsm(2,:)-meansm;
        
        CIsp = bootci(1000, @mean, splus_traces);
        meansp=mean(splus_traces,1);
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
        
        
        
        pct1=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],1);
        pct99=prctile([mean(sminus_traces,1)'; mean(splus_traces(:,1:szSp(2)),1)'],99);
        
        
        
        [hlsm, hpsm] = boundedline(time_to_eventSm',mean(sminus_traces,1)', CIsm', 'b');
        [hlsp, hpsp] = boundedline(time_to_eventSp',mean(splus_traces,1)', CIsp', 'r');
        
        %Odor on markers
        plot([0 0],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        odorhl=plot([0 mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct1-0.1*(pct99-pct1)],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[pct1-0.1*(pct99-pct1) pct99+0.1*(pct99-pct1)],'-r')
        
        
        
        title("Ca changes aligned to odor onset, all trials")
        legend([hlsp hlsm odorhl reinfhl],'S+','S-','Odor','Reinforcement')
        xlabel('Time (sec)')
        ylabel('dF/F')
        ylim([pct1-0.2*(pct99-pct1) pct99+0.2*(pct99-pct1)])
        xlim([-10 19.8])
    end
     
    %Initialize variables
    time_bins=length(handles_out.time_to_eventSp);
    time=time_to_eventSp;
    trNo=0;
    tr_trNo=0;
    
    %First and last sp trial numbers
    if (first_sp_sm_trial_no<handles_out.no_sp_trials)&(first_sp_sm_trial_no+no_sp_sm_trials_to_use-1<=handles_out.no_sp_trials)
        first_sp_trial=first_sp_sm_trial_no;
        last_sp_trial=first_sp_sm_trial_no+no_sp_sm_trials_to_use-1;
    else
        if first_sp_sm_trial_no+no_sp_sm_trials_to_use-1>handles_out.no_sp_trials
            if handles_out.no_sp_trials-no_sp_sm_trials_to_use+1>0
                first_sp_trial=handles_out.no_sp_trials-no_sp_sm_trials_to_use+1;
            else
                first_sp_trial=1;
            end
            last_sp_trial=handles_out.no_sp_trials;
        end
    end
    
    %First and last sm trial numbers
    if (first_sp_sm_trial_no<handles_out.no_sm_trials)&(first_sp_sm_trial_no+no_sp_sm_trials_to_use-1<handles_out.no_sm_trials)
        first_sm_trial=first_sp_sm_trial_no;
        last_sm_trial=first_sp_sm_trial_no+no_sp_sm_trials_to_use-1;
    else
        if first_sp_sm_trial_no+no_sp_sm_trials_to_use-1>handles_out.no_sm_trials
            if handles_out.no_sm_trials-no_sp_sm_trials_to_use+1>0
                first_sm_trial=handles_out.no_sm_trials-no_sp_sm_trials_to_use+1;
            else
                first_sm_trial=1;
            end
            last_sm_trial=handles_out.no_sm_trials;
        end
    end
    
    
    
    no_sp_trials_to_use=last_sp_trial-first_sp_trial+1;
    no_sm_trials_to_use=last_sm_trial-first_sm_trial+1;
    
    
    %Note: Training is done with no_sp_sm_trials_to_use and outcome is
    %calculated for all trials
    all_trials=handles_out.no_sm_trials+handles_out.no_sp_trials;
    decisions=zeros(1,all_trials);
    neural_recordings=zeros(all_trials,handles_out.no_components,time_bins);
    ii_all=zeros(1,no_sp_trials_to_use+no_sm_trials_to_use);
    training_decisions=zeros(1,no_sp_trials_to_use+no_sm_trials_to_use);
    training_neural_recordings=zeros(no_sp_trials_to_use+no_sm_trials_to_use,handles_out.no_components,time_bins);
    
    %Save S+
    %All trials, and training trials if adequate
    ii=1;
    for trialNo=1:handles_out.no_sp_trials
        trNo=trNo+1;
        if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
            tr_trNo=tr_trNo+1;
        end
        for traceNo=1:handles_out.no_components
            neural_recordings(trNo,traceNo,:)=splus_traces(ii,:);
            if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
                training_neural_recordings(tr_trNo,traceNo,:)=splus_traces(ii,:);
            end
            ii=ii+1;
        end
        decisions(trNo)=1;
        if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
            training_decisions(tr_trNo)=1;
            ii_all(tr_trNo)=trialNo;
        end
    end
    
    
    %Save S-
    %All trials
    ii=1;
    for trialNo=1:handles_out.no_sm_trials
        trNo=trNo+1;
        if (trialNo>=first_sm_trial)&(trialNo<=last_sm_trial)
            tr_trNo=tr_trNo+1;
        end
        for traceNo=1:handles_out.no_components
            neural_recordings(trialNo+handles_out.no_sp_trials,traceNo,:)=sminus_traces(ii,:);
            if (trialNo>=first_sm_trial)&(trialNo<=last_sm_trial)
                training_neural_recordings(tr_trNo,traceNo,:)=sminus_traces(ii,:);
            end
            ii=ii+1;
        end
        decisions(trNo)=0;
        if (trialNo>=first_sm_trial)&(trialNo<=last_sm_trial)
            training_decisions(tr_trNo)=0;
            ii_all(tr_trNo)=trialNo;
        end
    end
    
    %For the training set calculate p values for the difference per neuron
    no_neurons=handles_out.no_components;
    tr_no_trials=tr_trNo;
    no_timepoints=length(time_to_eventLDA);
    
    p=[];
    for ii_neuron=1:handles_out.no_components
        these_recordings_sp=zeros(no_sp_trials_to_use,no_timepoints);
        these_recordings_sp(:,:)=training_neural_recordings(logical(training_decisions),ii_neuron,:);
        
        these_recordings_sm=zeros(no_sm_trials_to_use,no_timepoints);
        these_recordings_sm(:,:)=training_neural_recordings(~logical(training_decisions),ii_neuron,:);
        
        mean_sp_odor=zeros(no_sp_trials_to_use,1);
        mean_sp_odor(:,1)=mean(these_recordings_sp(:,time_to_eventLDA>0),2);
        
        mean_sm_odor=zeros(no_sm_trials_to_use,1);
        mean_sm_odor(:,1)=mean(these_recordings_sm(:,time_to_eventLDA>0),2);
        
        [h,p(ii_neuron)]=ttest2(mean_sp_odor,mean_sm_odor);
    end
    
    handles_out2.p=p;
    
else
    handles_out2.decoding_processed=0;
end

if handles_out2.decoding_processed==1
    
    %Threshold out the spikes
    std_traces=std(training_neural_recordings(:));
    thr=5*std_traces;
    
    %These are the placeholders for the spike counts
    Sp_spikes_timecourse=zeros(no_timepoints,handles_out.no_components);
    Sm_spikes_timecourse=zeros(no_timepoints,handles_out.no_components);
     
    training_neural_spikes=zeros(no_sp_trials_to_use+no_sm_trials_to_use,handles_out.no_components,time_bins);

    for ii_comps=1:handles_out.no_components
        
        these_recordings=zeros(length(training_decisions),no_timepoints);
        these_recordings(:,:)=training_neural_recordings(:,ii_comps,:);
        baseline_dF=prctile(these_recordings(:),5);
        
        %S+
        for ii_t=1:no_timepoints
            these_dFFs=zeros(1,no_sp_trials_to_use);
            these_dFFs(1,:)=training_neural_recordings(logical(training_decisions),ii_comps,ii_t)-baseline_dF;
            Sp_spikes_timecourse(ii_t,ii_comps)=sum(these_dFFs>thr);
        end
        
        %S-
        for ii_t=1:no_timepoints
            these_dFFs=zeros(1,no_sm_trials_to_use);
            these_dFFs(1,:)=training_neural_recordings(~logical(training_decisions),ii_comps,ii_t)-baseline_dF;
            Sm_spikes_timecourse(ii_t,ii_comps)=sum(these_dFFs>thr);
        end
        
        for ii_t=1:no_timepoints
            for trNo=1:no_sp_trials_to_use+no_sm_trials_to_use
                if training_neural_recordings(trNo,ii_comps,ii_t)-baseline_dF>thr
                    training_neural_spikes(trNo,ii_comps,ii_t)=1;
                end
            end
        end
        
    end
    
    max_spikes=max([max(Sp_spikes_timecourse(:)) max(Sm_spikes_timecourse(:))]);
    min_spikes=min([min(Sp_spikes_timecourse(:)) min(Sm_spikes_timecourse(:))]);
    
    
    
    Sp_spikes=zeros(1,length(time_to_eventLDA));
    Sp_spikes(1,:)=sum(Sp_spikes_timecourse,2);
    
    
    Sm_spikes=zeros(1,length(time_to_eventLDA));
    Sm_spikes(1,:)=sum(Sm_spikes_timecourse,2);
    
    handles_out2.Sm_spikes=Sm_spikes;
    handles_out2.Sp_spikes=Sp_spikes;
    handles_out2.Sp_spikes_timecourse=Sp_spikes_timecourse;
    handles_out2.Sm_spikes_timecourse=Sm_spikes_timecourse;
    handles_out2.thr=thr;
    
    pffft=1;
    
    
    
    
    
    %Do KL divergence and ROC using dFF
    
    %These are used to calculate the z values
    Nall=length(training_decisions);
    
    
    %Calculate z from dFF
    mean_per_neuron=zeros(1,sum(p<=p_threshold));
    STD_per_neuron=zeros(1,sum(p<=p_threshold));
    z_training_neural_recordings=zeros(Nall,sum(p<=p_threshold),no_timepoints);
    jj_neurons=0;
    for ii_neurons=1:no_neurons
        if p(ii_neurons)<=p_threshold
            jj_neurons=jj_neurons+1;
            all_pre=[];
            for ii_trials=1:tr_no_trials
                these_pre=zeros(1,sum((time_to_eventLDA<0)&(time_to_eventLDA>=-5)));
                these_pre(1,:)=training_neural_recordings(ii_trials,ii_neurons,(time_to_eventLDA<0)&(time_to_eventLDA>=-5));
                all_pre=[all_pre these_pre];
            end
            STD_per_neuron(jj_neurons)=std(all_pre);
            for ii_trials=1:tr_no_trials
                these_dFFs=[];
                these_dFFs=training_neural_recordings(ii_trials,ii_neurons,(time_to_eventLDA<0)&(time_to_eventLDA>=-5));
                z_training_neural_recordings(ii_trials,jj_neurons,:)=(training_neural_recordings(ii_trials,ii_neurons,:)-mean(these_dFFs))/STD_per_neuron(jj_neurons);
            end
        end
    end
    
    
    handles_out2.time_to_eventLDA=time_to_eventLDA;
    
    %Now do the KL Divergence
    KLd=zeros(1,no_timepoints);
    auROC=zeros(1,no_timepoints);
    for time_point=1:no_timepoints
        
        these_z=zeros(Nall,sum(p<=p_threshold));
        these_z(:,:)=z_training_neural_recordings(:,:,time_point);
        
        odor1_mean=zeros(1,sum(p<=p_threshold));
        odor2_mean=zeros(1,sum(p<=p_threshold));
        odor1_mean(1,:)=mean(these_z(training_decisions==0,:),1);
        odor2_mean(1,:)=mean(these_z(training_decisions==1,:),1);
        
        %Use pdist to find all distances for odor 1
        iiod1=0;
        iiod2=0;
        distances_od1=zeros(1,sum(training_decisions==0));
        distances_od2=zeros(1,sum(training_decisions==1));
        for ii=1:Nall
            odor_queary=zeros(1,sum(p<=p_threshold));
            odor_queary(1,:)=these_z(ii,:);
            all_points=[odor1_mean; odor2_mean; odor_queary];
            all_distances=pdist(all_points);
            d12=all_distances(1);
            dq1=all_distances(2);
            dq2=all_distances(3);
            if training_decisions(ii)==0
                iiod1=iiod1+1;
                distances_od1(1,iiod1)=(d12^2 + dq1^2 -dq2^2)/(2*d12);
            else
                iiod2=iiod2+1;
                distances_od2(1,iiod2)=(d12^2 + dq1^2 -dq2^2)/(2*d12);
            end
        end
        
        
        %KL divergence
        num_bins=20;
        max_d=max([max(distances_od1) max(distances_od2)]);
        min_d=min([min(distances_od1) min(distances_od2)]);
        X=[min_d:(max_d-min_d)/(num_bins-1):max_d];
        
        p1=zeros(1,length(X));
        p2=zeros(1,length(X));
        
        for ii=1:length(distances_od1)
            [min_d min_ii]=min(abs(distances_od1(ii)-X));
            p1(min_ii)=p1(min_ii)+1;
        end
        p1=p1/sum(p1);
        p1(p1==0)=eps;
        
        for ii=1:length(distances_od2)
            [min_d min_ii]=min(abs(distances_od2(ii)-X));
            p2(min_ii)=p2(min_ii)+1;
        end
        p2=p2/sum(p2);
        p2(p2==0)=eps;
        
        KLd(time_point) = kldiv(X,p1,p2);
        
        %ROC
        roc_data=[];
        roc_data(1:length(distances_od1),1)=distances_od1;
        roc_data(1:length(distances_od1),2)=zeros(length(distances_od1),1);
        roc_data(length(distances_od1)+1:Nall,1)=distances_od2;
        roc_data(length(distances_od1)+1:Nall,2)=ones(length(distances_od2),1);
        
        roc=roc_calc(roc_data,0,0.05,0);
        auROC(time_point)=roc.AUC-0.5;
    end
    
    handles_out2.KLdivergence=KLd;
    handles_out2.auROC=auROC;
    
    
    
    
    if show_figures==1
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .2 .25 .7])
        
        
        time=repmat(time_to_eventLDA,handles_out.no_components,1);
        comps_index=1:handles_out.no_components;
        comps=repmat(comps_index,length(time_to_eventLDA),1);
        
        
        drg_pcolor(time',comps,Sp_spikes_timecourse)
        
        
        colormap fire
        % shading interp
        caxis([min_spikes max_spikes]);
        xlabel('Time (sec)')
        ylabel('ROI');
        title(['delta F/F spikes vs time for all ROIs for S+, p thr=' num2str(p_threshold)])
        
        
        
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        set(hFig, 'units','normalized','position',[.05 .2 .25 .7])
        
        
        time=repmat(time_to_eventLDA,handles_out.no_components,1);
        comps_index=1:handles_out.no_components;
        comps=repmat(comps_index,length(time_to_eventLDA),1);
        
        
        drg_pcolor(time',comps,Sm_spikes_timecourse)
        
        
        colormap fire
        % shading interp
        % caxis([mindFF maxdFF]);
        caxis([min_spikes max_spikes]);
        xlabel('Time (sec)')
        ylabel('ROI');
        title(['delta F/F spikes vs time for all ROIs for S-, p thr=' num2str(p_threshold)])
        
        
        %Plot the average Splus and Sminus
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        hold on
        
        plot(time_to_eventSm',Sm_spikes,'-b')
        plot(time_to_eventSp',Sp_spikes,'-r')
        
        
        title(['Thresholded dF/F , p thr=' num2str(p_threshold)])
        
        
        %First plot the S+ traces
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        hold on
        
        plot(splus_traces')
        
        
        title('All S+ traces')
        
        
        
        %First plot the S- traces
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        
        hold on
        
        plot(sminus_traces')
        
        
        title('All S- traces')
        
        
        
        %KL divergence
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        plot(time_to_eventLDA,KLd)
        title(['KL divergence, p thr=' num2str(p_threshold)])
        xlabel('Time(sec)')
        ylabel('KL divergence')
        
        %auROC
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        plot(time_to_eventLDA,auROC)
        title(['ROC, p thr=' num2str(p_threshold)])
        xlabel('Time(sec)')
        ylabel('ROC')
        
        
    end
    
end
pffft=1;

