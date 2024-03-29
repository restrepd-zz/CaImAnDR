function handles_out2=drgCaImAn_pre_per_to_PCA_fsdz_new(pre_perBatchPathName, pre_perFileName, figNo, show_figures,no_sp_sm_trials_to_use,first_sp_sm_trial_no)
%
% reads the pre_per file and saves .mat files to process with
% Kording's lab
%

% close all
% clearvars -except pre_perBatchPathName pre_perFileName show_figures no_sp_sm_trials_to_use first_sp_sm_trial_no figNo
warning('off')

% simulation=0;

classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';

handles_out2.classifier_names=classifier_names;

% if nargin==0
%     [pre_perFileName,pre_perBatchPathName] = uigetfile({'*pre_per.mat'},'Select the pre_per.mat file');
%     fprintf(1, ['\ndrgCaImAn_pre_per_to_pydec run for ' pre_perFileName '\n\n']);
%
%     p_threshold=0.05;
%
%     %Choose the machine learning algorithm
%     % 1 = linear discriminant analysis
%     % 2 = SVM
%     % 3 = Naive Bayes Classifier
%     % 4 = neural network
%     % 5 = decision tree
%     MLalgo=1;
%     show_figures=1;
% end

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
        
        
        
        title("Ca changes aligned to odor onset")
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
            first_sp_trial=handles_out.no_sp_trials-no_sp_sm_trials_to_use+1;
            last_sp_trial=handles_out.no_sp_trials;
        end
    end
    
    %First and last sm trial numbers
    if (first_sp_sm_trial_no<handles_out.no_sm_trials)&(first_sp_sm_trial_no+no_sp_sm_trials_to_use-1<handles_out.no_sm_trials)
        first_sm_trial=first_sp_sm_trial_no;
        last_sm_trial=first_sp_sm_trial_no+no_sp_sm_trials_to_use-1;
    else
        if first_sp_sm_trial_no+no_sp_sm_trials_to_use-1>handles_out.no_sm_trials
            first_sm_trial=handles_out.no_sm_trials-no_sp_sm_trials_to_use+1;
            last_sm_trial=handles_out.no_sm_trials;
        end
    end
    
    %Note: Training is done with no_sp_sm_trials_to_use and outcome is
    %calculated for all trials
    all_trials=handles_out.no_sm_trials+handles_out.no_sp_trials;
    decisions=zeros(1,all_trials);
    neural_recordings=zeros(all_trials,handles_out.no_components,time_bins);
    ii_all=zeros(1,2*no_sp_sm_trials_to_use);
    training_decisions=zeros(1,2*no_sp_sm_trials_to_use);
    training_neural_recordings=zeros(2*no_sp_sm_trials_to_use,handles_out.no_components,time_bins);
    
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
        if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
            tr_trNo=tr_trNo+1;
        end
        for traceNo=1:handles_out.no_components
            neural_recordings(trialNo+handles_out.no_sp_trials,traceNo,:)=sminus_traces(ii,:);
            if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
                training_neural_recordings(tr_trNo,traceNo,:)=sminus_traces(ii,:);
            end
            ii=ii+1;
        end
        decisions(trNo)=0;
        if (trialNo>=first_sp_trial)&(trialNo<=last_sp_trial)
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
        these_recordings_sp=zeros(no_sp_sm_trials_to_use,no_timepoints);
        these_recordings_sp(:,:)=training_neural_recordings(logical(training_decisions),ii_neuron,:);
        
        these_recordings_sm=zeros(no_sp_sm_trials_to_use,no_timepoints);
        these_recordings_sm(:,:)=training_neural_recordings(~logical(training_decisions),ii_neuron,:);
        
        mean_sp_odor=zeros(no_sp_sm_trials_to_use,1);
        mean_sp_odor(:,1)=mean(these_recordings_sp(:,time_to_eventLDA>0),2);
        
        mean_sm_odor=zeros(no_sp_sm_trials_to_use,1);
        mean_sm_odor(:,1)=mean(these_recordings_sm(:,time_to_eventLDA>0),2);
        
        [h,p(ii_neuron)]=ttest2(mean_sp_odor,mean_sm_odor);
    end
    
    handles_out2.p=p;
    
else
    handles_out2.decoding_processed=0;
end
% else
%     sp_trials_per_file=[];
%     sm_trials_per_file=[];
%     len_time_to_eventSp=[];
%     len_time_to_eventSm=[];
%     no_components=[];
%     for ii_file=1:length(pre_perFileName)
%         load([pre_perBatchPathName{ii_file} pre_perFileName{ii_file}])
%         sp_trials_per_file(ii_file)=handles_out.no_sp_trials;
%         sm_trials_per_file(ii_file)=handles_out.no_sm_trials;
%         len_time_to_eventSp(ii_file)=length(handles_out.time_to_eventSp);
%         len_time_to_eventSm(ii_file)=length(handles_out.time_to_eventSm);
%         no_components(ii_file)=handles_out.no_components;
%     end
%
%
%
%     if (max(sp_trials_per_file)<min_no_trials)||(max(sm_trials_per_file)<min_no_trials)
%         handles_out2.decoding_processed=0;
%     else
%         no_sp_trials=max([min_no_trials min(sp_trials_per_file)]);
%         no_sm_trials=max([min_no_trials min(sm_trials_per_file)]);
%
%
%         no_timepoints=min([min(len_time_to_eventSp) min(len_time_to_eventSm)]);
%         time=handles_out.time_to_eventSp(1:no_timepoints);
%         neural_recordings=zeros(no_sp_trials+no_sm_trials,sum(no_components),no_timepoints);
%         decisions=zeros(1,no_sp_trials+no_sm_trials);
%
%         no_neurons=0;
%         for ii_file=1:length(pre_perFileName)
%
%             load([pre_perBatchPathName{ii_file} pre_perFileName{ii_file}])
%
%             %Save S+
%
%             trNo=1;
%             ii=1;
%             for trialNo=1:no_sp_trials
%                 for traceNo=1:handles_out.no_components
%                     neural_recordings(trialNo,traceNo+no_neurons,1:no_timepoints)=splus_traces(ii,1:no_timepoints);
%                     ii=ii+1;
%                 end
%                 decisions(trNo)=1;
%                 trNo=trNo+1;
%             end
%
%             %Save S-
%             ii=1;
%             for trialNo=1:no_sm_trials
%                 for traceNo=1:handles_out.no_components
%                     neural_recordings(trialNo+no_sp_trials,traceNo+no_neurons,1:no_timepoints)=sminus_traces(ii,1:no_timepoints);
%                     ii=ii+1;
%                 end
%                 decisions(trNo)=0;
%                 trNo=trNo+1;
%             end
%
%             no_neurons=no_neurons+handles_out.no_components;
%
%         end
%
%         %Calculate p values for the difference per neuron
%         p=[];
%         for ii_neuron=1:no_neurons
%             these_recordings_sp=zeros(no_sp_trials,no_timepoints);
%             these_recordings_sp=neural_recordings(logical(decisions),ii_neuron,:);
%
%             these_recordings_sm=zeros(no_sm_trials,no_timepoints);
%             these_recordings_sm=neural_recordings(~logical(decisions),ii_neuron,:);
%
%             mean_sp_odor=zeros(no_sp_trials,1);
%             mean_sp_odor(:,1)=mean(these_recordings_sp(:,time>0),2);
%
%             mean_sm_odor=zeros(no_sm_trials,1);
%             mean_sm_odor(:,1)=mean(these_recordings_sm(:,time>0),2);
%
%             [h,p(ii_neuron)]=ttest2(mean_sp_odor,mean_sm_odor);
%         end
%
%         handles_out2.p=p;
%         handles_out2.time=time;
%         no_trials=trNo-1;
%     end
% end


if handles_out2.decoding_processed==1
    
    %Now do the PCA
    for time_point=1:no_timepoints
        these_neural_recs=zeros(tr_no_trials,no_neurons);
        these_neural_recs(:,:)=training_neural_recordings(:,:,time_point);
        [coeff,score,latent] = pca(these_neural_recs);
        if time_point==1
            pca_scores=zeros(no_timepoints,tr_no_trials,size(score,2));
        end
        pca_scores(time_point,:,:)=score;
    end
    
    show_figures=1;
    if show_figures==1
        
        %Plot trajectory in 3D space
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        for trNo=1:tr_no_trials
            x=zeros(1,no_timepoints);
            x(1,:)=pca_scores(:,trNo,1);
            y=zeros(1,no_timepoints);
            y(1,:)=pca_scores(:,trNo,2);
            z=zeros(1,no_timepoints);
            z(1,:)=pca_scores(:,trNo,3);
            if training_decisions(trNo)==1
                plot3(x,y,z,'-r')
            else
                plot3(x,y,z,'-b')
            end
        end
        
        %Plot average trajectory in 3D space
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        x=zeros(1,no_timepoints);
        x(1,:)=mean(pca_scores(:,training_decisions==1,1),2);
        y=zeros(1,no_timepoints);
        y(1,:)=mean(pca_scores(:,training_decisions==1,2),2);
        z=zeros(1,no_timepoints);
        z(1,:)=mean(pca_scores(:,training_decisions==1,3),2);
        plot3(x,y,z,'-r')
        
        x=zeros(1,no_timepoints);
        x(1,:)=mean(pca_scores(:,training_decisions==0,1),2);
        y=zeros(1,no_timepoints);
        y(1,:)=mean(pca_scores(:,training_decisions==0,2),2);
        z=zeros(1,no_timepoints);
        z(1,:)=mean(pca_scores(:,training_decisions==0,3),2);
        plot3(x,y,z,'-b')
        
        %Plot average single points per trial in 3D space
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        x=zeros(1,sum(training_decisions==1));
        x(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==1,1),1);
        y=zeros(1,sum(training_decisions==1));
        y(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==1,2),1);
        z=zeros(1,sum(training_decisions==1));
        z(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==1,3),1);
        plot3(x,y,z,'or')
        
        x=zeros(1,sum(training_decisions==0));
        x(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==0,1),1);
        y=zeros(1,sum(training_decisions==1));
        y(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==0,2),1);
        z=zeros(1,sum(training_decisions==1));
        z(1,:)=mean(pca_scores(time_to_eventLDA>0,training_decisions==0,3),1);
        plot3(x,y,z,'ob')
        
        %Now do the PCA for mean(time>0)
    
        these_neural_recs=zeros(tr_no_trials,no_neurons);
        these_neural_recs(:,:)=mean(training_neural_recordings(:,:,time_to_eventLDA>0),3);
        [coeff,pca_scores_past_zero,latent] = pca(these_neural_recs);
        
        %Plot average single points per trial in 3D space
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        x=zeros(1,sum(training_decisions==1));
        x(1,:)=pca_scores_past_zero(training_decisions==1,1);
        y=zeros(1,sum(training_decisions==1));
        y(1,:)=pca_scores_past_zero(training_decisions==1,2);
        z=zeros(1,sum(training_decisions==1));
        z(1,:)=pca_scores_past_zero(training_decisions==1,3);
        plot3(x,y,z,'or')
        
        x=zeros(1,sum(training_decisions==0));
        x(1,:)=pca_scores_past_zero(training_decisions==0,1);
        y=zeros(1,sum(training_decisions==0));
        y(1,:)=pca_scores_past_zero(training_decisions==0,2);
        z=zeros(1,sum(training_decisions==0));
        z(1,:)=pca_scores_past_zero(training_decisions==0,3);
        plot3(x,y,z,'ob')
        
        %Now do the PCA for all points time>0
        these_neural_recs=zeros(tr_no_trials,no_neurons*length(time_to_eventLDA>0));
        ii_time_first=find(time_to_eventLDA>0,1,'first');
        ii_comps=0;
        for ii_time=ii_time_first:length(time_to_eventLDA)
            these_neural_recs(:,ii_comps+1:ii_comps+no_neurons)=training_neural_recordings(:,:,ii_time);
            ii_comps=ii_comps+no_neurons;
        end
        [coeff,pca_scores_past_zero,latent] = pca(these_neural_recs);
        
        %Plot average single points per trial in 3D space
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        x=zeros(1,sum(training_decisions==1));
        x(1,:)=pca_scores_past_zero(training_decisions==1,1);
        y=zeros(1,sum(training_decisions==1));
        y(1,:)=pca_scores_past_zero(training_decisions==1,2);
        z=zeros(1,sum(training_decisions==1));
        z(1,:)=pca_scores_past_zero(training_decisions==1,3);
        plot3(x,y,z,'or')
        
        x=zeros(1,sum(training_decisions==0));
        x(1,:)=pca_scores_past_zero(training_decisions==0,1);
        y=zeros(1,sum(training_decisions==0));
        y(1,:)=pca_scores_past_zero(training_decisions==0,2);
        z=zeros(1,sum(training_decisions==0));
        z(1,:)=pca_scores_past_zero(training_decisions==0,3);
        plot3(x,y,z,'ob')
        
    end
    
    pffft=1;
    
else
    handles_out2.decoding_processed=0;
end

%     %Note that the accuracy output is only for the training trial set
%     handles_out2.accuracy=accuracy;
%     handles_out2.sh_accuracy=sh_accuracy;
%     handles_out2.mean_accuracy=mean(accuracy(time>=0));
%     handles_out2.mean_sh_accuracy=mean(sh_accuracy(time>=0));
%     handles_out2.delta_odor=mean(delta_odor);
%     handles_out2.delta_odor_on_reinf_on=mean(delta_odor_on_reinf_on);
%     handles_out2.delta_reinf=mean(delta_reinf);
%     handles_out2.no_sp_trials=handles_out.no_sp_trials;
%     handles_out2.no_sm_trials=handles_out.no_sm_trials;
%     
%     %Now use Bishop's (2006) majority rule
%     winning_label=zeros(1,Nall);
%     for this_tr_No=1:Nall
%         these_labels=zeros(1,sum(time>=0));
%         these_labels(1,:)=training_output_labels(time>=0,this_tr_No);
%         if sum(these_labels==1)>sum(these_labels==0)
%             winning_label(this_tr_No)=1;
%         else
%             if sum(these_labels==0)>sum(these_labels==1)
%                 winning_label(this_tr_No)=0;
%             else
%                 if rand>0.5
%                     winning_label(this_tr_No)=1;
%                 else
%                     winning_label(this_tr_No)=0;
%                 end
%             end
%         end
%     end
%     
%     handles_out2.bishop_accuracy=sum(winning_label==training_decisions)/Nall;
%     these_bishop_sh_accuracy=[];
%     for ii_no=1:5
%         shuffled_training_decisions=training_decisions(randperm(Nall));
%         these_bishop_sh_accuracy(ii_no)=sum(winning_label==shuffled_training_decisions)/Nall;
%     end
%     handles_out2.bishop_sh_accuracy=mean(these_bishop_sh_accuracy);
%     
%     handles_out2.winning_label=winning_label;
%     handles_out2.training_decisions=training_decisions;
%     
%     %Now decode for the rest of the trials using the entire training set for training
%     if sum(p<=p_threshold)>0
%         
%         
%         %Calculate the z values
%         z_neural_recordings=zeros(all_trials,handles_out.no_components,no_timepoints);
%         for ii_neurons=1:no_neurons
%             z_neural_recordings(:,ii_neurons,:)=(neural_recordings(:,ii_neurons,:)-mean_per_neuron(ii_neurons))/STD_per_neuron(ii_neurons);
%         end
%         
%         
%         
%         
%         for time_point=1:no_timepoints
%             
%             %dFF per trial per component
%             measurements=zeros(Nall,sum(p<p_threshold));
%             measurements(:,:)=z_training_neural_recordings(:,p<p_threshold,time_point);
%             
%             
%             
%             correct_predict=[];
%             correct_predict_shuffled=[];
%             
%             
%             
%             %Now train the decoder with all trials in the decoding data set
%             
%             
%             
%             %Store the training data in a table.
%             tblTrn=[];
%             tblTrn = array2table(measurements);
%             
%             %Store the decisions in Y
%             Y=training_decisions;
%             
%             %Train a discriminant analysis model using the training set and default options.
%             %By default this is a regularized linear discriminant analysis (LDA)
%             switch MLalgo
%                 case 1
%                     Mdl = fitcdiscr(tblTrn,Y);
%                 case 2
%                     Mdl = fitcsvm(tblTrn,Y);
%                 case 3
%                     Mdl = fitcnb(tblTrn,Y);
%                 case 4
%                     Mdl = fitcnet(tblTrn,Y);
%                 case 5
%                     Mdl = fitctree(tblTrn,Y);
%             end
%             
%             %Now predict the outcome for all trials excluding the training set
%             all_measurements=zeros(all_trials,sum(p<p_threshold));
%             all_measurements(:,:)=z_neural_recordings(:,p<p_threshold,time_point);
%             
%             %Predict labels for all trials excluding the training set
%             for jj=1:all_trials
%                 if sum(jj==ii_all)>0
%                     [label,score] = predict(Mdl,all_measurements(jj,:));
%                     
%                     %label is the predicted label, and score is the predicted class
%                     %posterior probability
%                     
%                     if label==decisions(jj)
%                         handles_out2.correct_predict(time_point,jj)=1;
%                     else
%                         handles_out2.correct_predict(time_point,jj)=0;
%                     end
%                     
%                     jj_shuffled=randperm(all_trials);
%                     
%                     if label==decisions(jj_shuffled(jj))
%                         handles_out2.correct_predict_shuffled(time_point,jj)=1;
%                     else
%                         handles_out2.correct_predict_shuffled(time_point,jj)=0;
%                     end
%                     
%                 end
%             end
%             
%             
%             if show_figures==1
%                 fprintf(1, ['For timepoint %d accuracy= %d and shuffled accuracy= %d\n'],time_point,accuracy(time_point),sh_accuracy(time_point));
%             end
%         end
%         
%         
%     end
%     
%     if show_figures==1
%         figNo=figNo+1;
%         try
%             close(figNo)
%         catch
%         end
%         
%         figure(figNo)
%         
%         %         subplot(1,2,1)
%         hold on
%         
%         
%         %Plot the bounded line for the 5 percentile for the shuffled trials
%         per95=prctile(sh_accuracy(1,:),95);
%         per5=prctile(sh_accuracy(1,:),5);
%         CIsh=[mean(sh_accuracy)-per5 per95-mean(sh_accuracy)]';
%         [hlCR, hpCR] = boundedline([time_to_eventLDA(1) time_to_eventLDA(end)],[mean(sh_accuracy) mean(sh_accuracy)], CIsh', 'r');
%         
%         %Plot the accuracy
%         plot(time',accuracy,'-k')
%         
%         %Now plot the bootstrapped confidence intervals for the original and
%         
%         %For the accuracy use odor 0 to the end
%         CI = bootci(1000, {@mean, accuracy(time>=0)},'type','cper');
%         bar_offset=max(time)+1;
%         plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
%         plot([bar_offset],mean(accuracy(time>=0)),'ok')
%         
%         %For the shuffled accuracy use all points
%         CI = bootci(1000, {@mean, sh_accuracy},'type','cper');
%         bar_offset=max(time)+1;
%         plot([bar_offset bar_offset],CI,'-r','LineWidth',3)
%         plot([bar_offset],mean(sh_accuracy),'ok')
%         
%         %Odor on markers
%         plot([0 0],[0 1.1],'-k')
%         odorhl=plot([0 mean(delta_odor)],[0.32 0.32],'-k','LineWidth',5);
%         plot([mean(delta_odor) mean(delta_odor)],[0 1.10],'-k')
%         
%         %Reinforcement markers
%         plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 1.10],'-r')
%         reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[.32 .32],'-r','LineWidth',5);
%         plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 1.10],'-r')
%         
%         
%         ylim([00 1.10])
%         
%         xlabel('Time (sec)')
%         ylabel('Accuracy')
%     end

pffft=1;


