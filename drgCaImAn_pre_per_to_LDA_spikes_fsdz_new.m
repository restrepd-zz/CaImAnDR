function handles_out2=drgCaImAn_pre_per_to_LDA_spikes_fsdz_new(pre_perBatchPathName, pre_perFileName, p_thr_less_than,p_thr_more_than, MLalgo, show_figures...
    ,no_sp_sm_trials_to_use,first_sp_sm_trial_no,figNo,fileNo,this_cost)
%
% reads the pre_per file and saves .mat files to process with
% Kording's lab
%Create input and target vectors leaving one trial out
%



close all
clearvars -except pre_perBatchPathName pre_perFileName p_thr_less_than p_thr_more_than MLalgo show_figures no_sp_sm_trials_to_use first_sp_sm_trial_no this_cost figNo this_cost
warning('off')

%Set the cost
handles_out2.cost=this_cost;
% this_cost=[0 5;5 0];


% What is the meaning of 1 and 2 in example ([0,1;2,0])?
%
% This matrix is of form
%
% 0 1
% 2 0
%
% thus it means:
%
% if point has class 1 and we assign it class 1, penalty is 0 (correct classification)
% if point has class 1 and we assign it class 2, penalty is 1
% if point has class 2 and we assign it class 1, penalty is 2
% if point has class 2 and we assign it class 2, penalty is 0 (correct classification)

% simulation=0;

classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';

handles_out2.classifier_names=classifier_names;


min_no_trials=10;
 

handles_out2.pre_perBatchPathName=pre_perBatchPathName;
handles_out2.pre_perFileName=pre_perFileName;
handles_out2.p_thr_less_than=p_thr_less_than;
handles_out2.p_thr_more_than=p_thr_more_than;
handles_out2.MLalgo=MLalgo;
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
    
    if (no_sp_sm_trials_to_use>handles_out.no_sp_trials)||(no_sp_sm_trials_to_use>handles_out.no_sm_trials)
        no_sp_sm_trials_to_use=min([handles_out.no_sp_trials handles_out.no_sm_trials]);
    end
    %First and last sp trial numbers
    first_sp_trial=first_sp_sm_trial_no;
    last_sp_trial=first_sp_sm_trial_no+no_sp_sm_trials_to_use-1;
   

    %First and last sm trial numbers
    first_sm_trial=first_sp_sm_trial_no;
    last_sm_trial=first_sp_sm_trial_no+no_sp_sm_trials_to_use-1;


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
                training_neural_recordings(tr_trNo,traceNo,:)=splus_traces(ii,1:time_bins);
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
            neural_recordings(trialNo+handles_out.no_sm_trials,traceNo,:)=sminus_traces(ii,1:time_bins);
            if (trialNo>=first_sm_trial)&(trialNo<=last_sm_trial)
                training_neural_recordings(tr_trNo,traceNo,:)=sminus_traces(ii,1:time_bins);
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



if handles_out2.decoding_processed==1

    %Threshold out the spikes
    std_traces=std(training_neural_recordings(:));
    thr=5*std_traces;
    
    %These is the placeholders for the spike counts     
    training_neural_spikes=zeros(2*no_sp_sm_trials_to_use,handles_out.no_components,time_bins);

    for ii_comps=1:handles_out.no_components
        
        these_recordings=zeros(length(training_decisions),no_timepoints);
        these_recordings(:,:)=training_neural_recordings(:,ii_comps,:);
        baseline_dF=prctile(these_recordings(:),5);
         
        for ii_t=1:no_timepoints
            for trNo=1:2*no_sp_sm_trials_to_use
                if training_neural_recordings(trNo,ii_comps,ii_t)-baseline_dF>thr
                    training_neural_spikes(trNo,ii_comps,ii_t)=1;
                end
            end
        end
        
    end

    %First decode for the training trials leaving one trial out
    Nall=length(training_decisions);
    accuracy=zeros(1,no_timepoints);
    sh_accuracy=zeros(1,no_timepoints);

    

    %Calculate z from dFF
    mean_per_neuron=zeros(1,sum((p<=p_thr_less_than)&(p>=p_thr_more_than)));
    STD_per_neuron=zeros(1,sum((p<=p_thr_less_than)&(p>=p_thr_more_than)));
    
    z_training_neural_spikes=zeros(Nall,sum((p<=p_thr_less_than)&(p>=p_thr_more_than)),no_timepoints);
    jj_neurons=0;
    all_pre=[];
    for ii_neurons=1:no_neurons
        if (p(ii_neurons)<=p_thr_less_than)&(p(ii_neurons)>=p_thr_more_than)
            jj_neurons=jj_neurons+1;
            
            for ii_trials=1:tr_no_trials
                these_pre=zeros(1,sum((time_to_eventLDA<0)&(time_to_eventLDA>=-5)));
                these_pre(1,:)=training_neural_spikes(ii_trials,ii_neurons,(time_to_eventLDA<0)&(time_to_eventLDA>=-5));
                all_pre=[all_pre these_pre];
            end
            
            for ii_trials=1:tr_no_trials
                these_spikes=[];
                these_spikes=training_neural_spikes(ii_trials,ii_neurons,(time_to_eventLDA<0)&(time_to_eventLDA>=-5));
                z_training_neural_spikes(ii_trials,jj_neurons,:)=(training_neural_spikes(ii_trials,ii_neurons,:)-mean(these_spikes));
            end
        end
    end
    
    STD_all_neurons=std(all_pre);
    z_training_neural_spikes=z_training_neural_spikes/STD_all_neurons;

    handles_out2.time_to_eventLDA=time_to_eventLDA;

    if sum((p<=p_thr_less_than)&(p>=p_thr_more_than))>0

        gcp;
        handles_out2.decoding_processed=1;
        handles_out2.decisions=decisions;
        handles_out2.correct_predict=zeros(no_timepoints,all_trials);
        handles_out2.correct_predict_shuffled=zeros(no_timepoints,all_trials);
        training_output_labels=zeros(no_timepoints,Nall);
        handles_out2.this_cost=this_cost;

        %         What is clear meaning of 1 and 2 in example ([0,1;2,0])?
        %
        % This matrix is of form
        %
        % 0 1
        % 2 0
        % thus it means:
        %
        % if point has class 1 and we assign it class 1, penalty is 0 (correct classification)
        % if point has class 1 and we assign it class 2, penalty is 1
        % if point has class 2 and we assign it class 1, penalty is 2
        % if point has class 2 and we assign it class 2, penalty is 0 (correct classification)

        for time_point=1:no_timepoints

            %dFF per trial per component
            measurements=zeros(Nall,sum((p<=p_thr_less_than)&(p>=p_thr_more_than)));
            measurements(:,:)=z_training_neural_spikes(:,:,time_point);


            labels=[];
            timepoint_processed=[];
            correct_predict=[];
            correct_predict_shuffled=[];



            parfor ii=1:Nall
%                 for ii=1:Nall

                %Partition the data into training and test sets.

                %Create input and target vectors leaving one trial out
                %For per_input each column has the dF/F for one trial
                %each row is a single time point for dF/F for one of the cells
                %For per_target the top row is 1 if the odor is S+ and 0 if it is
                %S-, and row 2 has 1 for S-
                idxTrn=ones(Nall,1);
                idxTrn(ii)=0;
                idxTest=zeros(Nall,1);
                idxTest(ii)=1;

                %Store the training data in a table.
                tblTrn=[];
                tblTrn = array2table(measurements(logical(idxTrn),:));

                %Store the decisions in Y
                Y=training_decisions(logical(idxTrn));

                %Train a discriminant analysis model using the training set and default options.
                %By default this is a regularized linear discriminant analysis (LDA)
                timepoint_processed(ii)=1;

                try
                    
                    switch MLalgo
                        case 1
                            Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        case 2
                            Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                        case 4
                            Mdl = fitcnet(tblTrn,Y,'Cost',this_cost);
                        case 5
                            Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                    end
                    
%                     if MLalgo==2
%                         pffft=1;
%                     end
%                     Mdl.Cost(1,2) = 10;

                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                    labels(ii)=label;

                    %label is the predicted label, and score is the predicted class
                    %posterior probability

                    if label==training_decisions(ii)
                        correct_predict(ii)=1;
                    else
                        correct_predict(ii)=0;
                    end

                    ii_shuffled=randperm(Nall);

                    if label==training_decisions(ii_shuffled(ii))
                        correct_predict_shuffled(ii)=1;
                    else
                        correct_predict_shuffled(ii)=0;
                    end
                catch
                    %fit did not work
                    timepoint_processed(ii)=0;
                    correct_predict(ii)=-1;
                    correct_predict_shuffled(ii)=-1;
                    labels(ii)=-1;
                end
            end
            if (p_thr_less_than==1.1)&(p_thr_more_than==0.5)
                pffft=1;
            end
            accuracy(time_point)=mean(correct_predict(correct_predict~=-1));
            sh_accuracy(time_point)=mean(correct_predict_shuffled(correct_predict_shuffled~=-1));
            training_output_labels(time_point,:)=labels;



            for jj=1:Nall
                handles_out2.timepoint_processed(time_point,ii_all(jj))=timepoint_processed(jj);
                handles_out2.correct_predict(time_point,ii_all(jj))=correct_predict(jj);
                handles_out2.correct_predict_shuffled(time_point,ii_all(jj))=correct_predict_shuffled(jj);
            end

            if show_figures==1
                fprintf(1, ['For timepoint %d accuracy= %d and shuffled accuracy= %d\n'],time_point,accuracy(time_point),sh_accuracy(time_point));
            end
        end



        %Note that the accuracy output is only for the training trial set
        handles_out2.accuracy=accuracy;
        handles_out2.sh_accuracy=sh_accuracy;
        handles_out2.mean_accuracy=mean(accuracy(time>=0));
        handles_out2.mean_sh_accuracy=mean(sh_accuracy(time>=0));
        handles_out2.delta_odor=mean(delta_odor);
        handles_out2.delta_odor_on_reinf_on=mean(delta_odor_on_reinf_on);
        handles_out2.delta_reinf=mean(delta_reinf);
        handles_out2.no_sp_trials=handles_out.no_sp_trials;
        handles_out2.no_sm_trials=handles_out.no_sm_trials;
 
        if sum(correct_predict==-1)>0
            pffft=1;
        end

        %Now use Bishop's (2006) majority rule
        if handles_out2.decoding_processed==1
            winning_label=zeros(1,Nall);
            for this_tr_No=1:Nall
                these_labels=zeros(1,sum(time>=0));
                these_labels(1,:)=training_output_labels(time>=0,this_tr_No);
                if sum(these_labels==1)>sum(these_labels==0)
                    winning_label(this_tr_No)=1;
                else
                    if sum(these_labels==0)>sum(these_labels==1)
                        winning_label(this_tr_No)=0;
                    else
                        if rand>0.5
                            winning_label(this_tr_No)=1;
                        else
                            winning_label(this_tr_No)=0;
                        end
                    end
                end
            end

            handles_out2.bishop_accuracy=sum(winning_label==training_decisions)/Nall;
            these_bishop_sh_accuracy=[];
            for ii_no=1:5
                shuffled_training_decisions=training_decisions(randperm(Nall));
                these_bishop_sh_accuracy(ii_no)=sum(winning_label==shuffled_training_decisions)/Nall;
            end
            handles_out2.bishop_sh_accuracy=mean(these_bishop_sh_accuracy);

            handles_out2.winning_label=winning_label;
            handles_out2.training_decisions=training_decisions;
        else
            handles_out2.bishop_accuracy=[];
            handles_out2.bishop_sh_accuracy=[];
            handles_out2.winning_label=[];
            handles_out2.training_decisions=[];
        end

    else
        handles_out2.decoding_processed=0;
    end
    
    if show_figures==1
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        
        %         subplot(1,2,1)
        hold on
        
        
        %Plot the bounded line for the 5 percentile for the shuffled trials
        per95=prctile(sh_accuracy(1,:),95);
        per5=prctile(sh_accuracy(1,:),5);
        CIsh=[mean(sh_accuracy)-per5 per95-mean(sh_accuracy)]';
        [hlCR, hpCR] = boundedline([time_to_eventLDA(1) time_to_eventLDA(end)],[mean(sh_accuracy) mean(sh_accuracy)], CIsh', 'r');
        
        %Plot the accuracy
        plot(time',accuracy,'-k')
        
        %Now plot the bootstrapped confidence intervals for the original and
        
        %For the accuracy use odor 0 to the end
        CI = bootci(1000, {@mean, accuracy(time>=0)},'type','cper');
        bar_offset=max(time)+1;
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot([bar_offset],mean(accuracy(time>=0)),'ok')
        
        %For the shuffled accuracy use all points
        CI = bootci(1000, {@mean, sh_accuracy},'type','cper');
        bar_offset=max(time)+1;
        plot([bar_offset bar_offset],CI,'-r','LineWidth',3)
        plot([bar_offset],mean(sh_accuracy),'ok')
        
        %Odor on markers
        plot([0 0],[0 1.1],'-k')
        odorhl=plot([0 mean(delta_odor)],[0.32 0.32],'-k','LineWidth',5);
        plot([mean(delta_odor) mean(delta_odor)],[0 1.10],'-k')
        
        %Reinforcement markers
        plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)],[0 1.10],'-r')
        reinfhl=plot([mean(delta_odor_on_reinf_on) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[.32 .32],'-r','LineWidth',5);
        plot([mean(delta_odor_on_reinf_on)+mean(delta_reinf) mean(delta_odor_on_reinf_on)+mean(delta_reinf)],[0 1.10],'-r')
        
        
        ylim([00 1.10])
        
        xlabel('Time (sec)')
        ylabel('Accuracy')
    end
end
pffft=1;


