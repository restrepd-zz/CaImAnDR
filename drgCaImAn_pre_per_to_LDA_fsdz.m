function handles_out2=drgCaImAn_pre_per_to_LDA_fsdz(pre_perBatchPathName, pre_perFileName, p_threshold, MLalgo, show_figures)
%
% reads the pre_per file and saves .mat files to process with
% Kording's lab
%

close all
clearvars -except pre_perBatchPathName pre_perFileName p_threshold MLalgo show_figures

simulation=0;

classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';

handles_out2.classifier_names=classifier_names;

if nargin==0
    [pre_perFileName,pre_perBatchPathName] = uigetfile({'*pre_per.mat'},'Select the pre_per.mat file');
    fprintf(1, ['\ndrgCaImAn_pre_per_to_pydec run for ' pre_perFileName '\n\n']);
    
    p_threshold=0.05;
    
    %Choose the machine learning algorithm
    % 1 = linear discriminant analysis
    % 2 = SVM
    % 3 = Naive Bayes Classifier
    % 4 = neural network
    % 5 = decision tree
    MLalgo=1;
    show_figures=1;
end

min_no_trials=10;
figNo=0;

handles_out2.pre_perBatchPathName=pre_perBatchPathName;
handles_out2.pre_perFileName=pre_perFileName;
handles_out2.p_threshold=p_threshold;
handles_out2.MLalgo=MLalgo;
handles_out2.decoding_processed=1;

if ~iscell(pre_perFileName)
    load([pre_perBatchPathName pre_perFileName])
    
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
        
        %Now save the data for python
        pffft=1;
        
        time_bins=length(handles_out.time_to_eventSp);
        neural_recordings=zeros(handles_out.no_sp_trials+handles_out.no_sm_trials,handles_out.no_components,time_bins);
        decisions=zeros(1,handles_out.no_sp_trials+handles_out.no_sm_trials);
        time=time_to_eventSp;
        
        %Save S+
        
        trNo=1;
        ii=1;
        for trialNo=1:handles_out.no_sp_trials
            for traceNo=1:handles_out.no_components
                neural_recordings(trialNo,traceNo,:)=splus_traces(ii,:);
                ii=ii+1;
            end
            decisions(trNo)=1;
            trNo=trNo+1;
        end
        
        %Save S-
        ii=1;
        for trialNo=1:handles_out.no_sm_trials
            for traceNo=1:handles_out.no_components
                neural_recordings(trialNo+handles_out.no_sp_trials,traceNo,:)=sminus_traces(ii,:);
                ii=ii+1;
            end
            decisions(trNo)=0;
            trNo=trNo+1;
        end
        
        %Calculate p values for the difference per neuron
        no_neurons=handles_out.no_components;
        no_trials=trNo-1;
        no_timepoints=length(time_to_eventLDA);
        
        p=[];
        for ii_neuron=1:handles_out.no_components
            these_recordings_sp=zeros(handles_out.no_sp_trials,no_timepoints);
            these_recordings_sp=neural_recordings(logical(decisions),ii_neuron,:);
            
            these_recordings_sm=zeros(handles_out.no_sm_trials,no_timepoints);
            these_recordings_sm=neural_recordings(~logical(decisions),ii_neuron,:);
            
            mean_sp_odor=zeros(handles_out.no_sp_trials,1);
            mean_sp_odor(:,1)=mean(these_recordings_sp(:,time_to_eventLDA>0),2);
            
            mean_sm_odor=zeros(handles_out.no_sm_trials,1);
            mean_sm_odor(:,1)=mean(these_recordings_sm(:,time_to_eventLDA>0),2);
            
            [h,p(ii_neuron)]=ttest2(mean_sp_odor,mean_sm_odor);
        end
        
        handles_out2.p=p;
        
    else
        handles_out2.decoding_processed=0;
    end
else
    sp_trials_per_file=[];
    sm_trials_per_file=[];
    len_time_to_eventSp=[];
    len_time_to_eventSm=[];
    no_components=[];
    for ii_file=1:length(pre_perFileName)
        load([pre_perBatchPathName{ii_file} pre_perFileName{ii_file}])
        sp_trials_per_file(ii_file)=handles_out.no_sp_trials;
        sm_trials_per_file(ii_file)=handles_out.no_sm_trials;
        len_time_to_eventSp(ii_file)=length(handles_out.time_to_eventSp);
        len_time_to_eventSm(ii_file)=length(handles_out.time_to_eventSm);
        no_components(ii_file)=handles_out.no_components;
    end
    
    
    
    if (max(sp_trials_per_file)<min_no_trials)||(max(sm_trials_per_file)<min_no_trials)
        handles_out2.decoding_processed=0;
    else
        no_sp_trials=max([min_no_trials min(sp_trials_per_file)]);
        no_sm_trials=max([min_no_trials min(sm_trials_per_file)]);
        
        
        no_timepoints=min([min(len_time_to_eventSp) min(len_time_to_eventSm)]);
        time=handles_out.time_to_eventSp(1:no_timepoints);
        neural_recordings=zeros(no_sp_trials+no_sm_trials,sum(no_components),no_timepoints);
        decisions=zeros(1,no_sp_trials+no_sm_trials);
        
        no_neurons=0;
        for ii_file=1:length(pre_perFileName)
            
            load([pre_perBatchPathName{ii_file} pre_perFileName{ii_file}])
            
            %Save S+
            
            trNo=1;
            ii=1;
            for trialNo=1:no_sp_trials
                for traceNo=1:handles_out.no_components
                    neural_recordings(trialNo,traceNo+no_neurons,1:no_timepoints)=splus_traces(ii,1:no_timepoints);
                    ii=ii+1;
                end
                decisions(trNo)=1;
                trNo=trNo+1;
            end
            
            %Save S-
            ii=1;
            for trialNo=1:no_sm_trials
                for traceNo=1:handles_out.no_components
                    neural_recordings(trialNo+no_sp_trials,traceNo+no_neurons,1:no_timepoints)=sminus_traces(ii,1:no_timepoints);
                    ii=ii+1;
                end
                decisions(trNo)=0;
                trNo=trNo+1;
            end
            
            no_neurons=no_neurons+handles_out.no_components;
            
        end
        
        %Calculate p values for the difference per neuron
        p=[];
        for ii_neuron=1:no_neurons
            these_recordings_sp=zeros(no_sp_trials,no_timepoints);
            these_recordings_sp=neural_recordings(logical(decisions),ii_neuron,:);
            
            these_recordings_sm=zeros(no_sm_trials,no_timepoints);
            these_recordings_sm=neural_recordings(~logical(decisions),ii_neuron,:);
            
            mean_sp_odor=zeros(no_sp_trials,1);
            mean_sp_odor(:,1)=mean(these_recordings_sp(:,time>0),2);
            
            mean_sm_odor=zeros(no_sm_trials,1);
            mean_sm_odor(:,1)=mean(these_recordings_sm(:,time>0),2);
            
            [h,p(ii_neuron)]=ttest2(mean_sp_odor,mean_sm_odor);
        end
        
        handles_out2.p=p;
        handles_out2.time=time;
        no_trials=trNo-1;
    end
end


if handles_out2.decoding_processed==1
    if simulation==1
        %This is here to troubleshoot the algorithm
        PathName='/Users/restrepd/Documents/Projects/SFTP/Artigo_HomeOdor/AmylAcetateAcetophenone/';
        load([PathName 'simulated_data.mat'])
    end
    
    Nall=length(decisions);
    accuracy=zeros(1,no_timepoints);
    sh_accuracy=zeros(1,no_timepoints);
    
    %These are used to calculate the z values
    %I use time<0
    mean_per_neuron=zeros(1,no_neurons);
    STD_per_neuron=zeros(1,no_neurons);
    for ii_neurons=1:no_neurons
        all_pre=[];
        for ii_trials=1:no_trials
            these_pre=zeros(1,sum(time_to_eventLDA<0));
            these_pre(1,:)=neural_recordings(ii_trials,ii_neuron,time_to_eventLDA<0);
            all_pre=[all_pre these_pre];
        end
        mean_per_neuron(ii_neurons)=mean(all_pre);
        STD_per_neuron(ii_neurons)=std(all_pre);
    end
    
    z_neural_recordings=zeros(handles_out.no_sp_trials+handles_out.no_sm_trials,handles_out.no_components,no_timepoints);
    for ii_neurons=1:no_neurons
        z_neural_recordings(:,ii_neurons,:)=(neural_recordings(:,ii_neurons,:)-mean_per_neuron(ii_neurons))/STD_per_neuron(ii_neurons);
    end
    
    if sum(p<p_threshold)>0
        
        gcp;
        handles_out2.decoding_processed=1;
        for time_point=1:no_timepoints
            
            %dFF per trial per component
            measurements=zeros(Nall,sum(p<p_threshold));
            measurements(:,:)=z_neural_recordings(:,p<p_threshold,time_point);
            
            %Now calculate the z value
            
            
            scores=[];
            trials_processed=[];
            correct_predict=[];
            correct_predict_shuffled=[];
            lda_no_trials=[];
            
            
            parfor ii=1:Nall
                
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
                Y=decisions(logical(idxTrn));
                
                %Train a discriminant analysis model using the training set and default options.
                %By default this is a regularized linear discriminant analysis (LDA)
                switch MLalgo
                    case 1
                        Mdl = fitcdiscr(tblTrn,Y);
                    case 2
                        Mdl = fitcsvm(tblTrn,Y);
                    case 3
                        Mdl = fitcnb(tblTrn,Y);
                    case 4
                        Mdl = fitcnet(tblTrn,Y);
                    case 5
                        Mdl = fitctree(tblTrn,Y);
                end
                
                
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                [label,score] = predict(Mdl,measurements(logical(idxTest),:));
                
                %label is the predicted label, and score is the predicted class
                %posterior probability
                
                if label==decisions(ii)
                    correct_predict(ii)=1;
                else
                    correct_predict(ii)=0;
                end
                
                ii_shuffled=randperm(Nall);
                
                if label==decisions(ii_shuffled(ii))
                    correct_predict_shuffled(ii)=1;
                else
                    correct_predict_shuffled(ii)=0;
                end
                
            end
            
            accuracy(time_point)=mean(correct_predict);
            sh_accuracy(time_point)=mean(correct_predict_shuffled);
            
            if show_figures==1
                fprintf(1, ['For timepoint %d accuracy= %d and shuffled accuracy= %d\n'],time_point,accuracy(time_point),sh_accuracy(time_point));
            end
        end
        
    else
        handles_out2.decoding_processed=0;
    end
    
    handles_out2.accuracy=accuracy;
    handles_out2.sh_accuracy=sh_accuracy;
    handles_out2.mean_accuracy=mean(accuracy(time>=0));
    handles_out2.mean_sh_accuracy=mean(sh_accuracy(time>=0));
    handles_out2.delta_odor=mean(delta_odor);
    handles_out2.delta_odor_on_reinf_on=mean(delta_odor_on_reinf_on);
    handles_out2.delta_reinf=mean(delta_reinf);
    
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


