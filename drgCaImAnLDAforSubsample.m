function handles_par=drgCaImAnLDAforSubsample(dFF_trial_mask,time_to_eventLDA,...
    all_lda_no_comp,min_trials,handles_par,caimanhandles,no_trial_windows,all_lda_fileNo...
    ,number_of_replicates,first_num_odor_trials,num_odor_trials_dFF,perCorr,all_lda_events,...
    all_lda_input_timecourse,figNo,supertitle_description,per_file_coms,no_of_components)


%Now find out what happens if we choose smaller subsets of components
%The components are chosen randomly and they are the same for all trials within each file

%First find the comps for each trial and save them in
% compsubset(components,replicates)
all_comps=[];

szwins=size(caimanhandles.caimandr_choices.wins);
% handles_par(no_trial_windows).number_of_components=[[1:20] [22:2:30] [35:5:45] [50:10:100] [120:20:160] [190:30:280]];

%Note that I am doing this only for percent correct > 80%, no_trial_windows=3
no_trial_windows=3;
handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
dFF_trial_mask=[];
jj=0;


% if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)

% fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
pct_windows=[45 65;65 80;80 100.1];

for ii=1:num_odor_trials_dFF
    if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
        dFF_trial_mask(ii)=1;
        jj=jj+1;
        handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
        events{jj,1}=all_lda_events{ii};
        if strcmp(events{jj,1},'S+')
            %S+
            per_targets(1,jj)=1;
            %S-
            per_targets(2,jj)=0;
        else
            %S+
            per_targets(1,jj)=0;
            %S-
            per_targets(2,jj)=1;
        end
    else
        dFF_trial_mask(ii)=0;
    end
end

handles_par(no_trial_windows).number_of_components=no_of_components;

for ii_no_comps=1:length(handles_par(no_trial_windows).number_of_components)+1
    for winNo=1:szwins(1)
        handles_par(no_trial_windows).timewin(winNo).lda_delta_N(ii_no_comps).no_ldas=0;
    end
    
    
    
    %dFF per trial per component
    which_file=[];
    which_file=all_lda_fileNo(logical(dFF_trial_mask));
    no_comps=[];
    no_comps=all_lda_no_comp(logical(dFF_trial_mask));
    
    
    for fileNo=min(which_file):max(which_file)
        
        %If there are enough trials process the LDA
        ii=find(which_file==fileNo,1,'first');
        if isempty(ii)
            N=0
        else
            N=sum(which_file(ii)==which_file);
        end
        
        
        if N>=min_trials
            
            if ii_no_comps==length(handles_par(no_trial_windows).number_of_components)+1
                %This will use all components
                do_lda=1;
                
            else
                %This will use a subset of components
                if  handles_par(no_trial_windows).number_of_components(ii_no_comps)<no_comps(ii)
                    do_lda=1;
                else
                    do_lda=0;
                end
            end
            
            if do_lda==1
                
                %Choose unique subsets of components
                if ii_no_comps<length(handles_par(no_trial_windows).number_of_components)+1
                    compsubset=[];
                    
                    no_replicates=number_of_replicates;
                    
                    while(size(compsubset,1)<handles_par(no_trial_windows).number_of_components(ii_no_comps))
                        compsubset=[compsubset;randi(no_comps(ii),[handles_par(no_trial_windows).number_of_components(ii_no_comps) no_replicates])];
                        compsubset=unique(compsubset,'rows');
                    end
                    compsubset=compsubset(1:handles_par(no_trial_windows).number_of_components(ii_no_comps),:);
                    this_no_comps=handles_par(no_trial_windows).number_of_components(ii_no_comps);
                    
                else
                    compsubset=[1:no_comps(ii)]';
                    this_no_comps=no_comps(ii);
                    no_replicates=1;
                end
                
                all_comps(fileNo).comps(ii_no_comps).compsubset=compsubset;
                
            end
            
        end
    end
    
end

% %Note that I am doing this only for percent correct > 80%, no_trial_windows=3
% no_trial_windows=3;
% handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
% dFF_trial_mask=[];
% jj=0;
% 
% 
% % if caimanhandles.caimandr_choices.start_reversal>length(first_num_odor_trials)
% 
% % fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
% pct_windows=[45 65;65 80;80 100.1];
% 
% for ii=1:num_odor_trials_dFF
%     if (perCorr(ii)>=pct_windows(no_trial_windows,1))&(perCorr(ii)<pct_windows(no_trial_windows,2))
%         dFF_trial_mask(ii)=1;
%         jj=jj+1;
%         handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
%         events{jj,1}=all_lda_events{ii};
%         if strcmp(events{jj,1},'S+')
%             %S+
%             per_targets(1,jj)=1;
%             %S-
%             per_targets(2,jj)=0;
%         else
%             %S+
%             per_targets(1,jj)=0;
%             %S-
%             per_targets(2,jj)=1;
%         end
%     else
%         dFF_trial_mask(ii)=0;
%     end
% end
% else
%     if no_trial_windows==1
%         %Forward trials
%         fprintf(1, '\n\nLDA processed for dF/F for trials before reversal \n');
%         
%         for ii=1:num_odor_trials_dFF
%             if (trial_dFF(ii)>=1)&(trial_dFF(ii)<=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal)-1)
%                 dFF_trial_mask(ii)=1;
%                 jj=jj+1;
%                 handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
%                 events{jj,1}=all_lda_events{ii};
%                 if strcmp(events{jj,1},'S+')
%                     %S+
%                     per_targets(1,jj)=1;
%                     %S-
%                     per_targets(2,jj)=0;
%                 else
%                     %S+
%                     per_targets(1,jj)=0;
%                     %S-
%                     per_targets(2,jj)=1;
%                 end
%             else
%                 dFF_trial_mask(ii)=0;
%             end
%         end
%         
%     else
%         %Trials at end of reversal
%         fprintf(1, '\n\nLDA processed for dF/F for trials after reversal \n');
%         for ii=1:num_odor_trials_dFF
%             if (trial_dFF(ii)>=max(trial_dFF)-100)
%                 dFF_trial_mask(ii)=1;
%                 jj=jj+1;
%                 handles_par(no_trial_windows).perCorr(jj)=perCorr(ii);
%                 events{jj,1}=all_lda_events{ii};
%                 if strcmp(events{jj,1},'S+')
%                     %S+
%                     per_targets(1,jj)=1;
%                     %S-
%                     per_targets(2,jj)=0;
%                 else
%                     %S+
%                     per_targets(1,jj)=0;
%                     %S-
%                     per_targets(2,jj)=1;
%                 end
%             else
%                 dFF_trial_mask(ii)=0;
%             end
%         end
%     end
%     
% end

Nall=sum(dFF_trial_mask);


fprintf(1,'For S+ vs S-\n')

%handles_par(no_trial_windows).number_of_components=[1 2 [5:5:30] [40:10:100] [120:20:160] [190:30:280]];

for ii_no_comps=1:length(handles_par(no_trial_windows).number_of_components)+1
    %Choose unique subsets of components
    
    
    for winNo=1:szwins(1)
        first_timepoint=find(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1),1,'first');
        last_timepoint=find(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,2),1,'first');
        no_timepoints=length(time_to_eventLDA(first_timepoint:last_timepoint));
        
        handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas=0;
        
        for time_point=first_timepoint:last_timepoint
            
            %dFF per trial per component
            measurements=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))));
            this_all_lda_input_timecourse=zeros(Nall,max(all_lda_no_comp(logical(dFF_trial_mask))))';
            this_all_lda_input_timecourse(:,:)=all_lda_input_timecourse(time_point,1:max(all_lda_no_comp(logical(dFF_trial_mask))),logical(dFF_trial_mask));
            measurements(:,:)=this_all_lda_input_timecourse';
            which_file=[];
            which_file=all_lda_fileNo(logical(dFF_trial_mask));
            no_comps=[];
            no_comps=all_lda_no_comp(logical(dFF_trial_mask));
            
            ii_predict=0;
            fileNo_pred=[];
            these_files=[];
            ii_files=0;
            correct_predict=[];
            correct_predict_shuffled=[];
            dimensionality=[];
            
            for ii=1:Nall
                
                %If there are enough trials process the LDA
                N=sum(which_file(ii)==which_file);
                
                
                if N>=min_trials
                    
                    
                    
                    if ii_no_comps==length(handles_par(no_trial_windows).number_of_components)+1
                        %This will use all components
                        do_lda=1;
                    else
                        %This will use a subset of components
                        if  handles_par(no_trial_windows).number_of_components(ii_no_comps)<no_comps(ii)
                            do_lda=1;
                        else
                            do_lda=0;
                        end
                    end
                    
                    if do_lda==1
                        
                        compsubset=[];
                        fileNo=which_file(ii);
                        if sum(these_files==fileNo)==0
                            ii_files=ii_files+1;
                            these_files(ii_files)=fileNo;
                        end
                        compsubset=all_comps(fileNo).comps(ii_no_comps).compsubset;
                        
                        ii_predict=ii_predict+1;
                        fileNo_pred(ii_predict)=fileNo;
                        
                        %Choose number of replicates
                        if ii_no_comps<length(handles_par(no_trial_windows).number_of_components)+1
                            no_replicates=number_of_replicates;
                        else
                            no_replicates=1;
                        end
                        
                        parfor ii_sub=1:no_replicates
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
                            tblTrn = array2table(measurements(logical(idxTrn),sort(compsubset(:,ii_sub))));
                            these_events=[];
                            all_these_events=[];
                            noallEvs=0;
                            noEvs=0;
                            
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
                            Mdl = fitcdiscr(tblTrn,'Y');
                            
                            
                            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                            [label,score] = predict(Mdl,measurements(logical(idxTest),sort(compsubset(:,ii_sub))));
                            
                            
                            correct_predict(ii_sub,ii_predict)=strcmp(events{ii},label);
                            ii_shuffled=randperm(N);
                            correct_predict_shuffled(ii_sub,ii_predict)=strcmp(all_these_events{ii_shuffled(ii_event)},label);
                            
                            %Get the data
                            %Columns: cells, Rows: dF/F
                            Signal=measurements(logical(idxTrn),sort(compsubset(:,ii_sub)));
                            dimensionality(ii_sub,ii_predict) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
                            
                        end
                        
                    end
                    
                end
            end
            
            if ii_predict>0
                handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).percent_correct...
                    (handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                    100*sum(correct_predict')/ii_predict;
                handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).percent_correct_shuffled...
                    (handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                    100*sum(correct_predict_shuffled')/ii_predict;
                handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).dimensionality...
                    (handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                    mean(dimensionality');
                
                for ii_f=1:length(these_files)
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).these_files=these_files;
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).percent_correct...
                        (:,handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                        handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                        100*sum(correct_predict(:,fileNo_pred==these_files(ii_f))')/sum((fileNo_pred==these_files(ii_f)));
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).percent_correct_shuffled...
                        (:,handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                        handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                        100*sum(correct_predict_shuffled(:,fileNo_pred==these_files(ii_f))')/sum((fileNo_pred==these_files(ii_f)));
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).dimensionality...
                        (handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+1:...
                        handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates)=...
                        mean(dimensionality(fileNo_pred==these_files(ii_f))');
                    handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).compsubset=...
                        all_comps(these_files(ii_f)).comps(ii_no_comps).compsubset;
                end
                
                handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas=handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas+no_replicates;
                
                %Now calculate percent correct per file
                
                szcompsub=size(compsubset);
                this_no_comps=szcompsub(1);
                fprintf(1, 'LDA percent correct computed in percent correct window No %d for %d components for window %d\n',no_trial_windows,this_no_comps,winNo);
                
            end
        end
        if (winNo==2)&(no_trial_windows==3)
            pffft=1;
        end
    end
end

%Calculate the mean and CI for percent correct and plot the percent correct vs number of components
min_ii_no_comps=zeros(1,szwins(1));
for winNo=1:szwins(1)
    
    percent_correct_exists=zeros(1,length(handles_par(no_trial_windows).number_of_components)+1);
    mean_percent_correct=zeros(1,length(handles_par(no_trial_windows).number_of_components)+1);
    CI_percent_correct=zeros(2,length(handles_par(no_trial_windows).number_of_components)+1);
    var_percent_correct=zeros(2,length(handles_par(no_trial_windows).number_of_components)+1);
    mean_percent_correct_shuffled=zeros(1,length(handles_par(no_trial_windows).number_of_components)+1);
    CI_percent_correct_shuffled=zeros(2,length(handles_par(no_trial_windows).number_of_components)+1);
    var_percent_correct_shuffled=zeros(2,length(handles_par(no_trial_windows).number_of_components)+1);
    mean_dimensionality=zeros(1,length(handles_par(no_trial_windows).number_of_components)+1);
    CI_dimensionality=zeros(2,length(handles_par(no_trial_windows).number_of_components)+1);
    
    
    for ii_no_comps=1:length(handles_par(no_trial_windows).number_of_components)+1
        if handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).no_ldas>0
            p_corrs=[];
            p_corrs_sh=[];
            dims=[];
            for ii_f=1:length(handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file)
                p_corrs=[p_corrs handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).percent_correct];
                p_corrs_sh=[p_corrs_sh handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).per_file(ii_f).percent_correct_shuffled];
                dims=[dims handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps).dimensionality];
            end
            mean_percent_correct(ii_no_comps)=mean(p_corrs);
            CI_percent_correct(:,ii_no_comps) = bootci(1000, @mean, p_corrs);
            var_percent_correct(:,ii_no_comps) = var(p_corrs);
            percent_correct_exists(ii_no_comps)=1;
            mean_percent_correct_shuffled(ii_no_comps)=mean(p_corrs_sh);
            CI_percent_correct_shuffled(:,ii_no_comps) = bootci(1000, @mean, p_corrs_sh);
            var_percent_correct_shuffled(:,ii_no_comps) = var(p_corrs_sh);
            mean_dimensionality(ii_no_comps)= mean(dims);
            CI_dimensionality(:,ii_no_comps) = bootci(1000, @mean, dims);
            
        end
    end
    
    if sum(percent_correct_exists)>2
        CI_percent_correct_shuffled(1,:)=mean_percent_correct_shuffled-CI_percent_correct_shuffled(1,:);
        CI_percent_correct_shuffled(2,:)=CI_percent_correct_shuffled(2,:)-mean_percent_correct_shuffled;
        
        CI_percent_correct(1,:)=mean_percent_correct-CI_percent_correct(1,:);
        CI_percent_correct(2,:)=CI_percent_correct(2,:)-mean_percent_correct;
        
        %Note that the last numcomps uses all components
        num_comps=handles_par(no_trial_windows).number_of_components;
        num_comps(end+1)=floor(mean(no_comps)); %This is all the components
        
        handles_par(no_trial_windows).timewin(winNo).percent_correct_exists=percent_correct_exists;
        handles_par(no_trial_windows).timewin(winNo).mean_percent_correct=mean_percent_correct;
        handles_par(no_trial_windows).timewin(winNo).CI_percent_correct=CI_percent_correct;
        handles_par(no_trial_windows).timewin(winNo).mean_percent_correct_shuffled=mean_percent_correct_shuffled;
        handles_par(no_trial_windows).timewin(winNo).CI_percent_correct_shuffled=CI_percent_correct_shuffled;
        handles_par(no_trial_windows).timewin(winNo).num_comps=num_comps;
        
        
        
        %plot the percent correct vs no components
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        boundedline(num_comps(1,logical(percent_correct_exists))', mean_percent_correct_shuffled(1,logical(percent_correct_exists))', CI_percent_correct_shuffled(:,logical(percent_correct_exists))', 'b');
        p1=plot(num_comps(1,logical(percent_correct_exists))', mean_percent_correct_shuffled(1,logical(percent_correct_exists))','ob','MarkerFace','b','MarkerSize',3);
        
        boundedline(num_comps(1,logical(percent_correct_exists))', mean_percent_correct(1,logical(percent_correct_exists))', CI_percent_correct(:,logical(percent_correct_exists))', 'r');
        p2=plot(num_comps(1,logical(percent_correct_exists))', mean_percent_correct(1,logical(percent_correct_exists))','or','MarkerFace','r','MarkerSize',3);
        
        ylim([0 110])
        xlim([min(num_comps(1,logical(percent_correct_exists))') max(num_comps(1,logical(percent_correct_exists))')])
        these_xlab=xticklabels;
        these_xlab{end}='all';
        xticklabels(these_xlab)
        xlabel('Number of components used for LDA')
        ylabel('Percent correct')
        legend([p1 p2],'shuffled','original')
        title(['Window  from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec for ' supertitle_description{no_trial_windows}])
        
        
        %Plot the cumulative histogram for representative
        %subsampling
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        %Find the minimum
        mpc=mean_percent_correct(logical(percent_correct_exists));
        ii_offset=3;
        [min_pc min_pc_ii]=min(mpc(ii_offset:end));
        ii_no_comps_min=min_pc_ii+ii_offset-1;
        min_ii_no_comps(winNo)=ii_no_comps_min;
        [f_pc,x_pc] = drg_ecdf(handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps_min).percent_correct);
        plot(x_pc,f_pc,'-b')
        
        %Find the maximum after the minimum
        [max_after_pc max_after_ii]=max(mpc(ii_no_comps_min:end-1));
        ii_no_comps_max_after=max_after_ii+ii_no_comps_min-1;
        [f_pc,x_pc] = drg_ecdf(handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps_max_after).percent_correct);
        plot(x_pc,f_pc,'-r')
        
        %Find the maximum before the minimum
        [max_before_pc ii_no_comps_max_before]=max(mpc(1:ii_no_comps_min));
        ii_no_comps_max=ii_no_comps_max_before+ii_offset-1;
        max_ii_no_comps_before(winNo)=ii_no_comps_max;
        [f_pc,x_pc] = drg_ecdf(handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(ii_no_comps_max_before).percent_correct);
        plot(x_pc,f_pc,'-m')
        
        %Plot cum histo for n=1
        [f_pc,x_pc] = drg_ecdf(handles_par(no_trial_windows).timewin(winNo).lda_delta_N1(1).percent_correct);
        plot(x_pc,f_pc,'-k')
        
        legend([num2str(num_comps(ii_no_comps_min)) ' components'],[num2str(num_comps(ii_no_comps_max_after)) ' components'],...
            [num2str(num_comps(ii_no_comps_max_before)) ' components'],['1 component'])
        title(['Window  from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec for ' supertitle_description{no_trial_windows}])
        ylabel('Probability')
        xlabel('LDA percent correct')
        
        %Plot the variance
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        plot(num_comps(1,logical(percent_correct_exists))', var_percent_correct(1,logical(percent_correct_exists))', '-or','MarkerFace','r');
        %plot(num_comps(1,logical(percent_correct_exists))', var_percent_correct_shuffled(1,logical(percent_correct_exists))', '-ob');
        
        xlim([min(num_comps(1,logical(percent_correct_exists))') max(num_comps(1,logical(percent_correct_exists))')])
        these_xlab=xticklabels;
        these_xlab{end}='all';
        xticklabels(these_xlab)
        xlabel('Number of components used for LDA')
        ylabel('Variance')
        legend([p1 p2],'shuffled','original')
        title(['Window  from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec for ' supertitle_description{no_trial_windows}])
        
        
        %plot the dimensionality vs no components
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        figure(figNo)
        hold on
        
        CI_dimensionality(1,:)=mean_dimensionality-CI_dimensionality(1,:);
        CI_dimensionality(2,:)=CI_dimensionality(2,:)-mean_dimensionality;
        
        boundedline(num_comps(1,logical(percent_correct_exists))', mean_dimensionality(1,logical(percent_correct_exists))', CI_dimensionality(:,logical(percent_correct_exists))', 'r');
        p2=plot(num_comps(1,logical(percent_correct_exists))', mean_dimensionality(1,logical(percent_correct_exists))','or','MarkerFace','r','MarkerSize',3);
        
        %             ylim([0 110])
        xlim([min(num_comps(1,logical(percent_correct_exists))') max(num_comps(1,logical(percent_correct_exists))')])
        these_xlab=xticklabels;
        these_xlab{end}='all';
        xticklabels(these_xlab)
        xlabel('Number of components used for LDA')
        ylabel('Dimensionality')
        
        title(['Window  from ' num2str(caimanhandles.caimandr_choices.wins(winNo,1)) ' to ' num2str(caimanhandles.caimandr_choices.wins(winNo,2)) ' sec for ' supertitle_description{no_trial_windows}])
        
        
        
        
    end
end
