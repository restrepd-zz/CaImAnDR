function drgCaImAn_batch_analysis_pre_per_to_SVZ_per_thr_fsdzv2pp

close all
clear all

[FileName,PathName] = uigetfile({'*.mat'},'Select the .mat file with drgCaImAn_batch_pre_per_to_LDA_fsdz output');
load([PathName FileName])

classifier_names{1}='LDA';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='NN';
classifier_names{5}='Tree';

group_names{1}='AAAP';
group_names{2}='hod Fem';
group_names{3}='hod Fem(ky)';
group_names{4}='hod Mal';
group_names{5}='hod hod';
group_names{6}='homeodor_male_pcdh21';

%Trial window for calculating habituation per_trial
trial_window=9; %This must be an odd number
% total_trials=2*handles_out.handles.no_sp_sm_trials_to_use;
% per_trial_trialNo=1+(trial_window-1)/2:total_trials-(trial_window-1)/2;


%Compare the different algorithms
accuracy=[];
sh_accuracy=[];
bishop_accuracy=[];
bishop_sh_accuracy=[];
per_trial=[];

if isfield(handles_out.handles,'p_threshold')
    no_thr=length(handles_out.handles.p_thresholds);
else
    no_thr=length(handles_out.handles.p_thr_more_than);
end

for ii_thr=1:no_thr
    for ii_MLalgo=1:5
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp=[];
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp=[];
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp=[];
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        for grNo=1:5
            per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(grNo).correct_decisions=[];
            per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(grNo).ii_decisions=0;
        end
       
    end
end

ii_05=0;
ii_MLalgo_05=0;

accuracy_1=[];
sh_accuracy_1=[];
ii_1=0;

min_timepoints=200000;
for ii_out=1:length(handles_out.ii_out)
    if handles_out.ii_out(ii_out).handles.decoding_processed==1
        min_timepoints=min([min_timepoints length(handles_out.ii_out(ii_out).handles.accuracy_pp)]);
    end
end

all_group_nos=[];
ii_gr=0;
for ii_out=1:length(handles_out.ii_out)
    if handles_out.ii_out(ii_out).handles.decoding_processed==1
        

        ii_thr=find((handles_out.handles.p_thr_less_than==handles_out.ii_out(ii_out).p_thr_less_than)&(handles_out.handles.p_thr_more_than==handles_out.ii_out(ii_out).p_thr_more_than));
        ii_MLalgo=handles_out.ii_out(ii_out).MLalgo;
        
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy+1;
        ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy;
        
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs(ii_accuracy)=sum((handles_out.ii_out(ii_out).handles.p>=handles_out.ii_out(ii_out).p_thr_more_than)&...
            (handles_out.ii_out(ii_out).handles.p<=handles_out.ii_out(ii_out).p_thr_less_than));
        
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_accuracy_pp;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_sh_accuracy_pp;
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no(ii_accuracy)=handles_out.ii_out(ii_out).grNo;
        ii_gr=ii_gr+1;
        all_group_nos(ii_gr)=handles_out.ii_out(ii_out).grNo;
        
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp(ii_accuracy)=handles_out.ii_out(ii_out).handles.bishop_accuracy_pp;
        bii_accuracy=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp(bii_accuracy+1:bii_accuracy+length(handles_out.ii_out(ii_out).handles.bishop_sh_accuracy_pp))=handles_out.ii_out(ii_out).handles.bishop_sh_accuracy_pp;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy+length(handles_out.ii_out(ii_out).handles.bishop_sh_accuracy_pp);
        
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.accuracy_pp(1:min_timepoints);
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.sh_accuracy_pp(1:min_timepoints);
         
        %Now calculate the correct decisions for the habituation per_trial
        trialNo_sp=handles_out.ii_out(ii_out).handles.trialNo_sp;
        trialNo_sm=handles_out.ii_out(ii_out).handles.trialNo_sm;
        if isfield(handles_out.ii_out(ii_out).handles,'no_sp_sm_trials_to_use')
            total_trials=2*handles_out.ii_out(ii_out).handles.no_sp_sm_trials_to_use;
        else
            total_trials=2*handles_out.handles.no_sp_sm_trials_to_use; 
        end
        per_trial_trialNo=1+(trial_window-1)/2:total_trials-(trial_window-1)/2;
        
        for ii_trial_no=1:total_trials-trial_window+1
            these_trials=ii_trial_no:ii_trial_no+trial_window-1;
            these_winning_labels=[];
            these_training_decisions=[];  
            if isfield(handles_out.ii_out(ii_out).handles,'winning_label')
                for ttii=1:length(these_trials)
                    if sum(trialNo_sp==these_trials(ttii))
                        this_ii_sp=find(trialNo_sp==these_trials(ttii));
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.winning_label_pp(this_ii_sp);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions_pp(this_ii_sp);
                    else
                        if isfield(handles_out.ii_out(ii_out).handles,'no_sp_sm_trials_to_use')
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.ii_out(ii_out).handles.no_sp_sm_trials_to_use;
                        else
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.handles.no_sp_sm_trials_to_use;
                        end 
                        if this_ii_sm<=length(handles_out.ii_out(ii_out).handles.winning_label_pp)
                            these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.winning_label_pp(this_ii_sm);
                            these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions_pp(this_ii_sm);
                        end
                    end
                end
            else
                for ttii=1:length(these_trials) 
                    if sum(trialNo_sp==these_trials(ttii))
                        this_ii_sp=find(trialNo_sp==these_trials(ttii));
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.bishop_choice(this_ii_sp);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions_pp(this_ii_sp);
                    else
                        if isfield(handles_out.ii_out(ii_out).handles,'no_sp_sm_trials_to_use')
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.ii_out(ii_out).handles.no_sp_sm_trials_to_use;
                        else
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.handles.no_sp_sm_trials_to_use;
                        end
                       
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.bishop_choice(this_ii_sm);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions_pp(this_ii_sm);
                    end
                end
            end
            
            these_correct_decisions=these_winning_labels==these_training_decisions;
            per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(handles_out.ii_out(ii_out).grNo).correct_decisions(ii_trial_no,per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(handles_out.ii_out(ii_out).grNo).ii_decisions+1:...
                per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(handles_out.ii_out(ii_out).grNo).ii_decisions+length(these_correct_decisions))=these_correct_decisions;
        end
        per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(handles_out.ii_out(ii_out).grNo).ii_decisions=per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(handles_out.ii_out(ii_out).grNo).ii_decisions+length(these_correct_decisions);
    end
end

all_group_nos=unique(all_group_nos);
MLalgos=handles_out.ii_out(ii_out).handles.MLalgo;

%Plot bar graphs for accuracy and do the glm analysis
edges=[0:0.03:1];
rand_offset=0.5;
figNo=0;
for ii_MLalgo=MLalgos
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
    
    bar_offset=0;
    
    glm_acc=[];
    glm_acc_ii=0;
    
    for grNo=all_group_nos
        
        
        %Shuffled accuracy
        %I am adding try catch end because when I increase the cost of
        %misclassifiaciton as subset of the algorithms crash
       
        
        
       
        

        
        %Shuffled data for glm
        
        all_these_accs=[];
        for ii_thr=1:no_thr
            these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
            these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
            these_accs_per_group=these_accs(these_groups==grNo);
            all_these_accs=[all_these_accs these_accs_per_group];
            
            glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs_per_group))=these_accs_per_group;
            glm_acc.thr(glm_acc_ii+1:glm_acc_ii+length(these_accs_per_group))=ii_thr*ones(1,length(these_accs_per_group));
            glm_acc.group(glm_acc_ii+1:glm_acc_ii+length(these_accs_per_group))=grNo*ones(1,length(these_accs_per_group));
            glm_acc.shuffled_vs_original(glm_acc_ii+1:glm_acc_ii+length(these_accs_per_group))=zeros(1,length(these_accs_per_group));
            glm_acc_ii=glm_acc_ii+length(these_accs_per_group);
        end
        
        bar(bar_offset,mean(all_these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
        plot(bar_offset*ones(1,length(all_these_accs)),all_these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        if length(all_these_accs)>2
            CI = bootci(1000, {@mean, all_these_accs},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        end
        
        bar_offset=bar_offset+1;
        
        for ii_thr=1:no_thr
            %Accuracy
                these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
                these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
                these_accs=these_accs(these_groups==grNo);
                switch ii_thr
                    case 1
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])
                    case 2
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
                    case 3
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                end

                plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                if length(these_accs)>2
                    CI = bootci(1000, {@mean, these_accs},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                end
                
             glm_acc.data(glm_acc_ii+1:glm_acc_ii+length(these_accs))=these_accs;
            glm_acc.thr(glm_acc_ii+1:glm_acc_ii+length(these_accs))=ii_thr*ones(1,length(these_accs));
            glm_acc.group(glm_acc_ii+1:glm_acc_ii+length(these_accs))=grNo*ones(1,length(these_accs));
            glm_acc.shuffled_vs_original(glm_acc_ii+1:glm_acc_ii+length(these_accs))=ones(1,length(these_accs));
            glm_acc_ii=glm_acc_ii+length(these_accs);
            
            bar_offset=bar_offset+1;
        end
        
        bar_offset=bar_offset+2;
        
    end
   
    
    ylim([0 1.1])
    
    
    title(['Decoding accuracy for  ' classifier_names{ii_MLalgo}])

    
    plot([-1 28],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 28])
    
    xticklabels({group_names{1},group_names{2},group_names{3},group_names{4},group_names{5}})
    xticks([2 8 14 19 25])
    ylabel('Accuracy')

    
%         %Perform the glm
%     fprintf(1, ['glm for accuracy\n'])
% %     fprintf(fileID, ['glm for zPRP per mouse per odor pair for '  bandwidth_names{pacii} '\n']);
%     
%     tbl = table(glm_acc.data',glm_acc.thr',glm_acc.group',glm_acc.shuffled_vs_original',...
%         'VariableNames',{'accuracy','threshold','group','shuffled_vs_original'});
%     mdl = fitglm(tbl,'accuracy~threshold+group+shuffled_vs_original+threshold*group*shuffled_vs_original'...
%         ,'CategoricalVars',[2,3,4])
%     
    
%     txt = evalc('mdl');
%     txt=regexp(txt,'<strong>','split');
%     txt=cell2mat(txt);
%     txt=regexp(txt,'</strong>','split');
%     txt=cell2mat(txt);
%     
%     fprintf(fileID,'%s\n', txt);
    
%     %Do the ranksum/t-test
%     fprintf(1, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' hippocampus\n'])
%     fprintf(fileID, ['\n\nRanksum or t-test p values for zPRPe per mouse per odor pair for ' bandwidth_names{pacii} ' hippocampus\n']);
%     
%     [output_data] = drgMutiRanksumorTtest(input_data, fileID);
    
end



%Plot bar graphs for winner takes all accuracy
 for ii_MLalgo=MLalgos

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
    
    bar_offset=0;
    
    
    for grNo=all_group_nos
        

        %Shuffled accuracy
        %I am adding try catch end because when I increase the cost of
        %misclassifiaciton as subset of the algorithms crash
        these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
        try
            these_accs=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
            
            these_accs=these_accs(these_groups==grNo);
            
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            if length(these_accs)>2
                CI = bootci(1000, {@mean, these_accs},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            end
            
        catch
        end
        
        bar_offset=bar_offset+1;
        
        for ii_thr=1:no_thr
            %Accuracy
            try
                these_accs=bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
                these_accs=these_accs(these_groups==grNo);
                switch ii_thr
                    case 1
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])
                    case 2
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])
                    case 3
                        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                end

                plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                if length(these_accs)>2
                    CI = bootci(1000, {@mean, these_accs},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                end
                
            catch
            end
            
            bar_offset=bar_offset+1;
        end
        
        bar_offset=bar_offset+2;
        
    end
   
    
    ylim([0 1.1])
    
    
    title(['Winner takes all decoding accuracy for  ' classifier_names{ii_MLalgo}])

    
    plot([-1 28],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 28])
    
    xticklabels({group_names{1},group_names{2},group_names{3},group_names{4},group_names{5}})
    xticks([2 8 14 19 25])
    ylabel('Accuracy')
%  
    
 end 
       

 
 %Plot timecourses
 time=handles_out.ii_out(ii_out).handles.time_to_eventSm(1:min_timepoints);
 for grNo=all_group_nos
     figNo=figNo+1;
     try
         close(figNo)
     catch
     end
     
     hFig = figure(figNo);
     
     ax=gca;ax.LineWidth=3;
     
     set(hFig, 'units','normalized','position',[.1 .5 .5 .3])
     
     hold on
     
     for ii_thr=1:no_thr
         subplot(1,3,ii_thr)
         hold on
         
         ii_MLalgo=2;
         
         %Shuffled accuracy
         all_these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp_timecourse;
         these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
         ii_acc=0;
         these_accs=[];
         for ii_tr=1:length(these_groups)
             if these_groups(ii_tr)==grNo
                 ii_acc=ii_acc+1;
                 these_accs(ii_acc,:)=all_these_accs(ii_tr,:);
             end
         end
         
         pruned_these_accs=[];
         ii_pruned=0;
         
         for ii=1:size(these_accs,1)
             if ~isnan(these_accs(ii,:))
                 ii_pruned=ii_pruned+1;
                 pruned_these_accs(ii_pruned,:)=these_accs(ii,:);
             end
         end
         
         this_mean_accs=mean(pruned_these_accs,1);
         
         if size(pruned_these_accs,1)>2
             CI=[];
             CI = bootci(1000, {@mean, pruned_these_accs})';
             CI(:,1)= this_mean_accs'-CI(:,1);
             CI(:,2)=CI(:,2)- this_mean_accs';
             
             
             [hlCR, hpCR] = boundedline(time',this_mean_accs', CI, 'k');
         else
             plot(time',this_mean_accs', 'k');
         end
         
         %Accuracy
         all_these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp_timecourse;
         these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
         ii_acc=0;
         these_accs=[];
         for ii_tr=1:length(these_groups)
             if these_groups(ii_tr)==grNo
                 ii_acc=ii_acc+1;
                 these_accs(ii_acc,:)=all_these_accs(ii_tr,:);
             end
         end
         
         
         pruned_these_accs=[];
         ii_pruned=0;
         
         for ii=1:size(these_accs,1)
             if ~isnan(these_accs(ii,:))
                 ii_pruned=ii_pruned+1;
                 pruned_these_accs(ii_pruned,:)=these_accs(ii,:);
             end
         end
         
         this_mean_accs=mean(pruned_these_accs,1);
         
         
         if size(pruned_these_accs,1)>2
             CI=[];
             CI = bootci(1000, {@mean, pruned_these_accs})';
             CI(:,1)= this_mean_accs'-CI(:,1);
             CI(:,2)=CI(:,2)- this_mean_accs';
             
             
             [hlCR, hpCR] = boundedline(time',this_mean_accs', CI, 'r');
         else
             plot(time',this_mean_accs', 'r');
         end
         
         ylim([0 1.1])
         
         
         
         
         title(['p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
             ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
         
         
         plot([min(time) max(time)],[0.5 0.5],'-k','LineWidth',2)
         
     end
     
       sgtitle(['Decoding accuracy for ' group_names{grNo}])
 end
 
 
%Now see whether the accuracy decreases as a function of trial number
C = {'k','b','r','g','m',[.5 .6 .7],[.8 .2 .6]};
for ii_thr=1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
    
    bar_offset=0;
    
    
        title(['Decoding accuracy for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
   
    
    p=[];
    
    for ii_MLalgo=MLalgos
        for grNo=all_group_nos
            this_per_trial_accuracy=zeros(1,length(per_trial_trialNo));
            for jj_trial=1:length(per_trial_trialNo)
                these_correct_decisions=per_trial.thr(ii_thr).MLalgo(ii_MLalgo).groupNo(grNo).correct_decisions(jj_trial,:);
                this_per_trial_accuracy(jj_trial)=sum(these_correct_decisions)/length(these_correct_decisions);
            end
            plot(per_trial_trialNo,this_per_trial_accuracy,'o','MarkerFaceColor',C{grNo})
        end
    end
    
    
    ylim([0.1 1.1])
    xlim([0 per_trial_trialNo(end)+1])
    
    
    
   
        title(['Decoding accuracy for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
   
    
    plot([0 per_trial_trialNo(end)],[0.5 0.5],'-k','LineWidth',2)
    
    ylabel('Accuracy')
    xlabel('Trial number')
    

    
end

 

%Plot bar graphs for accuracy per group
for ii_thr=1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
    
    bar_offset=0;
    
       if isfield(handles_out.handles,'p_threshold')
        fprintf(1,['t test wta p values for p<' num2str(handles_out.handles.p_thresholds(ii_thr)) '\n'])
    else
        fprintf(1,['t test wta p values for p from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
    end
  
    p=[];
    
    for ii_MLalgo=MLalgos
        
        ii_stats=0;
        accuracy_stats=[];
        
        %Shuffled accuracy all groups
        these_accs=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
        
        if length(these_accs)>2
            CI = bootci(1000, {@mean, these_accs},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        end
        
        plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        %         [mean_out, CIout]=drgViolinPoint(these_accs,edges,bar_offset,rand_offset,'k','k',3);
        these_accs_sh=these_accs;
        
        ii_stats=ii_stats+1;
        accuracy_stats(ii_stats).data=these_accs;
        accuracy_stats(ii_stats).description='Shuffled';
        
        
        %Accuracy per group
        these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
        these_accs=bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_pp;
        
        for groupNo=1:6
            bar_offset=bar_offset+1;
            this_group_accs=these_accs(these_groups==groupNo);
            
            if ~isempty(this_group_accs)
                bar(bar_offset,mean(this_group_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                
                if length(this_group_accs)>2
                    CI = bootci(1000, {@mean, this_group_accs},'type','cper');
                    plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
                end
                
                plot(bar_offset*ones(1,length(this_group_accs)),this_group_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
                
                ii_stats=ii_stats+1;
                accuracy_stats(ii_stats).data=this_group_accs;
                accuracy_stats(ii_stats).description=group_names{groupNo};
            end
            
            if (ii_MLalgo==1)
                if isfield(handles_out.handles,'p_threshold')
                    fprintf(1,['Number of experiments for ' group_names{groupNo} ' for p<' num2str(handles_out.handles.p_thresholds(ii_thr)) ' is ' num2str(length(this_group_accs)) '\n\n'])
                else
                    fprintf(1,['Number of experiments for ' group_names{groupNo} ' for p from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
                        ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr))...
                        ' is ' num2str(length(this_group_accs)) '\n\n'])
                end
                
                
            end
        end
        %         [h p(ii_MLalgo)]=ttest2(these_accs_sh,these_accs);
        %         fprintf(1,[classifier_names{ii_MLalgo} ' p=' num2str(p(ii_MLalgo)) '\n'])
        
        bar_offset=bar_offset+2;
        
        %Do ranksum/t test
        if isfield(handles_out.handles,'p_threshold')
            fprintf(1, ['\n\nRanksum or t-test p values for accuracy for ' classifier_names{ii_MLalgo} ' for p<' num2str(handles_out.handles.p_thresholds(ii_thr)) '\n'])
        else
            fprintf(1, ['\n\nRanksum or t-test p values for accuracy for ' classifier_names{ii_MLalgo} ' for  p from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
                ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) '\n'])
        end
         
        [output_data] = drgMutiRanksumorTtest(accuracy_stats);
        fprintf(1, '\n\n')
        
    end
    
    
    
    ylim([0 1.1])
    
    if isfield(handles_out.handles,'p_threshold')
        title(['Decoding accuracy wta for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    else
        title(['Decoding accuracy wta for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
    end
    
    
    plot([-1 6],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 6])
    
    xticks([0:5])
    expression='xticklabels({''Shuff'',';
    for grNo=1:5
        expression=[expression '''' group_names{grNo} '''' ];
        if grNo~=5
            expression=[expression ', '];
        end
    end
    expression=[expression '})'];
    eval(expression)

end



%Plot bar graphs for no of ROIs per group
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

ax=gca;ax.LineWidth=3;

set(hFig, 'units','normalized','position',[.1 .5 .5 .4])

hold on

bar_offset=0;



for ii_thr=1
    
    ii_MLalgo=2;
    
    
    fprintf(1, ['\n\nNumber of mice per group\n'])
       
    
    %Accuracy per group
    these_groups=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no;
    these_ROIs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs;
    
    for groupNo=1:5
        bar_offset=bar_offset+1;
        this_group_ROIs=these_ROIs(these_groups==groupNo);
        
        if ~isempty(this_group_ROIs)
            bar(bar_offset,mean(this_group_ROIs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            
            if length(this_group_ROIs)>2
                CI = bootci(1000, {@mean, this_group_ROIs},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            end
            plot(bar_offset*ones(1,length(this_group_ROIs)),this_group_ROIs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',7)
        end
         fprintf(1, [group_names{groupNo} ' ' num2str(sum(these_groups==groupNo)) '\n'])
       
    end
    %         [h p(ii_MLalgo)]=ttest2(these_accs_sh,these_accs);
    %         fprintf(1,[classifier_names{ii_MLalgo} ' p=' num2str(p(ii_MLalgo)) '\n'])
    bar_offset=bar_offset+2;
    
end

title(['Number of ROIs used for decoding'])


xlim([-1 6])

xticks([0:5])
expression='xticklabels({''Shuff'',';
for grNo=1:5
    expression=[expression '''' group_names{grNo} '''' ];
    if grNo~=5
        expression=[expression ', '];
    end
end
expression=[expression '})'];
eval(expression)

ylabel('No ROIs')



pffft=1
