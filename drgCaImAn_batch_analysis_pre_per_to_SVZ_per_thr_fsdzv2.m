function drgCaImAn_batch_analysis_pre_per_to_SVZ_per_thr_fsdzv2

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
trial_window=5; %This must be an odd number
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
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no=[];
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
%         per_trial.thr(ii_thr).MLalgo(ii_MLalgo).correct_decisions=zeros(total_trials-trial_window-1,handles_out.handles.no_files*trial_window);
        per_trial.thr(ii_thr).MLalgo(ii_MLalgo).correct_decisions=[];
        per_trial.thr(ii_thr).MLalgo(ii_MLalgo).ii_decisions=0;
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
        min_timepoints=min([min_timepoints length(handles_out.ii_out(ii_out).handles.accuracy)]);
    end
end

for ii_out=1:length(handles_out.ii_out)
    if handles_out.ii_out(ii_out).handles.decoding_processed==1
        
        
        if isfield(handles_out.handles,'p_threshold')
            ii_thr=find(handles_out.handles.p_thresholds==handles_out.ii_out(ii_out).p_threshold);
        else
            ii_thr_more=find(handles_out.handles.p_thr_more_than==handles_out.ii_out(ii_out).p_thr_more_than);
            ii_thr_less=find(handles_out.handles.p_thr_less_than==handles_out.ii_out(ii_out).p_thr_less_than);
            
            for ii_th_m=1:length(ii_thr_more)
                for ii_th_l=1:length(ii_thr_less)
                    if (ii_thr_more(ii_th_m)==ii_thr_less(ii_th_l))
                        ii_thr=ii_thr_more(ii_th_m);
                    end
                end
            end
        end
        
        
        ii_MLalgo=handles_out.ii_out(ii_out).MLalgo;

        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy+1;
        ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy;
        if ~isfield(handles_out.handles,'p_threshold')
            accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs(ii_accuracy)=sum((handles_out.ii_out(ii_out).handles.p>=handles_out.ii_out(ii_out).p_thr_more_than)&...
                (handles_out.ii_out(ii_out).handles.p<=handles_out.ii_out(ii_out).p_thr_less_than));
        else
            accuracy.thr(ii_thr).MLalgo(ii_MLalgo).no_ROIs(ii_accuracy)=sum(handles_out.ii_out(ii_out).handles.p<handles_out.ii_out(ii_out).handles.p_threshold);
        end
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_accuracy;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_sh_accuracy;
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).gr_no(ii_accuracy)=handles_out.ii_out(ii_out).grNo;
        
        bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(ii_accuracy)=handles_out.ii_out(ii_out).handles.bishop_accuracy;
        bii_accuracy=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(bii_accuracy+1:bii_accuracy+length(handles_out.ii_out(ii_out).handles.bishop_sh_accuracy))=handles_out.ii_out(ii_out).handles.bishop_sh_accuracy;
        bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy+length(handles_out.ii_out(ii_out).handles.bishop_sh_accuracy);
        
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.accuracy(1:min_timepoints);
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.sh_accuracy(1:min_timepoints);
        
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
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.winning_label(this_ii_sp);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions(this_ii_sp);
                    else
                        if isfield(handles_out.ii_out(ii_out).handles,'no_sp_sm_trials_to_use')
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.ii_out(ii_out).handles.no_sp_sm_trials_to_use;
                        else
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.handles.no_sp_sm_trials_to_use;
                        end
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.winning_label(this_ii_sm);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions(this_ii_sm);
                    end
                end
            else
                for ttii=1:length(these_trials) 
                    if sum(trialNo_sp==these_trials(ttii))
                        this_ii_sp=find(trialNo_sp==these_trials(ttii));
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.bishop_choice(this_ii_sp);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions(this_ii_sp);
                    else
                        if isfield(handles_out.ii_out(ii_out).handles,'no_sp_sm_trials_to_use')
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.ii_out(ii_out).handles.no_sp_sm_trials_to_use;
                        else
                            this_ii_sm=find(trialNo_sm==these_trials(ttii))+handles_out.handles.no_sp_sm_trials_to_use;
                        end
                       
                        these_winning_labels(ttii)=handles_out.ii_out(ii_out).handles.bishop_choice(this_ii_sm);
                        these_training_decisions(ttii)=handles_out.ii_out(ii_out).handles.training_decisions(this_ii_sm);
                    end
                end
            end
            
            these_correct_decisions=these_winning_labels==these_training_decisions;
            per_trial.thr(ii_thr).MLalgo(ii_MLalgo).correct_decisions(ii_trial_no,per_trial.thr(ii_thr).MLalgo(ii_MLalgo).ii_decisions+1:...
                per_trial.thr(ii_thr).MLalgo(ii_MLalgo).ii_decisions+trial_window)=these_correct_decisions;
        end
        per_trial.thr(ii_thr).MLalgo(ii_MLalgo).ii_decisions=per_trial.thr(ii_thr).MLalgo(ii_MLalgo).ii_decisions+trial_window;
    end
end

MLalgos=handles_out.ii_out(ii_out).handles.MLalgo;

%Plot bar graphs for accuracy
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
    
    
    fprintf(1,['t test p values for ' classifier_names{ii_MLalgo} '\n'])
    
    
    p=[];
    mean_acc_per_th=[];
    for ii_thr=1:no_thr
         
        %Shuffled accuracy
        %I am adding try catch end because when I increase the cost of
        %misclassifiaciton as subset of the algorithms crash
        try
            these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
            plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            if length(these_accs)>2
                CI = bootci(1000, {@mean, these_accs},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            end
            these_accs_sh=these_accs;
        catch
        end
        
        bar_offset=bar_offset+1;
        
        %Accuracy
        try
            these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
            bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
            plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
            if length(these_accs)>2
                CI = bootci(1000, {@mean, these_accs},'type','cper');
                plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
            end
           
        catch
        end
         
        bar_offset=bar_offset+2;
        
        [h p(ii_thr)]=ttest(these_accs_sh,these_accs);
        
        if isfield(handles_out.handles,'p_threshold')
            fprintf(1,['For p threshold p<' num2str(handles_out.handles.p_thresholds(ii_thr)) ' p=' num2str(p(ii_thr)) '\n'])
        else
            fprintf(1,['t test p values for p from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
                ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ' p=' num2str(p(ii_thr)) '\n' ])
        end
        
        
    end
    
    fprintf(1,['p FDR =' num2str(drsFDRpval(p)) '\n\n'])
    
    ylim([0 1.1])
    
    
    title(['Decoding accuracy for  ' classifier_names{ii_MLalgo}])

    
    plot([-1 no_thr*3],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 no_thr*3])
    
    xticks([1:3:(no_thr-1)*3+1])
    expression='xticklabels({';
    for ii_thr=1:no_thr
        expression=[expression '''' num2str(handles_out.handles.p_thr_less_than(ii_thr)) '''' ];
        if ii_thr~=no_thr
            expression=[expression ', '];
        end
    end
    expression=[expression '})'];
    eval(expression)
    
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
    
    
     fprintf(1,['t test wta p values for ' classifier_names{ii_MLalgo} '\n'])
     

    
    p=[];
    
    for ii_thr=1:no_thr
        
        %Shuffled accuracy
        these_accs=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
        plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        if length(these_accs)>2
            CI = bootci(1000, {@mean, these_accs},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        end
        %         [mean_out, CIout]=drgViolinPoint(these_accs,edges,bar_offset,rand_offset,'k','k',3);
        these_accs_sh=these_accs;
        
        bar_offset=bar_offset+1;
        
        %Accuracy
        these_accs=bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
        plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        if length(these_accs)>2
            CI = bootci(1000, {@mean, these_accs},'type','cper');
            plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        end
         mean_acc_per_th(ii_thr)=mean(these_accs);
         
        bar_offset=bar_offset+2;
        
        [h p(ii_thr)]=ttest2(these_accs_sh,these_accs);
        
        if isfield(handles_out.handles,'p_threshold')
            fprintf(1,['t test wta for p threshold p<' num2str(handles_out.handles.p_thresholds(ii_thr)) ' p=' num2str(p(ii_thr)) '\n'])
        else
            fprintf(1,['t test wta p values for p from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
                ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ' p=' num2str(p(ii_thr)) '\n' ])
        end
        
    end
    
    fprintf(1,['p FDR =' num2str(drsFDRpval(p)) '\n\n'])
    
    ylim([0 1.1])
    
    title(['Decoding accuracy for  ' classifier_names{ii_MLalgo}])

    
    plot([-1 no_thr*3],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 no_thr*3])
    
    xticks([1:3:(no_thr-1)*3+1])
    expression='xticklabels({';
    for ii_thr=1:no_thr
        expression=[expression '''' num2str(handles_out.handles.p_thr_less_than(ii_thr)) '''' ];
        if ii_thr~=no_thr
            expression=[expression ', '];
        end
    end
    expression=[expression '})'];
    eval(expression)
    
 end 
       
 [max_acc ii_max_acc]=max(mean_acc_per_th)
 
%Plot timecourses
time=handles_out.ii_out(ii_out).handles.time_to_eventSm(1:min_timepoints);
for ii_thr=ii_max_acc
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
    
     
    ii_MLalgo=2;
    
    %Shuffled accuracy
    these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse;
    

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
    these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse;

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
    
    
    
    if isfield(handles_out.handles,'p_threshold')
        title(['Decoding accuracy for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    else
        title(['Decoding accuracy for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
    end
    
    plot([min(time) max(time)],[0.5 0.5],'-k','LineWidth',2)
    
    
end

%Now see whether the accuracy decreases as a function of trial number
%Plot bar graphs for bishop
C = {'k','b','r','g','m',[.5 .6 .7],[.8 .2 .6]};
for ii_thr=ii_max_acc
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
        title(['Decoding accuracy for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    else
        title(['Decoding accuracy for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
    end
    
    p=[];
    
    for ii_MLalgo=MLalgos
        
        this_per_trial_accuracy=zeros(1,length(per_trial_trialNo));
        for jj_trial=1:length(per_trial_trialNo)
            these_correct_decisions=per_trial.thr(ii_thr).MLalgo(ii_MLalgo).correct_decisions(jj_trial,:);
            this_per_trial_accuracy(jj_trial)=sum(these_correct_decisions)/length(these_correct_decisions);
        end
        plot(per_trial_trialNo,this_per_trial_accuracy,'o','MarkerFaceColor',C{ii_MLalgo})
    end
    
    
    ylim([0.1 1.1])
    xlim([0 per_trial_trialNo(end)+1])
    
    
    
    if isfield(handles_out.handles,'p_threshold')
        title(['Decoding accuracy for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    else
        title(['Decoding accuracy for p threshold from ' num2str(handles_out.handles.p_thr_more_than(ii_thr))...
            ' to ' num2str(handles_out.handles.p_thr_less_than(ii_thr)) ])
    end
    
    plot([0 per_trial_trialNo(end)],[0.5 0.5],'-k','LineWidth',2)
    
    ylabel('Accuracy')
    xlabel('Trial number')
    
%     text(2,0.4,'LDA','Color',C{1})
%     text(2,0.37,'SVM','Color',C{2})
%     text(2,0.34,'Bayes','Color',C{3})
%     text(2,0.31,'NN','Color',C{4})
%     text(2,0.28,'Tree','Color',C{5})
    
end

 

%Plot bar graphs for accuracy per group
for ii_thr=ii_max_acc
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
        these_accs=bishop_sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
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
        these_accs=bishop_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
        
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



for ii_thr=ii_max_acc
    
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
