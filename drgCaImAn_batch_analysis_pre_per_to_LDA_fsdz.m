function drgCaImAn_batch_analysis_pre_per_to_LDA_fsdz

[FileName,PathName] = uigetfile({'*.mat'},'Select the .mat file with drgCaImAn_batch_pre_per_to_LDA_fsdz output');
load([PathName FileName])

classifier_names{1}='LDA';
classifier_names{2}='SVM';
classifier_names{3}='Bayes';
classifier_names{4}='NN';
classifier_names{5}='Tree';

%Compare the different algorithms
accuracy=[];
sh_accuracy=[];
  
for ii_thr=1:length(handles_out.handles.p_thresholds)
    for ii_MLalgo=1:5
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy=[];
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=0;
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
        min_timepoints=min([min_timepoints length(handles_out.ii_out(ii_out).handles.accuracy)])
    end
end

for ii_out=1:length(handles_out.ii_out)
    if handles_out.ii_out(ii_out).handles.decoding_processed==1
        ii_thr=find(handles_out.handles.p_thresholds==handles_out.ii_out(ii_out).p_threshold);
        ii_MLalgo=handles_out.ii_out(ii_out).MLalgo;
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy+1;
        ii_accuracy=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).ii_accuracy;
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_accuracy;
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy(ii_accuracy)=handles_out.ii_out(ii_out).handles.mean_sh_accuracy;
        accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.accuracy(1:min_timepoints);
        sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse(ii_accuracy,:)=handles_out.ii_out(ii_out).handles.sh_accuracy(1:min_timepoints);
    end
end

%Plot bar graphs
figNo=0;
for ii_thr=1:length(handles_out.handles.p_thresholds)
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
    
    for ii_MLalgo=1:5
        
        %Shuffled accuracy
        these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
        CI = bootci(1000, {@mean, these_accs},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        bar_offset=bar_offset+1;
        
        %Accuracy
        these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy;
        bar(bar_offset,mean(these_accs),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
        CI = bootci(1000, {@mean, these_accs},'type','cper');
        plot([bar_offset bar_offset],CI,'-k','LineWidth',3)
        plot(bar_offset*ones(1,length(these_accs)),these_accs,'o','MarkerFaceColor', [0.7 0.7 0.7],'MarkerEdgeColor',[0 0 0],'MarkerSize',5)
        
        bar_offset=bar_offset+2;
    end
    ylim([0.4 0.7])
    
    title(['Decoding accuracy for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    
    plot([-1 14],[0.5 0.5],'-k','LineWidth',2)
    xlim([-1 14])
    
    xticks([0.5 3.5 6.5 9.5 12.5])
    xticklabels({'LDA', 'SVM','Bayes', 'NN', 'Tree'})
    
end


%Plot timecourses
time=handles_out.ii_out(ii_out).handles.time_to_eventSm(1:min_timepoints);
for ii_thr=1:length(handles_out.handles.p_thresholds)
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    ax=gca;ax.LineWidth=3;
    
    set(hFig, 'units','normalized','position',[.1 .5 .5 .4])
    
    hold on
   
    
    ii_MLalgo=3;
    
    %Shuffled accuracy
    these_accs=sh_accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse;
    this_mean_accs=mean(these_accs,1);
    
    CI=[];
    CI = bootci(1000, {@mean, these_accs})';
    CI(:,1)= this_mean_accs'-CI(:,1);
    CI(:,2)=CI(:,2)- this_mean_accs';
    
    
    [hlCR, hpCR] = boundedline(time',this_mean_accs', CI, 'k');
    
      %Accuracy
    these_accs=accuracy.thr(ii_thr).MLalgo(ii_MLalgo).accuracy_timecourse;
    this_mean_accs=mean(these_accs,1);
    
    CI=[];
    CI = bootci(1000, {@mean, these_accs})';
    CI(:,1)= this_mean_accs'-CI(:,1);
    CI(:,2)=CI(:,2)- this_mean_accs';
    
    
    [hlCR, hpCR] = boundedline(time',this_mean_accs', CI, 'r');
        
    ylim([0.2 0.9])
    
    title(['Decoding accuracy for p threshold ' num2str(handles_out.handles.p_thresholds(ii_thr))])
    
    plot([min(time) max(time)],[0.5 0.5],'-k','LineWidth',2)

    
    
end

pffft=1
