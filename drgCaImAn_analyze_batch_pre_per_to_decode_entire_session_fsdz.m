%drgCaImAn_analyze_batch_pre_per_to_decode_entire_session_fsdz
close all
clear all

%Load file
[pre_perFileName,pre_perPathName] = uigetfile({'*.mat'},'Select the file with the entire session batch analysis');

load([pre_perPathName pre_perFileName])

figNo=0;

%Now do the analysis of whether the ROIs that contribute to memory are the same ROIs as those decoding the odor
MLalgo=6;  %This is glm

for grNo=1:unique(handles_out.handles.group)
    
    %Find the ROIs that contribute to the glm
    ROIs=[];
    ii_grno=0;
    for ii_out=1:length(handles_out.ii_out)
        if handles_out.ii_out(ii_out).grNo==grNo
            ii_grno=ii_grno+1;
            k_fold=handles_out.ii_out(ii_out).handles_out.k_fold;
            ii_post_shift=find(handles_out.handles.post_shift==handles_out.ii_out(ii_out).handles_out.post_shift);
            for ii_k_fold=1:k_fold
                thisMdl=handles_out.ii_out(ii_out).handles_out.MLalgo(MLalgo).models(ii_k_fold).Mdl;
                pValues = thisMdl.Coefficients.pValue;
                pValues=pValues(2:end); %The first p value is for the intercept
                ROIs.fileNo(handles_out.ii_out(ii_out).fileNo).post_shift(ii_post_shift).k_fold(ii_k_fold).included_ROIs=find(pValues<=0.05);
                ROIs.fileNo(handles_out.ii_out(ii_out).fileNo).post_shift(ii_post_shift).k_fold(ii_k_fold).noROIs=length(pValues);
            end
            ROIs.fileNo(handles_out.ii_out(ii_out).fileNo).post_shift(ii_post_shift).accuracy_wta=handles_out.ii_out(ii_out).handles_out.MLalgo(MLalgo).accuracy_tr_wta;
        end
    end
    
    if ii_grno>2
        %For each file find the ROIs in common with post_shift 1
        fraction_ROIs_included_glm=zeros(length(handles_out.handles.post_shift),handles_out.handles.no_files);
        accuracy_wta=zeros(length(handles_out.handles.post_shift),handles_out.handles.no_files);
        fraction_ROIs_in_common=zeros(length(handles_out.handles.post_shift),length(handles_out.handles.post_shift),handles_out.handles.no_files);
        fraction_ROIs_in_common_vs_dt=[];
        all_dts=handles_out.handles.post_shift; %Note that this assumes that these are spaced by the same dt
        fraction_ROIs_in_common_vs_dt_ii=zeros(1,length(handles_out.handles.post_shift));
        for fileNo=1:handles_out.handles.no_files
            for ii_post_shift1=1:length(handles_out.handles.post_shift)
                these_f_included=zeros(1,k_fold);
                for ii_k_fold=1:k_fold
                    these_f_included(k_fold)=length(ROIs.fileNo(fileNo).post_shift(ii_post_shift1).k_fold(ii_k_fold).included_ROIs)/ROIs.fileNo(fileNo).post_shift(ii_post_shift1).k_fold(ii_k_fold).noROIs;
                end
                fraction_ROIs_included_glm(ii_post_shift1,fileNo)=mean(these_f_included);
                accuracy_wta(ii_post_shift1,fileNo)=ROIs.fileNo(fileNo).post_shift(ii_post_shift1).accuracy_wta;
                for ii_post_shift2=1:length(handles_out.handles.post_shift)
                    these_f_in_common=zeros(1,k_fold);
                    for ii_k_fold=1:k_fold
                        [val,pos]=intersect(ROIs.fileNo(fileNo).post_shift(ii_post_shift1).k_fold(ii_k_fold).included_ROIs,ROIs.fileNo(fileNo).post_shift(ii_post_shift2).k_fold(ii_k_fold).included_ROIs);
                        these_f_in_common(ii_k_fold)=length(val)/length(ROIs.fileNo(fileNo).post_shift(ii_post_shift1).k_fold(ii_k_fold).included_ROIs);
                    end
                    fraction_ROIs_in_common(ii_post_shift1,ii_post_shift2,fileNo)=mean(these_f_in_common);
                    this_dt=abs(handles_out.handles.post_shift(ii_post_shift1)-handles_out.handles.post_shift(ii_post_shift2));
                    this_dt_ii=find(this_dt==handles_out.handles.post_shift);
                    fraction_ROIs_in_common_vs_dt_ii(1,this_dt_ii)=fraction_ROIs_in_common_vs_dt_ii(1,this_dt_ii)+1;
                    fraction_ROIs_in_common_vs_dt(this_dt_ii,fraction_ROIs_in_common_vs_dt_ii(1,this_dt_ii))=fraction_ROIs_in_common(ii_post_shift1,ii_post_shift2,fileNo);
                end
            end
        end
        
        %Now plot the fraction of ROIs in common vs distance
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        set(hFig, 'units','normalized','position',[.05 .1 .4 .3])
        
        hold on
        
        edges=[0:0.05:1.2];
        rand_offset=2.5;
        
        mean_f_common=ones(1,length(handles_out.handles.post_shift));
        
        plot(0,1,'ok')
        for ii=2:length(handles_out.handles.post_shift)
            
            this_f_common=fraction_ROIs_in_common_vs_dt(ii,1:fraction_ROIs_in_common_vs_dt_ii(ii));
            this_dt=handles_out.handles.post_shift(ii);
            mean_f_common(ii)=mean(this_f_common);
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(this_f_common...
                ,edges,this_dt,rand_offset,'k','k',4);
        end
        plot(handles_out.handles.post_shift,mean_f_common,'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2)
        title(['Fraction of common ROIs included in glm for ' handles_out.handles.group_names{grNo}])
        xlabel('Difference in time (sec)')
        ylabel('Fraction')
        
        %Now plot the fraction of ROIs included in glm vs distance
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        set(hFig, 'units','normalized','position',[.05 .1 .4 .3])
        
        hold on
        
        edges=[0:0.05:1.2];
        rand_offset=2.5;
        
        mean_f_included=zeros(1,length(handles_out.handles.post_shift));
        
        
        for ii=1:length(handles_out.handles.post_shift)
            this_f_included=fraction_ROIs_included_glm(ii,:);
            this_dt=handles_out.handles.post_shift(ii)+2.5;
            mean_f_included(ii)=mean(this_f_included);
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(this_f_included...
                ,edges,this_dt,rand_offset,'k','k',4);
        end
        plot(handles_out.handles.post_shift+2.5,mean_f_included,'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2)
        title(['Fraction of ROIs included in glm for ' handles_out.handles.group_names{grNo}])
        xlabel('Time after odor on (sec)')
        ylabel('Fraction')
        
        %Now plot wta
        %Now plot the fraction of ROIs included in glm vs distance
        figNo=figNo+1;
        try
            close(figNo)
        catch
        end
        
        hFig = figure(figNo);
        
        set(hFig, 'units','normalized','position',[.05 .1 .4 .3])
        
        hold on
        
        edges=[0:0.05:1.2];
        rand_offset=2.5;
        
        mean_accuracy_wta=zeros(1,length(handles_out.handles.post_shift));
        
        
        for ii=1:length(handles_out.handles.post_shift)
            this_accuracy_wta=accuracy_wta(ii,:);
            this_dt=handles_out.handles.post_shift(ii)+2.5;
            mean_accuracy_wta(ii)=mean(this_accuracy_wta);
            
            %Violin plot
            [mean_out, CIout]=drgViolinPoint(this_accuracy_wta...
                ,edges,this_dt,rand_offset,'k','k',4);
        end
        plot(handles_out.handles.post_shift+2.5,mean_accuracy_wta,'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2)
        title(['Accuracy wta glm for ' handles_out.handles.group_names{grNo}])
        xlabel('Time after odor on (sec)')
        ylabel('Accuracy')
         
        pffft=1;
    end
end