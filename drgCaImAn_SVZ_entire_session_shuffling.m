function handles_out=drgCaImAn_SVZ_entire_session_shuffling(handles_choices)
%This program trains several decoding algorithms with the post odorant and then determines what happens throughout the entire timecouse
%The user enters the choices entered under exist('handles_choices')==0
%
% processing_algorithm= 1 and 2 were used for troublehsooting and do not
% produce reliable results because of overtraining, use
% processing_algoritm=3, that was vetted for our manuscript
%
%
% the input is a pre_per file version 2
if exist('handles_choices')==0
    clear all
    close all
    
    %Load file
    [pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');
    
    processing_algorithm=3; %Use 3
    k_fold=5; %Only used for processing_algorithm=2,
    post_time=5; %The decoding model will be trained with all points in post_time sec interval starting post_shift secs after odor on
    post_shift=0; %Set to 0 if you want to train with odor on points
    pre_time=5; %Used to calculate the decoding accuracy pre_time sec before post_shift
    MLalgo_to_use=[6]; %Vector with the decoding algorithms you want to use
    
    %     switch MLalgo
    %         case 1
    %             fitcdiscr is LDA
    %         case 2
    %             fitcsvm support vector machine
    %         case 3
    %            fitcnb naive Bayes
    %         case 4
    %             fitcnet neural network
    %         case 5
    %             fitctree binary decision tree
    %         case 6
    %             fitglm binomial generalized linear model
    %     end
    
    p_threshold=0.3; %This limits the ROIs used in the decoding model to those whose p value in a ttest for the two odors in <=p_threshold
    dt_p_threshold=20; %Time to be used after the odor on for the p_threshold t_test
    show_figures=1; %Show the figures
else
    pre_perPathName=handles_choices.pre_per_PathName;
    pre_perFileName=handles_choices.pre_per_FileName;
    processing_algorithm=handles_choices.processing_algorithm;
    post_time=handles_choices.post_time;
    k_fold=handles_choices.k_fold;
    post_shift=handles_choices.post_shift;
    MLalgo_to_use=handles_choices.MLalgo_to_use;
    pre_time=handles_choices.pre_time;
    p_threshold=handles_choices.p_threshold;
    dt_p_threshold=handles_choices.dt_p_threshold;
    show_figures=handles_choices.show_figures;
end

warning('off')

dt_span=40; %Seconds shown before and after odor on in the figures
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=10; %Number of shuffles for per trial shuffling

load([pre_perPathName pre_perFileName])
fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
fprintf(1, 'post_time = %d, k_fold= %d, post_shift= %d\n',post_time,k_fold,post_shift);


classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';

handles_out.pre_perFileName=pre_perFileName;
handles_out.pre_perPathName=pre_perPathName;
handles_out.post_time=post_time;
handles_out.k_fold=k_fold;
handles_out.post_shift=post_shift;
handles_out.pre_time=pre_time;
handles_out.MLalgo_to_use=MLalgo_to_use;

figNo=0;

%time has the time for the dF/F traces(ROI,time)
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on
    
    % Determine the y spacing of the traces
    y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));
    
    %Plot the traces and do z normalization
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
        end
    end
    
    for trNo=1:no_traces
        % for trNo=1:20
        plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        traces(trNo,:)=traces(trNo,:)/std(traces(trNo,:));
    end
    
    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
    title(['All dFF timecourses ' num2str(size(traces,1)) ' ROIs'])
end

%epochs is a vector of the length of time that gives information on
%behavior
% 1=Final Valve
% 6=Hit (on for the duration of odor on)
% 7=Miss
% 8=FA
% 9=CR

%For example Hit||Miss shows S+ odor application times (red)
%and FA||CR gives S- (blue)

%Post points
Nall=size(traces,1);
dt=time(2)-time(1);
ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_post=floor(post_time/dt);
no_points_post_shift=floor(post_shift/dt);
no_points_pre=floor(pre_time/dt);
measurements_post=[];
measurements_pre=[];
epochs_sp_post=zeros(1,length(time));
epochs_sp_pre=zeros(1,length(time));
epochs_sm_post=zeros(1,length(time));
epochs_sm_pre=zeros(1,length(time));
which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));

%Do both S+ and S-
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;
trial_no=0;
ii_sp_post=0;
ii_sm_post=0;
ii_which_model=0;
dt_post_which_model=floor(20/dt); %Points that model will be used beyond the training period
ii_span=ceil(dt_span/dt);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI=[];
dFFs_sm_per_trial_per_ROI=[];

%training_decisions is 1 if S+ and 2 if S-
while (at_end==0)
    next_ii_sp=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    next_ii_sm=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    
    if (isempty(next_ii_sp))&(isempty(next_ii_sm))
        at_end=1;
    else
        
        if isempty(next_ii_sm)
            next_ii=next_ii_sp;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                ii_sp_post=ii_sp_post+1;
                dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=1;
                which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+no_points_post;
                ii_post=ii_post+no_points_post;
                ii_pre=ii_pre+no_points_pre;
            else
                at_end=1;
            end
        end
        
        if isempty(next_ii_sp)
            next_ii=next_ii_sm;
            if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
                measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                ii_sm_post=ii_sm_post+1;
                dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                trial_no=trial_no+1;
                decisions_per_trial(trial_no)=0;
                which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                this_ii=this_ii+next_ii+no_points_post;
                ii_post=ii_post+no_points_post;
                ii_pre=ii_pre+no_points_pre;
            else
                at_end=1;
            end
        end
        
        if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                next_ii=next_ii_sm;
                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=0;
                    ii_sm_post=ii_sm_post+1;
                    dFFs_sm_per_trial_per_ROI(ii_sm_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sm(ii_sm_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=0;
                    which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                    epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                    measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                    epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+no_points_post;
                    ii_post=ii_post+no_points_post;
                    ii_pre=ii_pre+no_points_pre;
                else
                    at_end=1;
                end
            else
                next_ii=next_ii_sp;
                if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
                    measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
                    training_decisions_post(ii_post+1:ii_post+no_points_post,:)=1;
                    ii_sp_post=ii_sp_post+1;
                    dFFs_sp_per_trial_per_ROI(ii_sp_post,:,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+ii_p_threshold);
                    dFFs_this_trial=traces(:,no_points_post_shift+this_ii+next_ii-ii_span:no_points_post_shift+this_ii+next_ii+ii_span);
                    dFF_per_trial_sp(ii_sp_post,:,:)=dFFs_this_trial;
                    trial_no=trial_no+1;
                    decisions_per_trial(trial_no)=1;
                    which_model_for_traces_loo(1,ii_which_model+1:no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model)=trial_no;
                    ii_which_model=no_points_post_shift+this_ii+next_ii+no_points_post-1+dt_post_which_model;
                    ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
                    epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
                    measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
                    epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
                    this_ii=this_ii+next_ii+no_points_post;
                    ii_post=ii_post+no_points_post;
                    ii_pre=ii_pre+no_points_pre;
                else
                    at_end=1;
                end
            end
        end
        
        
    end
    
end
%
% %Do S+
% at_end=0;
% this_ii=0;
% ii_post=0;
% ii_pre=0;
% ii=0;
% trial_no=0;
%
% while (at_end==0)
%     next_ii=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
%     if ~isempty(next_ii)
%         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
%             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
%             trial_no=trial_no+1;
%             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
%             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
%             epochs_sp_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
%             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
%             epochs_sp_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
%             this_ii=this_ii+next_ii+no_points_post;
%             ii_post=ii_post+no_points_post;
%             ii_pre=ii_pre+no_points_pre;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
%     end
% end
%
% training_decisions_post=ones(1,size(measurements_post,1));
% ii_sp_post=size(measurements_post,1);
%
%
% %Do S-
% at_end=0;
% this_ii=0;
% ii=0;
%
%
% while (at_end==0)
%     next_ii=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
%     if ~isempty(next_ii)
%         if (no_points_post_shift+this_ii+next_ii+no_points_post<length(epochs))&(no_points_post_shift+this_ii+next_ii-no_points_pre>0)
%             measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)';
%             trial_no=trial_no+1;
%             which_model_for_traces_loo(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=trial_no;
%             ii_pointer_to_td(ii_post+1:ii_post+no_points_post)=no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1;
%             epochs_sm_post(1,no_points_post_shift+this_ii+next_ii:no_points_post_shift+this_ii+next_ii+no_points_post-1)=1;
%             measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)';
%             epochs_sm_pre(1,no_points_post_shift+this_ii+next_ii-no_points_pre:no_points_post_shift+this_ii+next_ii-1)=1;
%             this_ii=this_ii+next_ii+no_points_post;
%             ii_post=ii_post+no_points_post;
%             ii_pre=ii_pre+no_points_pre;
%         else
%             at_end=1;
%         end
%     else
%         at_end=1;
%     end
% end
%
% ii_sm_post=ii_post-ii_sp_post;
%
% training_decisions_post=[training_decisions_post zeros(1,size(measurements_post,1)-ii_sp_post)];

%Now let's limit the ROIs to those below p_threshold

p_values=ones(1,size(dFFs_sm_per_trial_per_ROI,2));
for iiROI=1:size(dFFs_sm_per_trial_per_ROI,2)
    dFF_sm=zeros(size(dFFs_sm_per_trial_per_ROI,1),size(dFFs_sm_per_trial_per_ROI,3));
    dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI(:,iiROI,:);
    dFF_sp=zeros(size(dFFs_sp_per_trial_per_ROI,1),size(dFFs_sp_per_trial_per_ROI,3));
    dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI(:,iiROI,:);
    
    [h,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
end

p_value_mask=logical(p_values<=p_threshold);

%Trim the number of ROIs in all matrices
noROIs_before_trimming=size(measurements_post,2);
dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
measurements_post=measurements_post(:,p_value_mask);
measurements_pre=measurements_pre(:,p_value_mask);
traces=traces(p_value_mask,:);
no_traces=size(traces,1);

%Plot the trimmed traces
%time has the time for the dF/F traces(ROI,time)
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on
    
    % Determine the y spacing of the traces
    y_shift=1.2*(prctile(traces(:),95)-prctile(traces(:),5));
    
    %Plot the traces and do z normalization
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
        end
    end
    
    for trNo=1:no_traces
        % for trNo=1:20
        plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
        traces(trNo,:)=traces(trNo,:)/std(traces(trNo,:));
    end
    
    ylim([-y_shift*0.2 (no_traces+2)*y_shift])
    xlabel('time(sec)')
    title(['dFF timecourses after p value trimming ' num2str(size(measurements_post,2)) ' ROIs'])
end

%Plot the mean timecourse for S+ and S-
if show_figures==1
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    
    hFig = figure(figNo);
    
    set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
    
    
    hold on
    
    time_span=[0:dt:dt*size(dFF_per_trial_sp,3)]-dt_span+dt;
    time_span=time_span(1:end-1)+post_shift;
    
    meandFF_per_trial_sp=zeros(size(dFF_per_trial_sp,1),size(dFF_per_trial_sp,3));
    meandFF_per_trial_sp(:,:)=mean(dFF_per_trial_sp,2);
    
    meandFF_per_trial_sm=zeros(size(dFF_per_trial_sm,1),size(dFF_per_trial_sm,3));
    meandFF_per_trial_sm(:,:)=mean(dFF_per_trial_sm,2);
    
    %S+ trials
    try
        CIsp = bootci(1000, @mean, meandFF_per_trial_sp);
        meansp=mean(dFF_per_trial_sp,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
    catch
    end
    
    
    [hlsp, hpsp] = boundedline(time_span',mean(meandFF_per_trial_sp,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
    
    %S-
    try
        CIsp = bootci(1000, @mean, meandFF_per_trial_sm);
        meansp=mean(dFF_per_trial_sm,1);
        CIsp(1,:)=meansp-CIsp(1,:);
        CIsp(2,:)=CIsp(2,:)-meansp;
    catch
    end
    
    
    [hlsp, hpsp] = boundedline(time_span',mean(meandFF_per_trial_sm,1)', CIsp', 'cmap',[238/255 111/255 179/255]);
    
    plot(time_span',mean(meandFF_per_trial_sp,1)', 'Color',[80/255 194/255 255/255]);
    plot(time_span',mean(meandFF_per_trial_sm,1)', 'Color',[238/255 111/255 179/255]);
    
    text(30,0.75,'S-','Color',[80/255 194/255 255/255])
    text(30,0.65,'S+','Color',[0 114/255 178/255])
    
    ylim([0 1])
    
    title(['Mean dFF'])
    xlabel('Time(sec)')
    ylabel('dFF')
end


training_decisions_post_sh=zeros(1,length(training_decisions_post));
training_decisions_post_sh(1,:)=training_decisions_post(randperm(length(training_decisions_post)));

for ii=1:no_shuffles
    these_shuffled_decisions_per_trial=decisions_per_trial(randperm(length(decisions_per_trial)));
    ww=0;
    for jj=1:length(decisions_per_trial)
        training_decisions_post_sh2(ii,ww+1:ww+no_points_post)=these_shuffled_decisions_per_trial(jj)*ones(1,no_points_post);
        ww=ww+no_points_post;
    end
end

handles_out.Nall=Nall;
handles_out.dt=dt;
handles_out.no_points_post=no_points_post;
handles_out.no_points_post_shift=no_points_post_shift;
handles_out.no_points_pre=no_points_pre;
handles_out.measurements_post=measurements_post;
handles_out.measurements_pre=measurements_pre;
handles_out.training_decisions_post=training_decisions_post;
handles_out.epochs_sp_post=epochs_sp_post;
handles_out.epochs_sp_pre=epochs_sp_pre;
handles_out.epochs_sm_post=epochs_sm_post;
handles_out.epochs_sm_pre=epochs_sm_pre;


fprintf(1, ['Training post with %d ROIs (original no ROIs %d)...\n'],size(measurements_post,2),noROIs_before_trimming);



for MLalgo=MLalgo_to_use
    for ii_cost=1:1  %Note: cost between 1 to 3 does not make a difference
        
        this_cost=[0 ii_cost;ii_cost 0];
        labels=[];
        timepoint_processed=[];
        correct_predict=[];
        correct_predict_shuffled=[];
        
        
        
        Nall_post=size(measurements_post,1);
        
        
        
        switch processing_algorithm
            case 1
                %Train with all the post data
                %Store the training data in a table.
                tblTrn=[];
                tblTrn = array2table(measurements_post);
                
                %Store the decisions in Y
                Y=training_decisions_post;
                
                switch MLalgo
                    case 1
                        Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                    case 2
                        Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                    case 3
                        Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                    case 4
                        Mdl = fitcnet(tblTrn,Y);
                    case 5
                        Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                end
                
                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason tisdid not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post);
                %I had to resolrt to the for loop:
                switch MLalgo
                    case 3,4
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label_traces(ii),score] = predict(Mdl,this_time_point);
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post(ii),score] = predict(Mdl,this_time_point);
                        end
                    otherwise
                        [label_traces,score] = predict(Mdl,traces');
                        [label_post,score] = predict(Mdl,measurements_post);
                end
                
            case 2
                %k-fold cross validation evaluating performance with left-out data
                %Store the training data in a table.
                
                handles_out.MLalgo(MLalgo).models=[];
                handles_out.MLalgo(MLalgo).processed_succesfully=1;
                points_masked=floor(no_points_post/k_fold);
                which_model=ones(1,size(measurements_post,1));
                for kk=1:k_fold
                    
                    
                    training_mask=ones(size(measurements_post,1),1);
                    at_end=0;
                    ii=0;
                    
                    while at_end==0
                        if ii+no_points_post<=size(measurements_post,1)
                            %This is what I originally used, these samples are adjacent
                            %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                            %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                            %
                            %I then changed it to this, spaced by delta_ii
                            delta_ii=floor(no_points_post/points_masked);
                            for ii_pm=1:points_masked
                                training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                            end
                            
                            ii=ii+no_points_post;
                        else
                            at_end=1;
                        end
                    end
                    
                    these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                    these_training_decisions=zeros(1,sum(training_mask));
                    
                    
                    jj=0;
                    for ii=1:size(measurements_post,1)
                        if training_mask(ii)==1
                            jj=jj+1;
                            these_training_measurements(jj,:)=measurements_post(ii,:);
                            these_training_decisions(jj)=training_decisions_post(ii);
                        end
                    end
                    
                    tblTrn=[];
                    tblTrn = array2table(these_training_measurements);
                    
                    %Store the decisions in Y
                    Y=these_training_decisions;
                    
                    switch MLalgo
                        case 1
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        case 2
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            %The try catch was entered here because of this
                            %error
                            % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                            % A normal distribution cannot be fit for the combination of class 1 and predictor
                            % these_training_measurements107. The data has zero variance.
                            try
                                handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                            catch
                                handles_out.MLalgo(MLalgo).processed_succesfully=0;
                            end
                        case 4
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                        case 5
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                        case 6
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                    end
                end
                
                if handles_out.MLalgo(MLalgo).processed_succesfully==1
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason tisdid not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resolrt to the for loop:
                    
                    which_model_for_traces=randi(k_fold,1,size(traces,2));
                    
                    for ii=1:length(which_model)
                        which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                    end
                    
                    if MLalgo==6
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_traces(ii)=1;
                            else
                                label_traces(ii)=0;
                            end
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_post(ii)=1;
                            else
                                label_post(ii)=0;
                            end
                        end
                    else
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label_traces(ii),score] = predict(handles_out.MLalgo(MLalgo).models(which_model_for_traces(ii)).Mdl,this_time_point);
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post(ii),score] = predict(handles_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                        end
                    end
                    
                    
                    handles_out.MLalgo(MLalgo).label_post=label_post;
                    handles_out.MLalgo(MLalgo).label_traces=label_traces;
                    
                    %Now do prediction with shuffled training decisions
                    %k-fold cross validation evaluating performance with left-out data
                    %Store the training data in a table.
                    
                    handles_out.MLalgo(MLalgo).sh_models=[];
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post,1));
                    for kk=1:k_fold
                        
                        
                        training_mask=ones(size(measurements_post,1),1);
                        at_end=0;
                        ii=0;
                        
                        while at_end==0
                            if ii+no_points_post<=size(measurements_post,1)
                                %This is what I originally used, these samples are adjacent
                                %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                %
                                %I then changed it to this, spaced by delta_ii
                                delta_ii=floor(no_points_post/points_masked);
                                for ii_pm=1:points_masked
                                    training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                    which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                end
                                
                                ii=ii+no_points_post;
                            else
                                at_end=1;
                            end
                        end
                        
                        these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                        these_training_decisions=zeros(1,sum(training_mask));
                        
                        
                        jj=0;
                        for ii=1:size(measurements_post,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post(ii,:);
                                these_training_decisions(jj)=training_decisions_post_sh(ii);
                            end
                        end
                        
                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);
                        
                        %Store the decisions in Y
                        Y=these_training_decisions;
                        
                        handles_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;
                        
                        switch MLalgo
                            case 1
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                            case 2
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                                end
                            case 4
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason this did not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resort to the for loop:
                    
                    which_model_for_traces=randi(k_fold,1,size(traces,2));
                    
                    for ii=1:length(which_model)
                        which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                    end
                    
                    if MLalgo==6
                        
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_traces_sh(ii)=1;
                            else
                                label_traces_sh(ii)=0;
                            end
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_post_sh(ii)=1;
                            else
                                label_post_sh(ii)=0;
                            end
                        end
                    else
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label_traces_sh(ii),score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model_for_traces(ii)).Mdl,this_time_point);
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post_sh(ii),score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        end
                    end
                    
                    handles_out.MLalgo(MLalgo).label_post=label_post;
                    handles_out.MLalgo(MLalgo).label_traces=label_traces;
                    
                    handles_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                    handles_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
                    
                    %Now do shuffling on a per trial basis
                    for ii_sh=1:no_shuffles
                        
                        
                        %Now do prediction with shuffled training decisions
                        %k-fold cross validation evaluating performance with left-out data
                        %Store the training data in a table.
                        
                        handles_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                        points_masked=floor(no_points_post/k_fold);
                        which_model=ones(1,size(measurements_post,1));
                        for kk=1:k_fold
                            
                            
                            training_mask=ones(size(measurements_post,1),1);
                            at_end=0;
                            ii=0;
                            
                            while at_end==0
                                if ii+no_points_post<=size(measurements_post,1)
                                    %This is what I originally used, these samples are adjacent
                                    %                             training_mask(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=0;
                                    %                             which_model(ii+(kk-1)*points_masked+1:ii+(kk-1)*points_masked+points_masked)=kk;
                                    %
                                    %I then changed it to this, spaced by delta_ii
                                    delta_ii=floor(no_points_post/points_masked);
                                    for ii_pm=1:points_masked
                                        training_mask(ii+kk+(ii_pm-1)*delta_ii)=0;
                                        which_model(ii+kk+(ii_pm-1)*delta_ii)=kk;
                                    end
                                    
                                    ii=ii+no_points_post;
                                else
                                    at_end=1;
                                end
                            end
                            
                            these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                            these_training_decisions=zeros(1,sum(training_mask));
                            
                            
                            jj=0;
                            for ii=1:size(measurements_post,1)
                                if training_mask(ii)==1
                                    jj=jj+1;
                                    these_training_measurements(jj,:)=measurements_post(ii,:);
                                    these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                                end
                            end
                            
                            tblTrn=[];
                            tblTrn = array2table(these_training_measurements);
                            
                            %Store the decisions in Y
                            Y=these_training_decisions;
                            
                            handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;
                            
                            switch MLalgo
                                case 1
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                case 2
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                case 3
                                    %The try catch was entered here because of this
                                    %error
                                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                                    % these_training_measurements107. The data has zero variance.
                                    try
                                        handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                    end
                                case 4
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                case 5
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                case 6
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                            end
                        end
                        
                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason this did not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post);
                        %I had to resort to the for loop:
                        
                        which_model_for_traces=randi(k_fold,1,size(traces,2));
                        
                        for ii=1:length(which_model)
                            which_model_for_traces(ii_pointer_to_td(ii))=which_model(ii);
                        end
                        
                        if MLalgo==6
                            
                            
                            
                            for ii=1:size(measurements_post,1)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=measurements_post(ii,:);
                                [label,score] = predict(handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_post_sh2(ii_sh,ii)=1;
                                else
                                    label_post_sh2(ii_sh,ii)=0;
                                end
                            end
                        else
                            
                            
                            for ii=1:size(measurements_post,1)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=measurements_post(ii,:);
                                [label_post_sh2(ii_sh,ii),score] = predict(handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                            end
                        end
                        
                        handles_out.MLalgo(MLalgo).label_post=label_post_sh2;
                        
                    end
                    
                    
                end
            case 3
                %leave one trial out
                %Store the training data in a table.
                
                handles_out.MLalgo(MLalgo).models=[];
                handles_out.MLalgo(MLalgo).processed_succesfully=1;
                points_masked=floor(no_points_post/k_fold);
                which_model=ones(1,size(measurements_post,1));
                no_trials=ii_post/no_points_post;
                for kk=1:no_trials
                    
                    
                    training_mask=ones(size(measurements_post,1),1);
                    training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                    which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                    
                    
                    these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                    these_training_decisions=zeros(1,sum(training_mask));
                    
                    jj=0;
                    for ii=1:size(measurements_post,1)
                        if training_mask(ii)==1
                            jj=jj+1;
                            these_training_measurements(jj,:)=measurements_post(ii,:);
                            these_training_decisions(jj)=training_decisions_post(ii);
                        end
                    end
                    
                    tblTrn=[];
                    tblTrn = array2table(these_training_measurements);
                    
                    %Store the decisions in Y
                    Y=these_training_decisions;
                    
                    switch MLalgo
                        case 1
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        case 2
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            %The try catch was entered here because of this
                            %error
                            % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                            % A normal distribution cannot be fit for the combination of class 1 and predictor
                            % these_training_measurements107. The data has zero variance.
                            try
                                handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                            catch
                                handles_out.MLalgo(MLalgo).processed_succesfully=0;
                            end
                        case 4
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                        case 5
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                        case 6
                            handles_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                    end
                end
                
                if handles_out.MLalgo(MLalgo).processed_succesfully==1
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason this did not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resort to the for loop:
                    
                    
                    if MLalgo==6
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_traces(ii)=1;
                            else
                                label_traces(ii)=0;
                            end
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_post(ii)=1;
                            else
                                label_post(ii)=0;
                            end
                        end
                    else
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            try
                                [label_traces(ii),score] = predict(handles_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            catch
                                if rand(1)>0.5
                                    label_traces(ii)=1;
                                else
                                    label_traces(ii)=0;
                                end
                            end
                            
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            try
                                [label_post(ii),score] = predict(handles_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                            catch
                                if rand(1)>0.5
                                    label_post(ii)=1;
                                else
                                    label_post(ii)=0;
                                end
                            end
                            
                        end
                    end
                    
                    
                    handles_out.MLalgo(MLalgo).label_post=label_post;
                    handles_out.MLalgo(MLalgo).label_traces=label_traces;
                    
                    %Now do prediction with shuffled training decisions
                    %k-fold cross validation evaluating performance with left-out data
                    %Store the training data in a table.
                    
                    handles_out.MLalgo(MLalgo).sh_models=[];
                    points_masked=floor(no_points_post/k_fold);
                    which_model=ones(1,size(measurements_post,1));
                    for kk=1:no_trials
                        
                        
                        training_mask=ones(size(measurements_post,1),1);
                        training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                        which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                        
                        
                        these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                        these_training_decisions=zeros(1,sum(training_mask));
                        
                        
                        jj=0;
                        for ii=1:size(measurements_post,1)
                            if training_mask(ii)==1
                                jj=jj+1;
                                these_training_measurements(jj,:)=measurements_post(ii,:);
                                these_training_decisions(jj)=training_decisions_post_sh(ii);
                            end
                        end
                        
                        tblTrn=[];
                        tblTrn = array2table(these_training_measurements);
                        
                        %Store the decisions in Y
                        Y=these_training_decisions;
                        
                        handles_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;
                        
                        switch MLalgo
                            case 1
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                            case 2
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                            case 3
                                
                                %The try catch was entered here because of this
                                %error
                                % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                % A normal distribution cannot be fit for the combination of class 1 and predictor
                                % these_training_measurements107. The data has zero variance.
                                try
                                    handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                catch
                                    handles_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                                end
                            case 4
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                            case 5
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                            case 6
                                handles_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                        end
                    end
                    
                    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                    %Note: for some reason this did not work for net when I did this:
                    % [label_traces,score] = predict(Mdl,traces');
                    % [label_post,score] = predict(Mdl,measurements_post);
                    %I had to resort to the for loop:
                    
                    
                    
                    if MLalgo==6
                        
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_traces_sh(ii)=1;
                            else
                                label_traces_sh(ii)=0;
                            end
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label,score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                            if label>0.5
                                label_post_sh(ii)=1;
                            else
                                label_post_sh(ii)=0;
                            end
                        end
                    else
                        for ii=1:size(traces,2)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=traces(:,ii);
                            [label_traces_sh(ii),score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        end
                        
                        for ii=1:size(measurements_post,1)
                            this_time_point=zeros(1,size(traces,1));
                            this_time_point(1,:)=measurements_post(ii,:);
                            [label_post_sh(ii),score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        end
                    end
                    
                    handles_out.MLalgo(MLalgo).label_post=label_post;
                    handles_out.MLalgo(MLalgo).label_traces=label_traces;
                    
                    handles_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
                    handles_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
                    
                    %Now do shuffling on a per trial basis
                    for ii_sh=1:no_shuffles
                        
                        
                        %Now do prediction with shuffled training decisions
                        %k-fold cross validation evaluating performance with left-out data
                        %Store the training data in a table.
                        
                        handles_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                        points_masked=floor(no_points_post/k_fold);
                        which_model=ones(1,size(measurements_post,1));
                        for kk=1:no_trials
                            
                            
                            training_mask=ones(size(measurements_post,1),1);
                            training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                            which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;
                            
                            these_training_measurements=zeros(sum(training_mask),size(measurements_post,2));
                            these_training_decisions=zeros(1,sum(training_mask));
                            
                            
                            jj=0;
                            for ii=1:size(measurements_post,1)
                                if training_mask(ii)==1
                                    jj=jj+1;
                                    these_training_measurements(jj,:)=measurements_post(ii,:);
                                    these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                                end
                            end
                            
                            tblTrn=[];
                            tblTrn = array2table(these_training_measurements);
                            
                            %Store the decisions in Y
                            Y=these_training_decisions;
                            
                            handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;
                            
                            switch MLalgo
                                case 1
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                                case 2
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                                case 3
                                    %The try catch was entered here because of this
                                    %error
                                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                                    % these_training_measurements107. The data has zero variance.
                                    try
                                        handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                                    catch
                                        handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                                    end
                                case 4
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                                case 5
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                                case 6
                                    handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                            end
                        end
                        
                        %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                        %Note: for some reason this did not work for net when I did this:
                        % [label_traces,score] = predict(Mdl,traces');
                        % [label_post,score] = predict(Mdl,measurements_post);
                        %I had to resort to the for loop:
                        
                        
                        if MLalgo==6
                            for ii=1:size(traces,2)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=traces(:,ii);
                                [label,score] = predict(handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_traces_sh2(ii_sh,ii)=1;
                                else
                                    label_traces_sh2(ii_sh,ii)=0;
                                end
                            end
                            
                            
                            for ii=1:size(measurements_post,1)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=measurements_post(ii,:);
                                [label,score] = predict(handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                                if label>0.5
                                    label_post_sh2(ii_sh,ii)=1;
                                else
                                    label_post_sh2(ii_sh,ii)=0;
                                end
                            end
                        else
                            
                            for ii=1:size(traces,2)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=traces(:,ii);
                                [label_traces_sh2(ii_sh,ii),score] = predict(handles_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                            end
                            
                            for ii=1:size(measurements_post,1)
                                this_time_point=zeros(1,size(traces,1));
                                this_time_point(1,:)=measurements_post(ii,:);
                                [label_post_sh2(ii_sh,ii),score] = predict(handles_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                            end
                        end
                        
                        handles_out.MLalgo(MLalgo).label_post=label_post_sh2;
                        
                    end
                    
                    
                end
        end
        
        if handles_out.MLalgo(MLalgo).processed_succesfully==1
            %label is the predicted label, and score is the predicted class
            %posterior probability
            for ii=1:length(training_decisions_post)
                if label_post(ii)==training_decisions_post(ii)
                    correct_predict_tr(ii)=1;
                else
                    correct_predict_tr(ii)=0;
                end
            end
            
            
            %Calculate wta for windows of no_points_post
            correct_predict_tr_wta=zeros(1,length(training_decisions_post));
            for ii=1:length(training_decisions_post)-no_points_post
                this_correct=zeros(1,no_points_post);
                for jj=1:no_points_post
                    if label_post(ii+jj-1)==training_decisions_post(ii+jj-1)
                        this_correct(jj)=1;
                    end
                end
                if sum(this_correct)>(no_points_post/2)
                    correct_predict_tr_wta(ii+floor(no_points_post/2))=1;
                else
                    correct_predict_tr_wta(ii+floor(no_points_post/2))=0;
                end
            end
            
            %Now do shuffled
            %label is the predicted label, and score is the predicted class
            %posterior probability
            for ii=1:length(training_decisions_post_sh)
                if label_post_sh(ii)==training_decisions_post_sh(ii)
                    correct_predict_tr_sh(ii)=1;
                else
                    correct_predict_tr_sh(ii)=0;
                end
            end
            
            %Calculate wta for windows of no_points_post
            correct_predict_tr_wta_sh=zeros(1,length(training_decisions_post_sh));
            for ii=1:length(training_decisions_post_sh)-no_points_post
                this_correct=zeros(1,no_points_post);
                for jj=1:no_points_post
                    if label_post_sh(ii+jj-1)==training_decisions_post_sh(ii+jj-1)
                        this_correct(jj)=1;
                    end
                end
                if sum(this_correct)>(no_points_post/2)
                    correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=1;
                else
                    correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=0;
                end
            end
            
            %Now do shuffled per trial
            correct_predict_tr_wta_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
            correct_predict_tr_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
            for ii_sh=1:no_shuffles
                %label is the predicted label, and score is the predicted class
                %posterior probability
                for ii=1:size(training_decisions_post_sh2,2)
                    if label_post_sh2(ii_sh,ii)==training_decisions_post_sh2(ii_sh,ii)
                        correct_predict_tr_sh2(ii_sh,ii)=1;
                    else
                        correct_predict_tr_sh2(ii_sh,ii)=0;
                    end
                end
                
                %Calculate wta for windows of no_points_post
                
                for ii=1:size(training_decisions_post_sh2,2)-no_points_post
                    this_correct=zeros(1,no_points_post);
                    for jj=1:no_points_post
                        if label_post_sh2(ii_sh,ii+jj-1)==training_decisions_post_sh2(ii_sh,ii+jj-1)
                            this_correct(jj)=1;
                        end
                    end
                    if sum(this_correct)>(no_points_post/2)
                        correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=1;
                    else
                        correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=0;
                    end
                end
                
            end
            
            
            fprintf(1, ['Training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                ,ii_cost, sum(correct_predict_tr)/length(correct_predict_tr),sum(correct_predict_tr_wta)/length(correct_predict_tr_wta));
            fprintf(1, ['Shuffled training accuracy for ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                ,ii_cost, sum(correct_predict_tr_sh)/length(correct_predict_tr_sh),sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh));
            fprintf(1, ['Training accuracy for shuffled per trial ' classifier_names{MLalgo} ' and cost %d is %d, wta accuracy is %d\n']...
                ,ii_cost, sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:)),sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:)));
            fprintf(1, ['Mean label trace %d, variance %d\n'],mean(label_traces),var(label_traces));
            
            
            handles_out.MLalgo(MLalgo).correct_predict_tr=correct_predict_tr;
            handles_out.MLalgo(MLalgo).correct_predict_tr_wta=correct_predict_tr_wta;
            handles_out.MLalgo(MLalgo).accuracy_tr=sum(correct_predict_tr)/length(correct_predict_tr);
            handles_out.MLalgo(MLalgo).accuracy_tr_wta=sum(correct_predict_tr_wta)/length(correct_predict_tr_wta);
            
            handles_out.MLalgo(MLalgo).mean_label_traces=mean(label_traces);
            handles_out.MLalgo(MLalgo).var_label_traces=var(label_traces);
            
            moving_mean_label_traces = movmean(label_traces,moving_mean_n);
         

            %Now let's do the carpentry
%             moving_mean_label_traces_sh = movmean(mean(label_traces_sh2),moving_mean_n);
            
            moving_mean_label_traces_sh2 = movmean(label_traces_sh2,moving_mean_n);
            
            if show_figures==1
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                hFig = figure(figNo);
                
                set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
                
                hold on
                
                
%                 CIsm = bootci(1000, @mean, moving_mean_label_traces_sh2);
%                 meansm=mean(moving_mean_label_traces_sh2,1);
%                 CIsm(1,:)=meansm-CIsm(1,:);
%                 CIsm(2,:)=CIsm(2,:)-meansm;
%                 
%                 %S- Proficient
%                 [hlsm, hpsm] = boundedline(time',mean(moving_mean_label_traces_sh2,1)', CIsm', 'cmap',[80/255 194/255 255/255]);
%                 
                
                per95=prctile(moving_mean_label_traces_sh2(:),95);
                per5=prctile(moving_mean_label_traces_sh2(:),5);
                CIsh=[mean(moving_mean_label_traces_sh2(:))-per5 per95-mean(moving_mean_label_traces_sh2(:))]';
                [hlCR, hpCR] = boundedline([time(1) time(end)],[mean(moving_mean_label_traces_sh2(:)) mean(moving_mean_label_traces_sh2(:))], CIsh', 'cmap',[80/255 194/255 255/255]);
                
                
%                 plot(time,moving_mean_label_traces_sh,'-','Color',[80/255 194/255 255/255])
                plot(time,moving_mean_label_traces,'-k','LineWidth',1)
                plot(time,1.1*((epochs==8)+(epochs==9)),'-b')
                plot(time,1.1*((epochs==6)+(epochs==7)),'-r')
                
                ylim([-0.2 1.2])
                title(['Label prediction for entire session for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                
            end
            %Now let's do accounting and show it in a bar graph
            
            %post Splus
            post_label_sp=label_traces(logical(epochs_sp_post));
            points_per_cut=no_points_post;
            no_cuts=floor(length(post_label_sp)/points_per_cut);
            mean_post_label_sp=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_post_label_sp(ii)=mean(post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_post_label_sp=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_post_label_sp=sh_post_label_sp(logical(epochs_sp_post));
            %     mean_sh_post_label_sp=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_post_label_sp(ii)=mean(sh_post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end
            
            %post Sminus
            post_label_sm=label_traces(logical(epochs_sm_post));
            points_per_cut=no_points_post;
            no_cuts=floor(length(post_label_sm)/points_per_cut);
            mean_post_label_sm=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_post_label_sm(ii)=mean(post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_post_label_sm=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_post_label_sm=sh_post_label_sm(logical(epochs_sm_post));
            %     mean_sh_post_label_sm=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_post_label_sm(ii)=mean(sh_post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end
            
            
            %pre Sminus
            pre_label_sm=label_traces(logical(epochs_sm_pre));
            points_per_cut=no_points_pre;
            no_cuts=floor(length(pre_label_sm)/points_per_cut);
            mean_pre_label_sm=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_pre_label_sm(ii)=mean(pre_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            
            
            %pre Splus
            pre_label_sp=label_traces(logical(epochs_sp_pre));
            points_per_cut=no_points_pre;
            no_cuts=floor(length(pre_label_sp)/points_per_cut);
            mean_pre_label_sp=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_pre_label_sp(ii)=mean(pre_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            
            %     sh_pre_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_pre_label=sh_pre_label(logical(epochs_sm_pre+epochs_sp_pre));
            %     mean_sh_pre_label=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_pre_label(ii)=mean(sh_pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end
            
            %all
            all_label=label_traces;
            %points_per_cut=no_points_pre;
            points_per_cut=no_points_post;
            no_cuts=floor(length(all_label)/points_per_cut);
            mean_all_label=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_all_label(ii)=mean(all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_all_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     mean_sh_all_label=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_all_label(ii)=mean(sh_all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end
            
            if show_figures==1
                %Note that pre and post are refrenced to the start of the training period
                edges=[0:0.033:1.2];
                rand_offset=0.8;
                
                
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                hFig = figure(figNo);
                
                set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                
                hold on
                
                bar_offset=1;
                
                %                 %S- pre
                %                 bar_offset=1;
                %
                %                 bar(bar_offset,mean(mean_pre_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[238/255 111/255 179/255])
                %
                %                 %Violin plot
                %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sm...
                %                     ,edges,bar_offset,rand_offset,'k','k',3);
                %
                %                 bar_offset=bar_offset+1;
                %
                %                 %S+ pre
                %                 bar(bar_offset,mean(mean_pre_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[80/255 194/255 255/255])
                %
                %                 %Violin plot
                %                 [mean_out, CIout]=drgViolinPoint(mean_pre_label_sp...
                %                     ,edges,bar_offset,rand_offset,'k','k',3);
                %
                %                 bar_offset=bar_offset+2;
                
                %S- post
                bar(bar_offset,mean(mean_post_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[158/255 31/255 99/255])
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(mean_post_label_sm...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                bar_offset=bar_offset+1;
                
                %S+ post
                bar(bar_offset,mean(mean_post_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 114/255 178/255])
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(mean_post_label_sp...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                
                bar_offset=bar_offset+2;
                bar(bar_offset,mean(mean_all_label),'LineWidth', 3,'EdgeColor','none','FaceColor','m')
                
                %Violin plot
                [mean_out, CIout]=drgViolinPoint(mean_all_label...
                    ,edges,bar_offset,rand_offset,'k','k',3);
                
                xticks([1 2 4 5 7])
                xticklabels({'S- post', 'S+ post','Entire session'})
                title(['Label prediction for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                
            end
            
            %Now show average labels for S+, S-, etc for 30 sec
            handles_out.MLalgo(MLalgo).mean_all_label=mean_all_label;
            handles_out.MLalgo(MLalgo).mean_pre_label_sm=mean_pre_label_sm;
            handles_out.MLalgo(MLalgo).mean_pre_label_sp=mean_pre_label_sp;
            handles_out.MLalgo(MLalgo).mean_post_label_sm=mean_post_label_sm;
            handles_out.MLalgo(MLalgo).mean_post_label_sp=mean_post_label_sp;
            handles_out.MLalgo(MLalgo).mean_all_label=mean_all_label;
            
            at_end=0;
            ii=1;
            sp_ii=0;
            per_trial_sp_timecourse=[];
            epoch_before_sp=[];
            sm_ii=0;
            per_trial_sm_timecourse=[];
            epoch_before_sm=[];
            
            
            last_sp_sm=-1;
            
            
            while at_end==0
                next_ii=[];
                next_ii_sp=find(epochs_sp_post(ii:end)==1,1,'first');
                next_ii_sm=find(epochs_sm_post(ii:end)==1,1,'first');
                if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
                    if next_ii_sp<next_ii_sm
                        next_ii=next_ii_sp;
                        this_sp_sm=1;
                    else
                        next_ii=next_ii_sm;
                        this_sp_sm=0;
                    end
                else
                    if ~isempty(next_ii_sp)
                        next_ii=next_ii_sp;
                        this_sp_sm=1;
                    end
                    if ~isempty(next_ii_sm)
                        next_ii=next_ii_sm;
                        this_sp_sm=0;
                    end
                end
                
                if ~isempty(next_ii)
                    if ((ii+next_ii-ii_span)>0)&((ii+next_ii+ii_span<length(label_traces)))
                        if this_sp_sm==1
                            sp_ii=sp_ii+1;
                            per_trial_sp_timecourse(sp_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            epoch_before_sp(sp_ii)=last_sp_sm;
                            last_sp_sm=1;
                            ii_next_post=find(epochs_sp_post(ii+next_ii:end)==0,1,'first');
                            ii=ii+next_ii+ii_next_post;
                        else
                            sm_ii=sm_ii+1;
                            per_trial_sm_timecourse(sm_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            epoch_before_sm(sm_ii)=last_sp_sm;
                            last_sp_sm=0;
                            ii_next_post=find(epochs_sm_post(ii+next_ii:end)==0,1,'first');
                            ii=ii+next_ii+ii_next_post;
                        end
                    else
                        if  ((ii+next_ii+ii_span>length(label_traces)))
                            at_end=1;
                        else
                            ii=ii+next_ii;
                        end
                    end
                else
                    at_end=1;
                end
                
            end
            
            handles_out.MLalgo(MLalgo).per_trial_sp_timecourse=per_trial_sp_timecourse;
            handles_out.MLalgo(MLalgo).epoch_before_sp=epoch_before_sp;
            handles_out.MLalgo(MLalgo).per_trial_sm_timecourse=per_trial_sm_timecourse;
            handles_out.MLalgo(MLalgo).epoch_before_sm=epoch_before_sm;
            
            
            this_moving_mean_n=10;
            moving_mean_per_trial_sp_timecourse = movmean(per_trial_sp_timecourse',this_moving_mean_n)';
            moving_mean_per_trial_sm_timecourse = movmean(per_trial_sm_timecourse',this_moving_mean_n)';
            
            time_span=[0:dt:dt*size(per_trial_sp_timecourse,2)]-dt_span+dt;
            time_span=time_span(1:end-1)+post_shift;
            
            %Calculate correct predict
            for ii_tr=1:size(per_trial_sp_timecourse,1)
                for ii_time=1:size(per_trial_sp_timecourse,2)
                    if per_trial_sp_timecourse(ii_tr,ii_time)==1
                        this_correct_predict(ii_tr,ii_time)=1;
                    else
                        this_correct_predict(ii_tr,ii_time)=0;
                    end
                end
            end
            
            for ii_tr=1:size(per_trial_sm_timecourse,1)
                for ii_time=1:size(per_trial_sm_timecourse,2)
                    if per_trial_sm_timecourse(ii_tr,ii_time)==0
                        this_correct_predict(ii_tr+sp_ii,ii_time)=1;
                    else
                        this_correct_predict(ii_tr+sp_ii,ii_time)=0;
                    end
                end
            end
            
            
            %Calculate correct predict shuffled
            ii_plus=0;
            for ww=1:10
                for ii_tr=1:size(per_trial_sp_timecourse,1)
                    rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                    for ii_time=1:size(per_trial_sp_timecourse,2)
                        if per_trial_sp_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                        else
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                        end
                    end
                end
                
                ii_plus=ii_plus+sp_ii;
                
                for ii_tr=1:size(per_trial_sm_timecourse,1)
                    rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                    for ii_time=1:size(per_trial_sm_timecourse,2)
                        if per_trial_sm_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                        else
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                        end
                    end
                end
                ii_plus=ii_plus+sm_ii;
            end
            
            handles_out.MLalgo(MLalgo).this_correct_predict=this_correct_predict;
            handles_out.MLalgo(MLalgo).this_correct_predict_sh=this_correct_predict_sh;
            
            if show_figures==1
                %Plot the prediction for S+ and S- and the correct predictions
                figNo=figNo+1;
                try
                    close(figNo)
                catch
                end
                
                hFig = figure(figNo);
                
                set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                
                subplot(2,1,1)
                hold on
                
                CIsm = bootci(1000, @mean, moving_mean_per_trial_sm_timecourse);
                meansm=mean(moving_mean_per_trial_sm_timecourse,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                [hlsm, hpsm] = boundedline(time_span',mean(moving_mean_per_trial_sm_timecourse,1)', CIsm', 'cmap',[158/255 31/255 99/255]);
                
                CIsp = bootci(1000, @mean, moving_mean_per_trial_sp_timecourse);
                meansp=mean(moving_mean_per_trial_sp_timecourse,1);
                CIsp(1,:)=meansp-CIsp(1,:);
                CIsp(2,:)=CIsp(2,:)-meansp;
                
                
                [hlsp, hpsp] = boundedline(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
                
                plot(time_span',mean(moving_mean_per_trial_sm_timecourse,1)','-','Color',[158/255 31/255 99/255],'DisplayName','S-')
                plot(time_span',mean(moving_mean_per_trial_sp_timecourse,1)', '-','Color',[0 114/255 178/255],'DisplayName','S+');
                
                text(30,0.75,'S-','Color',[158/255 31/255 99/255])
                text(30,0.65,'S+','Color',[0 114/255 178/255])
                
                ylim([0 1])
                title(['Label prediction per trial for ' classifier_names{MLalgo} ' and p value threshold ' num2str(p_threshold)])
                xlabel('Time(sec)')
                ylabel('Label prediction, S+=1, S-=0')
                
                subplot(2,1,2)
                hold on
                
                CIsm = bootci(1000, @mean, this_correct_predict_sh);
                meansm=mean(this_correct_predict_sh,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict_sh,1)', CIsm', 'k');
                
                
                CIsm = bootci(1000, @mean, this_correct_predict);
                meansm=mean(this_correct_predict,1);
                CIsm(1,:)=meansm-CIsm(1,:);
                CIsm(2,:)=CIsm(2,:)-meansm;
                
                [hlsm, hpsm] = boundedline(time_span',mean(this_correct_predict,1)', CIsm', 'cmap',[0 114/255 178/255]);
                
                plot(time_span',mean(this_correct_predict_sh,1)','-k','DisplayName','Shuffled')
                plot(time_span',mean(this_correct_predict,1)', '-','Color',[0 114/255 178/255]);
                
                text(30,0.75,'Shuffled','Color','k')
                text(30,0.65,'S+ vs S-','Color',[0 114/255 178/255])
                
                ylim([0 1])
                title(['Percent correct accuracy per trial for ' classifier_names{MLalgo} ' and cost ' num2str(ii_cost)])
                xlabel('Time(sec)')
                ylabel('Accuracy')
                
                %                 %Plot the prediction for S+ and S- sorted by previous trial
                %                 figNo=figNo+1;
                %                 try
                %                     close(figNo)
                %                 catch
                %                 end
                %
                %                 hFig = figure(figNo);
                %
                %                 set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
                %
                %                 subplot(2,1,1)
                %                 hold on
                %
                %                 %S+ trials
                %                 %Previous trial was S-
                %                 this_sp_sm_beforetimecourse=moving_mean_per_trial_sp_timecourse(epoch_before_sp==0,:);
                %                 try
                %                     CIsp = bootci(1000, @mean, this_sp_sm_beforetimecourse);
                %                     meansp=mean(this_sp_sm_beforetimecourse,1);
                %                     CIsp(1,:)=meansp-CIsp(1,:);
                %                     CIsp(2,:)=CIsp(2,:)-meansp;
                %                 catch
                %                 end
                %
                %
                %                 [hlsp, hpsp] = boundedline(time_span',mean(this_sp_sm_beforetimecourse,1)', CIsp', 'cmap',[80/255 194/255 255/255]);
                %
                %                 %Previous trial was S+
                %                 this_sp_sp_beforetimecourse=moving_mean_per_trial_sp_timecourse(epoch_before_sp==1,:);
                %                 try
                %                     CIsp = bootci(1000, @mean, this_sp_sp_beforetimecourse);
                %                     meansp=mean(this_sp_sp_beforetimecourse,1);
                %                     CIsp(1,:)=meansp-CIsp(1,:);
                %                     CIsp(2,:)=CIsp(2,:)-meansp;
                %                 catch
                %                 end
                %
                %
                %                 [hlsp, hpsp] = boundedline(time_span',mean(this_sp_sp_beforetimecourse,1)', CIsp', 'cmap',[0 114/255 178/255]);
                %
                %                 plot(time_span',mean(this_sp_sm_beforetimecourse,1)','-','Color',[80/255 194/255 255/255])
                %                 plot(time_span',mean(this_sp_sp_beforetimecourse,1)', '-','Color',[0 114/255 178/255]);
                %                 text(30,0.75,'pre S-','Color',[80/255 194/255 255/255])
                %                 text(30,0.65,'pre S+','Color',[0 114/255 178/255])
                %
                %                 ylim([0 1])
                %
                %                 title(['S+ prediction sorted by previous S+ or S- for ' classifier_names{MLalgo} ' and cost ' num2str(ii_cost)])
                %                 xlabel('Time(sec)')
                %                 ylabel('Prediction')
                %
                %
                %                 %S- trials
                %                 subplot(2,1,2)
                %                 hold on
                %
                %                 %Previous trial was S+
                %                 this_sm_sp_beforetimecourse=moving_mean_per_trial_sm_timecourse(epoch_before_sm==1,:);
                %                 try
                %                     CIsp = bootci(1000, @mean, this_sm_sp_beforetimecourse);
                %                     meansp=mean(this_sm_sp_beforetimecourse,1);
                %                     CIsp(1,:)=meansp-CIsp(1,:);
                %                     CIsp(2,:)=CIsp(2,:)-meansp;
                %                 catch
                %                 end
                %
                %
                %                 [hlsp, hpsp] = boundedline(time_span',mean(this_sm_sp_beforetimecourse,1)', CIsp', 'cmap',[238/255 111/255 179/255]);
                %
                %                 %Previous trial was S-
                %                 this_sm_sm_beforetimecourse=moving_mean_per_trial_sm_timecourse(epoch_before_sm==0,:);
                %                 try
                %                     CIsp = bootci(1000, @mean, this_sm_sm_beforetimecourse);
                %                     meansp=mean(this_sm_sm_beforetimecourse,1);
                %                     CIsp(1,:)=meansp-CIsp(1,:);
                %                     CIsp(2,:)=CIsp(2,:)-meansp;
                %                 catch
                %                 end
                %
                %
                %                 [hlsp, hpsp] = boundedline(time_span',mean(this_sm_sm_beforetimecourse,1)', CIsp', 'cmap',[158/255 31/255 99/255]);
                %
                %                 plot(time_span',mean(this_sm_sp_beforetimecourse,1)','-','Color',[238/255 111/255 179/255])
                %                 plot(time_span',mean(this_sm_sm_beforetimecourse,1)', '-','Color',[158/255 31/255 99/255]);
                %                 text(30,0.75,'pre S-','Color',[158/255 31/255 99/255])
                %                 text(30,0.65,'pre S+','Color',[238/255 111/255 179/255])
                %
                %                 ylim([0 1])
                %
                %                 title(['S- prediction sorted by previous S+ or S- for ' classifier_names{MLalgo} ' and cost ' num2str(ii_cost)])
                %                 xlabel('Time(sec)')
                %                 ylabel('Prediction')
            end
        else
            fprintf(1, [classifier_names{MLalgo} ' was not processed succesfully\n']);
        end
        
    end
end
pffft=1;