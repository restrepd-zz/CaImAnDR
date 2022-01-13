%drgCaImAn_SVZ_entire_session
%This program trains the SVZ with the post odorant and then determines what happens throughout the entire timecouse 
clear all
close all

%Load file
[pre_perFileName,pre_perPathName] = uigetfile({'*pre_per.mat'},'Select the .m file with all the choices for analysis');
load([pre_perPathName pre_perFileName])

post_time=20; %The SVZ will be trained with all points 20 sec following odor on
pre_time=5;

moving_mean_n=20;
no_cuts=20;

this_cost=[0 1.5;1.5 0];

show_figures=1;

MLalgo=2;

%time has the time for the dF/F traces(ROI,time)
figNo=0;
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

%Plot the traces
for trNo=1:no_traces
    % for trNo=1:20
    plot(time,traces(trNo,:)+y_shift*trNo,'-k','LineWidth',1)
end

ylim([-y_shift*0.2 (no_traces+2)*y_shift])
xlabel('time(sec)')

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
no_points_post=length(time(1:1+post_time));
no_points_pre=length(time(1:1+pre_time));
measurements_post=[];
measurements_pre=[];
training_decisions_post=[];

%Do S+
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;

epochs_sp_post=zeros(1,length(time));
sp_ii=[];

epochs_sp_pre=zeros(1,length(time));
while (at_end==0)
    next_ii=find((epochs(this_ii+1:end)==7)|(epochs(this_ii+1:end)==6),1,'first');
    if ~isempty(next_ii)
        if (this_ii+next_ii+no_points_post<length(epochs))&(this_ii+next_ii-no_points_pre>0)
            measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_post-1)';
            epochs_sp_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
            measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
            epochs_sp_pre(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
            sp_ii(ii+1)=this_ii+next_ii;
            this_ii=this_ii+next_ii+no_points_post;
            ii_post=ii_post+no_points_post;
            ii_pre=ii_pre+no_points_pre;
            ii=ii+1;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

training_decisions_post=ones(1,size(measurements_post,1));
ii_sp_post=size(measurements_post,1);


%Do S-
at_end=0;
this_ii=0;
ii=0;
sm_ii=[];

epochs_sm_post=zeros(1,length(time));
epochs_sm_pre=zeros(1,length(time));
while (at_end==0)
    next_ii=find((epochs(this_ii+1:end)==8)|(epochs(this_ii+1:end)==9),1,'first');
    if ~isempty(next_ii)
        if (this_ii+next_ii+no_points_post<length(epochs))&(this_ii+next_ii-no_points_pre>0)
            measurements_post(ii_post+1:ii_post+no_points_post,:)=traces(:,this_ii+next_ii:this_ii+next_ii+no_points_post-1)';
            epochs_sm_post(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
            measurements_pre(ii_pre+1:ii_pre+no_points_pre,:)=traces(:,this_ii+next_ii-no_points_pre:this_ii+next_ii-1)';
            epochs_sm_pre(1,this_ii+next_ii:this_ii+next_ii+no_points_post-1)=1;
            sm_ii(ii+1)=this_ii+next_ii;
            this_ii=this_ii+next_ii+no_points_post;
            ii_post=ii_post+no_points_post;
            ii_pre=ii_pre+no_points_pre;
            ii=ii+1;
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

training_decisions_post=[training_decisions_post zeros(1,size(measurements_post,1)-ii_sp_post)];

%Do z normalization
mean_pre=mean(measurements_pre,1);
mean_pre=repmat(mean_pre,size(measurements_post,1),1);
SD_pre=std(measurements_pre,1);
SD_pre=repmat(SD_pre,size(measurements_post,1),1);
measurements_post=(measurements_post-mean_pre)./SD_pre;
%
labels=[];
timepoint_processed=[];
correct_predict=[];
correct_predict_shuffled=[];



Nall_post=size(measurements_post,1);

if show_figures==1
    fprintf(1, ['Training post with %d ROIs...\n'],size(measurements_post,2));
end


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

%label is the predicted label, and score is the predicted class
%posterior probability
for ii=1:length(training_decisions_post)
    if label_post(ii)==training_decisions_post(ii)
        correct_predict_tr(ii)=1;
    else
        correct_predict_tr(ii)=0;
    end
end

moving_mean_label_traces = movmean(label_traces,moving_mean_n);



if show_figures==1
    fprintf(1, ['Training accuracy %d\n'],sum(correct_predict)/length(correct_predict));
end

%Do shuffled SVZ
sh_label_traces=[];
sh_label_post=[];
sh_correct_predict_shuffled=[];

for jj=1:50
    ii_shuffled=randperm(length(training_decisions_post));
    
    for ii=1:length(training_decisions_post)
        sh_training_decisions_post(ii)=training_decisions_post(ii_shuffled(ii));
    end
    
    %Store the decisions in Y
    Y=sh_training_decisions_post;
    
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
    
    
    %                     if MLalgo==2
    %                         pffft=1;
    %                     end
    %                     Mdl.Cost(1,2) = 10;
    
    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
    
    [sh_label_traces(jj,:),score] = predict(Mdl,traces');
    [sh_label_post(jj,:),score] = predict(Mdl,measurements_post);
    
    %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
    %Note: for some reason this did not work for net:
    %    [sh_label_traces(jj,:),score] = predict(Mdl,traces');
    %     [sh_label_post(jj,:),score] = predict(Mdl,measurements_post);
    %I had to resolrt to the for loop:
    switch MLalgo
        case 3,4
            for ii=1:size(traces,2)
                this_time_point=zeros(1,size(traces,1));
                this_time_point(1,:)=traces(:,ii);
                [sh_label_traces(jj,ii),score] = predict(Mdl,this_time_point);
            end
            
            for ii=1:size(measurements_post,1)
                this_time_point=zeros(1,size(traces,1));
                this_time_point(1,:)=measurements_post(ii,:);
                [sh_label_post(jj,ii),score] = predict(Mdl,this_time_point);
            end
        otherwise
            [sh_label_traces(jj,:),score] = predict(Mdl,traces');
            [sh_label_post(jj,:),score] = predict(Mdl,measurements_post);
    end
    
    this_shuffled_post=zeros(1,size(sh_label_post,2));
    this_shuffled_post(:,:)=sh_label_post(jj,:);
    
    for ii=1:length(training_decisions_post)
        if this_shuffled_post(ii)==training_decisions_post(ii)
            sh_correct_predict_shuffled(jj,ii)=1;
        else
            sh_correct_predict_shuffled(ii)=0;
        end
    end
    
end

%Now let's do the carpentry

moving_mean_sh_label_traces = movmean(sh_label_traces,moving_mean_n);


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
% subplot(2,1,1)
% hold on
%
% plot(time,(epochs==8)+(epochs==9),'-b')
% plot(time,(epochs==6)+(epochs==7),'-r')
%
% xlabel('time(sec)')
% plot(time,label_traces,'ok','MarkerSize',2)
% ylim([-0.2 1.2])
% title('prediction timecourse')
%
% subplot(2,1,2)
hold on


CIsm = bootci(1000, @mean, moving_mean_sh_label_traces);
meansm=mean(moving_mean_sh_label_traces,1);
CIsm(1,:)=meansm-CIsm(1,:);
CIsm(2,:)=CIsm(2,:)-meansm;

%S- Proficient
[hlsm, hpsm] = boundedline(time',mean(moving_mean_sh_label_traces,1)', CIsm', 'cmap',[80/255 194/255 255/255]);

plot(time,(epochs==8)+(epochs==9),'-b')
plot(time,(epochs==6)+(epochs==7),'-r')
plot(time,moving_mean_label_traces,'-k')

ylim([-0.2 1.2])
title('moving window mean of prediction')

%Now let's do accounting and show it in a bar graph

%post Splus
post_label_sp=label_traces(logical(epochs_sp_post));
points_per_cut=floor(length(post_label_sp)/no_cuts);
mean_post_label_sp=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_post_label_sp(ii)=mean(post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end
sh_post_label_sp=sum(sh_label_traces,1)/size(sh_label_traces,1);
sh_post_label_sp=sh_post_label_sp(logical(epochs_sp_post));
mean_sh_post_label_sp=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_sh_post_label_sp(ii)=mean(sh_post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end

%post Sminus
post_label_sm=label_traces(logical(epochs_sm_post));
points_per_cut=floor(length(post_label_sm)/no_cuts);
mean_post_label_sm=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_post_label_sm(ii)=mean(post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end
sh_post_label_sm=sum(sh_label_traces,1)/size(sh_label_traces,1);
sh_post_label_sm=sh_post_label_sm(logical(epochs_sm_post));
mean_sh_post_label_sm=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_sh_post_label_sm(ii)=mean(sh_post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end


%pre
pre_label=label_traces(logical(epochs_sm_pre+epochs_sp_pre));
points_per_cut=floor(length(pre_label)/no_cuts);
mean_pre_label=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_pre_label(ii)=mean(pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end
sh_pre_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
sh_pre_label=sh_pre_label(logical(epochs_sm_pre+epochs_sp_pre));
mean_sh_pre_label=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_sh_pre_label(ii)=mean(sh_pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end

%all
all_label=label_traces;
points_per_cut=floor(length(all_label)/no_cuts);
mean_all_label=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_all_label(ii)=mean(all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end
sh_all_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
mean_sh_all_label=zeros(1,no_cuts);
for ii=1:no_cuts
    mean_sh_all_label(ii)=mean(sh_all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
end


edges=[0:0.033:1.2];
rand_offset=0.8;


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .85 .3])
% subplot(2,1,1)
% hold on
%
% plot(time,(epochs==8)+(epochs==9),'-b')
% plot(time,(epochs==6)+(epochs==7),'-r')
%
% xlabel('time(sec)')
% plot(time,label_traces,'ok','MarkerSize',2)
% ylim([-0.2 1.2])
% title('prediction timecourse')
%
% subplot(2,1,2)
hold on

%Splus post
bar_offset=1;
bar(bar_offset,mean(mean_sh_post_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_sh_post_label_sp...
    ,edges,bar_offset,rand_offset,'k','k',3);

bar_offset=bar_offset+1;
bar(bar_offset,mean(mean_post_label_sp),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_post_label_sp...
    ,edges,bar_offset,rand_offset,'k','k',3);

%Sminus post
bar_offset=bar_offset+2;
bar(bar_offset,mean(mean_sh_post_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_sh_post_label_sm...
    ,edges,bar_offset,rand_offset,'k','k',3);

bar_offset=bar_offset+1;
bar(bar_offset,mean(mean_post_label_sm),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_post_label_sm...
    ,edges,bar_offset,rand_offset,'k','k',3);

%pre
bar_offset=bar_offset+2;
bar(bar_offset,mean(mean_sh_pre_label),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_sh_pre_label...
    ,edges,bar_offset,rand_offset,'k','k',3);

bar_offset=bar_offset+1;
bar(bar_offset,mean(mean_pre_label),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_pre_label...
    ,edges,bar_offset,rand_offset,'k','k',3);

%all
bar_offset=bar_offset+2;
bar(bar_offset,mean(mean_sh_all_label),'LineWidth', 3,'EdgeColor','none','FaceColor',[213/255 94/255 0/255])

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_sh_all_label...
    ,edges,bar_offset,rand_offset,'k','k',3);

bar_offset=bar_offset+1;
bar(bar_offset,mean(mean_all_label),'LineWidth', 3,'EdgeColor','none','FaceColor','m')

%Violin plot
[mean_out, CIout]=drgViolinPoint(mean_all_label...
    ,edges,bar_offset,rand_offset,'k','k',3);

if show_figures==1
    
end
ii=0;
for time_point=1:ii_zero-1
    this_correct_predict=correct_predict(ii+1:ii+Nall);
    accuracy(time_point)=mean(this_correct_predict(this_correct_predict~=-1));
    this_correct_predict_shuffled=correct_predict_shuffled(ii+1:ii+Nall);
    sh_accuracy(time_point)=mean(this_correct_predict_shuffled(this_correct_predict_shuffled~=-1));
    training_output_labels(time_point,:)=labels(ii+1:ii+Nall);
    this_timepoint_processed=timepoint_processed(ii+1:ii+Nall);
    for jj=1:Nall
        handles_out2.timepoint_processed_pp(time_point,jj)=this_timepoint_processed(jj);
        handles_out2.correct_predict_pp(time_point,jj)=this_correct_predict(jj);
        handles_out2.correct_predict_shuffled_pp(time_point,jj)=this_correct_predict_shuffled(jj);
    end
    if show_figures==1
        fprintf(1, ['For timepoint %d accuracy= %d and shuffled accuracy= %d\n'],time_point,accuracy(time_point),sh_accuracy(time_point));
    end
    ii=ii+Nall;
end
