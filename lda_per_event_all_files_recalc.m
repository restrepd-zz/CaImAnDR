%lda_per_event_all_files_recalc
%This program generates a summary figure for the per event LDA analysis for Fig. 5


%LDA accuracy for forward go-no go experiments
close all
clear all

no_tests(3,4)=0;
shuffled_percent_correct=[];
perCorr_per_experiment_all=[];
%mmPVG04 20180917_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_LDA_events_events_lda.mat')
experimentNo=1;
for winNo=1:length(handles_events.win)
   shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG04 20180910_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_new_out_events_lda.mat')
experimentNo=2;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG05 20181017
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_mmPVG05_Cerebellum_out_events_lda.mat')
experimentNo=3;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG02 20180515_18
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_lda_new_events_lda.mat')
experimentNo=4;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmG7f09 20180702_05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_lda_sum_events_lda.mat')
experimentNo=5;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmG06  20180419
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_events_lda.mat')
experimentNo=6;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%Plot the figure
try
    close 5
catch
end

hFig5 = figure(5);
set(hFig5, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x=0;
offset=-50;
for winNo=1:length(handles_events.win)
    these_perCorrs=zeros(4,6);
    for evNo=1:4
        if evNo==3
            x=x+1;
        end
        switch evNo
            case 1
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'r')
            case 2
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'c')
            case 3
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'b')
            case 4
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'m')
        end
        CI = bootci(1000, @mean, perCorr_per_experiment(winNo,evNo,:)+offset);
        plot([x x],CI,'-k','Linewidth',3)
        pcorr=zeros(1,6);
        pcorr(1,:)=perCorr_per_experiment(winNo,evNo,:)+offset;
        these_perCorrs(evNo,:)=pcorr(1,:);
        plot(x*ones(1,6),pcorr,'ok','MarkerSize',10)
        x=x+1;
    end
    for ii=1:6
        plot([x-5 x-4 x-2 x-1],these_perCorrs(:,ii),'-k')
    end
    x=x+2;
end

bar(21,mean(shuffled_percent_correct)+offset,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k')
CI = bootci(1000, @mean, shuffled_percent_correct+offset);
 plot([x x],CI,'-k','Linewidth',3)
 
ylim([-60 60])
yticks([-50 0 50])
yticklabels({'0','50','100'})
ylabel('Percent correct')
xticks([2 9 14])
xticklabels({'3-4 sec','4.5-6 sec','shuffled'})
title('Percent correct odor classification for LDA analysis')
xlim([-2 16])

%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{1}='Odorant';
win_labels{2}='Reinforcement';


ii_rank=0;
decod_acc=[];
glm_lda_ii=0;
glm_lda=[];

ii_rank=ii_rank+1;
decode_acc(ii_rank).data=shuffled_percent_correct';
decode_acc(ii_rank).description=['shuffled'];



for winNo=1:2
    for evNo1=1:4
        ii_rank=ii_rank+1;
        this_perCorr=zeros(1,6);
        this_perCorr(1,:)=perCorr_per_experiment(winNo,evNo1,:);
        decode_acc(ii_rank).data=this_perCorr';
        decode_acc(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo1}];
        
        
        glm_lda.data(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=this_perCorr;
        glm_lda.event(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=evNo1;
        glm_lda.window(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=winNo;
        glm_lda_ii=glm_lda_ii+length(this_perCorr);
        
    end
end


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_lda.data',glm_lda.event',glm_lda.window',...
    'VariableNames',{'percent_correct_LDA','event','window'});
mdl = fitglm(tbl,'percent_correct_LDA~window+event+window*event'...
    ,'CategoricalVars',[2,3])

%  fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
%                 tbl = table(glm_averagewave.data',glm_averagewave.perCorr',glm_averagewave.event',...
%                     'VariableNames',{'Average_wave','perCorr','event'});
%                 mdl = fitglm(tbl,'Average_wave~perCorr+event+perCorr*event'...
%                     ,'CategoricalVars',[2,3])
                
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(decode_acc);
    fprintf(1, '\n\n')
catch
end

%Do vartest2
fprintf(1, ['\n\nTest of difference in variance for decoding accuracy \n\n'])
try
    [output_data] = drgMutiVartest2(decode_acc);
    fprintf(1, '\n\n')
catch
end

%Now do the lick frequency
all_licking_frequency=[];
all_dFFs=[];
no_exps=zeros(3,4);
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_slopes.mat')
experimentNo=1;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%mmPVG04 20180910_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_slopes.mat')
experimentNo=2;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%mmPVG05 20181017
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19_mmPVG05_Cerebellum_out_slopes.mat')
experimentNo=3;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%mmPVG02 20180515_18
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_slopes.mat')
experimentNo=4;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%mmG7f09 20180702_05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_out_slopes.mat')
experimentNo=5;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%mmG06  20180419
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellum_licks_slopes.mat')
experimentNo=6;
for winNo=1:3
    for evNo=1:4
        if handles_outs.no_lick_freq_choice_win(evNo)>0
            no_exps(winNo,evNo)=no_exps(winNo,evNo)+1;
            all_licking_frequency(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).lick_freq_choice_win);
            all_dFFs(winNo,evNo,no_exps(winNo,evNo))=mean(handles_outs.window(winNo).event(evNo).dFF_choice_win);
        end
    end
end

%Plot the figure for lick frequency
try
    close 6
catch
end

hFig6 = figure(6);
set(hFig6, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x=0;
for winNo=2:3
    these_lfs=zeros(4,6);
    for evNo=1:4
        if evNo==3
            x=x+1;
        end
        switch evNo
            case 1
                bar(x,mean(all_licking_frequency(winNo,evNo,1:no_exps(winNo,evNo))),'r')
            case 2
                bar(x,mean(all_licking_frequency(winNo,evNo,1:no_exps(winNo,evNo))),'c')
            case 3
                bar(x,mean(all_licking_frequency(winNo,evNo,1:no_exps(winNo,evNo))),'b')
            case 4
                bar(x,mean(all_licking_frequency(winNo,evNo,1:no_exps(winNo,evNo))),'m')
        end
        CI = bootci(1000, @mean, all_licking_frequency(winNo,evNo,:));
        plot([x x],CI,'-k','Linewidth',3)
        lick_f=zeros(1,6);
        lick_f(1,:)=all_licking_frequency(winNo,evNo,1:no_exps(winNo,evNo));
        these_lfs(evNo,:)=lick_f(1,:);
        plot(x*ones(1,6),lick_f,'ok','MarkerSize',10)
        x=x+1;
    end
    for ii=1:6
        plot([x-5 x-4 x-2 x-1],these_lfs(:,ii),'-k')
    end
    x=x+2;
end


 

ylabel('Lick frequency')
xticks([2 9])
xticklabels({'Odor','Reinforcement'})
title('Lick Frequency')
xlim([-2 12])

%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{1}='Pre-Odor';
win_labels{2}='Odorant';
win_labels{3}='Reinforcement';


ii_rank=0;
lick_freq=[];
glm_lickf_ii=0;
glm_lickf=[];


for winNo=2:3
    for evNo1=1:4
        ii_rank=ii_rank+1;
        this_lick_f=zeros(1,6);
        this_lick_f(1,:)=all_licking_frequency(winNo,evNo1,:);
        lick_freq(ii_rank).data=this_lick_f';
        lick_freq(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo1}];
        
        glm_lickf.data(glm_lickf_ii+1:glm_lickf_ii+length(this_lick_f))=this_lick_f;
        glm_lickf.event(glm_lickf_ii+1:glm_lickf_ii+length(this_lick_f))=evNo1;
        glm_lickf.window(glm_lickf_ii+1:glm_lickf_ii+length(this_lick_f))=winNo;
        glm_lickf_ii=glm_lickf_ii+length(this_lick_f);
        
    end
end


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for lick frequency\n'])
tbl = table(glm_lickf.data',glm_lickf.event',glm_lickf.window',...
    'VariableNames',{'lick_frequency','event','window'});
mdl = fitglm(tbl,'lick_frequency~window+event+window*event'...
    ,'CategoricalVars',[2,3])

%  fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
%                 tbl = table(glm_averagewave.data',glm_averagewave.lick_f',glm_averagewave.event',...
%                     'VariableNames',{'Average_wave','lick_f','event'});
%                 mdl = fitglm(tbl,'Average_wave~lick_f+event+lick_f*event'...
%                     ,'CategoricalVars',[2,3])
                
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for lick frequency \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(lick_freq);
    fprintf(1, '\n\n')
catch
end

%Do vartest2
fprintf(1, ['\n\nTest of difference in variance for lick frequency \n\n'])
try
    [output_data] = drgMutiVartest2(lick_freq);
    fprintf(1, '\n\n')
catch
end

%Plot the figure for dFF
try
    close 7
catch
end

hFig7 = figure(7);
set(hFig7, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x=0;
 
for winNo=2:3
    these_dFFs=zeros(4,6);
    for evNo=1:4
        if evNo==3
            x=x+1;
        end 
        switch evNo
            case 1
                bar(x,mean(all_dFFs(winNo,evNo,1:no_exps(winNo,evNo))),'r')
            case 2
                bar(x,mean(all_dFFs(winNo,evNo,1:no_exps(winNo,evNo))),'c')
            case 3
                bar(x,mean(all_dFFs(winNo,evNo,1:no_exps(winNo,evNo))),'b')
            case 4
                bar(x,mean(all_dFFs(winNo,evNo,1:no_exps(winNo,evNo))),'m')
        end
        CI = bootci(1000, @mean, all_dFFs(winNo,evNo,:));
        plot([x x],CI,'-k','Linewidth',3)
        this_dFF=zeros(1,6);
        this_dFF(1,:)=all_dFFs(winNo,evNo,1:no_exps(winNo,evNo));
        these_dFFs(evNo,:)=this_dFF(1,:);
        plot(x*ones(1,6),this_dFF,'ok','MarkerSize',10)
        x=x+1;
        if winNo==2
           pffft=1; 
        end
    end
    for ii=1:6
        plot([x-5 x-4 x-2 x-1],these_dFFs(:,ii),'-k')
    end
    x=x+2;
end


 

ylabel('dFF')
xticks([2 9])
xticklabels({'Pre-Odor','Odor','Reinforcement'})
title('dFF')
xlim([-2 12])

%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{2}='Odorant';
win_labels{3}='Reinforcement';


ii_rank=0;
dFF_stats=[];
glm_dFF_ii=0;
glm_dFF=[];


for winNo=2:3
    for evNo1=1:4
        ii_rank=ii_rank+1;
        this_lick_f=zeros(1,6);
        this_lick_f(1,:)=all_dFFs(winNo,evNo1,:);
        dFF_stats(ii_rank).data=this_lick_f';
        dFF_stats(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo1}];
        
        glm_dFF.data(glm_dFF_ii+1:glm_dFF_ii+length(this_lick_f))=this_lick_f;
        glm_dFF.event(glm_dFF_ii+1:glm_dFF_ii+length(this_lick_f))=evNo1;
        glm_dFF.window(glm_dFF_ii+1:glm_dFF_ii+length(this_lick_f))=winNo;
        glm_dFF_ii=glm_dFF_ii+length(this_lick_f);
        
    end
end


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for dFF\n'])
tbl = table(glm_dFF.data',glm_dFF.event',glm_dFF.window',...
    'VariableNames',{'dFF','event','window'});
mdl = fitglm(tbl,'dFF~window+event+window*event'...
    ,'CategoricalVars',[2,3])

%  fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
%                 tbl = table(glm_averagewave.data',glm_averagewave.lick_f',glm_averagewave.event',...
%                     'VariableNames',{'Average_wave','lick_f','event'});
%                 mdl = fitglm(tbl,'Average_wave~lick_f+event+lick_f*event'...
%                     ,'CategoricalVars',[2,3])
                
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for dFF \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(dFF_stats);
    fprintf(1, '\n\n')
catch
end

%Do vartest2
fprintf(1, ['\n\nTest of difference in variance for dFF \n\n'])
try
    [output_data] = drgMutiVartest2(dFF_stats);
    fprintf(1, '\n\n')
catch
end

%LDA accuracy for reversed go-no go experiments

no_tests(3,4)=0;
shuffled_percent_correct=[];
perCorr_per_experiment=[];
%mmPVG04 20180919_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180919_mmPVG04_Cerebellum new analysis/20180919_mmPVG04_Cerebellum_LDA_eventsrev_events_lda.mat')
experimentNo=1;
for winNo=1:length(handles_events.win)
   shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo+6)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end


%mmG7f09 20180711_15
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180611_mmG7f09_Cerebellum new analysis/20180611_15_mmG7f09_LDAeventsrev_events_lda.mat')
experimentNo=2;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo+6)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmG06  20180423
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180423_mmG06_cerebellum new analysis/20180423_mmG06_LDAeventsrev_events_lda.mat')
experimentNo=3;
for winNo=1:length(handles_events.win)
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       perCorr_per_experiment_all(winNo,evNo,experimentNo+6)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%Plot the figure
try
    close 8
catch
end

hFig8 = figure(8);
set(hFig8, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x=0;
offset=-50;
for winNo=1:length(handles_events.win)
    these_perCorrs=zeros(4,3);
    for evNo=1:4
        if evNo==3
            x=x+1;
        end
        switch evNo
            case 1
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'r')
            case 2
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'c')
            case 3
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'b')
            case 4
                bar(x,mean(perCorr_per_experiment(winNo,evNo,:))+offset,'m')
        end
        CI = bootci(1000, @mean, perCorr_per_experiment(winNo,evNo,:)+offset);
        plot([x x],CI,'-k','Linewidth',3)
        pcorr=zeros(1,3);
        pcorr(1,:)=perCorr_per_experiment(winNo,evNo,:)+offset;
        these_perCorrs(evNo,:)=pcorr(1,:);
        plot(x*ones(1,3),pcorr,'ok','MarkerSize',10)
        x=x+1;
    end
    for ii=1:3
        plot([x-5 x-4 x-2 x-1],these_perCorrs(:,ii),'-k')
    end
    x=x+2;
end

bar(21,mean(shuffled_percent_correct)+offset,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k')
CI = bootci(1000, @mean, shuffled_percent_correct+offset);
 plot([x x],CI,'-k','Linewidth',3)
 
ylim([-60 60])
yticks([-50 0 50])
yticklabels({'0','50','100'})
ylabel('Percent correct')
xticks([2 9 14])
xticklabels({'3-4 sec','4.5-6 sec','shuffled'})
title('Percent correct odor classification for LDA analysis, reversal')
xlim([-2 16])

%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{1}='Odorant';
win_labels{2}='Reinforcement';


ii_rank=0;
decod_acc=[];
glm_lda_ii=0;
glm_lda=[];

ii_rank=ii_rank+1;
decode_acc(ii_rank).data=shuffled_percent_correct';
decode_acc(ii_rank).description=['shuffled'];



for winNo=1:2
    for evNo1=1:4
        ii_rank=ii_rank+1;
        this_perCorr=zeros(1,3);
        this_perCorr(1,:)=perCorr_per_experiment(winNo,evNo1,:);
        decode_acc(ii_rank).data=this_perCorr';
        decode_acc(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo1}];
        
        
        glm_lda.data(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=this_perCorr;
        glm_lda.event(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=evNo1;
        glm_lda.window(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=winNo;
        glm_lda_ii=glm_lda_ii+length(this_perCorr);
        
    end
end


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for LDA percent correct, reversal\n'])
tbl = table(glm_lda.data',glm_lda.event',glm_lda.window',...
    'VariableNames',{'percent_correct_LDA','event','window'});
mdl = fitglm(tbl,'percent_correct_LDA~window+event+window*event'...
    ,'CategoricalVars',[2,3])

%  fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
%                 tbl = table(glm_averagewave.data',glm_averagewave.perCorr',glm_averagewave.event',...
%                     'VariableNames',{'Average_wave','perCorr','event'});
%                 mdl = fitglm(tbl,'Average_wave~perCorr+event+perCorr*event'...
%                     ,'CategoricalVars',[2,3])
                
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy, reversal \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(decode_acc);
    fprintf(1, '\n\n')
catch
end

%Do vartest2
fprintf(1, ['\n\nTest of difference in variance for decoding accuracy, reversal \n\n'])
try
    [output_data] = drgMutiVartest2(decode_acc);
    fprintf(1, '\n\n')
catch
end


%Plot the figure for forward+reversed
try
    close 9
catch
end

hFig9 = figure(9);
set(hFig9, 'units','normalized','position',[.07 .7 .7 .3])
hold on

x=0;
offset=-50;
for winNo=1:length(handles_events.win)
    these_perCorrs=zeros(4,9);
    for evNo=1:4
        if evNo==3
            x=x+1;
        end
        switch evNo
            case 1
                bar(x,mean(perCorr_per_experiment_all(winNo,evNo,:))+offset,'r')
            case 2
                bar(x,mean(perCorr_per_experiment_all(winNo,evNo,:))+offset,'c')
            case 3
                bar(x,mean(perCorr_per_experiment_all(winNo,evNo,:))+offset,'b')
            case 4
                bar(x,mean(perCorr_per_experiment_all(winNo,evNo,:))+offset,'m')
        end
        CI = bootci(1000, @mean, perCorr_per_experiment_all(winNo,evNo,:)+offset);
        plot([x x],CI,'-k','Linewidth',3)
        pcorr=zeros(1,9);
        pcorr(1,:)=perCorr_per_experiment_all(winNo,evNo,:)+offset;
        these_perCorrs(evNo,:)=pcorr(1,:);
        plot(x*ones(1,9),pcorr,'ok','MarkerSize',10)
        x=x+1;
    end
    for ii=1:9
        plot([x-5 x-4 x-2 x-1],these_perCorrs(:,ii),'-k')
    end
    x=x+2;
end

bar(21,mean(shuffled_percent_correct)+offset,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k')
CI = bootci(1000, @mean, shuffled_percent_correct+offset);
 plot([x x],CI,'-k','Linewidth',3)
 
ylim([-60 60])
yticks([-50 0 50])
yticklabels({'0','50','100'})
ylabel('Percent correct')
xticks([2 9 14])
xticklabels({'3-4 sec','4.5-6 sec','shuffled'})
title('Percent correct odor classification for LDA analysis, forward+reversal')
xlim([-2 16])

%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{1}='Odorant';
win_labels{2}='Reinforcement';


ii_rank=0;
decod_acc=[];
glm_lda_ii=0;
glm_lda=[];

ii_rank=ii_rank+1;
decode_acc(ii_rank).data=shuffled_percent_correct';
decode_acc(ii_rank).description=['shuffled'];



for winNo=1:2
    for evNo1=1:4
        ii_rank=ii_rank+1;
        this_perCorr=zeros(1,9);
        this_perCorr(1,:)=perCorr_per_experiment_all(winNo,evNo1,:);
        decode_acc(ii_rank).data=this_perCorr';
        decode_acc(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo1}];
        
        
        glm_lda.data(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=this_perCorr;
        glm_lda.event(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=evNo1;
        glm_lda.window(glm_lda_ii+1:glm_lda_ii+length(this_perCorr))=winNo;
        glm_lda_ii=glm_lda_ii+length(this_perCorr);
        
    end
end


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for LDA percent correct, forward+reversal\n'])
tbl = table(glm_lda.data',glm_lda.event',glm_lda.window',...
    'VariableNames',{'percent_correct_LDA','event','window'});
mdl = fitglm(tbl,'percent_correct_LDA~window+event+window*event'...
    ,'CategoricalVars',[2,3])

%  fprintf(1, ['\n\nglm for average wavelet power for Theta/' freq_names{pacii+1} '\n'])
%                 tbl = table(glm_averagewave.data',glm_averagewave.perCorr',glm_averagewave.event',...
%                     'VariableNames',{'Average_wave','perCorr','event'});
%                 mdl = fitglm(tbl,'Average_wave~perCorr+event+perCorr*event'...
%                     ,'CategoricalVars',[2,3])
                
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for decoding accuracy, forward+reversal \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(decode_acc);
    fprintf(1, '\n\n')
catch
end

%Do vartest2
fprintf(1, ['\n\nTest of difference in variance for decoding accuracy, forward+reversal \n\n'])
try
    [output_data] = drgMutiVartest2(decode_acc);
    fprintf(1, '\n\n')
catch
end

pffft=1;
