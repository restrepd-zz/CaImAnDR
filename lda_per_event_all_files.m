%lda_per_event_all_files
%This program generates a summary figure for the per event LDA analysis for Fig. 5
close all
clear all

no_tests(3,4)=0;
shuffled_percent_correct=[];
%mmPVG04 20180917_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/20180917_mmPVG04_Cerebellum_PCA_events_lda.mat')
experimentNo=1;
for winNo=1:3
   shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG04 20180910_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum/20180910_mmPVG04_Cerebellum_events_lda.mat')
experimentNo=2;
for winNo=1:3
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG05 20181017
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum/20181017_mmPVG05_Cerebellum_out_events_lda.mat')
experimentNo=3;
for winNo=1:3
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmPVG02 20180515_18
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum/20180515_18_mmPVG02_Cerebellum_out_events_lda.mat')
experimentNo=4;
for winNo=1:3
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmG7f09 20180702_05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09_Cerebellum/20180702_05_mmG7f09-Cerebellumbatch_events_lda.mat')
experimentNo=5;
for winNo=1:3
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
       no_tests(winNo,evNo)=no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo);
   end
end

%mmG06  20180419
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum/20180419_mmG06_cerebellumPCAevents_events_lda.mat')
experimentNo=6;
for winNo=1:3
    shuffled_percent_correct=[shuffled_percent_correct handles_events.win(winNo).shuffled_percent_correct];
   for evNo=1:4
       perCorr_per_experiment(winNo,evNo,experimentNo)=100*sum(handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo))/sum(handles_events.win(winNo).events_miss_FA==evNo);
       cum_corr_pred(winNo,evNo,no_tests(winNo,evNo)+1:no_tests(winNo,evNo)+sum(handles_events.win(winNo).events_miss_FA==evNo))=handles_events.win(winNo).correct_predict(handles_events.win(winNo).events_miss_FA==evNo);
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
for winNo=1:3
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
        plot(x*ones(1,6),pcorr,'ok')
        x=x+1;
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
xticks([2 9 16 21])
xticklabels({'3-4 sec','4.5-6 sec','8.5-10 sec','shuffled'})
title('Percent correct odor classification for LDA analysis')


%Now do the ranksum/ttests
event_labels{1}='Hit';
event_labels{2}='Miss';
event_labels{3}='CR';
event_labels{4}='FA';
win_labels{1}='3-4 sec';
win_labels{2}='4.5-6 sec';
win_labels{3}='8.5-10 sec';

for winNo=1:3
    fprintf(1, ['\nTests for significant differences for ' win_labels{winNo} '\n\n']);
    p_vals_corr=[];
    for evNo1=1:4
        pcorr1=zeros(1,6);
            pcorr1(1,:)=perCorr_per_experiment(winNo,evNo1,:);
            [p_val,r_or_t] = drg_ranksum_or_ttest(pcorr1,shuffled_percent_correct);
            p_vals_corr=[p_vals_corr p_val];
            if r_or_t==0
                fprintf(1, ['ranksum p value for pFDR for ' event_labels{evNo1} ' vs. shuffled =%d\n'],p_val);
            else
                fprintf(1, ['t test p value for pFDR for ' event_labels{evNo1} ' vs. shuffled  =%d\n'],p_val);
            end
        
        for evNo2=evNo1+1:4
            pcorr1=zeros(1,6);
            pcorr1(1,:)=perCorr_per_experiment(winNo,evNo1,:);
            pcorr2=zeros(1,6);
            pcorr2(1,:)=perCorr_per_experiment(winNo,evNo2,:);
            [p_val,r_or_t] = drg_ranksum_or_ttest(pcorr1,pcorr2);
            p_vals_corr=[p_vals_corr p_val];
            if r_or_t==0
                fprintf(1, ['ranksum p value for pFDR for ' event_labels{evNo1} ' vs. ' event_labels{evNo2} ' =%d\n'],p_val);
            else
                fprintf(1, ['t test p value for pFDR for ' event_labels{evNo1} ' vs. ' event_labels{evNo2} ' =%d\n'],p_val);
            end
        end
    end
    pFDRcorr=drsFDRpval(p_vals_corr);
    fprintf(1, ['\n\npFDR for significant difference percent correct  = %d\n\n'],pFDRcorr);
end

