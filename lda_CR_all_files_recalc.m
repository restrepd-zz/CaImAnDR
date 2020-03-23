%lda_CR_all_files_recalc
%This program generates a summary figure for the per event LDA analysis for Fig. 5
close all
clear all


%Now do the lick frequency
all_licking_frequency=[];
no_CR_trials=zeros(6,2);
lick_freq=zeros(6,2,3);
dFF=zeros(6,2,3);

event_labels{1}='CR licking';
event_labels{2}='CR no licks';
win_labels{1}='Pre-odor';
win_labels{2}='Odor';
win_labels{3}='Reinforcement';


load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_slopes.mat')
no_exps=zeros(length(handles_outs.window),length(handles_outs.window(1).event));
experimentNo=1;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end



%mmPVG04 20180910_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_slopes.mat')
experimentNo=2;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end

%mmPVG05 20181017
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19_mmPVG05_Cerebellum_out_slopes.mat')
experimentNo=3;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end

%mmPVG02 20180515_18
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_slopes.mat')
experimentNo=4;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end

%mmG7f09 20180702_05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_out_slopes.mat')
experimentNo=5;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end

%mmG06  20180419
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_slopes.mat')
experimentNo=6;
for evNo_CR=1:2
    no_CR_trials(experimentNo,evNo_CR)=handles_outs.window(1).no_lick_freq_choice_win(evNo_CR+4);
    for winNo=1:length(handles_outs.window)
        lick_freq(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).lick_freq_choice_win);
        dFF(experimentNo,evNo_CR,winNo)=mean(handles_outs.window(winNo).event(4+evNo_CR).dFF_choice_win);
    end
end

%Plot the percent of zero FR CRs
figNo=0;
figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)
hold on

percent_CR_trials=zeros(6,2);
percent_CR_trials(:,1)=100*no_CR_trials(:,1)./sum(no_CR_trials,2);
percent_CR_trials(:,2)=100*no_CR_trials(:,2)./sum(no_CR_trials,2);

x=0;

for evNo=1:2
    switch evNo
        case 1
            bar(x,mean(percent_CR_trials(:,evNo)),'r')
        case 2
            bar(x,mean(percent_CR_trials(:,evNo)),'b')
    end
    this_percent=zeros(6,1);
    this_percent(:,1)=percent_CR_trials(:,evNo);
    CI = bootci(1000, @mean, this_percent);
    plot([x x],CI,'-k','Linewidth',3)
    plot(x*ones(1,6),this_percent,'ok','MarkerSize',10)
    x=x+1;
    perCorr(evNo).data=this_percent;
end


ylabel('Percent of CR trials')
xticks([0 1])
xticklabels({'Zero Hz','> Zero Hz'})
title('Percent CR trials')
xlim([-1 2])

perCorr(1).description='Zero frequency';
perCorr(2).description='> Zero frequency';

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for percent of trials \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(perCorr);
    fprintf(1, '\n\n')
catch
end


%Plot the lick frequency for CR zero licks vs CR and do glm/drgMutiRanksumorTtest
glm_lickf_ii=0;
glm_lickf=[];
ii_rank=0;
lick_freq_stats=[];

figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)

hold on

ii_bar=0;

for winNo=1:3
    
    for evNo_CR=1:2
        
        
        these_lfs=zeros(6,1);
        these_lfs(:,1)=lick_freq(:,evNo_CR,winNo);
        this_mean_lf=mean(these_lfs);
        switch evNo_CR
            case 1
                bar(ii_bar,this_mean_lf, 'r');
            case 2
                bar(ii_bar,this_mean_lf, 'b');
        end
        
        
        CI=[];
        CI = bootci(1000, {@mean, these_lfs})';
        plot([ii_bar ii_bar],CI,'-k','LineWidth',2)
        plot(ii_bar*ones(1,6),these_lfs,'ok','MarkerSize',10)
        ii_bar=ii_bar+1;
        
        ii_rank=ii_rank+1;
        lick_freq_stats(ii_rank).data=these_lfs;
        lick_freq_stats(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo_CR}];
        
        glm_lickf.data(glm_lickf_ii+1:glm_lickf_ii+length(these_lfs))=these_lfs;
        glm_lickf.event(glm_lickf_ii+1:glm_lickf_ii+length(these_lfs))=evNo_CR;
        glm_lickf.window(glm_lickf_ii+1:glm_lickf_ii+length(these_lfs))=winNo;
        glm_lickf_ii=glm_lickf_ii+length(these_lfs);
        
    end
    ii_bar=ii_bar+1;
end
title('Lick frequency for CR no licks=blue, CR licks=red')


%Perform the glm for lick freequency 
fprintf(1, ['\n\nglm for lick frequency\n'])
tbl = table(glm_lickf.data',glm_lickf.event',glm_lickf.window',...
    'VariableNames',{'lick_frequency','event','window'});
mdl = fitglm(tbl,'lick_frequency~window+event+window*event'...
    ,'CategoricalVars',[2,3])
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for lick frequency \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(lick_freq_stats);
    fprintf(1, '\n\n')
catch
end

%Plot the dFF for CR zero licks vs CR and do glm/drgMutiRanksumorTtest
figNo=figNo+1;
try
    close(figNo)
catch
end
figure(figNo)

hold on

glm_dFF_ii=0;
glm_dFF=[];
ii_rank=0;
dFF_stats=[];

ii_bar=0;

for winNo=1:3
    
    for evNo_CR=1:2
        
        
        these_dFFs=zeros(6,1);
        these_dFFs(:,1)=dFF(:,evNo_CR,winNo);
        this_mean_dFFs=mean(these_dFFs);
        switch evNo_CR
            case 1
                bar(ii_bar,this_mean_dFFs, 'r');
            case 2
                bar(ii_bar,this_mean_dFFs, 'b');
        end
        
        
        CI=[];
        CI = bootci(1000, {@mean, these_dFFs})';
        plot([ii_bar ii_bar],CI,'-k','LineWidth',2)
        plot(ii_bar*ones(1,6),these_dFFs,'ok','MarkerSize',10)
        ii_bar=ii_bar+1;
        
        ii_rank=ii_rank+1;
        dFF_stats(ii_rank).data=these_dFFs;
        dFF_stats(ii_rank).description=[win_labels{winNo} ' ' event_labels{evNo_CR}];
        
        glm_dFF.data(glm_dFF_ii+1:glm_dFF_ii+length(these_dFFs))=these_dFFs;
        glm_dFF.event(glm_dFF_ii+1:glm_dFF_ii+length(these_dFFs))=evNo_CR;
        glm_dFF.window(glm_dFF_ii+1:glm_dFF_ii+length(these_dFFs))=winNo;
        glm_dFF_ii=glm_dFF_ii+length(these_dFFs);
    end
    ii_bar=ii_bar+1;
end
title('dFFfor CR no licks=blue, CR licks=red')          


%Perform the glm for lick freequency 
fprintf(1, ['\n\nglm for dFF\n'])
tbl = table(glm_dFF.data',glm_dFF.event',glm_dFF.window',...
    'VariableNames',{'lick_frequency','event','window'});
mdl = fitglm(tbl,'lick_frequency~window+event+window*event'...
    ,'CategoricalVars',[2,3])
                

%Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for dFF \n\n'])
try
    [output_data] = drgMutiRanksumorTtest(dFF_stats);
    fprintf(1, '\n\n')
catch
end


 

pffft=1;
