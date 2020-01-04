function drgCaImAnBatchPerSessionReversal

%Ask user for the m file that contains information on what the user wants the analysis to be
%This file has all the information on what the user wants done, which files
%to process, what groups they fall into, etc
%
% An example of this file: drgbChoicesDanielPrelim
%
%

close all
clear all

tic
 
[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAnChoices*.m'},'Select the .m file with all the choices for analysis');
fprintf(1, ['\ndrgCaImAnBatchPerSession run for ' choiceFileName '\n\n']);

addpath(choiceBatchPathName)
eval(['handles=' choiceFileName(1:end-2) ';'])
handles.choiceFileName=choiceFileName;
handles.choiceBatchPathName=choiceBatchPathName;

caimanhandles=handles;

%Read the files and calculate the dF/F in each window
num_odor_trials=0;
epochs_per_trial=[];
num_odor_trials_dFF=0;



for filNum=1:caimanhandles.caimandr_choices.no_files
        
    %Read the file
    if iscell(caimanhandles.caimandr_choices.PathName)==0
        load([caimanhandles.caimandr_choices.PathName caimanhandles.caimandr_choices.FileName{filNum}])
    else
        load([caimanhandles.caimandr_choices.PathName{filNum} caimanhandles.caimandr_choices.FileName{filNum}])
    end
    
    
    first_num_odor_trials(filNum)=num_odor_trials+1;
    
    for trNo=1:no_odor_trials
        
        %Save epoch
        num_odor_trials=num_odor_trials+1;
        
        
        if epoch_per_trial(trNo)==6
            %Hit
            epochs_per_trial(1,num_odor_trials)=1;
            epochs_per_trial(2:4,num_odor_trials)=0;
             
            %Was dF/F calculated?
            if sum(which_trial_Hit==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Hit_traces(which_trial_Hit==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Hit_traces(which_trial_Hit==trNo,win),2);
                    szhit=size(Hit_traces(which_trial_Hit==trNo,win));
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_Hitii_lick=which_Hitii_lick(find(which_trial_Hit==trNo,1));
                    these_Hitii_lick_times=[];
                    these_Hitii_lick_times=Hit_lick_times(this_Hitii_lick,1:Hit_no_lick_times(this_Hitii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Hitii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Hitii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=1;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
                
            end
        end
         
        if epoch_per_trial(trNo)==7
            %Miss
            epochs_per_trial(2,num_odor_trials)=1;
            epochs_per_trial(1,num_odor_trials)=0;
            epochs_per_trial(3:4,num_odor_trials)=0;
            
                        %Was dF/F calculated?
            if sum(which_trial_Miss==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(Miss_traces(which_trial_Miss==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(Miss_traces(which_trial_Miss==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                     
                    %Calculate lick frequency for this window
                    this_Missii_lick=which_Missii_lick(find(which_trial_Miss==trNo,1));
                    these_Missii_lick_times=[];
                    these_Missii_lick_times=Miss_lick_times(this_Missii_lick,1:Miss_no_lick_times(this_Missii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_Missii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_Missii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=2;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
        
        if epoch_per_trial(trNo)==8
            %FA
            epochs_per_trial(3,num_odor_trials)=1;
            epochs_per_trial(1:2,num_odor_trials)=0;
            epochs_per_trial(4,num_odor_trials)=0;
            
            %Was dF/F calculated?
            if sum(which_trial_FA==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(FA_traces(which_trial_FA==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                for winNo=1:szwins(1)
                    %Calculate dFFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(FA_traces(which_trial_FA==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_FAii_lick=which_FAii_lick(find(which_trial_FA==trNo,1));
                    these_FAii_lick_times=[];
                    these_FAii_lick_times=FA_lick_times(this_FAii_lick,1:FA_no_lick_times(this_FAii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_FAii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_FAii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=3;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
        
        if epoch_per_trial(trNo)==9
            %CR
            epochs_per_trial(4,num_odor_trials)=1;
            epochs_per_trial(1:3,num_odor_trials)=0;
            
            %Was dF/F calculated?
            if sum(which_trial_CR==trNo)>0
                %Calculate norm dFF
                ref_win=(time_to_event>=caimanhandles.caimandr_choices.ref_win(1))&(time_to_event<=caimanhandles.caimandr_choices.ref_win(2));
                ref_dFF=[];
                ref_dFF=mean(CR_traces(which_trial_CR==trNo,ref_win),2);
                num_odor_trials_dFF=num_odor_trials_dFF+1;
                this_num_odor_trial(num_odor_trials_dFF)=num_odor_trials;
                szwins=size(caimanhandles.caimandr_choices.wins);
                for winNo=1:szwins(1)
                    %Calculate dFF
                    win=(time_to_event>=caimanhandles.caimandr_choices.wins(winNo,1))&(time_to_event<=caimanhandles.caimandr_choices.wins(winNo,2));
                    win_dFF=[];
                    win_dFF=mean(CR_traces(which_trial_CR==trNo,win),2);
%                     norm_dFF=win_dFF./ref_dFF;
                    all_win_dFF(winNo,num_odor_trials_dFF,1:length(win_dFF))=win_dFF;
                    no_traces_win_dFF(winNo,num_odor_trials_dFF)=length(win_dFF);
                    SD_win_dFF(winNo,num_odor_trials_dFF)=std(win_dFF);
                    mean_win_dFF(winNo,num_odor_trials_dFF)=mean(win_dFF);
                    CI_win_dFF(winNo,num_odor_trials_dFF,:) = bootci(1000, @mean, win_dFF);
                    
                    %Calculate lick frequency for this window
                    this_CRii_lick=which_CRii_lick(find(which_trial_CR==trNo,1));
                    these_CRii_lick_times=[];
                    these_CRii_lick_times=CR_lick_times(this_CRii_lick,1:CR_no_lick_times(this_CRii_lick));
                    lick_freq(winNo,num_odor_trials_dFF)=sum( (these_CRii_lick_times>=caimanhandles.caimandr_choices.wins(winNo,1))&...
                        (these_CRii_lick_times<=caimanhandles.caimandr_choices.wins(winNo,2)))/(caimanhandles.caimandr_choices.wins(winNo,2)-...
                        caimanhandles.caimandr_choices.wins(winNo,1));
                end
                epochs_per_trial_dFF(num_odor_trials_dFF)=4;
                trial_dFF(num_odor_trials_dFF)=num_odor_trials;
            end
        end
  
    end
    noROIs(filNum)=szhit(1);
end

%Calculate percent correct
sliding_window=20; %Trials for determination of behavioral performance
min_precent_high_beh=80; %Minimum percent correct for good behavior blocks
max_percent_low_beh=65;

perCorr=[];

%Note: Because this is a reversal I am moving the window for calculation of perCorr to the right by nine points
for ii=1:num_odor_trials-sliding_window+1  
    no_Hits=sum(epochs_per_trial(1,ii:ii+sliding_window-1)==1);
    no_CRs=sum(epochs_per_trial(4,ii:ii+sliding_window-1)==1);
    perCorr(ii+sliding_window-1)=100*(no_Hits+no_CRs)/sliding_window;
end


perCorr(1:sliding_window)=perCorr(2*sliding_window+1);

%Note, this is here so that perCorr=0 is included in the 0-10 % bin.
perCorr(perCorr==0)=0.00001;




%Plot percent correct vs trial
try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.25 .65 .5 .25])

jj_low=find(perCorr<max_percent_low_beh);
plot(jj_low,perCorr(jj_low),'ob')
hold on
jj_high=find(perCorr>min_precent_high_beh);
plot(jj_high,perCorr(jj_high),'or')

jj_mid=find((perCorr<=min_precent_high_beh)&(perCorr>=max_percent_low_beh));
plot(jj_mid,perCorr(jj_mid),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7])
hold on
plot([0 num_odor_trials],[50 50],'-k')

%Draw the boundaries of each file
for filNum=2:caimanhandles.caimandr_choices.no_files
%     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
    if isfield(caimanhandles.caimandr_choices,'start_reversal')
        if caimanhandles.caimandr_choices.start_reversal==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_gogo')
        if caimanhandles.caimandr_choices.start_gogo==filNum
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,80,'Go-go','Color','k','FontSize',18)
        end
    end
    if isfield(caimanhandles.caimandr_choices,'start_session')
        if sum(caimanhandles.caimandr_choices.start_session==filNum)>0
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[0 100],'-k')
%             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
        end
    end
end

title(['Percent correct vs. trial number ' choiceFileName(18:end-2)])
xlabel('Trial number')
ylabel('Percent correct')
ylim([0 100])

%Plot norm dFF for each window
for winNo=1:szwins(1)
    figNo=1+winNo;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .05+0.2*(winNo-1) .5 .2])
    hold on
    
    for dFF_trNo=1:num_odor_trials_dFF
        %Confidence interval
        this_CI=zeros(1,2);
        this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
        plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
        
        if epochs_per_trial_dFF(dFF_trNo)==1
            %Hit
            plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==4
            %CR
            plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==2
            %Miss
            plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==3
            %FA
            plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
        end
        
    end
     
    if isfield(caimanhandles.caimandr_choices,'start_reversal')
        filNum=caimanhandles.caimandr_choices.start_reversal;
        if filNum>0
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Reversal','Color','k','FontSize',18)
        end
    end
    
     if isfield(caimanhandles.caimandr_choices,'start_gogo')
        filNum=caimanhandles.caimandr_choices.start_gogo;
        if filNum>0
            plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
            text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Go-go','Color','k','FontSize',18)
        end
     end
    
     if isfield(caimanhandles.caimandr_choices,'start_session')
         if length(caimanhandles.caimandr_choices.start_session)>=2
             for sessionNo=2:length(caimanhandles.caimandr_choices.start_session)
                 filNum=caimanhandles.caimandr_choices.start_session(sessionNo);
                 plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k')
                 %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
             end
         end
     end
    
    plot([1 num_odor_trials_dFF],[0 0], '-k')
    title(['Normalized dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])
    xlabel('Trial number')
    ylabel('Normalized dF/F')
    ylim([prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))])
 
end

%Now do the analysis of normdFF vs percent correct
% unique_perCorr=25:10:95;
unique_perCorr=5:10:95;
dFF_per_perCorr=[];
lick_freq_per_perCorr=[];
no_per_perCorr=zeros(szwins(1),4,length(unique_perCorr));
no_freqs_per_perCorr=zeros(szwins(1),4,length(unique_perCorr));

dFFspm_per_perCorr=[];
lickspm_freq_per_perCorr=[];
nospm_per_perCorr=zeros(szwins(1),2,length(unique_perCorr));
nospm_freqs_per_perCorr=zeros(szwins(1),2,length(unique_perCorr));

perCorrsubsampled=zeros(szwins(1),floor(length(unique_perCorr)),4,2);
minpc=20000;
maxpc=-20000;

all_mean_ndFF=[];
all_ndFF_CIs=[];
all_ndFF_vaild=[];

for winNo=1:szwins(1)
    for epochs=1:4
        for ii=1:floor(length(unique_perCorr))
            handles_out.mean_exists(winNo,ii,epochs)=0;
        end
    end
    
     for epochs=1:2
        for ii=1:floor(length(unique_perCorr))
            handles_out.mean_existsspm(winNo,ii,epochs)=0;
        end
    end
end


for winNo=1:szwins(1)

    for dFF_trNo=1:num_odor_trials_dFF
        
        if epochs_per_trial_dFF(dFF_trNo)==1
            %Hit
            these_dFFs=zeros(1,1,no_traces_win_dFF(winNo,dFF_trNo));
            these_dFFs(1,1,:)=all_win_dFF(winNo,dFF_trNo,1:no_traces_win_dFF(winNo,dFF_trNo));
%             ii_perCorr=find(unique_perCorr==perCorr(trial_dFF(dFF_trNo)),1);
            ii_perCorr=find( ((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)<=5)&((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)>-5),1);
            
            %Save for Hit
            dFF_per_perCorr(winNo,1,ii_perCorr,no_per_perCorr(winNo,1,ii_perCorr)+1:no_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            no_per_perCorr(winNo,1,ii_perCorr)=no_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs);
           
            no_freqs_per_perCorr(winNo,1,ii_perCorr)=no_freqs_per_perCorr(winNo,1,ii_perCorr)+1;
            lick_freq_per_perCorr(winNo,1,ii_perCorr,no_freqs_per_perCorr(winNo,1,ii_perCorr))=lick_freq(winNo,dFF_trNo);
            
            %Save for S+
            dFFspm_per_perCorr(winNo,1,ii_perCorr,nospm_per_perCorr(winNo,1,ii_perCorr)+1:nospm_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            nospm_per_perCorr(winNo,1,ii_perCorr)=nospm_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs);
            nospm_freqs_per_perCorr(winNo,1,ii_perCorr)=nospm_freqs_per_perCorr(winNo,1,ii_perCorr)+1;
            lickspm_freq_per_perCorr(winNo,1,ii_perCorr,no_freqs_per_perCorr(winNo,1,ii_perCorr))=lick_freq(winNo,dFF_trNo);

        end
        
        if epochs_per_trial_dFF(dFF_trNo)==4
            %CR, note that CR is saved as 2
            these_dFFs=zeros(1,1,no_traces_win_dFF(winNo,dFF_trNo));
            these_dFFs(1,1,:)=all_win_dFF(winNo,dFF_trNo,1:no_traces_win_dFF(winNo,dFF_trNo));
%             ii_perCorr=find(unique_perCorr==perCorr(trial_dFF(dFF_trNo)),1);
            ii_perCorr=find( ((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)<=5)&((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)>-5),1);
            
            %Save CR
            dFF_per_perCorr(winNo,2,ii_perCorr,no_per_perCorr(winNo,2,ii_perCorr)+1:no_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            no_per_perCorr(winNo,2,ii_perCorr)=no_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs);
            no_freqs_per_perCorr(winNo,2,ii_perCorr)=no_freqs_per_perCorr(winNo,2,ii_perCorr)+1;
            lick_freq_per_perCorr(winNo,2,ii_perCorr,no_freqs_per_perCorr(winNo,2,ii_perCorr))=lick_freq(winNo,dFF_trNo);
            
            %Save S-
            dFFspm_per_perCorr(winNo,2,ii_perCorr,nospm_per_perCorr(winNo,2,ii_perCorr)+1:nospm_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            nospm_per_perCorr(winNo,2,ii_perCorr)=nospm_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs);
            nospm_freqs_per_perCorr(winNo,2,ii_perCorr)=nospm_freqs_per_perCorr(winNo,2,ii_perCorr)+1;
            lickspm_freq_per_perCorr(winNo,2,ii_perCorr,nospm_freqs_per_perCorr(winNo,2,ii_perCorr))=lick_freq(winNo,dFF_trNo);
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==2
            %Miss
            these_dFFs=zeros(1,1,no_traces_win_dFF(winNo,dFF_trNo));
            these_dFFs(1,1,:)=all_win_dFF(winNo,dFF_trNo,1:no_traces_win_dFF(winNo,dFF_trNo));
%             ii_perCorr=find(unique_perCorr==perCorr(trial_dFF(dFF_trNo)),1);
            ii_perCorr=find( ((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)<=5)&((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)>-5),1);
            
            %Save Miss
            dFF_per_perCorr(winNo,3,ii_perCorr,no_per_perCorr(winNo,3,ii_perCorr)+1:no_per_perCorr(winNo,3,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            no_per_perCorr(winNo,3,ii_perCorr)=no_per_perCorr(winNo,3,ii_perCorr)+length(these_dFFs);
            no_freqs_per_perCorr(winNo,3,ii_perCorr)=no_freqs_per_perCorr(winNo,3,ii_perCorr)+1;
            lick_freq_per_perCorr(winNo,3,ii_perCorr,no_freqs_per_perCorr(winNo,3,ii_perCorr))=lick_freq(winNo,dFF_trNo);
            
            %Save S+
            dFFspm_per_perCorr(winNo,1,ii_perCorr,nospm_per_perCorr(winNo,1,ii_perCorr)+1:nospm_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            nospm_per_perCorr(winNo,1,ii_perCorr)=nospm_per_perCorr(winNo,1,ii_perCorr)+length(these_dFFs);
            nospm_freqs_per_perCorr(winNo,1,ii_perCorr)=nospm_freqs_per_perCorr(winNo,1,ii_perCorr)+1;
            lickspm_freq_per_perCorr(winNo,1,ii_perCorr,nospm_freqs_per_perCorr(winNo,1,ii_perCorr))=lick_freq(winNo,dFF_trNo);
        end
        
        if epochs_per_trial_dFF(dFF_trNo)==3
            %FA
            these_dFFs=zeros(1,1,no_traces_win_dFF(winNo,dFF_trNo));
            these_dFFs(1,1,:)=all_win_dFF(winNo,dFF_trNo,1:no_traces_win_dFF(winNo,dFF_trNo));
%             ii_perCorr=find(unique_perCorr==perCorr(trial_dFF(dFF_trNo)),1);
            ii_perCorr=find( ((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)<=5)&((perCorr(trial_dFF(dFF_trNo))-unique_perCorr)>-5),1);
            
            %Save FA
            dFF_per_perCorr(winNo,4,ii_perCorr,no_per_perCorr(winNo,4,ii_perCorr)+1:no_per_perCorr(winNo,4,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            no_per_perCorr(winNo,4,ii_perCorr)=no_per_perCorr(winNo,4,ii_perCorr)+length(these_dFFs);
            no_freqs_per_perCorr(winNo,4,ii_perCorr)=no_freqs_per_perCorr(winNo,4,ii_perCorr)+1;
            lick_freq_per_perCorr(winNo,4,ii_perCorr,no_freqs_per_perCorr(winNo,4,ii_perCorr))=lick_freq(winNo,dFF_trNo);
            
            %Save S-
            dFFspm_per_perCorr(winNo,2,ii_perCorr,nospm_per_perCorr(winNo,2,ii_perCorr)+1:nospm_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs))=...
                these_dFFs;
            nospm_per_perCorr(winNo,2,ii_perCorr)=nospm_per_perCorr(winNo,2,ii_perCorr)+length(these_dFFs);
            nospm_freqs_per_perCorr(winNo,2,ii_perCorr)=nospm_freqs_per_perCorr(winNo,2,ii_perCorr)+1;
            lickspm_freq_per_perCorr(winNo,2,ii_perCorr,nospm_freqs_per_perCorr(winNo,2,ii_perCorr))=lick_freq(winNo,dFF_trNo);
        end
        
    end
    
    %For Hit, CR, Miss and FA plot each point and save data for the anova
        figNo=4+winNo;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
    hold on
    
    data_dFF=[];
    perCorr_dFF=[];
    epoch_dFF=[];
    mindFF=prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1));
    maxdFF=prctile(mean_win_dFF(:),99)+0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1));
    edges=[mindFF:(maxdFF-mindFF)/100:maxdFF];
    rand_offset=5;
    for epochs=1:4
        subplot(2,2,epochs)
        hold on
        kk=0;
        line_perCorr=[];
        line_dFF=[];
        for ii_perCorrs=1:floor(length(unique_perCorr))
            
            perCorrsubsampled(winNo,ii_perCorrs,epochs)=sum(unique_perCorr(ii_perCorrs));
            
            this_dFF=zeros(1,sum(no_per_perCorr(winNo,epochs,ii_perCorrs)));
            tt=0;
            for pp=1:length(ii_perCorrs)
                for jj=1:no_per_perCorr(winNo,epochs,ii_perCorrs(pp))
                    tt=tt+1;
                    this_dFF(1,tt)=dFF_per_perCorr(winNo,epochs,ii_perCorrs(pp),jj);
                end
            end
            
            if ~isempty(this_dFF)
                
                h=histogram(this_dFF(:),edges,'Visible','off');
                normval=h.Values/max(h.Values);
                ii_this_dFF=floor(length(normval)*(this_dFF-mindFF)/(1.2*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1))))+1.0;
                plot_this_dFF=this_dFF((this_dFF>mindFF)&(this_dFF<maxdFF));
                if ~isempty(plot_this_dFF)
                    data_dFF=[data_dFF this_dFF(~isnan(this_dFF))];
                    perCorr_dFF=[perCorr_dFF perCorrsubsampled(winNo,ii_perCorrs,epochs)*ones(1,length(this_dFF(~isnan(this_dFF))))];
                    epoch_dFF=[epoch_dFF epochs*ones(1,length(this_dFF(~isnan(this_dFF))))];
                    pcorr_this_dFF=perCorrsubsampled(winNo,ii_perCorrs,epochs)*ones(1,sum((this_dFF>mindFF)&(this_dFF<maxdFF)))...
                        +rand_offset*normval(ii_this_dFF((this_dFF>mindFF)&(this_dFF<maxdFF)))...
                        .*(rand(1,sum((this_dFF>mindFF)&(this_dFF<maxdFF)))-0.5);
                    switch epochs
                        case 1
                            plot(pcorr_this_dFF,plot_this_dFF,'.r','MarkerSize',3)
                        case 2
                            plot(pcorr_this_dFF,plot_this_dFF,'.b','MarkerSize',3)
                        case 3
                            plot(pcorr_this_dFF,plot_this_dFF,'.c','MarkerSize',3)
                        case 4
                            plot(pcorr_this_dFF,plot_this_dFF,'.m','MarkerSize',3)
                    end
                    plot(perCorrsubsampled(winNo,ii_perCorrs,epochs),mean(plot_this_dFF),'ok','MarkerFaceColor','k','MarkerEdgeColor','k')
                    errorbar(perCorrsubsampled(winNo,ii_perCorrs,epochs),mean(plot_this_dFF),std(plot_this_dFF),'-k','LineWidth',2)
                    handles_out.mean_dFF(winNo,ii_perCorrs,epochs)=mean(plot_this_dFF);
                    handles_out.mean_exists(winNo,ii_perCorrs,epochs)=1;
                    kk=kk+1;
                    line_perCorr(kk)=perCorrsubsampled(winNo,ii_perCorrs,epochs);
                    line_dFF(kk)=mean(plot_this_dFF);
                    all_mean_ndFF(winNo,epochs,ii_perCorrs)=mean(plot_this_dFF);
                    thisCI=[];
                    if length(plot_this_dFF)>1
                        thisCI = bootci(1000, @mean, plot_this_dFF);
                    else
                        thisCI = mean(plot_this_dFF)*ones(2,1);
                    end
                    all_ndFF_CIs(winNo,epochs,ii_perCorrs,1:2)=thisCI;
                    all_ndFF_vaild(winNo,epochs,ii_perCorrs)=1;
                end
                
            end
        end
        plot(line_perCorr,line_dFF,'-k','LineWidth',2)
        minpc=min([minpc min(line_perCorr)]);
        maxpc=max([maxpc max(line_perCorr)]);
    end
 
    subplot(2,2,1)
    hold on
    title(['Hit'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_name{1}='Hit';
    
    
    subplot(2,2,2)
    hold on
    title(['CR'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_name{2}='CR';
    
    subplot(2,2,3)
    hold on
    title(['Miss'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_name{3}='Miss';
    
    subplot(2,2,4)
    hold on
    title(['FA'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_name{4}='FA';
    
    suptitle(['Normalized dF/F for window ' num2str(winNo)])
    
    %Now do the anovan
    %Calculate anovan for inteaction
    [p,tbl,stats]=anovan(data_dFF,{perCorr_dFF epoch_dFF},'varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['\n\nAnovan results for window No %d\n'],winNo)
    fprintf(1, ['p value for anovan normalized dFF per session for percent correct= %d \n'],  p(1));
    fprintf(1, ['p value for anovan normalized dFF per session for epoch= %d \n'],  p(2));
    
    %I am doing this because when I ask for intearaction anovan sometimes
    %gives NaN for percent correct when it should give p=0
    [p,tbl,stats]=anovan(data_dFF,{perCorr_dFF epoch_dFF},'model','interaction','varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['p value for anovan normalized dFF per session for percent correctxepoch= %d \n'],  p(3));

    %For S+ and S- plot each point and save data for the anova
    figNo=7+winNo;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
    hold on
    
    data_dFF=[];
    perCorr_dFF=[];
    epoch_dFF=[];
    mindFF=prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1));
    maxdFF=prctile(mean_win_dFF(:),99)+0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1));
    edges=[mindFF:(maxdFF-mindFF)/100:maxdFF];
    rand_offset=5;
    for epochs=1:2
        
        subplot(2,1,epochs)
        hold on
        kk=0;
        line_perCorr=[];
        line_dFF=[];
        for ii_perCorrs=1:floor(length(unique_perCorr))
            
            perCorrsubsampled(winNo,ii,epochs)=sum(unique_perCorr(ii_perCorrs));
            
            this_dFF=zeros(1,sum(nospm_per_perCorr(winNo,epochs,ii_perCorrs)));
            tt=0;
            
                for jj=1:nospm_per_perCorr(winNo,epochs,ii_perCorrs)
                    tt=tt+1;
                    this_dFF(1,tt)=dFFspm_per_perCorr(winNo,epochs,ii_perCorrs,jj);
                end
           
            
            if ~isempty(this_dFF)
                
                h=histogram(this_dFF(:),edges,'Visible','off');
                normval=h.Values/max(h.Values);
                ii_this_dFF=floor(length(normval)*(this_dFF-mindFF)/(1.2*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1))))+1.0;
                plot_this_dFF=this_dFF((this_dFF>mindFF)&(this_dFF<maxdFF));
                if ~isempty(plot_this_dFF)
                    data_dFF=[data_dFF this_dFF(~isnan(this_dFF))];
                    perCorr_dFF=[perCorr_dFF perCorrsubsampled(winNo,ii_perCorrs,epochs)*ones(1,length(this_dFF(~isnan(this_dFF))))];
                    epoch_dFF=[epoch_dFF epochs*ones(1,length(this_dFF(~isnan(this_dFF))))];
                    pcorr_this_dFF=perCorrsubsampled(winNo,ii_perCorrs,epochs)*ones(1,sum((this_dFF>mindFF)&(this_dFF<maxdFF)))...
                        +rand_offset*normval(ii_this_dFF((this_dFF>mindFF)&(this_dFF<maxdFF)))...
                        .*(rand(1,sum((this_dFF>mindFF)&(this_dFF<maxdFF)))-0.5);
                    switch epochs
                        case 1
                            plot(pcorr_this_dFF,plot_this_dFF,'.r','MarkerSize',3)
                        case 2
                            plot(pcorr_this_dFF,plot_this_dFF,'.b','MarkerSize',3)
                    end
                    plot(perCorrsubsampled(winNo,ii_perCorrs,epochs),mean(plot_this_dFF),'ok','MarkerFaceColor','k','MarkerEdgeColor','k')
                    errorbar(perCorrsubsampled(winNo,ii_perCorrs,epochs),mean(plot_this_dFF),std(plot_this_dFF),'-k','LineWidth',2)
                    handles_out.mean_dFFspm(winNo,ii_perCorrs,epochs)=mean(plot_this_dFF);
                    handles_out.mean_existsspm(winNo,ii_perCorrs,epochs)=1;
                    kk=kk+1;
                    line_perCorr(kk)=perCorrsubsampled(winNo,ii_perCorrs,epochs);
                    line_dFF(kk)=mean(plot_this_dFF);
                    all_mean_ndFF(winNo,epochs,ii_perCorrs)=mean(plot_this_dFF);
                    thisCI=[];
                    if length(plot_this_dFF)>1
                        thisCI = bootci(1000, @mean, plot_this_dFF);
                    else
                        thisCI = mean(plot_this_dFF)*ones(2,1);
                    end
                    all_ndFF_CIs(winNo,epochs,ii_perCorrs,1:2)=thisCI;
                    all_ndFF_vaild(winNo,epochs,ii_perCorrs)=1;
                end
                
            end
        end
        plot(line_perCorr,line_dFF,'-k','LineWidth',2)
        minpc=min([minpc min(line_perCorr)]);
        maxpc=max([maxpc max(line_perCorr)]);
    end
 
    subplot(2,1,1)
    hold on
    title(['S+'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_namespm{1}='S+';
    
    
    subplot(2,1,2)
    hold on
    title(['S-'])
    xlabel('Percent correct')
    ylabel('Normalized dF/F')
    ylim([mindFF maxdFF])
    xlim([minpc-0.1*(maxpc-minpc) maxpc+0.1*(maxpc-minpc)])
    plot([40 95],[0 0],'-k')
    handles_out.epoch_namepsm{2}='S-';

    suptitle(['Normalized dF/F for window ' num2str(winNo)])
    
    %Now do the anovan
    %Calculate anovan for inteaction
    [p,tbl,stats]=anovan(data_dFF,{perCorr_dFF epoch_dFF},'varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['\n\nAnovan results for spm for window No %d\n'],winNo)
    fprintf(1, ['p value for anovan normalized dFF per session for percent correct= %d \n'],  p(1));
    fprintf(1, ['p value for anovan normalized dFF per session for epoch= %d \n'],  p(2));
    
    %I am doing this because when I ask for intearaction anovan sometimes
    %gives NaN for percent correct when it should give p=0
    [p,tbl,stats]=anovan(data_dFF,{perCorr_dFF epoch_dFF},'model','interaction','varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['p value for anovan normalized dFF per session for percent correctxepoch= %d \n'],  p(3));
end

handles_out.perCorr=perCorrsubsampled;
  
%Analysis of lick frequency vs percent correct
max_freq=12;
min_freq=0;
freq_edges=[-0.0001+min_freq:0.333:max_freq];
shift_p=1;
data_lickf=[];
perCorr_lickf=[];
epoch_lickf=[];
mean_freq=zeros(szwins(1),floor(length(unique_perCorr)),4);
CI_freq=zeros(szwins(1),floor(length(unique_perCorr)),4,2);
perCorrsubsampled=zeros(szwins(1),floor(length(unique_perCorr)),4,2);

all_mean_freqs=[];
all_freq_CIs=[];
all_freq_vaild=[];

for winNo=1:szwins(1)
    
    %Plot the licks for Hit, CR, Miss, FA
    figNo=10+winNo;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
    hold on
     
    for epochs=4:-1:1

        for ii_perCorrs=1:floor(length(unique_perCorr))
            
            
            perCorrsubsampled(winNo,ii_perCorrs,epochs)=unique_perCorr(ii_perCorrs);
            this_lickf=zeros(1,sum(no_freqs_per_perCorr(winNo,epochs,ii_perCorrs)));
            tt=0;
           
            for jj=1:no_freqs_per_perCorr(winNo,epochs,ii_perCorrs)
                tt=tt+1;
                this_lickf(1,tt)=lick_freq_per_perCorr(winNo,epochs,ii_perCorrs,jj);
            end
            

            if length(this_lickf)>=1
                valid_freq(winNo,ii_perCorrs,epochs)=1;
                
                mean_freq(winNo,ii_perCorrs,epochs)=mean(this_lickf);
                if length(this_lickf)>=2
                    CI_freq(winNo,ii_perCorrs,epochs,1:2) = bootci(1000, @mean, this_lickf);
                else
                    CI_freq(winNo,ii_perCorrs,epochs,1:2) = [mean(this_lickf) mean(this_lickf)];
                end
            else
                valid_freq(winNo,ii_perCorrs,epochs)=0;
            end
 
        end
        
        these_meanf=zeros(1,floor(length(unique_perCorr)));
        these_meanf(1,:)=mean_freq(winNo,:,epochs);
        these_CIs=zeros(floor(length(unique_perCorr)),2);
        these_CIs(:,1:2)=CI_freq(winNo,:,epochs,1:2);
        valid_points=zeros(1,floor(length(unique_perCorr)));
        valid_points(1,:)=valid_freq(winNo,:,epochs);
        valid_points=logical(valid_points);
        pcsub=zeros(1,floor(length(unique_perCorr)));
        pcsub(1,:)=perCorrsubsampled(winNo,:,epochs);
        
        switch epochs
            case 1
                plot(pcsub(valid_points),these_meanf(valid_points),'-or','LineWidth',3,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
                for ii=1:length(pcsub)
                    if valid_points(ii)==1
                        this_CI=zeros(1,2);
                        this_CI(1,:)=CI_freq(winNo,ii,epochs,1:2);
                        plot([pcsub(ii) pcsub(ii)],this_CI,'-r','LineWidth',3)
                        all_mean_freqs(winNo,epochs,ii)=these_meanf(ii);
                        all_freq_CIs(winNo,epochs,ii,1:2)=this_CI;
                        all_freq_vaild(winNo,epochs,ii)=1;
                    else
                        all_mean_freqs(winNo,epochs,ii)=-1;
                        all_freq_CIs(winNo,epochs,ii,1:2)=[-1 -1];
                        all_freq_vaild(winNo,epochs,ii)=0;
                    end
                    
                end
            case 2
                plot(pcsub(valid_points),these_meanf(valid_points),'-ob','LineWidth',3,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)
                for ii=1:length(pcsub)
                    if valid_points(ii)==1
                        this_CI=zeros(1,2);
                        this_CI(1,:)=CI_freq(winNo,ii,epochs,1:2);
                        plot([pcsub(ii) pcsub(ii)],this_CI,'-b','LineWidth',3)
                        all_mean_freqs(winNo,epochs,ii)=these_meanf(ii);
                        all_freq_CIs(winNo,epochs,ii,1:2)=this_CI;
                        all_freq_vaild(winNo,epochs,ii)=1;
                    else
                        all_mean_freqs(winNo,epochs,ii)=-1;
                        all_freq_CIs(winNo,epochs,ii,1:2)=[-1 -1];
                        all_freq_vaild(winNo,epochs,ii)=0;
                    end
                    
                end
            case 3
                plot(pcsub(valid_points),these_meanf(valid_points),'-oc','LineWidth',3,'MarkerFaceColor','c','MarkerEdgeColor','c','MarkerSize',10)
                for ii=1:length(pcsub)
                    if valid_points(ii)==1
                        this_CI=zeros(1,2);
                        this_CI(1,:)=CI_freq(winNo,ii,epochs,1:2);
                        plot([pcsub(ii) pcsub(ii)],this_CI,'-c','LineWidth',3)
                        all_mean_freqs(winNo,epochs,ii)=these_meanf(ii);
                        all_freq_CIs(winNo,epochs,ii,1:2)=this_CI;
                        all_freq_vaild(winNo,epochs,ii)=1;
                    else
                        all_mean_freqs(winNo,epochs,ii)=-1;
                        all_freq_CIs(winNo,epochs,ii,1:2)=[-1 -1];
                        all_freq_vaild(winNo,epochs,ii)=0;
                    end
                    
                end
            case 4
                plot(pcsub(valid_points),these_meanf(valid_points),'-om','LineWidth',3,'MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10)
                for ii=1:length(pcsub)
                    if valid_points(ii)==1
                        this_CI=zeros(1,2);
                        this_CI(1,:)=CI_freq(winNo,ii,epochs,1:2);
                        plot([pcsub(ii) pcsub(ii)],this_CI,'-m','LineWidth',3)
                        all_mean_freqs(winNo,epochs,ii)=these_meanf(ii);
                        all_freq_CIs(winNo,epochs,ii,1:2)=this_CI;
                        all_freq_vaild(winNo,epochs,ii)=1;
                    else
                        all_mean_freqs(winNo,epochs,ii)=-1;
                        all_freq_CIs(winNo,epochs,ii,1:2)=[-1 -1];
                        all_freq_vaild(winNo,epochs,ii)=0;
                    end
                end
        end

                     %anovan data
        for ii_perCorr=1:length(unique_perCorr)

            this_lickf=zeros(1,no_freqs_per_perCorr(winNo,epochs,ii_perCorr));
            if no_freqs_per_perCorr(winNo,epochs,ii_perCorr)>0
                this_lickf(1,:)=lick_freq_per_perCorr(winNo,epochs,ii_perCorr,1:no_freqs_per_perCorr(winNo,epochs,ii_perCorr));
                
                if length(this_lickf)>=1
                    data_lickf=[data_lickf this_lickf];
                    perCorr_lickf=[perCorr_lickf unique_perCorr(ii_perCorr)*ones(1,length(this_lickf))];
                    epoch_lickf=[epoch_lickf epochs*ones(1,length(this_lickf))];
                end
            end
        end
 
    end
    
    title(['Lick frequency vs percent correct for window No ' num2str(winNo)])
    xlabel('Percent correct')
    ylabel('Lick frequency(Hz)')
    ylim([0 1.2*max(CI_freq(:))])
    xlim([pcsub(1)-(pcsub(2)-pcsub(1))/2 pcsub(end)+(pcsub(2)-pcsub(1))/2])
        
    %Now do the anovan
    %Calculate anovan for inteaction
    [p,tbl,stats]=anovan(data_lickf,{perCorr_lickf epoch_lickf},'varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['\n\nAnovan results for window No %d\n'],winNo)
    fprintf(1, ['p value for anovan lick rate per session for percent correct= %d \n'],  p(1));
    fprintf(1, ['p value for anovan lick rate per session for epoch= %d \n'],  p(2));
    
    %I am doing this because when I ask for intearaction anovan sometimes
    %gives NaN for percent correct when it should give p=0
    %This anovan is a first attempt, but it is likely not granted, we hav
    %to use the same ROIs in the entire experiment
    [p,tbl,stats]=anovan(data_lickf,{perCorr_lickf epoch_lickf},'model','interaction','varnames',{'percent_correct','epoch'},'display','off');
    fprintf(1, ['p value for anovan lick rate per session for percent correctxepoch= %d \n'],  p(3));
    
    
end

pfff=1;

figNo=14;
try
    close figNo
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
hold on

afreqs=[];
adFFs=[];
numpts=0;

for winNo=1:szwins(1)
    
    
    for epochs=1:4
        for ii=1:floor(length(unique_perCorr))
            
            if (all_freq_vaild(winNo,epochs,ii)==1)&(all_ndFF_vaild(winNo,epochs,ii)==1)
                numpts=numpts+1;
                afreqs(numpts)=all_mean_freqs(winNo,epochs,ii);
                adFFs(numpts)=all_mean_ndFF(winNo,epochs,ii);
                switch epochs
                    case 1
                        plot(all_mean_freqs(winNo,epochs,ii),all_mean_ndFF(winNo,epochs,ii),'or')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_freq_CIs(winNo,epochs,ii,:);
                        plot(theseCIs,[all_mean_ndFF(winNo,epochs,ii) all_mean_ndFF(winNo,epochs,ii)],'-r')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_ndFF_CIs(winNo,epochs,ii,:);
                        plot([all_mean_freqs(winNo,epochs,ii) all_mean_freqs(winNo,epochs,ii)],theseCIs,'-r')
                    case 2
                        plot(all_mean_freqs(winNo,epochs,ii),all_mean_ndFF(winNo,epochs,ii),'ob')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_freq_CIs(winNo,epochs,ii,:);
                        plot(theseCIs,[all_mean_ndFF(winNo,epochs,ii) all_mean_ndFF(winNo,epochs,ii)],'-b')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_ndFF_CIs(winNo,epochs,ii,:);
                        plot([all_mean_freqs(winNo,epochs,ii) all_mean_freqs(winNo,epochs,ii)],theseCIs,'-b')
                    case 3
                        plot(all_mean_freqs(winNo,epochs,ii),all_mean_ndFF(winNo,epochs,ii),'oc')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_freq_CIs(winNo,epochs,ii,:);
                        plot(theseCIs,[all_mean_ndFF(winNo,epochs,ii) all_mean_ndFF(winNo,epochs,ii)],'-c')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_ndFF_CIs(winNo,epochs,ii,:);
                        plot([all_mean_freqs(winNo,epochs,ii) all_mean_freqs(winNo,epochs,ii)],theseCIs,'-c')
                    case 4
                        plot(all_mean_freqs(winNo,epochs,ii),all_mean_ndFF(winNo,epochs,ii),'om')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_freq_CIs(winNo,epochs,ii,:);
                        plot(theseCIs,[all_mean_ndFF(winNo,epochs,ii) all_mean_ndFF(winNo,epochs,ii)],'-m')
                        theseCIs=zeros(1,2);
                        theseCIs(1,1:2)=all_ndFF_CIs(winNo,epochs,ii,:);
                        plot([all_mean_freqs(winNo,epochs,ii) all_mean_freqs(winNo,epochs,ii)],theseCIs,'-m')
                end
            end
        end
        
    end
    title('Normalized dFF vs. lick frequency')
    xlabel('Lick frequency (Hz)')
    ylabel('Normalized dFF')
end

[rho,pval]=corr(afreqs',adFFs');

fprintf(1, ['\ncorrelation coefficient %d, p value %d\n'],  rho,pval);

%plot number of ROIs
figNo=15;
try
    close figNo
catch
end

hFig = figure(figNo);
set(hFig, 'units','normalized','position',[.25 .1 .7 .7])
hold on
plot([1:caimanhandles.caimandr_choices.no_files],noROIs,'-ok')
ylim([0 1.2*max(noROIs)])
title('Number of ROIs per file')
xlabel('File no')
ylabel('No ROIs')

%Plot norm dFF for each window
if caimanhandles.caimandr_choices.start_reversal>0
    tr_reversal=first_num_odor_trials(caimanhandles.caimandr_choices.start_reversal);
else
    tr_reversal=0;
end

for winNo=1:szwins(1)
    figNo=15+winNo;
    try
        close figNo
    catch
    end
    
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.25 .05+0.3*(winNo-1) .5 .25])
    hold on
    
    for dFF_trNo=1:num_odor_trials_dFF
        
        %Plot odor 1 (S+ forward)
        subplot(2,1,1)
        hold on
        if dFF_trNo<tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
             if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
        
        
        %Plot odor 2 (S- forward)
        subplot(2,1,2)
        hold on
        if dFF_trNo>=tr_reversal
            %If Hit or Miss
            if (epochs_per_trial_dFF(dFF_trNo)==1)||(epochs_per_trial_dFF(dFF_trNo)==2)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
        else
            %If CR or FA
             if (epochs_per_trial_dFF(dFF_trNo)==3)||(epochs_per_trial_dFF(dFF_trNo)==4)
                %Confidence interval
                this_CI=zeros(1,2);
                this_CI(1,1:2)=CI_win_dFF(winNo,dFF_trNo,:);
                plot([trial_dFF(dFF_trNo) trial_dFF(dFF_trNo)],this_CI,'-k')
                
                if epochs_per_trial_dFF(dFF_trNo)==1
                    %Hit
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'or')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==4
                    %CR
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'ob')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==2
                    %Miss
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'oc')
                end
                
                if epochs_per_trial_dFF(dFF_trNo)==3
                    %FA
                    plot(trial_dFF(dFF_trNo), mean_win_dFF(winNo,dFF_trNo),'om')
                end
            end
            
        end
    end
    
    for spno=1:2
        subplot(2,1,spno)
        if isfield(caimanhandles.caimandr_choices,'start_reversal')
            filNum=caimanhandles.caimandr_choices.start_reversal;
            if filNum>0
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Reversal','Color','k','FontSize',18)
            end
        end
        
         if isfield(caimanhandles.caimandr_choices,'start_gogo')
            filNum=caimanhandles.caimandr_choices.start_gogo;
            if filNum>0
                plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k','LineWidth',4)
                text(first_num_odor_trials(filNum)+2,prctile(mean_win_dFF(:),1)+0.9*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)),'Go-go','Color','k','FontSize',18)
            end
         end
         
         if isfield(caimanhandles.caimandr_choices,'start_session')
             if length(caimanhandles.caimandr_choices.start_session)>=2
                 for sessionNo=2:length(caimanhandles.caimandr_choices.start_session)
                     filNum=caimanhandles.caimandr_choices.start_session(sessionNo);
                     plot([first_num_odor_trials(filNum) first_num_odor_trials(filNum)],[prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))],'-k')
                     %             text(first_num_odor_trials(filNum)+2,80,'Reversal','Color','k','FontSize',18)
                 end
             end
         end
         
         plot([1 num_odor_trials_dFF],[0 0], '-k')
         
         if spno==1
        title(['Odor 1 (S+ forward)'])
        else
          title(['Odor 2 (S- forward)'])  
        end
        xlabel('Trial number')
        ylabel('Normalized dF/F')
        ylim([prctile(mean_win_dFF(:),1)-0.1*(prctile(mean_win_dFF(:),99)-prctile(mean_win_dFF(:),1)) 1.2*prctile(mean_win_dFF(:),99)-0.1*(prctile(mean_win_dFF(:),99)+prctile(mean_win_dFF(:),1))])
    end
    suptitle(['Normalized dF/F for window No ' num2str(winNo) ' Hit(red) Miss(cyan) FA(magenta) CR(blue)'])
end

%Show all responses for each odor in a pseudocolor for two ranges of trials

save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName],'handles_out')
 
pffft=1
