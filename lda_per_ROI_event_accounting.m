%lda_per_event_all_files
%This program generates a summary figure for the per event LDA analysis for Fig. 5
close all
clear all

handles.no_ts=0;

%mmPVG04 20180917_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_LDA_events_account_events.mat')

for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=1;
end
 
%mmPVG04 20180910_mmPVG04_Cerebellum
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_new_out_account_events.mat')
for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=2;
end

%mmPVG05 20181017
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_mmPVG05_Cerebellum_out_account_events.mat')
for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=3;
end

%mmPVG02 20180515_18
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_lda_new_account_events.mat')
for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=4;
end

%mmG7f09 20180702_05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_lda_sum_account_events.mat')
for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=5;
end

%mmG06  20180419
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_account_events.mat')
for nts=1:handles_per_comp.no_tseries_included
    handles.no_ts=handles.no_ts+1;
    handles.tseries(handles.no_ts)=handles_per_comp.tseries(nts);
    handles.sessionNo(handles.no_ts)=6;
end
 
%Do the calculation of responses
percent_hits=[];
no_Hit=0;
percent_CRs=[];
no_CR=0;
percent_FAs=[];
no_FA=0;
percent_Miss=[];
no_Miss=0;
percent_Hit_per_Miss=[];
no_Hit_per_Miss=0;
precent_Hit_per_FA=[];
no_Hit_per_FA=0;
no_miss_responses=[];
no_FA_responses=[];
no_ROIs_Miss=[];
no_ROIs_FA=[];
no_ROIs_Hit=[];
sessionNo=[];
fileNo=[];
percent_FA=[];
percent_Miss=[];
session_FA=[];
session_Miss=[];
 
for nts=1:handles.no_ts
    if (~isempty(handles.tseries(nts).noHit_responses))&(~isempty(handles.tseries(nts).noCR_responses))
        
        no_Hit=no_Hit+1;
        percent_hits(no_Hit)=100*handles.tseries(nts).noHit_responses/handles.tseries(nts).noROIs;
        
        no_CR=no_CR+1;
        percent_CRs(no_CR)=100*handles.tseries(nts).noCR_responses/handles.tseries(nts).noROIs;
        
        sessionNo(no_Hit)=handles.sessionNo(nts);
        fileNo(no_Hit)=handles.tseries(nts).fileNo;
        no_ROIs_Hit(no_Hit)=handles.tseries(nts).noROIs;
        
        %If there are Miss trials
        if (handles.tseries(nts).noMiss_trials>0)&(handles.tseries(nts).noMissHit_responses>0)
            no_Miss=no_Miss+1;
            percent_Miss(no_Miss)=100*handles.tseries(nts).noMiss_responses/handles.tseries(nts).noROIs;
            percent_Hit_per_Miss(no_Miss)=100*handles.tseries(nts).noMissHit_responses/handles.tseries(nts).noMiss_responses;
            no_miss_responses(no_Miss)=handles.tseries(nts).noMiss_responses;
            session_Miss(no_Miss)=handles.sessionNo(nts);
            no_ROIs_Miss(no_Miss)=handles.tseries(nts).noROIs;
        end
        
        %If there are FA trials
        if (handles.tseries(nts).noFA_trials>0)&(handles.tseries(nts).noFAHit_responses>0)
            no_FA=no_FA+1;
            percent_FA(no_FA)=100*handles.tseries(nts).noFA_responses/handles.tseries(nts).noROIs;
            percent_Hit_per_FA(no_FA)=100*handles.tseries(nts).noFAHit_responses/handles.tseries(nts).noFA_responses;
            no_FA_responses(no_FA)=handles.tseries(nts).noFA_responses;
            session_FA(no_FA)=handles.sessionNo(nts);
            no_ROIs_FA(no_FA)=handles.tseries(nts).noROIs;
        end
    end
end

fprintf(1,'\n\nPercent of ROIs for Hits mean = %d, SD = %d, n= %d\n',mean(percent_hits),std(percent_hits),length(percent_hits))
fprintf(1,'Percent of ROIs for CRs mean = %d, SD = %d, n= %d\n',mean(percent_CRs),std(percent_CRs),length(percent_CRs))
fprintf(1,'Percent of ROIs for Misss mean = %d, SD = %d, n= %d\n',mean(percent_Miss),std(percent_Miss),length(percent_Miss))
fprintf(1,'Percent of ROIs for FAs mean = %d, SD = %d, n= %d\n',mean(percent_FA),std(percent_FA),length(percent_FA))


fprintf(1,'\nPercent of Hits per Misss mean = %d, SD = %d, n= %d\n',mean(percent_Hit_per_Miss),std(percent_Hit_per_Miss),length(percent_Hit_per_Miss))
fprintf(1,'Percent of Hits per FAs mean = %d, SD = %d, n= %d\n',mean(percent_Hit_per_FA),std(percent_Hit_per_FA),length(percent_Hit_per_FA))

pfft=1;