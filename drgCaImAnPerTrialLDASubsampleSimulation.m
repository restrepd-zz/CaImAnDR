%drgCaImAnPerTrialLDASubsampleSimulation
%This code simulates activity of multiple components with different
%covariant activity to analyze the dependence of LDA prediction
%on the number of components
clear all
close all

load_lda_temp=1;

if load_lda_temp==1
    %This is here to load actual data to troubleshoot how to simulate the data
    load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/20180917_mmPVG04_Cerebellum_PCA_ldatemp.mat')
    pffft=1;
    
    fileNo=1;
    
    %Note: I will only use the data in file 1
    last_trNo=find(all_lda_fileNo==1,1,'last');
    
    %Trim variables
    all_lda_fileNo=all_lda_fileNo(1,1:last_trNo);
    old_all_lda_events=all_lda_events;
    all_lda_events=[];
    for trNo=1:last_trNo
        all_lda_events{trNo}=old_all_lda_events{trNo};
        if strcmp(all_lda_events{trNo},'S+')
            mean_dFF(trNo)=1;
        else
            mean_dFF(trNo)=0;
        end
    end
    all_lda_no_comp=all_lda_no_comp(1,1:last_trNo);
    dFF_trial_mask=dFF_trial_mask(1,1:last_trNo);
    perCorr=perCorr(1,1:last_trNo);
    
    winNo=1;
    first_timepoint=find(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,1),1,'first');
    last_timepoint=find(time_to_eventLDA>=caimanhandles.caimandr_choices.wins(winNo,2),1,'first');
    no_timepoints=length(time_to_eventLDA(first_timepoint:last_timepoint));
    
    all_lda_input_timecourse=all_lda_input_timecourse(first_timepoint:last_timepoint,1:all_lda_no_comp(1),1:last_trNo);
    
    time_to_eventLDA=0:0.1:0.6; %Using the window defined above this results in using six time points
    
    %Output file name
    caimanhandles.caimandr_choices.outFileName='simulation_pointtwocov.mat';
    caimanhandles.caimandr_choices.outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/Simulations/';
    caimanhandles.caimandr_choices.start_reversal=200;
    
    %Enter the windows for evaluation of dF/F
    caimanhandles.caimandr_choices.wins=[0 0.5];
    
    figNo=0;
    first_num_odor_trials=[1 41];
    number_of_replicates=10;
    no_trial_windows=3;
    num_odor_trials_dFF=last_trNo;
    
    supertitle_description{3}='dF/F LDA analysis for with file 1 of 20180917_mmPVG04';
    
    handles_par=[];
    handles_par.time_to_eventLDA=time_to_eventLDA;
    handles_par.perCorr=perCorr;
    min_trials=20;
    
    
else
    no_trial_windows=3;
    
    %If you want to input COV make this variable =1, and enter input pathname and filename
    % import_cov=0;
    
    % if import_cov==1
    %     COV_FileName='1-Registered_CalmAn_COV.mat';
    %     COV_PathName='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/';
    %     load([COV_PathName COV_FileName])
    % else
    %Set covariance
    handles_par(no_trial_windows).COV=0.2; %Covariance
    % end
    
    %Output file name
    caimanhandles.caimandr_choices.outFileName='simulation_pointtwocov.mat';
    caimanhandles.caimandr_choices.outPathName='/Users/restrepd/Documents/Projects/MOM slidebook/Simulations/';
    caimanhandles.caimandr_choices.start_reversal=200;
    
    %Enter the windows for evaluation of dF/F
    caimanhandles.caimandr_choices.wins=[0 0.5];
    
    %Data for drgCaImAnLDAforSubsample
    %Please note that I am using all of these input values so that I can use
    %the same function for this code and for drgCaImAnBatchPerSessionPerTrialLDASubsample
    no_trials=80;
    handles_par(no_trial_windows).no_trials=no_trials;
    % if import_cov==1
    %     no_comps=size(COV,1);
    % else
    no_comps=250;
    % end
    handles_par(no_trial_windows).no_comps=no_comps;
    min_trials=20;
    dFF_trial_mask=ones(1,no_trials);
    for ii=1:no_trials
        if rem(ii,2)==0
            mean_dFF(ii)=1;
            all_lda_events{ii}='S+';
        else
            mean_dFF(ii)=0;
            all_lda_events{ii}='S-';
        end
    end
    handles_par(no_trial_windows).mean_dFF=mean_dFF;
    
    num_odor_trials_dFF=no_trials;
    time_to_eventLDA=0:0.1:0.6; %Using the window defined above this results in using six time points
    all_lda_no_comp=no_comps*ones(1,no_trials);
    perCorr=90*ones(1,no_trials);
    handles_par(no_trial_windows).time_to_eventLDA=time_to_eventLDA;
    handles_par(no_trial_windows).perCorr=perCorr;
    all_lda_fileNo=ones(1,no_trials);
    number_of_replicates=10;
    first_num_odor_trials=[1 41];
    figNo=0;
    
    %all_lda_input_timecourse(timepoints,no_comps,no_trials)
    
    all_lda_input_timecourse=zeros(length(time_to_eventLDA),no_comps,no_trials);
    
    % diff=1;
    diff=1*(1+0.2*randn(1,no_comps));
    for ii_t=1:length(time_to_eventLDA)
        
        %     if import_cov==1
        %         Signal = mvnrnd(zeros(1,no_comps),COV,no_trials); % Make a set of signals with a known cov.
        %     else
        if handles_par(no_trial_windows).COV~=0
            SIGMA = ones(no_comps,no_comps)*handles_par(no_trial_windows).COV; % Make a covariance matrix.
            SIGMA(diag(ones(1,no_comps))==1)= 1; % Set the self-covariance to 1.
            Signal = mvnrnd(ones(1,no_comps),SIGMA,no_trials); % Make a set of signals with a known cov.
        else
            Signal=randn(no_trials,no_comps);
        end
        %     end
        
        %     if import_cov==1
        %         for trNo=1:no_trials
        %             if mean_dFF(trNo)==1
        %                 Signal(trNo,:)=Signal(trNo,:)+dFF_change;
        %             end
        %         end
        %     else
        %     for trNo=1:no_trials
        %         Signal(trNo,:)=Signal(trNo,:)+diff*mean_dFF(trNo);
        %     end
        for trNo=1:no_trials
            for ncomps=1:no_comps
                Signal(trNo,ncomps)=Signal(trNo,ncomps)+diff(ncomps)*mean_dFF(trNo);
            end
        end
        %     end
        
        all_lda_input_timecourse(ii_t,:,:)=Signal';
        Dim_out(ii_t) = nansum(eig(cov(Signal)))^2/nansum(eig(cov(Signal)).^2);
    end
    
    % if import_cov==1
    %     fprintf(1, 'Covariance imported, Dimensionality = %d\n', mean(Dim_out));
    % else
    fprintf(1, 'COV= %d Dimensionality = %d\n',handles_par(no_trial_windows).COV, mean(Dim_out));
    % end
    
    %Note: number 3 will be used
    supertitle_description{1}='dF/F LDA analysis for trials with percent correct <65';
    supertitle_description{2}='dF/F LDA analysis for trials with percent correct >=65&<80';
    supertitle_description{3}='dF/F LDA analysis for simulation with independent components';
end

%Show the data for two components
for figNo=1:16
    try
        close(figNo)
    catch
    end
    
    figure(figNo)
    hold on
    
    %Show 16 relationships
    ii_comp1=randperm(all_lda_no_comp(1),16);
    
    comps_differ=0;
    while comps_differ==0
        ii_comp2=randperm(all_lda_no_comp(1),16);
        comps_differ=1;
        for ii=1:16
            if ii_comp1(11)==ii_comp2(ii)
                comps_differ=1;
            end
            for jj=1:16
                if ii~=jj
                    if (ii_comp1(ii)==ii_comp1(jj))&(ii_comp2(ii)==ii_comp2(jj))
                        comps_differ=1;
                    end
                end
            end
        end
    end
    
    no_t=size(all_lda_input_timecourse,1);
    for sub_ii=1:16
        subplot(4,4,sub_ii)
        hold on
        for ii_t=1:no_t
            s_minus1=ones(1,length(mean_dFF)-sum(mean_dFF));
            s_minus1(1,:)=all_lda_input_timecourse(ii_t,ii_comp1(sub_ii),~logical(mean_dFF));
            s_minus2=ones(1,length(mean_dFF)-sum(mean_dFF));
            s_minus2(1,:)=all_lda_input_timecourse(ii_t,ii_comp2(sub_ii),~logical(mean_dFF));
            plot(s_minus1,s_minus2,'ob','MarkerFaceColor','b','MarkerSize',2)
            
            s_plus1=ones(1,sum(mean_dFF));
            s_plus1(1,:)=all_lda_input_timecourse(ii_t,ii_comp1(sub_ii),logical(mean_dFF));
            s_plus2=ones(1,sum(mean_dFF));
            s_plus2(1,:)=all_lda_input_timecourse(ii_t,ii_comp2(sub_ii),logical(mean_dFF));
            plot(s_plus1,s_plus2,'or','MarkerFaceColor','r','MarkerSize',2)
        end
        title([num2str(ii_comp1(sub_ii)) ' x ' num2str(ii_comp2(sub_ii)) ])
    end
end

% if import_cov==1
%     title(['Covariance imported, Dimensionality = ' num2str(mean(Dim_out))])
% else

% end




handles_par=drgCaImAnLDAforSubsample(dFF_trial_mask,time_to_eventLDA,...
    all_lda_no_comp,min_trials,handles_par,caimanhandles,no_trial_windows,all_lda_fileNo...
    ,number_of_replicates,first_num_odor_trials,num_odor_trials_dFF,perCorr,all_lda_events,...
    all_lda_input_timecourse,figNo,supertitle_description)

save([caimanhandles.caimandr_choices.outPathName caimanhandles.caimandr_choices.outFileName(1:end-4) '_ldasubsim.mat'],'handles_par')

pffft=1;
