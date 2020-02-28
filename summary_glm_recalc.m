%summary_glm_recalc
%This program generates a summary figure for the rho of DtdFF vs Dtlick
%rate for Fig. 6
close all
clear all

kinematics=zeros(3,6);
no_kinematics=0;
licks=zeros(3,6);
no_licks=0;
odors=zeros(3,6);
no_odor=0;
reinforcement_history=zeros(3,6);
no_rh=0;
accuracy=zeros(3,6);
no_accuracy=0;
perCorr=zeros(3,6);
no_pc=0;
perVar=zeros(3,6);
no_pvar=0;
discriminant_correct_shuffled=[];
no_included_files=0;

for fileNo=1:6
    
    switch fileNo
        case 1
            %mmPVG04
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum/20180910_mmPVG04_Cerebellumopflow_slopes.mat')
            load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellumopflow_slopestest_glmtest.mat')
        case 2
            %mmPVG02
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180518_mmPVG02_Cerebellum/2018018_mmPVG02_Cerebellum_OpFlowout_slopes.mat')
            load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180518_mmPVG02_Cerebellum new analysis/2018018_mmPVG02_Cerebellum_OpFlowout_slopestest_glmtest.mat')
        case 3
            %mmG7f09
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09_Cerebellum/20180702_05_mmG7f09-Cerebellumbatch_slopes_glm.mat')
            load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellumbatch_slopestest_glmtest.mat')
        case 4
            %mmPVG05
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum/20181017_mmPVG05_Cerebellum_out_slopes_glm.mat')
                load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_mmPVG05_Cerebellum_out_slopestest_glmtest.mat')
        case 5
            %mmPVG04
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum/20180917_mmPVG04_Cerebellum_PCA_slopes.mat')
             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180917_mmPVG04_Cerebellum new analysis/20180917_mmPVG04_Cerebellum_PCA_slopestest_glmtest.mat')
        case 6
            %mmG06
%             load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum/20180419_mmG06_cerebellumPCAevents_slopes_glm.mat')
            load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG06/20180419_mmG06_cerebellum new analysis/20180419_mmG06_cerebellumPCAevents_slopestest_glmtest.mat')
    end
    
    if isfield(handles_outs,'models')
        %Note: if "models" does not exist the file was not processed with
        %drgCaImAn_glm
        fprintf(1, ['\n\nFile No %d did have the model field\n'],fileNo)
        no_included_files=no_included_files+1;
        for ii=2:4
            %Note: In drgCaImAn_glm I saved the three windows in ii+1 (i.e. 2 to 4)
            
            %Do kinematics
            %Find accelearation
            in_found=0;
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'acceleration')
                    found_jj=jj;
                    kin_found=1;
                end
            end
            if found_jj>0
                kinematics(ii-1,no_included_files)=kinematics(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            %Find velocity
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'velocities')
                    found_jj=jj;
                    kin_found=1;
                end
            end
            if found_jj>0
                kinematics(ii-1,no_included_files)=kinematics(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            if (kin_found==1)&(ii==2)
                no_kinematics=no_kinematics+1;
            end
            
            %Do licks
            %Find lick_freq
            lick_found=0;
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'lick_freq')
                    found_jj=jj;
                    lick_found=1;
                end
            end
            if found_jj>0
                licks(ii-1,no_included_files)=licks(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            %Find lick_derivatives
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'lick_derivatives')
                    found_jj=jj;
                    lick_found=1;
                end
            end
            if found_jj>0
                licks(ii-1,no_included_files)=licks(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            if (lick_found==1)&(ii==2)
                no_licks=no_licks+1;
            end
            
            if (ii==3)&(fileNo==5)
                pffft=1;
            end
            
            %Do odors
            %Find splus_odor
            found_jj=-1;
            odor_found=0;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'splus_odor')
                    found_jj=jj;
                    odor_found=1;
                end
            end
            if found_jj>0
                odors(ii-1,no_included_files)=odors(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
%             %Find sminus_odor
%             found_jj=-1;
%             for jj=1:length(handles_outs.models(1).predictor_names)
%                 if strcmp(handles_outs.models(1).predictor_names{jj},'sminus_odor')
%                     found_jj=jj;
%                     odor_found=1;
%                 end
%             end
%             if found_jj>0
%                 odors(ii-1,no_included_files)=odors(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
%             end
            
            if (odor_found==1)&(ii==2)
                no_odor=no_odor+1;
            end
            
            %Find reinforcement_history
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'reinforcement_history')
                    found_jj=jj;
                end
            end
            if found_jj>0
                if ii==2
                    no_rh=no_rh+1;
                end
                reinforcement_history(ii-1,no_included_files)=reinforcement_history(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            %Find accuracy
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'accuracy')
                    found_jj=jj;
                end
            end
            if found_jj>0
                if ii==2
                    no_accuracy=no_accuracy+1;
                end
                accuracy(ii-1,no_included_files)=accuracy(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            
            %Find perCorr
            found_jj=-1;
            for jj=1:length(handles_outs.models(1).predictor_names)
                if strcmp(handles_outs.models(1).predictor_names{jj},'perCorr')
                    found_jj=jj;
                end
            end
            if found_jj>0
                if ii==2
                    no_pc=no_pc+1;
                end
                perCorr(ii-1,no_included_files)=perCorr(ii-1,no_included_files)+handles_outs.models(ii).contributions(found_jj);
            end
            fprintf(1, ['For file number %d the sum of contributions is %d for window %d\n'],fileNo,sum(handles_outs.models(ii).contributions),ii-1)
            
            %Find the percent of the variance that is explained by the model
            
            if ii==2
                no_pvar=no_pvar+1;
            end
            perVar(ii-1,no_included_files)=100*handles_outs.models(ii-1).mdl.SSE/handles_outs.models(ii-1).mdl.SST;
            
        end
        
%         for ii=1:5
%             discriminant_correct(ii,no_included_files)=handles_outs.LDA(ii).discriminant_correct;
%             discriminant_correct_shuffled=[discriminant_correct_shuffled handles_outs.LDA(ii).discriminant_correct_shuffled];
%         end
        
    else
        fprintf(1, ['\nFile No %d did not have the model field\n'],fileNo)
    end
    
end


%Plot a bar graph showing the number of significant inclusions as a
%predictive variable
figure(1)
hold on

%kinematics
bar(1,100*no_kinematics/6,'b')
%licks
bar(2,100*no_licks/6,'b')

%Odors
bar(3,100*no_odor/6,'b')

%Reinforcement history
bar(4,100*no_rh/6,'b')

%Accuracy
bar(5,100*no_accuracy/6,'b')

%perCorr
bar(6,100*no_pc/6,'b')

ylabel('% included in glm')
xticks([1 2 3 4 5 6])
xticklabels({'Kinematics','Licks','Odors','Past','Accuracy','% correct'})
ylim([0 120])

%Now plot the contributions
hFig=figure(2)
set(hFig, 'units','normalized','position',[.05 .05 .5 .8])

%kinematics
subplot(3,2,1)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=kinematics(ii,:);
    switch ii
        case 1
            bar(ii,mean(kinematics(ii,:)),'m')
            contributions(ii_cont).description='kinematics pre-odor';
        case 2
            bar(ii,mean(kinematics(ii,:)),'b')
            contributions(ii_cont).description='kinematics odor';
        case 3
            bar(ii,mean(kinematics(ii,:)),'r')
            contributions(ii_cont).description='kinematics outcome';
    end
            
    CIkin = bootci(1000, {@mean, kinematics(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot(ii*ones(1,6),kinematics(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);

end
title('Kinematics')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for kinematics\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

%licks
subplot(3,2,2)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=licks(ii,:);
    switch ii
        case 1
            bar(ii,mean(licks(ii,:)),'m')
            contributions(ii_cont).description='licks pre-odor';
        case 2
            bar(ii,mean(licks(ii,:)),'b')
            contributions(ii_cont).description='licks odor';
        case 3
            bar(ii,mean(licks(ii,:)),'r')
            contributions(ii_cont).description='licks outcome';
    end
            
    CIkin = bootci(1000, {@mean, licks(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot((ii)*ones(1,6),licks(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
end
title('Licks')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for licks\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

%odors
subplot(3,2,3)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=odors(ii,:);
    switch ii
        case 1
            bar(ii,mean(odors(ii,:)),'m')
            contributions(ii_cont).description='odors pre-odor';
        case 2
            bar(ii,mean(odors(ii,:)),'b')
            contributions(ii_cont).description='odors odor';
        case 3
            bar(ii,mean(odors(ii,:)),'r')
            contributions(ii_cont).description='odors outcome';
    end
            
    CIkin = bootci(1000, {@mean, odors(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot((ii)*ones(1,6),odors(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
end
title('Odors')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for odors\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

%Reinforcement history
subplot(3,2,4)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=reinforcement_history(ii,:);
    switch ii
        case 1
            bar(ii,mean(reinforcement_history(ii,:)),'m')
            contributions(ii_cont).description='reinforcement_history pre-odor';
        case 2
            bar(ii,mean(reinforcement_history(ii,:)),'b')
            contributions(ii_cont).description='reinforcement_history odor';
        case 3
            bar(ii,mean(reinforcement_history(ii,:)),'r')
            contributions(ii_cont).description='reinforcement_history outcome';
    end
            
    CIkin = bootci(1000, {@mean, reinforcement_history(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot((ii)*ones(1,6),reinforcement_history(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
end
title('Reinforcment history')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for reinforcement history\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

%Accuracy
subplot(3,2,5)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=accuracy(ii,:);
    switch ii
        case 1
            bar(ii,mean(accuracy(ii,:)),'m')
            contributions(ii_cont).description='accuracy pre-odor';
        case 2
            bar(ii,mean(accuracy(ii,:)),'b')
            contributions(ii_cont).description='accuracy odor';
        case 3
            bar(ii,mean(accuracy(ii,:)),'r')
            contributions(ii_cont).description='accuracy outcome';
    end
            
    CIkin = bootci(1000, {@mean, accuracy(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot((ii)*ones(1,6),accuracy(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
end
title('Accuracy')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for acuracy\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

%perCorr
subplot(3,2,6)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=perCorr(ii,:);
    switch ii
        case 1
            bar(ii,mean(perCorr(ii,:)),'m')
            contributions(ii_cont).description='perCorr pre-odor';
        case 2
            bar(ii,mean(perCorr(ii,:)),'b')
            contributions(ii_cont).description='perCorr odor';
        case 3
            bar(ii,mean(perCorr(ii,:)),'r')
            contributions(ii_cont).description='perCorr outcome';
    end
            
    CIkin = bootci(1000, {@mean, perCorr(ii,:)})';
    plot([ii ii],CIkin,'-k')
    plot((ii)*ones(1,6),perCorr(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
end
title('Percent correct')
ylabel('% Contribution to glm')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])

fprintf(1, ['\n\nttest or rnaksum test for percent correct\n'])
[output_data] = drgMutiRanksumorTtest(contributions);

figure(3)
hold on
contributions=[];
ii_cont=0;
for ii=1:3  
    ii_cont=ii_cont+1;
    contributions(ii_cont).data=perVar(ii,:);
    switch ii
        case 1
            bar(ii,mean(perVar(ii,:)),'m')
            contributions(ii_cont).description='perCorr pre-odor';
        case 2
            bar(ii,mean(perVar(ii,:)),'b')
            contributions(ii_cont).description='perCorr odor';
        case 3
            bar(ii,mean(perVar(ii,:)),'r')
            contributions(ii_cont).description='perCorr outcome';
    end
            
    CIkin = bootci(1000, {@mean, perVar(ii,:)})';
    plot([ii ii],CIkin,'-k','LineWidth',2)
    plot((ii)*ones(1,6),perVar(ii,:),'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',6);
end
title('Percent variance explained by glm')
ylabel('% variance explained')
xticks([1 2 3])
xticklabels({'Pre-odor','Odor','Outcome'})
ylim([0 120])
pffft=1;

%Plot results of LDA analysis for glm


% 
% %Plot a bar graph showing the number of significant inclusions as a
% %predictive variable
% figure(3)
% hold on
% 
% disc_corr=[];
% 
% %Plot discriminant_correct
% for ii=1:5
%     these_dc=zeros(1,6);
%     these_dc=discriminant_correct(ii,:);
%     disc_corr(ii).data=these_dc;
%     bar(ii,mean(these_dc),'b')
%     thisCI = bootci(1000, {@mean, these_dc})';
%     plot([ii ii],thisCI,'-k')
%     plot((ii)*ones(1,6),these_dc,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
% end
% 
% these_dc=[];
% these_dc=discriminant_correct_shuffled;
% disc_corr(6).data=these_dc;
% bar(6,mean(these_dc),'b')
% thisCI = bootci(1000, {@mean, these_dc})';
% plot([6 6],thisCI,'-k')
% plot((6)*ones(1,length(these_dc)),these_dc,'o','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',4);
% 
% ylim([0 120])
% xticks([1 2 3 4 5 6])
% xticklabels({'dF/F','pdF/F all','pdF/F kinematics','pdF/F licks','pdF/F odors','dF/F shuffled'})
% title('Linear discriminant analysis of predicted dF/F')
% 
% disc_corr(1).description='dF/F';
% disc_corr(2).description='pdF/F all';
% disc_corr(3).description='pdF/F kinematics';
% disc_corr(4).description='pdF/F licks';
% disc_corr(5).description='pdF/F odors';
% disc_corr(6).description='dF/F shuffled';
% 
% 
% fprintf(1, ['\n\nttest or rnaksum test for inclusions\n'])
% [output_data] = drgMutiRanksumorTtest(disc_corr);
% 
% 
% 
% %kinematics
% bar(1,100*no_kinematics/6,'b')
% 
% %licks
% bar(2,100*no_licks/6,'b')
% 
% %Odors
% bar(4,100*no_odor/6,'b')
% 
% %Reinforcement history
% bar(6,100*no_rh/6,'b')
% 
% %Accuracy
% bar(7,100*no_accuracy/6,'b')
% 
% %perCorr
% bar(8,100*no_pc/6,'b')
% 
% ylabel('% included in glm')
% xticks([1 2 4 6 7 8])
% xticklabels({'Kinematics','Licks','Odors','Past','Accuracy','% correct'})
% ylim([0 120])