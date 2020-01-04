%This program generates a summary figure for the LDA analysis for Fig. 3
%This uses the output produced by drgCaImAnBatchPerSessionReversalPerTrialLDA
close all
clear all

ldaPerCorr=[];
glm_lda=[];
glm_lda_ii=0;
pcorr_lda=[];
pcorr_lda_ii=0;


%mmPVG04
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum/20180910_mmPVG04_Cerebellum_lda_sum.mat')
%20180910_mmPVG04_Cerebellum_lda

%We use <65 (low) and >=80 (high) for mmPVG04
%Note, drgCaImAnBatchPerSessionReversalPerTrialLDA saves for each percent
%window that worked the percent correct for the three windows and the
%percent correct for the shuffled trials

%<65
ii_l=0;
twins=[2 3 4 1];
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,1)=mean(handles_sig.win(ii).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,1);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=1;
    pcorr_lda(twins(ii)).data(1)=ldaPerCorr(ii_l,1);
    pcorr_lda(twins(ii)).description=handles_sig.win(ii).description;
end
 
%>=80
ii_l=3;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,1)=mean(handles_sig.win(ii+4).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,1);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(1)=ldaPerCorr(ii_l,1);
    pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

%mmPVG02
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum/20180515_18_mmPVG02_Cerebellum_out_lda_sum.mat')

%We use >=65&<80 (low) and >80 (high) for mmPVG04

%<65
ii_l=0;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,2)=mean(handles_sig.win(ii).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,2);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=1;
    pcorr_lda(twins(ii)).data(2)=ldaPerCorr(ii_l,2);
%     pcorr_lda(twins(ii)).description=handles_sig.win(ii).description;
end

%>=80
ii_l=3;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,2)=mean(handles_sig.win(ii+4).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,2);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(2)=ldaPerCorr(ii_l,2);
%     pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

%mmG7f09
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09_Cerebellum/20180702_05_mmG7f09-Cerebellum_lda_sum_lda.mat')

%<65
ii_l=0;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,3)=mean(handles_sig.win(ii).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,3);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=1;
    pcorr_lda(twins(ii)).data(3)=ldaPerCorr(ii_l,3);
%     pcorr_lda(twins(ii)).description=handles_sig.win(ii).description;
end

%>=80
ii_l=3;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,3)=mean(handles_sig.win(ii+8).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,3);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(3)=ldaPerCorr(ii_l,3);
%     pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

%mmPVG05
load('/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum/20181017_19mmPVG05_Cerebellum_out_lda_sum_lda.mat')

%<65
ii_l=0;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,4)=mean(handles_sig.win(ii).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,4);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=1;
    pcorr_lda(twins(ii)).data(4)=ldaPerCorr(ii_l,4);
%     pcorr_lda(twins(ii)).description=handles_sig.win(ii).description;
end

%>=80
ii_l=3;
for ii=[1 2 4]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,4)=mean(handles_sig.win(ii+8).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,4);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(4)=ldaPerCorr(ii_l,4);
%     pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

CIldaPerCorr=zeros(2,8);

for ii=1:6
     CIldaPerCorr(:,ii) = bootci(1000, {@mean, ldaPerCorr(ii,:)})';
end
    

%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_lda.data',glm_lda.time_win',glm_lda.pcorr_win',...
    'VariableNames',{'percent_correct_LDA','time_window','naive_proficient'});
mdl = fitglm(tbl,'percent_correct_LDA~time_window+naive_proficient+time_window*naive_proficient'...
    ,'CategoricalVars',[2,3])


% Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for for LDA percent correct\n'])
try
    [output_data] = drgMutiRanksumorTtest(pcorr_lda);
    fprintf(1, '\n\n')
catch
end


% 
% description{1}='Learning, 3-4 sec';
% description{2}='Learning, 4.5-6 sec';
% description{3}='Learning, 8.5-10 sec';
% description{4}='Learning, shuffled';
% description{5}='Proficient, 3-4 sec';
% description{6}='Proficient, 4.5-6 sec';
% description{7}='Proficient, 8.5-10 sec';
% description{8}='Proficient, shuffled';
% 
% fprintf(1, 'Tests of significance for difference in percent correct LDA\n')
% fprintf(1, ['Note: For shuffled trials we use all time points\n\n'])
% p_vals_LDA=[];
% no_LDA_pvals=0;
% for ii=1:8
%     for jj=ii+1:8
%         
%         no_LDA_pvals=no_LDA_pvals+1;
%         if (length(ldaPerCorr(ii,:))<4)||(length(ldaPerCorr(jj,:))<4)
%             %adtest does not work with n<4. In that case go the
%             %safe way ranksum
%             p_vals_LDA(no_LDA_pvals)=ranksum(ldaPerCorr(ii,:),ldaPerCorr(jj,:));
%             fprintf(1, ['p values ranksum for ' description{ii} ' vs. ' description{jj} ' =%d\n'],p_vals_LDA(no_LDA_pvals));
%         else
%             if (adtest(ldaPerCorr(ii,:))==1)||(adtest(ldaPerCorr(jj,:))==1)
%                 p_vals_LDA(no_LDA_pvals)=ranksum(ldaPerCorr(ii,:),ldaPerCorr(jj,:));
%                 fprintf(1, ['p values ranksum for ' description{ii} ' vs. ' description{jj} ' =%d\n'],p_vals_LDA(no_LDA_pvals));
%             else
%                 [h p_vals_LDA(no_LDA_pvals)]=ttest2(ldaPerCorr(ii,:),ldaPerCorr(jj,:));
%                 fprintf(1, ['p values t test for ' description{ii} ' vs. ' description{jj} ' =%d\n'],p_vals_LDA(no_LDA_pvals));
%             end
%         end
%         
%         
%     end
% end
%   
% pFDRLDA=drsFDRpval(p_vals_LDA);
% fprintf(1, ['\npFDR for significant difference percent correct  = %d\n\n'],pFDRLDA);

%Compute the CIs for each measurement


try
    close 5
catch
end

hFig5 = figure(5);
set(hFig5, 'units','normalized','position',[.07 .7 .7 .3])
hold on

%Learning
bar(1,mean(ldaPerCorr(3,:)),'b')
plot([1 1],CIldaPerCorr(:,3),'-k','LineWidth',3)

for ii=2:3
    bar(ii,mean(ldaPerCorr(ii-1,:)),'b')
    plot([ii ii],CIldaPerCorr(:,ii-1),'-k','LineWidth',3)
end

%Proficient
bar(5,mean(ldaPerCorr(6,:)),'r')
plot([5 5],CIldaPerCorr(:,6),'-k','LineWidth',3)

for ii=6:7
    bar(ii,mean(ldaPerCorr(ii-2,:)),'r')
    plot([ii ii],CIldaPerCorr(:,ii-2),'-k','LineWidth',3)
end

xticks([1 2 3 5 6 7])
xticklabels({'shuffled','odor','reinforcement','shuffled','odor','reinforcement'})
ylim([40 110])
text(1.5,105,'Learning')
text(5.5,105,'Proficient')
title('Percent correct predicted by linear discriminant analysis')

pffft=1;


