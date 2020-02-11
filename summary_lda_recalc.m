%This program generates a summary figure for the LDA analysis for Fig. 3d
%This uses the output produced by drgCaImAnBatchPerSessionReversalPerTrialLDA
close all
clear all

ldaPerCorr=[];
glm_lda=[];
glm_lda_ii=0;
pcorr_lda=[];
pcorr_lda_ii=0;


%mmPVG04
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_new_out_lda.mat';
load(file_name)

fprintf(1,['Reading ' file_name '\n\n'])
for jj=1:length(handles_sig.win)
    fprintf(1,[handles_sig.win(jj).description '\n'])
end
fprintf(1,'\n\n')


%We use <65 (low) and >=80 (high) for mmPVG04
%Note, drgCaImAnBatchPerSessionReversalPerTrialLDA saves for each percent
%window that worked the percent correct for the three windows and the
%percent correct for the shuffled trials

%<65
ii_l=0;
twins=[2 3 1]; %Shuffled is stored in 1, odor in 2 and reinforcement in 3
for ii=[1 2 3] %These are the windows for <65%
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
for ii=[1 2 3]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,1)=mean(handles_sig.win(ii+3).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,1);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(1)=ldaPerCorr(ii_l,1);
    pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+3).description;
end



%mmPVG02
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_lda_new_lda.mat';
load(file_name)

fprintf(1,['Reading ' file_name '\n\n'])
for jj=1:length(handles_sig.win)
    fprintf(1,[handles_sig.win(jj).description '\n'])
end
fprintf(1,'\n\n')

%We use >=65&<80 (low) and >80 (high) for mmPVG04
 
%<65
ii_l=0;
for ii=[1 2 3]
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
for ii=[1 2 3]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,2)=mean(handles_sig.win(ii+3).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,2);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(2)=ldaPerCorr(ii_l,2);
%     pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

%mmG7f09
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_lda_sum_lda.mat';
load(file_name)
 
fprintf(1,['Reading ' file_name '\n\n'])
for jj=1:length(handles_sig.win)
    fprintf(1,[handles_sig.win(jj).description '\n'])
end
fprintf(1,'\n\n')

%<65
ii_l=0;
for ii=[1 2 3]
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
for ii=[1 2 3]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,3)=mean(handles_sig.win(ii+6).discriminant_correct);
    glm_lda_ii=glm_lda_ii+1;
    glm_lda.data(glm_lda_ii)=ldaPerCorr(ii_l,3);
    glm_lda.time_win(glm_lda_ii)=twins(ii);
    glm_lda.pcorr_win(glm_lda_ii)=3;
    pcorr_lda(twins(ii)+3).data(3)=ldaPerCorr(ii_l,3);
%     pcorr_lda(twins(ii)+3).description=handles_sig.win(ii+4).description;
end

%mmPVG05
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19mmPVG05_Cerebellum_out_lda_sum_lda.mat';
load(file_name)


fprintf(1,['Reading ' file_name '\n\n'])
for jj=1:length(handles_sig.win)
    fprintf(1,[handles_sig.win(jj).description '\n'])
end
fprintf(1,'\n\n')

%<65
ii_l=0;
for ii=[1 2 3]
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
for ii=[1 2 3]
    ii_l=ii_l+1;
    ldaPerCorr(ii_l,4)=mean(handles_sig.win(ii+6).discriminant_correct);
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

%Enter the points

%mmPVG04 is GCaMP6f
plot([1 2 3],[ldaPerCorr(3,1) ldaPerCorr(1,1) ldaPerCorr(2,1)],'-ok','MarkerSize',10)
plot([5 6 7],[ldaPerCorr(6,1) ldaPerCorr(4,1) ldaPerCorr(5,1)],'-ok','MarkerSize',10)

%mmPVG02 is GCaMP6s
plot([1 2 3],[ldaPerCorr(3,2) ldaPerCorr(1,2) ldaPerCorr(2,2)],'-sk','MarkerSize',10)
plot([5 6 7],[ldaPerCorr(6,2) ldaPerCorr(4,2) ldaPerCorr(5,2)],'-sk','MarkerSize',10)

%mmG7f09 is GCaMP7f
plot([1 2 3],[ldaPerCorr(3,3) ldaPerCorr(1,3) ldaPerCorr(2,3)],'-xk','MarkerSize',10)
plot([5 6 7],[ldaPerCorr(6,3) ldaPerCorr(4,3) ldaPerCorr(5,3)],'-xk','MarkerSize',10)

%mmPVG05 is GCaMP6f
plot([1 2 3],[ldaPerCorr(3,4) ldaPerCorr(1,4) ldaPerCorr(2,4)],'-ok','MarkerSize',10)
plot([5 6 7],[ldaPerCorr(6,4) ldaPerCorr(4,4) ldaPerCorr(5,4)],'-ok','MarkerSize',10)


%CIs

%Learning
plot([1 1],CIldaPerCorr(:,3),'-k','LineWidth',3)

for ii=2:3
    plot([ii ii],CIldaPerCorr(:,ii-1),'-k','LineWidth',3)
end

%Proficient
plot([5 5],CIldaPerCorr(:,6),'-k','LineWidth',3)

for ii=6:7
    plot([ii ii],CIldaPerCorr(:,ii-2),'-k','LineWidth',3)
end

xticks([1 2 3 5 6 7])
xticklabels({'shuffled','odor','reinforcement','shuffled','odor','reinforcement'})
ylim([40 110])
text(1.5,105,'Learning')
text(5.5,105,'Proficient')
title('Percent correct predicted by linear discriminant analysis')

pffft=1;


