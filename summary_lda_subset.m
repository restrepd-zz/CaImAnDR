%This program generates a summary figure for the LDA analysis for Fig. 3d
%This uses the output produced by drgCaImAnBatchPerSessionReversalPerTrialLDA
close all
clear all

mean_percent_correct_shuffled=[];
mean_percent_correct=[];
ii=0;
glm_lda=[];
glm_ii=0;
pcorr_lda=[];
pcorr_lda_ii=0;


%mmPVG04
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG04/20180910_mmPVG04_Cerebellum_new_analysis/20180910_mmPVG04_Cerebellum_new_out_ldasub.mat';
load(file_name)

this_num_comps=handles_par(3).timewin(1).num_comps;
this_mean_pc_shuffled=handles_par(3).timewin(1).mean_percent_correct_shuffled;
this_mean_pc=handles_par(3).timewin(1).mean_percent_correct;

ii=ii+1;

mean_percent_correct_shuffled(ii,1)=this_mean_pc_shuffled(this_num_comps==1);
mean_percent_correct_shuffled(ii,2)=this_mean_pc_shuffled(this_num_comps==5);
mean_percent_correct_shuffled(ii,3)=this_mean_pc_shuffled(this_num_comps==100);
mean_percent_correct_shuffled(ii,4)=this_mean_pc_shuffled(end);

mean_percent_correct(ii,1)=this_mean_pc(this_num_comps==1);
mean_percent_correct(ii,2)=this_mean_pc(this_num_comps==5);
mean_percent_correct(ii,3)=this_mean_pc(this_num_comps==100);
mean_percent_correct(ii,4)=this_mean_pc(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=this_mean_pc_shuffled;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=[1:length(this_mean_pc_shuffled)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=ones(1,length(this_mean_pc_shuffled));
glm_ii=glm_ii+length(this_mean_pc_shuffled);

pcorr_lda(1).data(ii)=this_mean_pc_shuffled(this_num_comps==1);
pcorr_lda(1).description='1 ROI shuffled';
pcorr_lda(2).data(ii)=this_mean_pc_shuffled(this_num_comps==5);
pcorr_lda(2).description='5 ROIs shuffled';
pcorr_lda(3).data(ii)=this_mean_pc_shuffled(this_num_comps==100);
pcorr_lda(3).description='100 ROIs shuffled';
pcorr_lda(4).data(ii)=this_mean_pc_shuffled(end);
pcorr_lda(4).description='All ROIs shuffled';

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc))=this_mean_pc;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc))=[1:length(this_mean_pc)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc))=zeros(1,length(this_mean_pc));
glm_ii=glm_ii+length(this_mean_pc);

pcorr_lda(5).data(ii)=this_mean_pc(this_num_comps==1);
pcorr_lda(5).description='1 ROI';
pcorr_lda(6).data(ii)=this_mean_pc(this_num_comps==5);
pcorr_lda(6).description='5 ROIs';
pcorr_lda(7).data(ii)=this_mean_pc(this_num_comps==100);
pcorr_lda(7).description='100 ROIs';
pcorr_lda(8).data(ii)=this_mean_pc(end);
pcorr_lda(8).description='All ROIs';

%mmPVG02
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG02/20180515_mmPVG02_Cerebellum_new_analysis/20180515_18_mmPVG02_Cerebellum_out_lda_new_ldasub.mat';
load(file_name)

this_num_comps=handles_par(3).timewin(1).num_comps;
this_mean_pc_shuffled=handles_par(3).timewin(1).mean_percent_correct_shuffled;
this_mean_pc=handles_par(3).timewin(1).mean_percent_correct;

ii=ii+1;

mean_percent_correct_shuffled(ii,1)=this_mean_pc_shuffled(this_num_comps==1);
mean_percent_correct_shuffled(ii,2)=this_mean_pc_shuffled(this_num_comps==5);
mean_percent_correct_shuffled(ii,3)=this_mean_pc_shuffled(this_num_comps==100);
mean_percent_correct_shuffled(ii,4)=this_mean_pc_shuffled(end);

mean_percent_correct(ii,1)=this_mean_pc(this_num_comps==1);
mean_percent_correct(ii,2)=this_mean_pc(this_num_comps==5);
mean_percent_correct(ii,3)=this_mean_pc(this_num_comps==100);
mean_percent_correct(ii,4)=this_mean_pc(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=this_mean_pc_shuffled;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=[1:length(this_mean_pc_shuffled)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=ones(1,length(this_mean_pc_shuffled));
glm_ii=glm_ii+length(this_mean_pc_shuffled);

pcorr_lda(1).data(ii)=this_mean_pc_shuffled(this_num_comps==1);
pcorr_lda(2).data(ii)=this_mean_pc_shuffled(this_num_comps==5);
pcorr_lda(3).data(ii)=this_mean_pc_shuffled(this_num_comps==100);
pcorr_lda(4).data(ii)=this_mean_pc_shuffled(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc))=this_mean_pc;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc))=[1:length(this_mean_pc)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc))=zeros(1,length(this_mean_pc));
glm_ii=glm_ii+length(this_mean_pc);

pcorr_lda(5).data(ii)=this_mean_pc(this_num_comps==1);
pcorr_lda(6).data(ii)=this_mean_pc(this_num_comps==5);
pcorr_lda(7).data(ii)=this_mean_pc(this_num_comps==100);
pcorr_lda(8).data(ii)=this_mean_pc(end);

%mmG7f09
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmG7f09/20180702_mmG7f09-Cerebellum new analysis/20180702_05_mmG7f09-Cerebellum_lda_sum_ldasub.mat';
load(file_name)
 
this_num_comps=handles_par(3).timewin(1).num_comps;
this_mean_pc_shuffled=handles_par(3).timewin(1).mean_percent_correct_shuffled;
this_mean_pc=handles_par(3).timewin(1).mean_percent_correct;

ii=ii+1;

mean_percent_correct_shuffled(ii,1)=this_mean_pc_shuffled(this_num_comps==1);
mean_percent_correct_shuffled(ii,2)=this_mean_pc_shuffled(this_num_comps==5);
mean_percent_correct_shuffled(ii,3)=this_mean_pc_shuffled(this_num_comps==100);
mean_percent_correct_shuffled(ii,4)=this_mean_pc_shuffled(end);

mean_percent_correct(ii,1)=this_mean_pc(this_num_comps==1);
mean_percent_correct(ii,2)=this_mean_pc(this_num_comps==5);
mean_percent_correct(ii,3)=this_mean_pc(this_num_comps==100);
mean_percent_correct(ii,4)=this_mean_pc(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=this_mean_pc_shuffled;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=[1:length(this_mean_pc_shuffled)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=ones(1,length(this_mean_pc_shuffled));
glm_ii=glm_ii+length(this_mean_pc_shuffled);

pcorr_lda(1).data(ii)=this_mean_pc_shuffled(this_num_comps==1);
pcorr_lda(2).data(ii)=this_mean_pc_shuffled(this_num_comps==5);
pcorr_lda(3).data(ii)=this_mean_pc_shuffled(this_num_comps==100);
pcorr_lda(4).data(ii)=this_mean_pc_shuffled(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc))=this_mean_pc;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc))=[1:length(this_mean_pc)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc))=zeros(1,length(this_mean_pc));
glm_ii=glm_ii+length(this_mean_pc);

pcorr_lda(5).data(ii)=this_mean_pc(this_num_comps==1);
pcorr_lda(6).data(ii)=this_mean_pc(this_num_comps==5);
pcorr_lda(7).data(ii)=this_mean_pc(this_num_comps==100);
pcorr_lda(8).data(ii)=this_mean_pc(end);


%mmPVG05
file_name='/Users/restrepd/Documents/Projects/MOM slidebook/mmPVG05/20181017_mmPVG05_Cerebellum new analysis/20181017_19mmPVG05_Cerebellum_out_lda_sum_ldasub.mat';
load(file_name)


this_num_comps=handles_par(3).timewin(1).num_comps;
this_mean_pc_shuffled=handles_par(3).timewin(1).mean_percent_correct_shuffled;
this_mean_pc=handles_par(3).timewin(1).mean_percent_correct;

ii=ii+1;

mean_percent_correct_shuffled(ii,1)=this_mean_pc_shuffled(this_num_comps==1);
mean_percent_correct_shuffled(ii,2)=this_mean_pc_shuffled(this_num_comps==5);
mean_percent_correct_shuffled(ii,3)=this_mean_pc_shuffled(this_num_comps==100);
mean_percent_correct_shuffled(ii,4)=this_mean_pc_shuffled(end);

mean_percent_correct(ii,1)=this_mean_pc(this_num_comps==1);
mean_percent_correct(ii,2)=this_mean_pc(this_num_comps==5);
mean_percent_correct(ii,3)=this_mean_pc(this_num_comps==100);
mean_percent_correct(ii,4)=this_mean_pc(end);

glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=this_mean_pc_shuffled;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=[1:length(this_mean_pc_shuffled)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc_shuffled))=ones(1,length(this_mean_pc_shuffled));
glm_ii=glm_ii+length(this_mean_pc_shuffled);

pcorr_lda(1).data(ii)=this_mean_pc_shuffled(this_num_comps==1);
pcorr_lda(2).data(ii)=this_mean_pc_shuffled(this_num_comps==5);
pcorr_lda(3).data(ii)=this_mean_pc_shuffled(this_num_comps==100);
pcorr_lda(4).data(ii)=this_mean_pc_shuffled(end);


glm_lda.data(glm_ii+1:glm_ii+length(this_mean_pc))=this_mean_pc;
glm_lda.subset(glm_ii+1:glm_ii+length(this_mean_pc))=[1:length(this_mean_pc)];
glm_lda.shuffled(glm_ii+1:glm_ii+length(this_mean_pc))=zeros(1,length(this_mean_pc));
glm_ii=glm_ii+length(this_mean_pc);

pcorr_lda(5).data(ii)=this_mean_pc(this_num_comps==1);
pcorr_lda(6).data(ii)=this_mean_pc(this_num_comps==5);
pcorr_lda(7).data(ii)=this_mean_pc(this_num_comps==100);
pcorr_lda(8).data(ii)=this_mean_pc(end);


%Perform the glm for LDA percent correct 
fprintf(1, ['\n\nglm for LDA percent correct\n'])
tbl = table(glm_lda.data',glm_lda.subset',glm_lda.shuffled',...
    'VariableNames',{'percent_correct_LDA','subset','shuffled'});
mdl = fitglm(tbl,'percent_correct_LDA~subset+shuffled+subset*shuffled'...
    ,'CategoricalVars',[2,3])

 
% Do ranksum/t test
fprintf(1, ['\n\nRanksum or t-test p values for for LDA percent correct\n'])
try
    [output_data] = drgMutiRanksumorTtest(pcorr_lda);
    fprintf(1, '\n\n')
catch
end

CIldaPerCorr=zeros(2,4);
CIldaPerCorrSh=zeros(2,4);

for ii=1:4
     CIldaPerCorr(:,ii) = bootci(1000, {@mean, mean_percent_correct(:,ii)})';
     CIldaPerCorrSh(:,ii) = bootci(1000, {@mean, mean_percent_correct_shuffled(:,ii)})';
end
    

try
    close 5
catch
end

hFig5 = figure(5);
set(hFig5, 'units','normalized','position',[.07 .7 .7 .3])
hold on

%Draw the bars
for ii=1:4
    bar(1+(ii-1)*3,mean(mean_percent_correct_shuffled(:,ii)),'b')
    bar(2+(ii-1)*3,mean(mean_percent_correct(:,ii)),'r')
end

%Plot the points
for ii=1:4
    %mmPVG04 is GCaMP6f
    plot([1+(ii-1)*3 2+(ii-1)*3],[mean_percent_correct_shuffled(1,ii) mean_percent_correct(1,ii)],'-ok','MarkerSize',10)
      %mmPVG02 is GCaMP6s
    plot([1+(ii-1)*3 2+(ii-1)*3],[mean_percent_correct_shuffled(2,ii) mean_percent_correct(2,ii)],'-sk','MarkerSize',10)
     %mmG7f09 is GCaMP7f
    plot([1+(ii-1)*3 2+(ii-1)*3],[mean_percent_correct_shuffled(3,ii) mean_percent_correct(3,ii)],'-xk','MarkerSize',10)
    %mmPVG05 is GCaMP6f
    plot([1+(ii-1)*3 2+(ii-1)*3],[mean_percent_correct_shuffled(4,ii) mean_percent_correct(4,ii)],'-ok','MarkerSize',10)
end


%CIs
for ii=1:4
    plot([1+(ii-1)*3 1+(ii-1)*3],CIldaPerCorrSh(:,ii),'-k','LineWidth',3)
    plot([2+(ii-1)*3 2+(ii-1)*3],CIldaPerCorr(:,ii),'-k','LineWidth',3)
end


xticks([1.5 4.5 7.5 10.5])
xticklabels({'1','5','100','All'})
ylim([40 110])
title('Percent correct predicted by linear discriminant analysis')

pffft=1;


