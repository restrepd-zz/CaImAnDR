%This scrript calculates the statistics for Figs. 7aii,bii in Ming's
%manuscript
glm_flt_ii=0;
glm_flt=[];
rst_ii=0;
rst_data=[];

figure(1)
hold on
 

% %Control, no CNO
flts=[68.8
    82.9
    86
    81.9
    90.8
    87];
flts=(flts-50)/50;
mean_prob=mean(flts);
CIprob = bootci(1000, {@mean, flts})';;
bar(1,mean_prob,'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1])
plot([1 1],CIprob,'-k','LineWidth',2)
plot(ones(1,length(flts)),flts,'ok')

glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
glm_flt_ii=glm_flt_ii+length(flts);

rst_ii=rst_ii+1;
rst_data(rst_ii).data=flts;
rst_data(rst_ii).description='Control, no CNO';

%Control, CNO
flts=[69.1
    73.9
    79.6
    71.2
    89.6
    84.7];
flts=(flts-50)/50;
mean_prob=mean(flts);
CIprob = bootci(1000, {@mean, flts})';
bar(2,mean_prob,'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7])
plot([2 2],CIprob,'-k','LineWidth',2)
plot(2*ones(1,length(flts)),flts,'ok')

glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
glm_flt_ii=glm_flt_ii+length(flts);

rst_ii=rst_ii+1;
rst_data(rst_ii).data=flts;
rst_data(rst_ii).description='Control, CNO';

% %hM4D, no CNO

flts=[77.6
    72.4
    88.5
    80.8
    86.1
    83.8];
flts=(flts-50)/50;
mean_prob=mean(flts);
CIprob = bootci(1000, {@mean, flts})';
bar(4,mean_prob,'b')
plot([4 4],CIprob,'-k','LineWidth',2)
plot(4*ones(1,length(flts)),flts,'ok')

glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=0;
glm_flt_ii=glm_flt_ii+length(flts);

rst_ii=rst_ii+1;
rst_data(rst_ii).data=flts;
rst_data(rst_ii).description='hM4Di, no CNO';

%hM4D, CNO

flts=[60.4
    56.4
    61.1
    56.6
    66.3
    66.1];
flts=(flts-50)/50;
mean_prob=mean(flts);
CIprob = bootci(1000, {@mean, flts})';
bar(5,mean_prob,'r')
plot([5 5],CIprob,'-k','LineWidth',2)
plot(5*ones(1,length(flts)),flts,'ok')

glm_flt.data(glm_flt_ii+1:glm_flt_ii+length(flts))=flts;
glm_flt.genotype(glm_flt_ii+1:glm_flt_ii+length(flts))=2;
glm_flt.CNO(glm_flt_ii+1:glm_flt_ii+length(flts))=1;
glm_flt_ii=glm_flt_ii+length(flts);

rst_ii=rst_ii+1;
rst_data(rst_ii).data=flts;
rst_data(rst_ii).description='HM4Di, CNO';

xticks([1 2 4 5])
ylim([0 100])
% ylim([0 1])
xticklabels({'Ctrl, noCNO','Ctrl, CNO','hM4d, noCNO','hM4d, CNO'})
title('S-')
ylabel('Frequency (Hz)')

suptitle('Average percent correct')

%Perform the glm
fprintf(1, '\n\nglm for percent correct\n')
tbl = table(glm_flt.data',glm_flt.genotype',glm_flt.CNO',...
    'VariableNames',{'frequency','genotype','CNO'});
mdl = fitglm(tbl,'frequency~genotype+CNO+genotype*CNO'...
    ,'CategoricalVars',[2,3])

%Do the ranksum/t-test
fprintf(1, '\n\nRanksum or t-test p values for  percent correct\n')
[output_data] = drgMutiRanksumorTtest(rst_data);
