%%supfigDREADDs is used for Supplementary Fig. 13a
glm_pcorr=[];
glm_ii=0

%4d11 forward
glm_pcorr.data(1:4)=[75.26 56.10 56.60 70.83];
glm_pcorr.genotype(1:4)=[2 2 2 2];
glm_pcorr.fwd_rev(1:4)=[1 1 1 1];


%4d11 reversed
glm_pcorr.data(5:9)=[59.55 51.85 73.39 71.15 53.76];
glm_pcorr.genotype(5:9)=[2 2 2 2 2];
glm_pcorr.fwd_rev(5:9)=[2 2 2 2 2];

%4d12 forward
glm_pcorr.data(10:13)=[62.32 58 63.53 65.25];
glm_pcorr.genotype(10:13)=[2 2 2 2];
glm_pcorr.fwd_rev(10:13)=[1 1 1 1];

%4d12 reversed
glm_pcorr.data(14:18)=[56.92 63.83 70.37 68.63 54.29];
glm_pcorr.genotype(14:18)=[2 2 2 2 2];
glm_pcorr.fwd_rev(14:18)=[2 2 2 2 2];

pct_data(1).data=[75.26 56.10 56.60 70.83 62.32 58 63.53 65.25];
pct_data(1).description='h4DMi forward';

pct_data(2).data=[59.55 51.85 73.39 71.15 53.76 56.92 63.83 70.37 68.63 54.29];
pct_data(2).description='h4DMi reversed';

%ctrl12 forward
glm_pcorr.data(19:22)=[89.77 86.14 96.33 81.91];
glm_pcorr.genotype(19:22)=[1 1 1 1];
glm_pcorr.fwd_rev(19:22)=[1 1 1 1];

%ctrl12 reversed
glm_pcorr.data(23:27)=[51.85 93.20 92.31 98.04 99.22];
glm_pcorr.genotype(23:27)=[1 1 1 1 1];
glm_pcorr.fwd_rev(23:27)=[2 2 2 2 2];

%ctrl13 forward
glm_pcorr.data(28:31)=[89.10 89.91 80.20 85.12];
glm_pcorr.genotype(28:31)=[1 1 1 1];
glm_pcorr.fwd_rev(28:31)=[1 1 1 1];

%ctrl12 reversed
glm_pcorr.data(32:36)=[43.12 70.19 83.59 76.22 80.53];
glm_pcorr.genotype(32:36)=[1 1 1 1 1 ];
glm_pcorr.fwd_rev(32:36)=[2 2 2 2 2];

pct_data(3).data=[89.77 86.14 96.33 81.91 89.10 89.91 80.20 85.12];
pct_data(3).description='Control forward';

pct_data(4).data=[51.85 93.20 92.31 98.04 99.22 43.12 70.19 83.59 76.22 80.53];
pct_data(4).description='Control reversed';

fprintf(1, ['\n\nglm for percent correct for supp fig 13\n'])
tbl = table(glm_pcorr.data',glm_pcorr.genotype',glm_pcorr.fwd_rev',...
    'VariableNames',{'frequency','genotype','fwdrev'});
mdl = fitglm(tbl,'frequency~genotype+fwdrev+genotype*fwdrev'...
    ,'CategoricalVars',[2,3])

%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for supp fig 11\n'])
[output_data] = drgMutiRanksumorTtest(pct_data);